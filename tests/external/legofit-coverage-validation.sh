#!/usr/bin/env bash
# Tier 3 honest bootstrap COVERAGE for the LEGOFIT bridge (master plan S-2).
#
# This replaces the retired self-bracketing assertion, which was
# vacuous when it passed and flaky when it failed. Here we draw K INDEPENDENT
# datasets from a known truth, build a parametric bootstrap CI for each, and ask
# how often the KNOWN TRUTH falls inside the CI, aggregated over K runs and the
# identifiable parameters. A real coverage check, not a tautology.
#
# Caps: singletons on; multi-admix uses fixed twoN
# (one=1) so there is no time/twoN scale degeneracy; the structurally
# non-identifiable admix-event time T_admix_m is EXCLUDED. Identifiable
# parameters scored: T_R, T_A, T_B (tree split times) and m_m (mixFrac).
#
# Slow (K*(B+1) fits); run manually. Set $LEGOFIT to the src/ dir.
set -euo pipefail

: "${LEGOFIT:?must point to legofit src dir (legofit, legosim, bootci.py)}"
for b in legosim legofit bootci.py flatfile.py; do
  [[ -e "$LEGOFIT/$b" ]] || { echo "no $b in $LEGOFIT"; exit 1; }
done
REPO="${REPO:-$(cd "$(dirname "$0")/../.." && pwd)}"
FIX="$REPO/tests/testthat/fixtures"
WORK="$(cd "$(dirname "$0")" && pwd)/work-coverage"
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"

K="${COVERAGE_K:-12}"        # independent datasets
B="${COVERAGE_B:-20}"        # parametric-bootstrap replicates per dataset
LGO="$FIX/multi-admix.lgo"
echo "==> S-2 coverage: $K datasets x $B bootstrap reps on multi-admix (singletons on)"

# Truth = the model's own free values (legosim draws from it).
Rscript -e "
lgo <- readLines('$LGO')
re  <- '^(?:time|twoN|mixFrac) free[[:space:]]+([A-Za-z0-9_.]+)[[:space:]]*=[[:space:]]*([-0-9.eE+]+)'
m   <- Filter(function(x) length(x)==3, regmatches(lgo, regexec(re, lgo)))
tr  <- setNames(as.numeric(vapply(m,'[',character(1),3)), vapply(m,'[',character(1),2))
writeLines(sprintf('%s %s', names(tr), tr), '$WORK/truth.txt')
"

# One dataset: draw a 'main' + B bootstrap replicates, fit all, bootci -> CI file.
one_run () {
  local k="$1" d="$WORK/run_$1"; mkdir -p "$d"
  local j
  for j in $(seq 0 "$B"); do
    "$LEGOFIT/legosim" -1 -i 4000 "$LGO" > "$d/rep_$j.opf" 2>/dev/null
  done
  for j in $(seq 0 "$B"); do
    "$LEGOFIT/legofit" -1 -d 1e-3 -p 8 "$LGO" "$d/rep_$j.opf" > "$d/rep_$j.legofit" 2>/dev/null &
    if (( j % 8 == 7 )); then wait; fi
  done
  wait
  # flatfile: main fit (rep_0) FIRST, then the bootstrap reps.
  python3 "$LEGOFIT/flatfile.py" $(for j in $(seq 0 "$B"); do echo "$d/rep_$j.legofit"; done) \
    > "$d/flat" 2>/dev/null
  python3 "$LEGOFIT/bootci.py" "$d/flat" > "$d/bootci" 2>/dev/null
}

for k in $(seq 1 "$K"); do
  one_run "$k"
  echo "  run $k/$K done"
done

# Aggregate coverage over runs x identifiable params.
Rscript -e "
suppressMessages(pkgload::load_all('$REPO', quiet=TRUE))
truth <- read.table('$WORK/truth.txt'); truth <- setNames(truth\$V2, truth\$V1)
ident <- c('T_R','T_A','T_B','m_m')              # exclude non-identifiable T_admix_m
hits <- 0L; total <- 0L; bad <- 0L; relhw <- numeric(0)
for (k in 1:$K) {
  f <- sprintf('$WORK/run_%d/bootci', k)
  r <- tryCatch(read_legofit_bootstrap(f), error=function(e) NULL)
  if (is.null(r)) { bad <- bad + 1L; next }
  v <- r[match(ident, r\$parameter), ]
  if (any(r\$family[r\$parameter %in% ident] == 'unknown')) stop('family unknown leaked')
  in_ci <- truth[ident] >= v\$lo & truth[ident] <= v\$hi
  hits  <- hits + sum(in_ci, na.rm=TRUE); total <- total + sum(!is.na(in_ci))
  relhw <- c(relhw, ((v\$hi - v\$lo)/2) / abs(v\$point_estimate))   # informativeness
}
cov <- hits/total; med_relhw <- median(relhw, na.rm=TRUE)
cat(sprintf('  coverage = %d/%d = %.1f%% over %d runs x %d identifiable params (%d unreadable)\n',
            hits, total, 100*cov, $K, length(ident), bad))
cat(sprintf('  median relative CI half-width = %.2f (informative if well below 1)\n', med_relhw))
# Two guards together make this non-vacuous, unlike the retired self-bracketing
# check. (1) High coverage of KNOWN TRUTH over INDEPENDENT draws: small-B
# percentile CIs over-cover, so we require >=0.80, not an exact 0.95. (2) The
# intervals must be INFORMATIVE (not trivially wide), else high coverage is
# meaningless: the median relative half-width must be well below 1, i.e. the CI
# is tighter than the estimate itself.
stopifnot(bad == 0, cov >= 0.80, med_relhw < 0.75)
cat('  S-2 PASS: known truth reliably bracketed by informative CIs over independent draws\n')
" && echo "==> S-2 COVERAGE PASSED" || { echo "==> S-2 COVERAGE FAILED"; exit 1; }
