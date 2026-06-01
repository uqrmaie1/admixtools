#!/usr/bin/env bash
# Tier 4 REAL third-party data integration (master plan B series), per research
# Fits the bundled LEGOFIT archaic .daf data (Altai Neanderthal,
# Denisovan, 3 modern genomes -> one "modern" pop) through a graph_to_lgo-
# generated model and reads the result back.
#
#   B-1: graph_to_lgo 3-leaf tree (m,(n,d)) -> legofit -1 -> read_legofit_output,
#        with REAL parameter families (not "unknown") and a populated convergence
#        attribute.
#   B-2: moving-block bootstrap (tabpat -1 -b -r) -> legofit per rep ->
#        flatfile.py -> bootci.py -> read_legofit_bootstrap, real CI tibble.
#
# This is an INTEGRATION check on real data, NOT a paper reproduction: the
# bundled .daf is a tiny chr22 unit-test region (~991 aligned sites), so the
# fitted numbers are not biologically meaningful (that is R-1/R-2). The model
# MUST be graph_to_lgo-generated so the T_<node>/twoN_one names classify into
# real families (a hand-written .lgo with arbitrary names classifies as
# "unknown").
#
# Requires LEGOFIT built ($LEGOFIT) and the bundled test/*.daf alongside it.
set -euo pipefail

: "${LEGOFIT:?must point to legofit src dir}"
for b in tabpat legofit bootci.py flatfile.py; do
  [[ -e "$LEGOFIT/$b" ]] || { echo "no $b in $LEGOFIT"; exit 1; }
done
DAF="$(cd "$LEGOFIT/../test" && pwd)"
for f in altai denisova Mgenomes3; do
  [[ -f "$DAF/$f.daf" ]] || { echo "missing $DAF/$f.daf (bundled archaic data)"; exit 1; }
done
REPO="${REPO:-$(cd "$(dirname "$0")/../.." && pwd)}"
WORK="$(cd "$(dirname "$0")" && pwd)/work-archaic"
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"

fails=0
pass () { echo "  PASS $1"; }
fail () { echo "  FAIL $1"; fails=$((fails+1)); }

echo "==> $("$LEGOFIT/legofit" --version 2>&1 | grep -m1 version | tr -s ' ')  real archaic .daf (n=Altai, d=Denisova, m=3 moderns)"

# graph_to_lgo model: ((n,d),m), moderns outgroup. Leaf segment names n/d/m must
# match the tabpat population labels. Pure tree, so time_handling="free".
Rscript -e "
suppressMessages(pkgload::load_all('$REPO', quiet=TRUE))
g <- tibble::tribble(
  ~from, ~to, ~type,     ~weight,
  'R',   'm', 'normal',  NA_real_,
  'R',   'nd','normal',  NA_real_,
  'nd',  'n', 'normal',  NA_real_,
  'nd',  'd', 'normal',  NA_real_)
writeLines(graph_to_lgo(g, time_handling='free', validate=FALSE), '$WORK/archaic.lgo')
" || fail "B-1: graph_to_lgo export"

# ---- B-1: fit the real data, read it back with real families ----------------
"$LEGOFIT/tabpat" -1 n="$DAF/altai.daf" d="$DAF/denisova.daf" m="$DAF/Mgenomes3.daf" \
  > archaic.opf 2>/dev/null || fail "B-1: tabpat"
"$LEGOFIT/legofit" -1 -d 1e-3 -p 8 archaic.lgo archaic.opf > archaic.legofit 2>/dev/null \
  || fail "B-1: legofit"
Rscript -e "
suppressMessages(pkgload::load_all('$REPO', quiet=TRUE))
r <- read_legofit_output('$WORK/archaic.legofit')
stopifnot(nrow(r) > 0)
stopifnot(!any(r\$family == 'unknown'))                 # real families, not unknown
stopifnot('time' %in% r\$family, 'twoN_one' %in% r\$family)
conv <- attr(r, 'fit_convergence')\$status
cat(sprintf('  B-1 read %d params; families={%s}; convergence=%s\n',
            nrow(r), paste(unique(r\$family), collapse=','), conv))
stopifnot(!is.na(conv))
" && pass "B-1 real archaic fit reads back with real families" || fail "B-1 read-back"

# ---- B-2: moving-block bootstrap -> bootci -> read_legofit_bootstrap ----------
NB=20
"$LEGOFIT/tabpat" -1 -b 50 -r "$NB" -f boot \
  n="$DAF/altai.daf" d="$DAF/denisova.daf" m="$DAF/Mgenomes3.daf" \
  > b2_main.opf 2>/dev/null || fail "B-2: tabpat bootstrap"
"$LEGOFIT/legofit" -1 -d 1e-3 -p 8 archaic.lgo b2_main.opf > b2_main.legofit 2>/dev/null \
  || fail "B-2: legofit main"
for j in $(seq 0 $((NB-1))); do
  "$LEGOFIT/legofit" -1 -d 1e-3 -p 8 archaic.lgo "boot${j}.opf" > "b2_boot${j}.legofit" 2>/dev/null &
  if (( j % 8 == 7 )); then wait; fi
done
wait
python3 "$LEGOFIT/flatfile.py" b2_main.legofit $(for j in $(seq 0 $((NB-1))); do echo "b2_boot${j}.legofit"; done) \
  > b2.flat 2>/dev/null || fail "B-2: flatfile.py"
python3 "$LEGOFIT/bootci.py" b2.flat > b2.bootci 2>/dev/null || fail "B-2: bootci.py"
Rscript -e "
suppressMessages(pkgload::load_all('$REPO', quiet=TRUE))
r <- read_legofit_bootstrap('$WORK/b2.bootci')
stopifnot(c('parameter','family','point_estimate','lo','hi') %in% names(r))
stopifnot(!any(r\$family == 'unknown'))
stopifnot(all(r\$lo <= r\$hi))                            # well-ordered intervals
cat(sprintf('  B-2 CI tibble: %d params, families={%s}\n',
            nrow(r), paste(unique(r\$family), collapse=',')))
" && pass "B-2 real bootstrap reads back as a CI tibble with real families" || fail "B-2 bootstrap"

echo
if [[ "$fails" -eq 0 ]]; then echo "==> ALL ARCHAIC CHECKS PASSED (B-1, B-2)"; exit 0
else echo "==> $fails ARCHAIC CHECK(S) FAILED"; exit 1; fi
