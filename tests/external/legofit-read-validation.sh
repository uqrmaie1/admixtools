#!/usr/bin/env bash
# Tier 2 functional validation for PR delta (read direction).
# Per at2-legofit-integration-pr-delta-test-plan.md.
#
# Requires LEGOFIT 1.87+ installed locally. Set $LEGOFIT to the src/ dir.
# bootci.py ships with LEGOFIT.
#
# Outputs vendored fixtures into tests/testthat/fixtures/ that can then be
# read by tests/testthat/test-read_legofit_tier1.R.
#
# Run from the package root:
#   LEGOFIT=/path/to/legofit/src tests/external/legofit-read-validation.sh

set -euo pipefail

: "${LEGOFIT:?must point to legofit src dir containing legofit, legosim, bootci.py}"
[[ -x "$LEGOFIT/legofit" ]]   || { echo "no legofit binary in $LEGOFIT"; exit 1; }
[[ -x "$LEGOFIT/legosim" ]]   || { echo "no legosim binary in $LEGOFIT"; exit 1; }
[[ -f "$LEGOFIT/bootci.py" ]] || { echo "no bootci.py in $LEGOFIT"; exit 1; }

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIXTURES="$REPO_ROOT/tests/testthat/fixtures"
WORK="$REPO_ROOT/tests/external/work"
mkdir -p "$WORK"

pass () { echo "  PASS $1"; }
fail () { echo "  FAIL $1"; exit 1; }

echo "==> T2.1a: end to end round trip on ourex1 (no admix, 2 free params)"

cd "$WORK"
# Use the bundled ex1.opf since it pairs with the ourex1 topology
"$LEGOFIT/legofit" --threads 1 -1 -d 0 \
  "$FIXTURES/ourex1.lgo" "$LEGOFIT/../exmpl/ex1/ex1.opf" \
  > t21-ourex1.legofit 2>&1
grep -q "DiffEv reached_goal" t21-ourex1.legofit \
  || fail "T2.1a: ourex1 did not converge"
pass "T2.1a: ourex1 fit converged"

Rscript -e "
pkgload::load_all('$REPO_ROOT', quiet = TRUE)
result <- suppressMessages(read_legofit_output(
  '$WORK/t21-ourex1.legofit',
  graph = make_ourex1_graph()
))
xy  <- result[result\$to == 'xy',  ]
xyz <- attr(result, 'node_times')[['xyz']]
# Ground truth from \$LEGOFIT/exmpl/ex1/true.lgo: Txy=1.0, Txyz=2.0
stopifnot(nrow(xy) > 0)
stopifnot(abs(xy\$time[1] - 1.0) < 0.2)   # T_xy fitted close to true 1.0
stopifnot(abs(xyz       - 2.0) < 0.2)     # T_xyz fitted close to true 2.0
cat('ourex1 round trip ok: T_xy=', xy\$time[1], 'T_xyz=', xyz, '\n')
" || fail "T2.1a: ourex1 read back did not match ground truth"
pass "T2.1a: ourex1 reader recovered fitted times within tolerance"

echo
echo "==> T2.1b: end to end round trip on ooa-admix (4 leaves, 1 admix, 5 free params)"
echo "  Real-world OOA archaic introgression scenario. Replacement for the"
echo "  PR gamma minimal.lgo fixture which cannot be fit (one-leaf graph)."

cd "$WORK"
# graph_to_lgo -> legosim -> legofit -> read_legofit_output round trip
Rscript -e "
pkgload::load_all('$REPO_ROOT', quiet = TRUE)
g <- make_ooa_admix_graph()
graph_to_lgo(g, file = '$WORK/t21b-ooa-admix.lgo',
             time_handling = 'free', validate = FALSE)
cat('graph_to_lgo wrote t21b-ooa-admix.lgo\n')
" || fail "T2.1b: graph_to_lgo failed"
pass "T2.1b: graph_to_lgo emitted .lgo from edge tibble"

"$LEGOFIT/legosim" -i 5000 t21b-ooa-admix.lgo > t21b-ooa-admix.opf 2>&1
[[ "$(grep -v '^#' t21b-ooa-admix.opf | grep -v '^$' | wc -l | tr -d ' ')" -ge "5" ]] \
  || fail "T2.1b: legosim produced too few site patterns"
pass "T2.1b: legosim produced site patterns from emitted .lgo"

"$LEGOFIT/legofit" -t 2 -d 1e-3 t21b-ooa-admix.lgo t21b-ooa-admix.opf \
  > t21b-ooa-admix.legofit 2>&1
grep -q "DiffEv reached_goal" t21b-ooa-admix.legofit \
  || fail "T2.1b: ooa-admix did not converge"
pass "T2.1b: ooa-admix fit converged to reached_goal"

Rscript -e "
pkgload::load_all('$REPO_ROOT', quiet = TRUE)
result <- suppressMessages(read_legofit_output(
  '$WORK/t21b-ooa-admix.legofit', graph = make_ooa_admix_graph()
))
stopifnot(nrow(result) == 8L)
admix <- result[result\$type == 'admix', ]
stopifnot(nrow(admix) == 2L)
stopifnot(all(!is.na(admix\$admix_event_time)))
nt <- attr(result, 'node_times')
stopifnot(setequal(names(nt),
                   c('R','out','anc','arch_anc','hum_anc','arch','afr','eur')))
cat('ooa-admix round trip ok:', nrow(result), 'edges,',
    length(nt), 'node times,',
    'admix proportion fitted to', round(admix\$weight[1], 3), '\n')
" || fail "T2.1b: reader did not recover ooa-admix structure"
pass "T2.1b: reader recovered 8-edge OOA admix structure from fitted output"

echo
echo "==> T2.2: end to end bootstrap pipeline (multi-admix, 10 reps)"
echo "  bootci.py's true workflow needs tabpat (or booma) for moving-block"
echo "  bootstrap on genotype data. With only legosim-simulated patterns"
echo "  available we fake the bootstrap by treating N independent legosim"
echo "  replicates as bootstrap samples. This still validates the format"
echo "  contract end to end: legofit fit, flatfile.py combine, bootci.py."

cd "$WORK"
N=10
mkdir -p boot
rm -f boot/rep_*.opf boot/rep_*.legofit
for i in $(seq 1 $N); do
  "$LEGOFIT/legosim" -i 1000 "$FIXTURES/multi-admix.lgo" > boot/rep_${i}.opf 2>&1
done
[[ "$(ls boot/rep_*.opf 2>/dev/null | wc -l)" -eq "$N" ]] \
  || fail "T2.2: legosim failed to produce $N replicate patterns"
pass "T2.2: $N legosim replicates produced"

for i in $(seq 1 $N); do
  "$LEGOFIT/legofit" -t 2 -d 1e-3 \
    "$FIXTURES/multi-admix.lgo" boot/rep_${i}.opf > boot/rep_${i}.legofit 2>&1 &
done
wait
converged="$(grep -l reached_goal boot/rep_*.legofit | wc -l | tr -d ' ')"
[[ "$converged" -eq "$N" ]] || fail "T2.2: only $converged/$N reps reached_goal"
pass "T2.2: all $N reps converged"

python3 "$LEGOFIT/flatfile.py" boot/rep_1.legofit boot/rep_2.legofit boot/rep_3.legofit \
  boot/rep_4.legofit boot/rep_5.legofit boot/rep_6.legofit boot/rep_7.legofit \
  boot/rep_8.legofit boot/rep_9.legofit boot/rep_10.legofit > t22-multi-admix.flat \
  || fail "T2.2: flatfile.py failed"
pass "T2.2: flatfile.py combined fits"

python3 "$LEGOFIT/bootci.py" t22-multi-admix.flat > t22-multi-admix.bootci \
  || fail "T2.2: bootci.py failed"
pass "T2.2: bootci.py produced CI table"

Rscript -e "
pkgload::load_all('$REPO_ROOT', quiet = TRUE)
result <- read_legofit_bootstrap('$WORK/t22-multi-admix.bootci',
                                  graph = make_multi_admix_graph())
stopifnot(c('T_R','T_A','T_B','T_admix_m','m_m') %in% result\$parameter)
stopifnot(c('parameter','family','point_estimate','lo','hi') %in% names(result))
admix <- result[result\$parameter == 'T_admix_m', ]
stopifnot(admix\$lo <= admix\$point_estimate && admix\$hi >= admix\$point_estimate)
cat('multi-admix bootstrap ok: T_admix_m CI [', admix\$lo, ',', admix\$hi, ']\n')
" || fail "T2.2: reader could not parse multi-admix bootci output"
pass "T2.2: bootstrap reader recovered CIs with correct column shape"

echo
echo "==> T2.3: partial fit on rha20 (17 free params, finished_iterations)"

cd "$WORK"
# Step 1: generate patterns from rha20.lgo via legosim
"$LEGOFIT/legosim" -i 1000 "$FIXTURES/rha20.lgo" > t23-rha20.opf 2>&1
[[ -s t23-rha20.opf ]] || fail "T2.3: legosim produced no patterns"

# Step 2: bounded fit (1 generation, 100 reps, 3 points per dim)
"$LEGOFIT/legofit" -t 2 -d 1e-4 -S 1@100 -p 3 \
  "$FIXTURES/rha20.lgo" t23-rha20.opf > t23-rha20.legofit 2>&1
grep -q "DiffEv finished_iterations" t23-rha20.legofit \
  || fail "T2.3: bounded fit did not report finished_iterations"
pass "T2.3: rha20 partial fit produced finished_iterations"

Rscript -e "
pkgload::load_all('$REPO_ROOT', quiet = TRUE)
res <- suppressMessages(read_legofit_output('$WORK/t23-rha20.legofit', graph = NULL))
stopifnot(nrow(res) >= 17)
stopifnot(attr(res, 'fit_convergence')\$status == 'finished_iterations')
cat('rha20 partial fit ok:', nrow(res), 'params parsed\n')
" || fail "T2.3: reader could not parse rha20 partial fit"
pass "T2.3: reader parsed rha20 partial fit with finished_iterations attr"

echo
echo "==> All Tier 2 tests passed"
echo "  Optional: vendor t23-rha20.legofit as tests/testthat/fixtures/rha20.legofit"
echo "  (already vendored; re-copy via: cp $WORK/t23-rha20.legofit $FIXTURES/rha20.legofit)"
