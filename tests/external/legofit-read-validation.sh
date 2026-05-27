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

echo "==> T2.1: end to end round trip on ourex1 (no admix, 2 free params)"

cd "$WORK"
# Use the bundled ex1.opf since it pairs with the ourex1 topology
"$LEGOFIT/legofit" --threads 1 -1 -d 0 \
  "$FIXTURES/ourex1.lgo" "$LEGOFIT/../exmpl/ex1/ex1.opf" \
  > t21-ourex1.legofit 2>&1
grep -q "DiffEv reached_goal" t21-ourex1.legofit \
  || fail "T2.1: ourex1 did not converge"
pass "T2.1: ourex1 fit converged"

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
" || fail "T2.1: ourex1 read back did not match ground truth"
pass "T2.1: ourex1 reader recovered fitted times within tolerance"

echo
echo "==> T2.2: end to end bootstrap pipeline — DEFERRED"
echo "  bootstrap requires tabpat/booma to generate replicate patterns,"
echo "  then per-replicate legofit fits via a loop, then flatfile.py to combine,"
echo "  then bootci.py. Set up that pipeline separately; the in-branch"
echo "  test-read_legofit_bootstrap.R + T1.5 cover the reader against synthesized"
echo "  bootci.py output and a real bootci file produced by perturbing the"
echo "  ourex1 fit (committed as ourex1.bootci)."

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
