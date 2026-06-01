#!/usr/bin/env bash
# UNIFIED external validation harness for the LEGOFIT bridge.
# Runs the write-side acceptance checks (I-1, I-7) AND the existing read-side
# checks behind one entry point. Exits non-zero on ANY failure so a nightly job
# gets a clean signal.
#
# Requires LEGOFIT 1.87+ built locally. Set $LEGOFIT to the src/ dir
# (legofit, legosim, bootci.py, flatfile.py). LEGOFIT is not a package
# dependency, so this runs manually and in the nightly job, not in R CMD check.
#
# Run from the package root:
#   LEGOFIT=/path/to/legofit/src tests/external/legofit-validation.sh
set -euo pipefail

: "${LEGOFIT:?must point to legofit src dir (legofit, legosim, bootci.py)}"
[[ -x "$LEGOFIT/legosim" ]] || { echo "no legosim in $LEGOFIT"; exit 1; }
[[ -x "$LEGOFIT/legofit" ]] || { echo "no legofit in $LEGOFIT"; exit 1; }

# Derive the package root from this script's location (tests/external/..),
# matching legofit-read-validation.sh. $REPO overrides if set.
REPO="${REPO:-$(cd "$(dirname "$0")/../.." && pwd)}"
FIX="$REPO/tests/testthat/fixtures"
WORK="$(cd "$(dirname "$0")" && pwd)/work-unified"
rm -rf "$WORK"; mkdir -p "$WORK"   # clean: a stale work dir caused spurious first-run failures

# version pinning, per the master plan provenance rule
VER="$("$LEGOFIT/legosim" --version 2>&1 | grep -m1 version | tr -s ' ')"
echo "==> LEGOFIT $VER"

fails=0
pass () { echo "  PASS $1"; }
fail () { echo "  FAIL $1"; fails=$((fails+1)); }

echo "==> WRITE SIDE I-1: legosim --network accepts every Tier 0 write golden"
for g in minimal minimal-init minimal-free minimal-twoN-scalar minimal-twoN-named; do
  if "$LEGOFIT/legosim" --network "$FIX/$g.lgo" >/dev/null 2>"$WORK/$g.err"; then
    pass "I-1 $g.lgo"
  else
    fail "I-1 $g.lgo: $(head -1 "$WORK/$g.err")"
  fi
done

echo "==> WRITE SIDE I-7: legosim accepts graph_to_lgo-generated .lgo"
echo "    (within the validity envelope: trees + small 1-admix; complex"
echo "     admix free-exports can be legosim-invalid and are excluded by design)"
Rscript -e "
  suppressMessages(pkgload::load_all('$REPO', quiet=TRUE))
  specs <- list(c(3,0),c(4,0),c(5,0),c(6,0),c(4,1),c(5,1))
  for (i in seq_along(specs)) {
    s <- specs[[i]]; set.seed(100+i)
    g <- random_sim(nleaf=s[1], nadmix=s[2])
    l <- graph_to_lgo(g\$edges, time_handling='free', validate=FALSE)
    writeLines(l, file.path('$WORK', sprintf('gen_n%da%d.lgo', s[1], s[2])))
  }
" || { echo "  FAIL I-7: generation step"; fails=$((fails+1)); }
for f in "$WORK"/gen_*.lgo; do
  [[ -e "$f" ]] || continue
  b="$(basename "$f")"
  if "$LEGOFIT/legosim" --network "$f" >/dev/null 2>"$WORK/$b.err"; then
    pass "I-7 $b"
  else
    fail "I-7 $b: $(head -1 "$WORK/$b.err")"
  fi
done

echo "==> READ SIDE: delegate to the in-branch read-validation harness"
if LEGOFIT="$LEGOFIT" bash "$REPO/tests/external/legofit-read-validation.sh"; then
  pass "read-side harness green"
else
  fail "read-side harness"
fi

echo "==> STATISTICAL FIDELITY: delegate to the S-series harness (S-1,S-3,S-4,I-4)"
if LEGOFIT="$LEGOFIT" bash "$REPO/tests/external/legofit-statistical-validation.sh"; then
  pass "statistical harness green"
else
  fail "statistical harness"
fi

echo "==> CROSS SIMULATOR: delegate to the msprime->simpat harness (S-5,S-6)"
if LEGOFIT="$LEGOFIT" bash "$REPO/tests/external/legofit-crosssim-validation.sh"; then
  pass "cross-simulator harness green (or skipped if msprime absent)"
else
  fail "cross-simulator harness"
fi

echo "==> REAL ARCHAIC DATA: delegate to the bundled-.daf harness (B-1,B-2)"
if LEGOFIT="$LEGOFIT" bash "$REPO/tests/external/legofit-archaic-validation.sh"; then
  pass "archaic real-data harness green"
else
  fail "archaic real-data harness"
fi

# S-2 bootstrap coverage is the slow tier (K x B fits). Opt in with RUN_COVERAGE=1.
if [[ "${RUN_COVERAGE:-0}" == "1" ]]; then
  echo "==> COVERAGE: delegate to the S-2 known-truth bootstrap coverage harness"
  if LEGOFIT="$LEGOFIT" bash "$REPO/tests/external/legofit-coverage-validation.sh"; then
    pass "S-2 coverage harness green"
  else
    fail "S-2 coverage harness"
  fi
else
  echo "==> COVERAGE: S-2 skipped (set RUN_COVERAGE=1 to run the slow coverage tier)"
fi

echo
if [[ "$fails" -eq 0 ]]; then
  echo "==> ALL CHECKS PASSED (write side I-1,I-7 + read side)"
  exit 0
else
  echo "==> $fails CHECK(S) FAILED"
  exit 1
fi
