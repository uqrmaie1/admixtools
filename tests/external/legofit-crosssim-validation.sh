#!/usr/bin/env bash
# Tier 3 cross-simulator fidelity for the LEGOFIT bridge.
# S series (cross simulator).
#
# The strong statistical test: generate data with an INDEPENDENT coalescent
# engine (msprime) rather than legosim, so the data-generating process differs
# from LEGOFIT's inference model. Validated separately (Path 1).
#   S-6: the msprime->simpat spectrum agrees with a legosim draw of the same
#        model within sampling noise (the ingest is not distorting the spectrum).
#   S-5: fitting the msprime patterns recovers the truth in coalescent units.
#
# Requires LEGOFIT built ($LEGOFIT) and python3 with msprime. If msprime is
# absent the harness SKIPS (exit 0) rather than failing, since it is an optional
# external dependency.
set -euo pipefail

: "${LEGOFIT:?must point to legofit src dir (legofit, legosim, simpat)}"
[[ -x "$LEGOFIT/simpat" ]]  || { echo "no simpat in $LEGOFIT"; exit 1; }

REPO="${REPO:-$(cd "$(dirname "$0")/../.." && pwd)}"
WORK="$(cd "$(dirname "$0")" && pwd)/work-crosssim"
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"

if ! python3 -c "import msprime" 2>/dev/null; then
  echo "==> SKIP cross-simulator harness: python3 msprime not available"
  exit 0
fi
echo "==> $("$LEGOFIT/legosim" --version 2>&1 | grep -m1 version | tr -s ' ')  msprime $(python3 -c 'import msprime;print(msprime.__version__)')"

fails=0
pass () { echo "  PASS $1"; }
fail () { echo "  FAIL $1"; fails=$((fails+1)); }

# Ground truth in msprime's OWN units: diploid Ne=5000 (2N=10000), splits at
# 10000 and 20000 generations => coalescent-unit truth Txy=1.0, Txyz=2.0.
cat > msp.py <<'PY'
import sys, msprime
NE, T1, T2, L, MU, RHO, SEED = 5000, 10000, 20000, int(1e8), 1.25e-8, 1e-8, 42
dem = msprime.Demography()
for p in ("x","y","z","xy","xyz"): dem.add_population(name=p, initial_size=NE)
dem.add_population_split(time=T1, derived=["x","y"], ancestral="xy")
dem.add_population_split(time=T2, derived=["xy","z"], ancestral="xyz")
ts = msprime.sim_ancestry(samples={"x":1,"y":1,"z":1}, demography=dem,
        sequence_length=L, recombination_rate=RHO, ploidy=2, random_seed=SEED)
ts = msprime.sim_mutations(ts, rate=MU, random_seed=SEED,
        model=msprime.BinaryMutationModel())
cols = [0,2,4]
out = open(sys.argv[1],"w")
out.write("npops = 3\npop sampsize\nx 1\ny 1\nz 1\n")
for var in ts.variants():
    g = var.genotypes[cols]
    if g.max() > 1: continue
    s = int(g.sum())
    if s == 0 or s == 3: continue
    out.write("1 " + " ".join(str(int(a)) for a in g) + "\n")
out.close()
PY
# simpat includes singletons by default (no -1 flag) and takes the file as an arg.
if python3 msp.py msp.simin 2>/dev/null && "$LEGOFIT/simpat" msp.simin > msp.simpat 2>/dev/null; then
  pass "msprime -> simpat ingest produced site patterns"
else
  fail "msprime -> simpat ingest"
fi

# truth model (for the legosim spectrum draw) and a fit model (free, seeded off truth)
cat > truth.lgo <<'LGO'
time fixed zero=0
twoN fixed one=1
time free  Txy=1.0
time free  Txyz=2.0
segment x   t=zero twoN=one samples=1
segment y   t=zero twoN=one samples=1
segment z   t=zero twoN=one samples=1
segment xy  t=Txy  twoN=one
segment xyz t=Txyz twoN=one
derive x  from xy
derive y  from xy
derive xy from xyz
derive z  from xyz
LGO
sed -E 's/Txy=1.0/Txy=0.7/; s/Txyz=2.0/Txyz=1.5/' truth.lgo > fit.lgo

# S-6: spectrum agreement -----------------------------------------------------
"$LEGOFIT/legosim" -1 -i 200000 truth.lgo > truth.opf 2>/dev/null
Rscript -e "
parse_opf <- function(p){ ln<-readLines(p); ln<-ln[!grepl('^#',ln)&nzchar(trimws(ln))]
  pr<-strsplit(trimws(ln),'[[:space:]]+'); v<-setNames(as.numeric(sapply(pr,'[',2)),sapply(pr,'[',1)); v/sum(v) }
a <- parse_opf('$WORK/msp.simpat'); b <- parse_opf('$WORK/truth.opf')
k <- intersect(names(a),names(b)); d <- max(abs(a[k]-b[k]))
cat(sprintf('  S-6 spectra: %d shared patterns, max|freq diff| = %.4f\n', length(k), d))
stopifnot(length(k)==6, d < 0.01)
" && pass "S-6 msprime spectrum agrees with legosim draw" || fail "S-6 spectrum"

# S-5: cross-simulator recovery -----------------------------------------------
"$LEGOFIT/legofit" -1 -d 1e-4 -p 12 -T 1e-6 fit.lgo msp.simpat > fit.legofit 2>/dev/null
Rscript -e "
suppressMessages(pkgload::load_all('$REPO', quiet=TRUE))
r <- read_legofit_output('$WORK/fit.legofit'); v <- setNames(r\$value, r\$name)
e_xy <- abs(v[['Txy']]-1.0); e_xyz <- abs(v[['Txyz']]-2.0)/2.0
cat(sprintf('  S-5 recovered Txy=%.3f (truth 1.0)  Txyz=%.3f (truth 2.0)\n', v[['Txy']], v[['Txyz']]))
stopifnot(e_xy < 0.20, e_xyz < 0.20)
" && pass "S-5 cross-simulator recovery within 20%" || fail "S-5 recovery"

echo
if [[ "$fails" -eq 0 ]]; then
  echo "==> ALL CROSS-SIMULATOR CHECKS PASSED (S-5, S-6)"
  exit 0
else
  echo "==> $fails CROSS-SIMULATOR CHECK(S) FAILED"
  exit 1
fi
