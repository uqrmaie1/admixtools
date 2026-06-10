# LEGOFIT test fixtures

Small LEGOFIT inputs and outputs used by the `graph_to_lgo()`,
`read_lgo()`, `read_legofit_output()`, and `read_legofit_bootstrap()`
tests. They let the reader and writer tests run without a LEGOFIT
binary installed.

## Write-direction golden files

Byte-exact expected output of `graph_to_lgo()`, compared verbatim by
the writer tests.

- `minimal.lgo` is the default `graph_to_lgo()` output for the minimal
  one-admixture graph.
- `minimal-init.lgo` and `minimal-free.lgo` are the same graph under
  `time_handling = "init"` and `time_handling = "free"`.
- `minimal-twoN-scalar.lgo` and `minimal-twoN-named.lgo` cover the
  scalar and named-vector forms of the `twoN=` argument.

## Read-direction pipeline fixtures (ourex1)

A three-leaf, no-admixture topology (the xyz/xy split) carried end to
end through the LEGOFIT pipeline, used by the `read_legofit_*` tests.

- `ourex1.lgo` is the `graph_to_lgo()` topology.
- `ourex1.opf` is the site-pattern file from `legosim`, the simulated
  "real data" input to `legofit`.
- `ourex1.legofit` is the converged fit.
- `ourex1.flat` is the combined flat file from `flatfile.py`.
- `ourex1.bootci` is the confidence-interval table from `bootci.py`.

## Special-case read fixtures

- `ourex1-twoN-scalar.lgo` and `ourex1-twoN-scalar.legofit` exercise
  the scalar `shared=N` twoN mode, which must read back without a
  spurious parameter-mismatch message.
- `ourex1-underconverged.legofit` is a fit capped at one generation
  (`legofit -S 1`) so it ends in `finished_iterations`. It drives the
  `legofit_fit_incomplete` message path.
- `multi-admix.lgo`, `multi-admix.legofit`, and `multi-admix.bootci`
  are a three-leaf, one-admixture topology, the smallest fixture that
  carries an admix-event time and a mixFrac through the reader.
- `ooa-admix.lgo` and `ooa-admix.legofit` are an out-of-Africa-style
  scenario with four sampled populations (out, arch, afr, eur) and one
  admixture event (archaic ancestry into eur). It is the multi-leaf,
  fittable counterpart to `minimal.lgo`, which has a single leaf and so
  cannot be fit.

## Published reference

- `rha20.lgo` and `rha20.legofit` are the Rogers, Harris & Achenbach
  2020 (Science Advances, doi:10.1126/sciadv.aay5483)
  superarchaic-introgression model. `rha20.lgo` ships with LEGOFIT
  under the ISC license (`src/rha20.lgo`, citation `Li:N-505-43-S88`).
  They confirm the reader handles real third-party output, not only
  fixtures generated here.

## Regenerating the fixtures

The tools come from a built LEGOFIT source tree. Set `LEGOFIT_SRC` to
its `src/` directory. The committed files were generated with LEGOFIT
git `3675373` (2026-05-18) and system Python 3.

```bash
LEGOFIT_SRC=/path/to/legofit/src

# ourex1 (three-leaf, no admixture)
"$LEGOFIT_SRC"/legosim -1 -d 0 ourex1.lgo > ourex1.opf
"$LEGOFIT_SRC"/legofit --threads 1 -1 -d 0 ourex1.lgo ourex1.opf > ourex1.legofit

# Bootstrap replicates are intermediates and are NOT committed. A real
# bootstrap repeats the fit on resampled data; here the fitted values
# are perturbed deterministically, which is enough for format testing.
for i in 0 1 2 3 4; do
  PERTURB_XYZ=$(awk "BEGIN { print 2 + ($i - 2) * 0.05 }")
  PERTURB_XY=$(awk "BEGIN { print 0.5 + ($i - 2) * 0.02 }")
  sed -e "s/T_xyz = 2$/T_xyz = $PERTURB_XYZ/" \
      -e "s/T_xy = 0.5$/T_xy = $PERTURB_XY/" \
      ourex1.legofit > ourex1_boot${i}.legofit
done
python3 "$LEGOFIT_SRC"/flatfile.py ourex1.legofit ourex1_boot*.legofit > ourex1.flat
python3 "$LEGOFIT_SRC"/bootci.py ourex1.flat > ourex1.bootci

# multi-admix (three-leaf, one admixture); multi-admix.opf is an
# intermediate and is NOT committed.
"$LEGOFIT_SRC"/legosim -i 5000 multi-admix.lgo > multi-admix.opf
"$LEGOFIT_SRC"/legofit -t 2 -d 1e-3 multi-admix.lgo multi-admix.opf > multi-admix.legofit
```
