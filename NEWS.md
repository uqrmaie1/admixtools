# admixtools 2.1.0

This is the first release since 2.0.10. It adds a full round trip LEGOFIT
integration, support for admixture graphs with sampled internal nodes, native
PLINK 2 (PFILE) support, OpenMP and multicore parallelism, large memory
reductions for big datasets, new tools for model exploration at scale and for
pipeline integration, and a number of correctness, compatibility, and
documentation fixes.

Note: 18 new functions are exported in this release (see the lists below).

## New features

* **LEGOFIT integration (round trip).** A new `R/legofit.R` module connects
  admixtools graphs to the LEGOFIT demographic history estimator:
    * `graph_to_lgo()` writes an admixtools or qpgraph admixture graph (igraph
      or edge tibble) to a LEGOFIT `.lgo` file, with three `time_handling` modes
      (`fix_admix`, `init`, `free`), flexible `twoN` parameterization, and a
      `fix_times` mode that removes the degeneracy in overall scale so absolute
      effective sizes are recoverable (#126, #158).
    * `read_lgo()` parses a `.lgo` file back into an edge tibble or igraph,
      preserving real intermediate nodes and each segment's `samples`, `twoN`,
      and `t` parameters, with grammar support for line continuation and bounded
      free declarations (#126, #145, #156).
    * `read_legofit_output()` pulls a fitted `.legofit` file onto graph edges
      (node times, admixture proportions, a classification of parameter
      identifiability), and `read_legofit_bootstrap()` loads a `bootci.py` file
      into a tidy confidence interval tibble, one row per parameter (#145).
    * `run_legofit()` drives the external `legofit` binary from start to finish:
      a graph goes in, a fitted edge tibble comes out, with typed error
      conditions when the binary is missing or the inputs are invalid (#160).
    * Validated by a faithful round trip of the published `rha20.lgo` (Rogers et
      al. 2020) superarchaic introgression model. A new `legofit.Rmd` vignette
      and bundled fixtures document the full workflow.
    * New exports: `graph_to_lgo()`, `read_lgo()`, `read_legofit_output()`,
      `read_legofit_bootstrap()`, `run_legofit()`, `default_drift_to_time()`.

* **Admixture graphs with sampled internal nodes.** Graphs can now carry a
  `nodes` tibble (`attr(graph, "nodes")`) so internal or ancestral nodes hold
  sample counts, effective sizes, and times that previously only leaf
  populations could. `add_sampled_tips()` rewrites such a graph into the tip and
  sibling form qpgraph fits, and `qpgraph()` now aborts with a clear, typed
  error (`admixtools_internal_samples_need_tips`) instead of silently producing
  a degenerate fit when handed a graph of this kind that has not been
  transformed. New exports: `add_sampled_tips()`, `graph_nodes()`,
  `set_node_attrs()`, `node_times()`, `as_edge_tibble()` (#157).

* **PLINK 2 PFILE support.** The new `pfile_to_afs()` reads PLINK 2
  `.pgen`/`.pvar`/`.psam` data directly into the f-statistic pipeline (via the
  new `pgenlibr` dependency), so `extract_f2()`, `qpadm()`, and `qpgraph()` now
  run on native PFILE panels with no conversion to BED. `anygeno_to_aftable()`
  detects a PFILE prefix automatically (including `.pvar.zst`) and dispatches to
  it, with PFILE taking precedence over BED at the same prefix. A `multiallelic`
  argument sets the multiallelic policy and an optional `cm_file` argument grafts
  on the centimorgan data that `.pvar` files lack (#109). The `qpfstats` path
  also gained a PFILE backend, so `qpfstats()` and `extract_f2(qpfstats = TRUE)`
  work on PFILE inputs, with results identical to the BED backend (#111).

* **`qpadm_sweep()`** fits qpAdm for every combination of targets, named source
  sets, and named right sets in a single call, loading the f2 cache only once
  and returning a flat tibble ready to filter (one row per model), with optional
  parallelism over models via `furrr` (#118). It also surfaces the rank
  deficiency diagnostics `f4_var_rcond` and `f4_var_singular_loadings` (and
  `popdrop`) for each model, and `qpadm_multi()` now validates forwarded `...`
  arguments with pointed errors (#144).

* **`result_to_json()`** serializes any admixtools fitter result
  (the list results from `qpadm`, `qpwave`, `qpgraph`, and `qpadm_sweep`, the 3D
  array from `qpfstats`, or a single tibble from `f4`, `f3`, `f2`, and
  `qpdstat`) into a stable JSON envelope with a versioned
  schema, making it easy to drive analyses from Python, Rust, Go, or other
  languages. Adds `jsonlite` to Imports (#119).

* **f2 cache fingerprint and metadata.** `compute_f2_cache_id()` returns a
  stable `sha256:<hex>` fingerprint of the inputs that produced an f2 directory,
  and `extract_f2()` now writes both a `.f2_cache_id` sidecar and a
  `cache_metadata.json` sidecar with a versioned schema (populations, SNP and
  block counts, filter arguments, build timestamp, qpfstats flag) on success.
  The new `read_f2_cache_metadata()` reads that metadata back as an R list.
  Together they let pipelines validate cached qpAdm or qpGraph results and
  recover build telemetry cheaply, without rehashing the block files on disk
  (#117, #134).

* **`discard_from_aftable()`** is now exported and documented. This is the SNP
  filter (`maxmiss`, `minmaf`, `maxmaf`, `minac2`, `outpop`, `auto_only`,
  `poly_only`, `transitions`, `transversions`, `keepsnps`) that `extract_f2()`
  applies between reading
  allele frequencies and computing f2 blocks; exporting it lets callers assemble
  AFS pipelines by hand without reaching into the namespace with `:::` (#110).

* **`find_graphs_old()`** is exported again. It was fully documented but had been
  dropped from NAMESPACE by accident (back in 2021), so calls failed with "could
  not find function"; the export is restored (#131).

## Minor improvements

* `qpadm()` and `qpwave()` gained a `singular_threshold` argument and now report
  the reciprocal condition number of the f4 variance matrix as `f4_var_rcond`.
  When that matrix is nearly singular they also compute
  `f4_var_singular_loadings`, flagging collinear right populations (for example
  sister clades) whose linear dependence causes a silent pseudoinverse fallback
  and misleadingly tight weight standard errors. The default
  (`singular_threshold = NA`) preserves prior behavior and existing return
  fields are unchanged (#121). The loadings diagnostic was later reworked so it
  uses only the right populations, does not depend on the chosen basis, and is
  reproducible across BLAS and LAPACK backends (#152).

* `write_dot()` was rewritten with new arguments: `fontsize`, `color`,
  `hide_weights`, `highlight_unidentifiable`, `nodesep`, `ranksep`, and
  `fix_names`. With `color = TRUE` (the new default) the Graphviz output now
  matches the colors of `plot_graph()`; `color = FALSE` skips the ggplot
  pipeline and is about 11x faster. The new arguments are appended after the
  existing positional ones, so prior positional calls are unaffected (#122).

* Long running loops now use a progress emitter that adapts to its output. In an
  interactive terminal a single progress line overwrites in place, with no
  leftover characters; on a redirected or piped stderr it emits one line per
  iteration so orchestrators and CI logs see live progress. This covers
  `extract_f2()` and `extract_afs()`, `read_f2()`, `afs_to_f2_blocks()`,
  `plink_to_afs()` and `eigenstrat_to_afs()`, the qpgraph topology search, and
  more, including a rollup for each chromosome on the `qpfstats = TRUE` path
  (#116, #130).

## Performance

* **f-statistic kernels parallelized with OpenMP.** The four hot C++ kernels
  (the accumulator that turns genotypes into allele frequencies, and the f4
  numerator and denominator kernels) now use OpenMP, cutting wall time on dense
  `qpdstat` runs by about 51% (and 5 to 10% on `qpadm`) with identical output.
  The thread budget is resolved via `parallelly::availableCores()` so it
  composes with `future::plan()` without oversubscription, and it falls back to
  serial where the toolchain lacks OpenMP (for example macOS Apple clang without
  libomp) (#143).

* **Parallel f-statistic extraction.** The loop over blocks in
  `f4blockdat_from_geno()` (behind `extract_f2()` and f4 or D-statistic
  extraction) now runs through `furrr::future_map()`, so setting a `future` plan
  (for example `plan(multisession, workers = 8)`) spreads it across cores, with
  speedups approaching the worker count; the sequential default is identical
  (#106). `read_f2()` likewise parallelizes its per pair `readRDS()` loads under
  a `future` plan (#114).

* **Lower memory.** Genotype readers (PLINK BED, PACKEDANCESTRYMAP, EIGENSTRAT)
  now return integer dosage matrices in {0, 1, 2, NA} instead of doubles,
  halving the largest allocation per block (about 26% lower peak R memory on a
  qpadm run straight from genotypes) with identical output (#161). Separately,
  `f4blockdat_from_geno()` uses a streaming C++ kernel that avoids building the
  full (npopcomb by nsnp) matrix, cutting peak memory per block from
  O(npopcomb by nsnp) to O(npopcomb) (for example about 45 GB down to about 12
  MB at 1.57M population combinations) and running up to about 2.6x faster in
  parallel (#107).

* **Faster, lighter `qpfstats()`.** `qpfstats()` now consumes the f4 block
  matrices directly through a jackknife computed over those matrices, instead of
  building a tibble in long format and running a dplyr `group_by` jackknife over
  each population combination. This removes an intermediate frame of hundreds of
  GB and the worse than linear slowdown on dense f4 or large population sets.
  `f4blockdat_from_geno()` gained a `return_matrices` argument to support this
  path (#108).

* The jackknife over pairs in `f2()` and `fst()` now uses a base R loop instead
  of dplyr `group_by`, removing a blowup in dplyr's data mask that hung for many
  minutes at large pair counts (now seconds; about 3x faster at smaller scales),
  with identical results (#115).

* Allele frequency accumulation per block (`gmat_to_aftable()`) now delegates to
  a single pass C++ kernel (`cpp_gmat_to_aftable`), about 2.5 to 3x faster per
  call with identical output (#113).

* The qpgraph scoring path (run once per scored graph during `find_graphs()`)
  was optimized: vertex ids are passed to `igraph::get_edge_ids()` directly
  rather than resolved through names on every path, in and out degrees are
  derived from the raw edge list, the graph root and admixture edges are
  computed once and shared, and the index tables are built with base R vector
  operations rather than a dplyr pipeline. Output is identical; `find_graphs()`
  wall time drops by roughly 30% with about 23% less allocation (#163, #164,
  #165, #167, #168).

## Bug fixes

* **Upward bias in f2, f3, and f4 with `maxmiss > 0`.** Partly missing genotypes
  leaked through a `tidyr::replace_na()` call that does nothing on matrices,
  silently inflating f2, f3, and f4 statistics by roughly N/(N-k) for population
  pairs with k missing SNPs per block. Missing cells are now zeroed one element
  at a time. Runs with the default `maxmiss = 0` are identical, but **f2 caches
  built with `maxmiss > 0` under the old code should be rebuilt** (#155).

* `qpfstats()` no longer silently returns a result that is entirely NaN on data
  with missingness: the regression is now solved block by block with a weighted
  least squares that drops the population combinations whose values are not
  finite, instead of letting a single NaN cell poison every statistic for that
  block through BLAS propagation. Results on clean data are unchanged to machine
  precision (#112).

* `qpadm()` and `qpwave()` now return finite weight standard errors on inputs
  with missing blocks (f2 block arrays, a supplied `f4blocks`, or
  `remove_na = FALSE`), where they previously errored or returned NaN; the block
  jackknife no longer leaks NA through a `tidyr::replace_na()` call that does
  nothing on arrays (#154).

* `qpadm()` no longer crashes on a precomputed f2 directory. The internal
  `get_f2()` call now drops arguments that only apply to genotype reading
  (`auto_only`, `blgsize`, `poly_only`, and so on) instead of forwarding them to
  `f2_from_precomp()`, restoring the usual `extract_f2()` then `qpadm()` workflow
  (#120, #92).

* Duplicate individual IDs across populations are now allowed: the genotype
  readers reject only duplicated (individual, population) pairs, so the same
  sample can be assigned to more than one population label, and the matching of
  samples to populations was corrected so requesting a subset of populations no
  longer misaligns the sample selection (thanks @floutt, #92).

* Fixed a crash in `qpadm(target = NULL, return_f4 = TRUE)` on a genotype prefix
  path ("Column 'pop1' doesn't exist"), including with `allsnps = TRUE` (#124).

* Fixed a regression where a SNP block containing exactly one kept SNP crashed
  `f4blockdat_from_geno()` with "Not a matrix" (a lost `drop = FALSE` guard),
  affecting blocks with few SNPs such as chromosome tails, a tight `blgsize`,
  and ancient DNA panels (#127).

* `geno_to_treemix()` no longer writes leftover `"NA,NA"` cells for partly
  missing SNP rows (another `tidyr::replace_na()` call that does nothing on a
  matrix); NA cells are zeroed one element at a time before output (#153).

* Fixed a `find_graphs(numgraphs = 1)` crash under dplyr 1.1+, caused by
  `ifelse()` on an empty test producing a `<logical>` column; it now uses
  `dplyr::if_else()` (#105, #104).

* Fixed six latent "could not find function" crashes from renamed or unimported
  helpers (in `joint_sfs()`, `f4blockdat_from_geno(allsnps = "qpfs")`,
  `proxypred`, `consistent_with_qpadm`, the path and leaf distance
  reconstruction, and `evaluate_moreadmix`), and made `R CMD check` pass with
  `--as-cran` (#159, and earlier #125).

* `run_shiny_admixtools()` now checks for the required GUI packages (DT, plotly,
  shinydashboard, shinyFiles, and so on) before launching and reports the exact
  `install.packages()` command, instead of loading a UI that dies silently on
  the first missing dependency (#133, #90).

* `plotly_comparison()` no longer clips a long title on the horizontal axis in
  the interactive plot (its bottom margin was increased) (#0f75732).

* **igraph 2.1.x compatibility.** `qpgraph()` broke under igraph 2.1.4 and
  later, which tightened `get.edge.ids()`. The qpgraph helpers now call the
  modern `igraph::get_edge_ids()`, DESCRIPTION requires `igraph (>= 2.1.0)`, and
  the bundled `example_igraph` and `example_opt` datasets were upgraded with
  `igraph::upgrade_graph()` to silence the warning about an older igraph version
  (#103, #149, #141, #94).

* **readr 2.0+ compatibility:** replaced the removed `readr::read_table2()` with
  `read_table()` across genotype I/O, fixing `read_plink()`,
  `read_packedancestrymap()`, `extract_f2()`, and related functions (#101).

* **RcppArmadillo 15+ build:** raised the C++ standard from C++11 to C++17
  (`SystemRequirements: C++17`), which RcppArmadillo 15 requires (#102, #100).
  Also silenced the harmless RcppArmadillo "system is singular" messages on
  stderr on modern Armadillo by setting `ARMA_WARN_LEVEL 1` (#147).

* Fixed an invalid `break` outside a loop in the msprime script writer
  `msprime_genome()`, and moved `qpgraph_resample_snps2()` off the removed
  `furrr::future_invoke_map()` (#125).

## Documentation

* New `legofit.Rmd` vignette documenting the full LEGOFIT round trip (#160).
* Rewrote the parallelization vignette for the current APIs (parallelism over
  blocks, over models, and over resamples) and the
  `plan(multisession/multicore, workers = N)` idiom (#129).
* Expanded the qpAdm (#132), f-statistics (#135), resampling (#136), plotting
  (#137), admixture graphs (#140, #162), main tutorial (#138), and I/O (#128)
  vignettes, covering PFILE, `qpadm_sweep`, the rcond diagnostics, the cache
  sidecars, and `plot_graph_map`.
* Added the missing help page for `discard_from_aftable()` (#123).

## Dependencies and build

* New Imports: `jsonlite` (JSON output) and `pgenlibr` (PLINK 2 PFILE reading).
* `igraph` now requires `>= 2.1.0`.
* `SystemRequirements`: `C++17, OpenMP (optional, for the parallel f-stat
  kernels)`.
