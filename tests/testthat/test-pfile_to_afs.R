# Tests for pfile_to_afs() — covers the regressions called out in PR #109's
# review (error handling, poly_only filter, ploidy probe corner case, dup
# SNP id rejection, FID-less psam) plus an equivalence check against
# plink_to_afs() on the matching BED fixture.
#
# Every test that builds a fixture calls build_pfile_fixture() (defined in
# helper-pfile.R), which skips the test if the plink2 binary isn't on PATH.

test_that("biallelic PFILE -> AFS matches BED -> AFS (plink_to_afs equivalence)", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  if(is.na(fx$bed_pref)) skip("BED conversion failed; cannot compare")

  pfile_afs = pfile_to_afs(fx$pfile_pref, verbose = FALSE)
  bed_afs   = plink_to_afs(fx$bed_pref,   verbose = FALSE)

  # Same SNP-row set + same population columns (popA, popB).
  expect_setequal(rownames(pfile_afs$afs), rownames(bed_afs$afs))
  expect_setequal(colnames(pfile_afs$afs), colnames(bed_afs$afs))

  # Align rows and compare numerically. The two readers may emit pops in
  # different orderings, so subset by name on both sides.
  common = intersect(rownames(pfile_afs$afs), rownames(bed_afs$afs))
  pf = pfile_afs$afs[common, colnames(bed_afs$afs)]
  bd = bed_afs$afs[common, ]
  na_mask = is.na(pf) & is.na(bd)
  expect_true(all(na_mask | abs(pf - bd) < 1e-10),
              info = paste("max abs diff =", max(abs(pf - bd), na.rm = TRUE)))
})

test_that("multiallelic = 'error' (default) lists the offending site", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = TRUE)
  expect_error(
    pfile_to_afs(fx$pfile_pref, verbose = FALSE),
    "multiallelic"
  )
  # Error message names the multi-allelic row(s) so the user can fix the input.
  err = tryCatch(pfile_to_afs(fx$pfile_pref, verbose = FALSE),
                 error = conditionMessage)
  expect_match(err, "row 10")
  expect_match(err, "--max-alleles 2")
})

test_that("multiallelic = 'skip' drops only the multi-allelic site", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = TRUE)
  with_multi = pfile_to_afs(fx$pfile_pref, multiallelic = "skip", verbose = FALSE)
  no_multi = build_pfile_fixture(withr::local_tempdir(), with_multi = FALSE)
  base = pfile_to_afs(no_multi$pfile_pref, verbose = FALSE)
  # Skipped fixture should have one less SNP than the all-biallelic equivalent
  # (the tri-allelic row was the only difference). Identity of SNP ids should
  # otherwise match (modulo the missing rs_multi).
  expect_equal(nrow(with_multi$afs), nrow(base$afs) - 1L)
  expect_false("rs_multi" %in% rownames(with_multi$afs))
})

test_that("FID-less .psam reads with a single 'unknown' population", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE, with_fid = FALSE)
  res = pfile_to_afs(fx$pfile_pref, verbose = FALSE)
  # When FID is absent, .read_psam falls back to FID = "0" for every sample,
  # so match_samples collapses all 6 samples into a single population.
  expect_equal(ncol(res$afs), 1L)
  expect_equal(colnames(res$afs), "0")
})

test_that("no .pvar CM column triggers cm = 0 + warning", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE, with_cm_col = FALSE)
  # alert_warning() is `cat(crayon::yellow(...))` — written to stdout, not
  # through R's signalling system — so grab the printed output and grep.
  out = capture.output(
    res <- pfile_to_afs(fx$pfile_pref, verbose = TRUE)
  )
  expect_true(any(grepl("All variants have cm = 0", out)),
              info = "expected cm = 0 warning in stdout")
  expect_true(all(res$snpfile$cm == 0))
})

test_that("cm_file argument populates cm via SNP-id lookup", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  cm_path = file.path(dir, "cm_map.tsv")
  # Build a cm map covering 19/20 SNPs (rs1..rs19) with cm = i/100, leaving
  # rs20 unmatched so the NA-fallback path is exercised.
  writeLines(c("variant_id\tcm",
               paste0("rs", 1:19, "\t", sprintf("%.6f", (1:19) / 100))),
             cm_path)
  res = pfile_to_afs(fx$pfile_pref, cm_file = cm_path, verbose = FALSE)
  # Expect cm == row position / 100 for SNPs whose id is in the map.
  by_id = setNames(res$snpfile$cm, res$snpfile$SNP)
  expect_equal(by_id[["rs1"]],  0.01)
  expect_equal(by_id[["rs10"]], 0.10)
  expect_equal(by_id[["rs19"]], 0.19)
  expect_true(is.na(by_id[["rs20"]]))
})

test_that("first > last is rejected with a clear error (issue #2)", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  expect_error(
    pfile_to_afs(fx$pfile_pref, first = 10, last = 5, verbose = FALSE),
    "'first' must be <= 'last'"
  )
})

test_that("adjust_pseudohaploid = 0 behaves as FALSE (diploid assumed; issue #3)", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  with_zero  = pfile_to_afs(fx$pfile_pref, adjust_pseudohaploid = 0, verbose = FALSE)
  with_false = pfile_to_afs(fx$pfile_pref, adjust_pseudohaploid = FALSE, verbose = FALSE)
  expect_equal(with_zero$afs,    with_false$afs)
  expect_equal(with_zero$counts, with_false$counts)
  # Counts should be 2 * n_samples in each populated pop = 6 per row. rs7
  # carries hand-injected missing data in popB so its count is < 6 — exclude
  # it from the all-rows check, then verify it independently.
  expect_true(all(with_zero$counts[setdiff(rownames(with_zero$counts), "rs7"), "popA"] == 6))
  expect_true(all(with_zero$counts[setdiff(rownames(with_zero$counts), "rs7"), "popB"] == 6))
  expect_lt(with_zero$counts["rs7", "popB"], 6)
})

test_that("poly_only keeps locally-polymorphic SNPs and drops globally-fixed ones (issue #4)", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  # SNP rs5 is fixed REF in popA + fixed ALT in popB — locally polymorphic
  # (cpp_is_polymorphic == TRUE), row-mean = 0.5 so the legacy rowMeans
  # check would keep it too. SNP rs6 is globally fixed (all REF) — must be
  # dropped under both filters.
  res_default = pfile_to_afs(fx$pfile_pref, poly_only = FALSE, verbose = FALSE)
  res_poly    = pfile_to_afs(fx$pfile_pref, poly_only = TRUE,  verbose = FALSE)

  expect_true("rs5" %in% rownames(res_default$afs))
  expect_true("rs6" %in% rownames(res_default$afs))
  expect_true("rs5" %in% rownames(res_poly$afs))
  expect_false("rs6" %in% rownames(res_poly$afs))
})

test_that("locally-polymorphic SNP (mean = 0.5) is kept under poly_only in both readers", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  if(is.na(fx$bed_pref)) skip("BED conversion failed; cannot compare")
  pf = pfile_to_afs(fx$pfile_pref, poly_only = TRUE, verbose = FALSE)
  bd = plink_to_afs(fx$bed_pref,   poly_only = TRUE, verbose = FALSE)
  # rs5 is fixed REF in popA + fixed ALT in popB — row-mean = 0.5, so the
  # legacy rowMeans filter (plink_to_afs) keeps it, and cpp_is_polymorphic
  # (pfile_to_afs after #109 review) also keeps it because the AFs differ
  # across populations.
  expect_true("rs5" %in% rownames(pf$afs))
  expect_true("rs5" %in% rownames(bd$afs))
  # pfile_to_afs uses the stricter cpp_is_polymorphic filter on purpose
  # (matches discard_from_aftable's filter), so its row set is a subset of
  # plink_to_afs's rowMeans-based filter.
  expect_true(all(rownames(pf$afs) %in% rownames(bd$afs)))
})

test_that("duplicate SNP ids in .pvar are rejected with a useful message (issue #5)", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  pvar_path = paste0(fx$pfile_pref, ".pvar")
  txt = readLines(pvar_path)
  # Force a duplicate by copying rs2's id over rs3.
  header_lineno = which(grepl("^#[^#]", txt))[1]
  header = strsplit(sub("^#", "", txt[header_lineno]), "\t", fixed = TRUE)[[1]]
  id_col = match("ID", header)
  body_idx = seq.int(header_lineno + 1L, length(txt))
  rows = strsplit(txt[body_idx], "\t", fixed = TRUE)
  rows[[3]][id_col] = rows[[2]][id_col]   # rs3 -> rs2
  txt[body_idx] = sapply(rows, paste, collapse = "\t")
  writeLines(txt, pvar_path)
  expect_error(
    pfile_to_afs(fx$pfile_pref, verbose = FALSE),
    "duplicate SNP id"
  )
})

test_that("ploidy divisor produces correct AFs for a fully-diploid fixture (issue #10)", {
  # With every sample diploid, ploidy_kept = 2 and (3 - ploidy_kept) = 1, so
  # geno values pass through unchanged into the rowsum; the divisor in
  # counts_block is 2 per sample. For a SNP with 2/2 ALT homozygotes in popB
  # (3 samples) the AF must be exactly 1.0 — anything off-by-one in the
  # ploidy math would show up as 0.5 or NaN.
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  res = pfile_to_afs(fx$pfile_pref, verbose = FALSE)
  # rs5 is hand-set to all-REF in popA, all-ALT in popB.
  expect_equal(res$afs["rs5", "popA"], 0)
  expect_equal(res$afs["rs5", "popB"], 1)
  # counts: 3 diploid samples per pop -> 6 per row
  expect_equal(res$counts["rs5", "popA"], 6)
  expect_equal(res$counts["rs5", "popB"], 6)
})

test_that(".read_pvar finds the header line within the 200-line scan window (issue #7)", {
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  # Augment the .pvar with ~50 extra metadata lines to confirm the single-pass
  # header detection still locates the #CHROM line. (Real .pvar files emitted
  # by plink2 have only a handful of ## lines, so 50 is comfortably within
  # the 200-line scan window.)
  pvar_path = paste0(fx$pfile_pref, ".pvar")
  txt = readLines(pvar_path)
  header_lineno = which(grepl("^#[^#]", txt))[1]
  extra = paste0("##bogus", seq_len(50))
  txt = c(txt[seq_len(header_lineno - 1L)],
          extra,
          txt[seq.int(header_lineno, length(txt))])
  writeLines(txt, pvar_path)

  res = pfile_to_afs(fx$pfile_pref, verbose = FALSE)
  expect_equal(nrow(res$afs), 20L)
})

test_that(".read_pvar falls back to a full read when the header is past line 200", {
  # >200 metadata lines forces the full-file fallback path.
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  pvar_path = paste0(fx$pfile_pref, ".pvar")
  txt = readLines(pvar_path)
  header_lineno = which(grepl("^#[^#]", txt))[1]
  extra = paste0("##bogus", seq_len(400))
  txt = c(txt[seq_len(header_lineno - 1L)],
          extra,
          txt[seq.int(header_lineno, length(txt))])
  writeLines(txt, pvar_path)
  res = pfile_to_afs(fx$pfile_pref, verbose = FALSE)
  expect_equal(nrow(res$afs), 20L)
})

test_that("anygeno_to_aftable auto-dispatches PFILE when .pvar.zst is at the prefix", {
  # Regression test for the dispatch fix: previously the PFILE branch in
  # anygeno_to_aftable() required plaintext .pvar, so a prefix with .pvar.zst
  # fell through to `stop('Genotype files not found!')`.
  if(Sys.which("zstd") == "") skip("zstd CLI not on PATH")
  dir = withr::local_tempdir()
  fx = build_pfile_fixture(dir, with_multi = FALSE)
  pvar_path = paste0(fx$pfile_pref, ".pvar")
  rc = system2("zstd", args = c("-q", "--rm", "-f", pvar_path), stdout = NULL, stderr = NULL)
  expect_equal(rc, 0L)
  expect_true(file.exists(paste0(pvar_path, ".zst")))
  expect_false(file.exists(pvar_path))
  # anygeno_to_aftable is internal — call via admixtools:::.
  res = admixtools:::anygeno_to_aftable(fx$pfile_pref, verbose = FALSE)
  expect_equal(nrow(res$afs), 20L)
})

test_that("pfile_to_afs reads .pvar.zst end-to-end and matches the plaintext path", {
  if(Sys.which("zstd") == "") skip("zstd CLI not on PATH")
  dir_plain = withr::local_tempdir()
  dir_zst   = withr::local_tempdir()
  fx_plain = build_pfile_fixture(dir_plain, with_multi = FALSE)
  fx_zst   = build_pfile_fixture(dir_zst,   with_multi = FALSE)
  pvar_zst = paste0(fx_zst$pfile_pref, ".pvar")
  rc = system2("zstd", args = c("-q", "--rm", "-f", pvar_zst), stdout = NULL, stderr = NULL)
  expect_equal(rc, 0L)

  res_plain = pfile_to_afs(fx_plain$pfile_pref, verbose = FALSE)
  res_zst   = pfile_to_afs(fx_zst$pfile_pref,   verbose = FALSE)
  expect_equal(res_plain$afs,    res_zst$afs)
  expect_equal(res_plain$counts, res_zst$counts)
})
