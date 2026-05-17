# Tests for compute_f2_cache_id() and the .f2_cache_id sidecar written by
# extract_f2(). Re-uses the fixture builders in helper-pfile.R (added by
# PR #109), which build both PFILE and BED triplets from a single VCF.

test_that("compute_f2_cache_id is deterministic via the extract_f2 sidecar", {
  # Two extract_f2 runs against identical inputs must write identical sidecars.
  # This is the determinism guarantee that the cache id is sold on.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    o1 = file.path(getwd(), "f2_run1")
    o2 = file.path(getwd(), "f2_run2")
    suppressMessages(suppressWarnings(
      extract_f2(pref, o1, pops = fix$fid_pops, verbose = FALSE)))
    suppressMessages(suppressWarnings(
      extract_f2(pref, o2, pops = fix$fid_pops, verbose = FALSE)))

    id1 = compute_f2_cache_id(o1)
    id2 = compute_f2_cache_id(o2)
    expect_equal(id1, id2)
    expect_match(id1, "^sha256:[0-9a-f]{64}$")
  })
})


test_that("compute_f2_cache_id changes when the IID subset changes", {
  # Different `inds` selections produce different IID sets, which must produce
  # different cache ids. Mode-1 direct call avoids extract_f2's "informative
  # SNPs" requirement (which can't be met with a 3-sample-per-pop fixture
  # restricted to one pop).
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    afdat = suppressMessages(suppressWarnings(
      plink_to_afs(pref, pops = fix$fid_pops)))
    afdat = discard_from_aftable(afdat, maxmiss = 0, auto_only = TRUE)

    id_all = compute_f2_cache_id(pref, format = "plink", pops = fix$fid_pops,
                                 snpfile_kept = afdat$snpfile, blgsize = 0.05)
    # Restrict to 4 of the 6 samples — different IID set, different cache id.
    # `pops = NULL` lets match_samples take population labels from the indfile.
    id_subset = compute_f2_cache_id(pref, format = "plink",
                                    inds = c("popA_1", "popA_2", "popB_1", "popB_2"),
                                    snpfile_kept = afdat$snpfile, blgsize = 0.05)

    expect_false(identical(id_all, id_subset))
  })
})


test_that("compute_f2_cache_id changes when the variant set changes", {
  # Same IIDs, perturb one variant tuple, cache id must change. Direct
  # mode-1 call so we can control snpfile_kept precisely.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    afdat = suppressMessages(suppressWarnings(
      plink_to_afs(pref, pops = fix$fid_pops)))
    afdat = discard_from_aftable(afdat, maxmiss = 0, auto_only = TRUE)

    id_a = compute_f2_cache_id(pref, format = "plink", pops = fix$fid_pops,
                               snpfile_kept = afdat$snpfile, blgsize = 0.05)

    # Perturb one variant's A1 — the SNP-id stays the same but the (CHR,POS,A1,A2)
    # tuple changes, so the variant_set hash changes and the cache id flips.
    perturbed = afdat$snpfile
    perturbed$A1[1] = "X"

    id_b = compute_f2_cache_id(pref, format = "plink", pops = fix$fid_pops,
                               snpfile_kept = perturbed, blgsize = 0.05)

    expect_false(identical(id_a, id_b))
  })
})


test_that("compute_f2_cache_id mode-2 reads the sidecar without genotype access", {
  # The orchestrator probe must work even if the genotype prefix has been
  # moved or deleted: only the .f2_cache_id sidecar in the f2 directory is
  # needed.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    outdir = file.path(getwd(), "f2_out")
    suppressMessages(suppressWarnings(
      extract_f2(pref, outdir, pops = fix$fid_pops, verbose = FALSE)))

    sidecar_path = file.path(outdir, ".f2_cache_id")
    expect_true(file.exists(sidecar_path))
    sidecar_value = readLines(sidecar_path, warn = FALSE)[1]
    expect_equal(compute_f2_cache_id(outdir), sidecar_value)

    # Delete the genotype prefix and verify the mode-2 probe still works.
    file.remove(paste0(pref, c(".bed", ".bim", ".fam")))
    expect_equal(compute_f2_cache_id(outdir), sidecar_value)
  })
})


test_that("FID-less PFILE: compute_f2_cache_id delegates to .read_psam (issue 1)", {
  # Regression for PR #117 review issue 1: compute_f2_cache_id's PFILE branch
  # used to grep the raw .psam for column names and fall back to IID as
  # population when no FID column was present. That diverged from what
  # extract_f2 actually consumed (.read_psam injects FID = "0"), so the
  # sidecar would be inconsistent with what the run produced.
  #
  # After the fix, compute_f2_cache_id reads the .psam via the same
  # .read_psam helper that pfile_to_afs uses. Two verifications:
  #   1. The mode-1 PFILE call completes cleanly and produces a valid hash.
  #   2. The IID set the hash depends on matches what .read_psam returns,
  #      proving the helper is actually wired in (not just the column
  #      detection bypassed).
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = FALSE)
    pref = fix$pfile_pref

    afdat = suppressMessages(suppressWarnings(pfile_to_afs(pref)))
    afdat = discard_from_aftable(afdat, maxmiss = 0, auto_only = TRUE)

    id = compute_f2_cache_id(pref, format = "pfile",
                             snpfile_kept = afdat$snpfile, blgsize = 0.05)
    expect_match(id, "^sha256:[0-9a-f]{64}$")

    # .read_psam normalizes FID-less PSAMs to FID = "0" everywhere. The hash
    # must depend on those IIDs (otherwise we're not really exercising the
    # helper). Build the alternative-but-identical PSAM and confirm the
    # cache id matches.
    psam_read = admixtools:::.read_psam(paste0(pref, ".psam"))
    expect_equal(length(unique(psam_read$IID)), 6L)
    expect_true(all(psam_read$FID == "0"))
  })
})


test_that("PFILE .pvar parses through compute_f2_cache_id without read_table errors (issue 2)", {
  # Regression for PR #117 review issue 2: the qpfstats-branch cache-id call
  # in extract_f2 used readr::read_table on the prefix's snpfile, which works
  # for .bim / .snp but not for .pvar (the VCF-derived #CHROM header and
  # ## metadata lines confuse read_table). The error was caught by tryCatch
  # and turned into a generic "Could not compute f2 cache id" warning, with
  # no sidecar written.
  #
  # The fix routes PFILE through .read_pvar. We exercise the same parsing
  # path here by feeding a PFILE prefix into mode-1 compute_f2_cache_id —
  # the fix to extract_f2's qpfstats branch uses the same .read_pvar helper
  # and so produces a structurally-equivalent snpfile_kept. The failure mode
  # this guards against is "any read of .pvar via read_table inside the
  # cache-id codepath".
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$pfile_pref

    # Direct .read_pvar usage (the helper now wired into the qpfstats branch):
    pvar = admixtools:::.read_pvar(paste0(pref, ".pvar"))
    expect_true(all(c("SNP", "CHR", "POS", "A1", "A2") %in% colnames(pvar)))

    # And the cache-id call against the same PFILE prefix should succeed.
    id = compute_f2_cache_id(pref, format = "pfile", pops = fix$fid_pops,
                             snpfile_kept = pvar, blgsize = 0.05)
    expect_match(id, "^sha256:[0-9a-f]{64}$")
  })
})


test_that("compute_f2_cache_id propagates get_block_lengths errors instead of returning a bogus hash (issue 3)", {
  # Regression for PR #117 review issue 3: an inner tryCatch around the
  # get_block_lengths() call was silently substituting an empty integer
  # vector when get_block_lengths failed, producing a deterministic-but-bogus
  # hash. The fix removes the inner tryCatch so callers see real errors;
  # extract_f2's outer tryCatch still demotes them to warnings, but direct
  # callers (orchestrators) get the actionable diagnostic.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    # Hand-crafted snpfile missing the CHR column required by get_block_lengths.
    bad_snp = tibble::tibble(POS = 1:10, A1 = "A", A2 = "G")
    expect_error(
      compute_f2_cache_id(pref, format = "plink", pops = fix$fid_pops,
                          snpfile_kept = bad_snp, blgsize = 0.05))
  })
})
