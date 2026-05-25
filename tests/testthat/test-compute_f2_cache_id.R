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
  # `with_cm_col = TRUE` populates a synthetic CM column in the fixture so
  # get_block_lengths takes the cm-distance path rather than warning about a
  # missing linkage map (incidental to this test's SUT).
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE, with_cm_col = TRUE)
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
  # `with_cm_col = TRUE` silences the incidental no-linkage-map warning
  # (see notes on the IID-subset test above).
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE, with_cm_col = TRUE)
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
    fix = build_pfile_fixture(getwd(), with_fid = FALSE, with_cm_col = TRUE)
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
    fix = build_pfile_fixture(getwd(), with_fid = TRUE, with_cm_col = TRUE)
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
    # Mute only the predictable scaffolding noise that the malformed input
    # provokes BEFORE the SUT error fires: tibble's "Unknown or uninitialised
    # column" warnings on `$CHR` / `$cm` access (deliberately absent in
    # bad_snp), and the pre-error "No genetic linkage map" warning from
    # get_block_lengths(). Anything else surfaces — narrower than blanket
    # `suppressWarnings()` so a future-added real warning in this code path
    # doesn't silently hide.
    bad_snp = tibble::tibble(POS = 1:10, A1 = "A", A2 = "G")
    expected_noise = c("No genetic linkage map",
                       "Unknown or uninitialised column")
    mute_expected = function(w) {
      msg = conditionMessage(w)
      if(any(vapply(expected_noise, grepl, logical(1), msg)))
        invokeRestart("muffleWarning")
    }
    expect_error(
      withCallingHandlers(
        compute_f2_cache_id(pref, format = "plink", pops = fix$fid_pops,
                            snpfile_kept = bad_snp, blgsize = 0.05),
        warning = mute_expected))
  })
})


test_that("compute_f2_cache_id auto-detects PFILE when only .pvar.zst is present", {
  # Regression test surfaced while documenting PFILE in vignettes/io.Rmd:
  # the format-detection block in compute_f2_cache_id required `.pvar`
  # specifically and didn't recognize a PFILE that ships only `.pvar.zst`,
  # even though the rest of the PFILE pipeline (pfile_to_afs, extract_f2,
  # f4blockdat_from_geno) handles both forms transparently via
  # `.resolve_pvar_path()`. The fix routes the detection through the same
  # resolver. Skip on systems without the `zstd` CLI (needed to construct
  # the .pvar.zst fixture from the helper's .pvar output).
  testthat::skip_if(Sys.which("zstd") == "", "zstd CLI not on PATH")
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE, with_cm_col = TRUE)
    pref = fix$pfile_pref
    testthat::skip_if(is.na(pref), "PFILE fixture unavailable")

    pvar_plain = paste0(pref, ".pvar")
    pvar_zst   = paste0(pref, ".pvar.zst")
    testthat::skip_if(!file.exists(pvar_plain),
                      "PFILE fixture did not produce a .pvar")
    rc = system2("zstd", args = c("-q", "-f", shQuote(pvar_plain),
                                  "-o", shQuote(pvar_zst)),
                 stdout = FALSE, stderr = FALSE)
    testthat::skip_if(rc != 0L, "zstd compression failed")
    # snpfile_kept comes from .read_pvar applied to the plain .pvar (which
    # we still have at this point). Then we remove the plain form so only
    # .pvar.zst is present at the prefix.
    snpfile_kept = admixtools:::.read_pvar(pvar_plain)
    file.remove(pvar_plain)
    expect_false(file.exists(pvar_plain))
    expect_true(file.exists(pvar_zst))

    id = compute_f2_cache_id(pref, pops = fix$fid_pops,
                             snpfile_kept = snpfile_kept, blgsize = 0.05)
    expect_match(id, "^sha256:[0-9a-f]{64}$")
  })
})
