# Tests for the PFILE backend in f4blockdat_from_geno (PR #111).
#
# The correctness claim of the PR: f4blockdat_from_geno produces bytewise-
# identical output when given a PFILE prefix vs the equivalent BED prefix.
# build_pfile_fixture() (added by PR #109) generates matched BED+PFILE
# triplets from a single VCF, which lets us exercise both backends on
# byte-equivalent input.

test_that("format_info detects PFILE prefix", {
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$bed_pref), "fixture build skipped (plink2 not on PATH)")

    l = admixtools:::format_info(fix$pfile_pref)
    expect_equal(l$format, "pfile")
    expect_true(isTRUE(l$is_pfile))
    expect_equal(l$snpend, ".pvar")
    expect_equal(l$indend, ".psam")
    expect_equal(l$genoend, ".pgen")
    # Cpp* fields are deliberately omitted -- callers that try to use them
    # on a PFILE prefix must route through .make_block_reader instead.
    expect_null(l$cpp_read_geno)
    expect_null(l$cpp_geno_to_afs)
    expect_null(l$cpp_geno_ploidy)
  })
})

test_that(".make_block_reader returns the same shape for BED and PFILE inputs", {
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$bed_pref), "fixture build skipped (plink2 not on PATH)")

    # nind_kept = all 6 samples
    indvec = rep(1L, fix$nsam)

    l_bed = admixtools:::format_info(fix$bed_pref)
    br_bed = admixtools:::.make_block_reader(l_bed, fix$bed_pref,
                                              nsnpall = fix$nsnp, nindall = fix$nsam,
                                              indvec = indvec)
    on.exit(br_bed$finalizer(), add = TRUE)

    l_pf = admixtools:::format_info(fix$pfile_pref)
    br_pf = admixtools:::.make_block_reader(l_pf, fix$pfile_pref,
                                             nsnpall = fix$nsnp, nindall = fix$nsam,
                                             indvec = indvec)
    on.exit(br_pf$finalizer(), add = TRUE)

    # Read the same block from each
    bed_block = br_bed$reader(0L, fix$nsnp)
    pf_block  = br_pf$reader(0L, fix$nsnp)

    expect_equal(dim(bed_block), dim(pf_block))
    expect_equal(dim(bed_block), c(fix$nsam, fix$nsnp))
  })
})

test_that("f4blockdat_from_geno: BED and PFILE backends produce bytewise-identical numer + cnt", {
  # The headline correctness claim of PR #111: both backends compute the
  # same per-block f4 statistics on byte-equivalent genotype input. The
  # PR validated this at 0.000e+00 over 10.5M cells on real ancient-DNA
  # data; here we pin the same equivalence on a 20-variant fixture as a
  # regression guard for the abstraction.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$bed_pref), "fixture build skipped (plink2 not on PATH)")

    # 2-left / 2-right popcombs (per-individual pops). Rewrite the .fam
    # so each sample is its own population; same trick test-qpadm-target-null.R
    # uses.
    fam_path = paste0(fix$bed_pref, ".fam")
    fam = read.table(fam_path, stringsAsFactors = FALSE)
    fam$V1 = fam$V2
    write.table(fam, fam_path, sep = " ", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    # Same for .psam (PFILE)
    psam_path = paste0(fix$pfile_pref, ".psam")
    psam_lines = readLines(psam_path)
    header_lineno = which(grepl("^#[^#]", psam_lines))[1]
    header = strsplit(sub("^#", "", psam_lines[header_lineno]), "\t", fixed = TRUE)[[1]]
    fid_col = match("FID", header)
    iid_col = match("IID", header)
    if(!is.na(fid_col) && !is.na(iid_col)) {
      body_idx = seq.int(header_lineno + 1L, length(psam_lines))
      body = strsplit(psam_lines[body_idx], "\t", fixed = TRUE)
      for(i in seq_along(body)) body[[i]][fid_col] = body[[i]][iid_col]
      psam_lines[body_idx] = sapply(body, paste, collapse = "\t")
      writeLines(psam_lines, psam_path)
    }

    popcombs = data.frame(
      pop1 = "popA_1", pop2 = "popA_2",
      pop3 = "popB_1", pop4 = "popB_2",
      stringsAsFactors = FALSE
    )

    # Use blgsize=100 (bp distance) since the fixture has no cm data.
    res_bed = suppressMessages(suppressWarnings(
      f4blockdat_from_geno(fix$bed_pref, popcombs = popcombs,
                           allsnps = TRUE, blgsize = 100, verbose = FALSE)
    ))
    res_pf = suppressMessages(suppressWarnings(
      f4blockdat_from_geno(fix$pfile_pref, popcombs = popcombs,
                           allsnps = TRUE, blgsize = 100, verbose = FALSE)
    ))

    # Bytewise-identical numer (the per-block f4 estimate) and n (count of
    # informative SNPs per block). If this ever fails, the BED-vs-PFILE
    # abstraction has leaked a difference somewhere downstream of
    # .make_block_reader.
    expect_identical(res_bed$est, res_pf$est)
    expect_identical(res_bed$n,   res_pf$n)
  })
})

test_that("f4blockdat_from_geno on PFILE emits the cm = 0 warning text when no cm data is available", {
  # PR #111 deliberately surfaces the cm = 0 / blgsize < 100 collapse
  # condition because get_block_lengths(blgsize < 100) on cm = 0 produces
  # a single block, which would give a nonsense jackknife. This warning is
  # the user-visible signal that they need to either supply cm_file or
  # switch to bp-distance blgsize. alert_warning routes through cat() to
  # stdout (it's not an R-level warning), so we capture stdout to inspect.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE, with_cm_col = FALSE)
    testthat::skip_if(is.na(fix$pfile_pref), "fixture build skipped")

    popcombs = data.frame(
      pop1 = "popA", pop2 = "popB",
      pop3 = "popA", pop4 = "popB",
      stringsAsFactors = FALSE
    )

    out = utils::capture.output(
      suppressWarnings(
        f4blockdat_from_geno(fix$pfile_pref, popcombs = popcombs,
                             allsnps = TRUE, blgsize = 100, verbose = TRUE)
      )
    )
    combined = paste(out, collapse = "\n")
    expect_match(combined, "blgsize < 100 will collapse",       fixed = TRUE)
    expect_match(combined, "PLINK 2 .pvar does not carry cm",   fixed = TRUE)
    expect_match(combined, "cm_file = <path>",                  fixed = TRUE)
  })
})


test_that("f4blockdat_from_geno on PFILE: cm_file argument grafts cm into block partitioning", {
  # The cm_file argument lets users supply per-variant centimorgan values
  # to PFILE inputs, since .pvar doesn't carry cm in the standard format.
  # Without cm, blgsize = 0.05 (Morgan) collapses every variant into a
  # single block; with cm spanning > 0.05 Morgan, we get multiple blocks.
  # This test verifies the cm grafting end-to-end by comparing the
  # number of distinct blocks in the output with vs. without cm_file.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE, with_cm_col = FALSE)
    testthat::skip_if(is.na(fix$pfile_pref), "fixture build skipped")

    # Write a cm_file with cm values that span > 0.05 Morgan (so block
    # partitioning at blgsize=0.05 produces > 1 block). The fixture has
    # 20 SNPs at positions i*100 for i in 1..20; assign cm = POS/100
    # so the 20 SNPs span 1 to 20 cM (0.01 to 0.20 Morgan).
    pvar = admixtools:::.read_pvar(paste0(fix$pfile_pref, ".pvar"))
    cm_file = file.path(getwd(), "cm.tsv")
    cm_dat = data.frame(SNP = pvar$SNP, cm = as.numeric(pvar$POS) / 100,
                        stringsAsFactors = FALSE)
    write.table(cm_dat, cm_file, sep = "\t", row.names = FALSE, quote = FALSE)

    popcombs = data.frame(
      pop1 = "popA", pop2 = "popB",
      pop3 = "popA", pop4 = "popB",
      stringsAsFactors = FALSE
    )

    res_with_cm = suppressMessages(suppressWarnings(
      f4blockdat_from_geno(fix$pfile_pref, popcombs = popcombs,
                           allsnps = TRUE, blgsize = 0.05,
                           cm_file = cm_file, verbose = FALSE)
    ))
    res_without_cm = suppressMessages(suppressWarnings(
      f4blockdat_from_geno(fix$pfile_pref, popcombs = popcombs,
                           allsnps = TRUE, blgsize = 0.05,
                           verbose = FALSE)
    ))

    # With cm grafted: 20 SNPs span 0.20 Morgan, blgsize = 0.05 -> ~4 blocks.
    # Without cm (cm = 0 everywhere): get_block_lengths collapses to 1 block.
    expect_gt(nrow(res_with_cm),    nrow(res_without_cm))
    expect_equal(nrow(res_without_cm), 1L)
    expect_gte(nrow(res_with_cm),  3L)
  })
})


test_that("extract_f2(pfile_pref, qpfstats = TRUE) end-to-end smoke test", {
  # The user-facing entry point that PR #111 enables. extract_f2's `...`
  # forwarding routes through qpfstats() to f4blockdat_from_geno; if any
  # link in that chain is broken on PFILE, this test catches it.
  #
  # qpfstats defaults to computing f2 + f3 + f4, which needs 3+ / 4+
  # populations respectively. The fixture has 2 pops (popA, popB), so we
  # rewrite .psam to split popA's 3 samples into popA1/popA2 and popB's 3
  # into popB1/popB2 -- giving 4 pops, enough for the full qpfstats set.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    testthat::skip_if(is.na(fix$pfile_pref), "fixture build skipped")

    # Rewrite .psam to give 4 populations (popA1, popA2, popB1, popB2),
    # one or two samples each, so qpfstats has enough pops for f3 and f4.
    psam_path = paste0(fix$pfile_pref, ".psam")
    psam_lines = readLines(psam_path)
    header_lineno = which(grepl("^#[^#]", psam_lines))[1]
    header = strsplit(sub("^#", "", psam_lines[header_lineno]), "\t",
                      fixed = TRUE)[[1]]
    fid_col = match("FID", header)
    iid_col = match("IID", header)
    body_idx = seq.int(header_lineno + 1L, length(psam_lines))
    body = strsplit(psam_lines[body_idx], "\t", fixed = TRUE)
    new_fids = c("popA1", "popA2", "popA2", "popB1", "popB2", "popB2")
    for(i in seq_along(body)) body[[i]][fid_col] = new_fids[i]
    psam_lines[body_idx] = sapply(body, paste, collapse = "\t")
    writeLines(psam_lines, psam_path)

    outdir = file.path(getwd(), "f2_qpfstats")
    expect_no_error(
      suppressMessages(suppressWarnings(
        extract_f2(fix$pfile_pref, outdir,
                   pops = c("popA1", "popA2", "popB1", "popB2"),
                   qpfstats = TRUE,
                   blgsize = 100,    # bp distance; fixture has no cm
                   verbose = FALSE)
      ))
    )
    # qpfstats path writes block_lengths and f2/ap arrays.
    expect_true(file.exists(file.path(outdir, "block_lengths_f2.rds")))
    bl = readRDS(file.path(outdir, "block_lengths_f2.rds"))
    expect_true(length(bl) >= 1L)
    expect_true(sum(bl) > 0L)
  })
})
