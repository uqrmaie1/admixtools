# Helpers for building small PFILE / BED fixtures for pfile_to_afs tests.
#
# Tests that need a fixture call `build_pfile_fixture()` (or its companion
# `build_bed_fixture()`) inside a `withr::with_tempdir()` (or just inside a
# tempdir of the test's own creation) and receive back a prefix path. The
# fixture is intentionally tiny (20 SNPs, 6 samples, 2 pops) so generation
# is sub-second.
#
# If `plink2` isn't on PATH, every fixture builder calls `testthat::skip()`
# so the rest of the suite still runs.

.plink2_bin = function() {
  bin = Sys.which("plink2")
  if(bin == "") testthat::skip("plink2 binary not on PATH; PFILE fixture build skipped")
  unname(bin)
}

# Write a small biallelic VCF and convert it to PFILE + BED so equivalence
# tests can compare. Returns a list with `pfile_pref` and `bed_pref`.
#
# Args:
#   dir          Output directory.
#   with_multi   If TRUE, includes one tri-allelic site (SNP_multi) so the
#                multiallelic policy paths can be exercised.
#   with_fid     If TRUE, the .psam carries an FID column with two
#                population labels (popA / popB, 3 samples each). If FALSE,
#                the .psam is written IID-only so the FID-less code path
#                is exercised.
#   with_cm_col  If TRUE, the .bim and (via plink2) the .pvar end up with
#                non-zero centimorgan values so the in-pvar CM branch is
#                exercised. If FALSE (default), cm columns stay at 0.
build_pfile_fixture = function(dir, with_multi = FALSE, with_fid = TRUE,
                               with_cm_col = FALSE) {
  plink2 = .plink2_bin()

  # Construct the genotype matrix: rows = 20 SNPs, cols = 6 samples.
  # popA samples are columns 1..3, popB samples are columns 4..6.
  # Encoding for the VCF: 0/0, 0/1, 1/1. withr::with_seed scopes the RNG
  # state to this block so tests elsewhere see an unperturbed global seed.
  nsnp = 20L
  nsam = 6L
  geno = withr::with_seed(20260517,
    matrix(sample(0:2, nsnp * nsam, replace = TRUE,
                  prob = c(0.5, 0.25, 0.25)),
           nrow = nsnp, ncol = nsam))

  # Force one across-population polymorphism with a row-mean of 0.5:
  #   popA homozygous REF, popB homozygous ALT. cpp_is_polymorphic must
  #   keep this; the old rowMeans %in% c(0,1) test also keeps it (mean is
  #   0.5), but it remains a useful regression check for issue #4.
  geno[5, ] = c(0, 0, 0, 2, 2, 2)

  # Add a globally-fixed SNP (all REF) so poly_only filters it out.
  geno[6, ] = c(0, 0, 0, 0, 0, 0)

  # One SNP with missingness in popB so the NA-handling path is touched.
  geno_missing_marker = 7L
  miss_mask = matrix(FALSE, nsnp, nsam)
  miss_mask[geno_missing_marker, 4:5] = TRUE

  vcf = file.path(dir, "fixture.vcf")
  con = file(vcf, "w")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  writeLines("##fileformat=VCFv4.2", con)
  writeLines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", con)
  for(chrom in 1) writeLines(sprintf("##contig=<ID=%d>", chrom), con)
  ids = if(with_fid) c("popA_1", "popA_2", "popA_3", "popB_1", "popB_2", "popB_3")
        else sprintf("S%d", seq_len(nsam))
  writeLines(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                     "INFO", "FORMAT", ids), collapse = "\t"), con)

  for(i in seq_len(nsnp)) {
    pos = i * 100L
    snp_id = sprintf("rs%d", i)
    if(with_multi && i == 10L) {
      # Tri-allelic site: REF=A, ALT=C,T
      snp_id = "rs_multi"
      ref = "A"; alt = "C,T"
      gts = c("0/0", "1/1", "2/2", "0/1", "0/2", "1/2")
    } else {
      ref = "A"; alt = "G"
      gts = sapply(seq_len(nsam), function(j) {
        if(miss_mask[i, j]) return("./.")
        c("0/0", "0/1", "1/1")[geno[i, j] + 1L]
      })
    }
    info = sprintf("PR")
    row = c("1", as.character(pos), snp_id, ref, alt, ".", "PASS",
            info, "GT", gts)
    writeLines(paste(row, collapse = "\t"), con)
  }
  close(con)

  pfile_pref = file.path(dir, "fixture")
  bed_pref   = file.path(dir, "fixture_bed")

  # Convert VCF -> PFILE
  args_pfile = c("--vcf", shQuote(vcf), "--make-pgen",
                 "--double-id",
                 "--out", shQuote(pfile_pref),
                 "--allow-extra-chr", "--silent")
  rc = system2(plink2, args_pfile, stdout = NULL, stderr = NULL)
  if(rc != 0) {
    testthat::skip(paste("plink2 VCF -> PFILE conversion failed (rc =", rc, ")"))
  }

  # If with_fid: rewrite .psam to put population labels in the FID column.
  # plink2 --double-id gives FID == IID; we want FID = popA / popB so that
  # the existing match_samples convention treats FID as the population.
  if(with_fid) {
    psam_path = paste0(pfile_pref, ".psam")
    txt = readLines(psam_path)
    header_lineno = which(grepl("^#[^#]", txt))[1]
    header = strsplit(sub("^#", "", txt[header_lineno]), "\t", fixed = TRUE)[[1]]
    body = strsplit(txt[-seq_len(header_lineno)], "\t", fixed = TRUE)
    fid_col = match("FID", header)
    iid_col = match("IID", header)
    pop_for = function(iid) if(startsWith(iid, "popA")) "popA" else "popB"
    body = lapply(body, function(r) { r[fid_col] = pop_for(r[iid_col]); r })
    txt[-seq_len(header_lineno)] = sapply(body, paste, collapse = "\t")
    writeLines(txt, psam_path)
  } else {
    # IID-only: strip FID from the .psam to exercise the FID-less branch.
    psam_path = paste0(pfile_pref, ".psam")
    txt = readLines(psam_path)
    header_lineno = which(grepl("^#[^#]", txt))[1]
    header = strsplit(sub("^#", "", txt[header_lineno]), "\t", fixed = TRUE)[[1]]
    body = strsplit(txt[-seq_len(header_lineno)], "\t", fixed = TRUE)
    fid_col = match("FID", header)
    iid_col = match("IID", header)
    if(!is.na(fid_col)) {
      keep_cols = setdiff(seq_along(header), fid_col)
      header_kept = header[keep_cols]
      body = lapply(body, function(r) r[keep_cols])
      txt = c(txt[seq_len(header_lineno - 1L)],
              paste0("#", paste(header_kept, collapse = "\t")),
              sapply(body, paste, collapse = "\t"))
      writeLines(txt, psam_path)
    }
  }

  # Also write a BED triplet for the equivalence test. Skip multi-allelic
  # for BED conversion (PLINK1 is biallelic-only).
  args_bed = c("--vcf", shQuote(vcf), "--make-bed",
               "--double-id",
               "--max-alleles", "2",
               "--out", shQuote(bed_pref),
               "--allow-extra-chr", "--silent")
  rc = system2(plink2, args_bed, stdout = NULL, stderr = NULL)
  if(rc != 0) {
    bed_pref = NA_character_
  } else if(with_fid) {
    # Rewrite .fam FIDs the same way as .psam.
    fam_path = paste0(bed_pref, ".fam")
    df = read.table(fam_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if(ncol(df) == 1) df = read.table(fam_path, header = FALSE,
                                      stringsAsFactors = FALSE)
    df$V1 = ifelse(startsWith(df$V2, "popA"), "popA", "popB")
    write.table(df, fam_path, sep = " ", quote = FALSE, row.names = FALSE,
                col.names = FALSE)
  }

  if(with_cm_col) {
    # Inject non-zero CM values into the .pvar header by rewriting it. plink2
    # does not synthesize CM, but downstream code paths only care that CM is
    # numeric and non-zero — we set CM = POS / 1e6 (rough genome-average rate).
    pvar_path = paste0(pfile_pref, ".pvar")
    txt = readLines(pvar_path)
    header_lineno = which(grepl("^#[^#]", txt))[1]
    header = strsplit(sub("^#", "", txt[header_lineno]), "\t", fixed = TRUE)[[1]]
    if(!"CM" %in% header) {
      pos_col = match("POS", header)
      header_new = c(header, "CM")
      txt[header_lineno] = paste0("#", paste(header_new, collapse = "\t"))
      body_idx = seq.int(header_lineno + 1L, length(txt))
      body = strsplit(txt[body_idx], "\t", fixed = TRUE)
      body = lapply(body, function(r) {
        cm = sprintf("%.6f", as.numeric(r[pos_col]) / 1e6)
        c(r, cm)
      })
      txt[body_idx] = sapply(body, paste, collapse = "\t")
      writeLines(txt, pvar_path)
    }
  }

  list(pfile_pref = pfile_pref,
       bed_pref = bed_pref,
       nsnp = nsnp,
       nsam = nsam,
       fid_pops = if(with_fid) c("popA", "popB") else NULL)
}
