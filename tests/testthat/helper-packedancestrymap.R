# Helpers for writing small packedancestrymap fixtures in either on-disk
# layout, in pure R (no plink2 / convertf dependency).
#
# Layouts follow DReichLab/AdmixTools transpose_info/transpose_format.txt:
#
#   header  "<GENO|TGENO> <nind> <nsnp> <indhash> <snphash>", NUL-padded to
#           the header size
#   GENO    header size = row size = max(48, ceil(nind/4)); one row per SNP,
#           2 bits per individual
#   TGENO   header size = 48; one row per individual of ceil(nsnp/4) bytes
#           (the reference writer) or max(48, ceil(nsnp/4)) bytes (plink2
#           --export eigt), 2 bits per SNP
#   codes   0/1/2 as in the EIGENSTRAT genotype value, 3 = missing (NA);
#           first code of each byte in the two high bits
#
# `geno` is an nsnp x nind matrix of 0/1/2/NA. Returns the prefix.

.pack_2bit_rows = function(codes, rowbytes) {
  # codes: nrow x ncol integer matrix in 0..3. Returns a raw vector holding
  # nrow rows of `rowbytes` bytes each, 4 codes per byte, zero-padded.
  n = nrow(codes)
  width = rowbytes * 4L
  padded = matrix(0L, n, width)
  padded[, seq_len(ncol(codes))] = codes
  bytes = padded[, seq(1L, width, 4L), drop = FALSE] * 64L +
          padded[, seq(2L, width, 4L), drop = FALSE] * 16L +
          padded[, seq(3L, width, 4L), drop = FALSE] * 4L +
          padded[, seq(4L, width, 4L), drop = FALSE]
  as.raw(t(bytes))
}

# Expected file size of a well-formed fixture, from the spec.
packedancestrymap_expected_size = function(nsnp, nind, transposed, min_rowbytes = 1L) {
  if(transposed) 48 + nind * max(min_rowbytes, ceiling(nsnp / 4))
  else (nsnp + 1) * max(48, ceiling(nind / 4))
}

write_packedancestrymap_fixture = function(pref, geno, ind, snp, transposed = FALSE, min_rowbytes = 1L) {
  stopifnot(is.matrix(geno), nrow(geno) == nrow(snp), ncol(geno) == nrow(ind))
  nsnp = nrow(geno)
  nind = ncol(geno)
  codes = geno
  codes[is.na(codes)] = 3L
  codes = matrix(as.integer(codes), nsnp, nind)

  if(transposed) {
    tag = 'TGENO'
    headerbytes = 48L
    rowbytes = as.integer(max(min_rowbytes, ceiling(nsnp / 4)))
    body = .pack_2bit_rows(t(codes), rowbytes)
  } else {
    tag = 'GENO'
    rowbytes = as.integer(max(48L, ceiling(nind / 4)))
    headerbytes = rowbytes
    body = .pack_2bit_rows(codes, rowbytes)
  }
  # Hash fields are not checked by admixtools; any hex string will do.
  header_txt = sprintf('%s %7d %7d %s %s', tag, as.integer(nind), as.integer(nsnp), 'deadbeef', 'beef')
  header = c(charToRaw(header_txt), raw(headerbytes - nchar(header_txt)))

  con = file(paste0(pref, '.geno'), 'wb')
  writeBin(c(header, body), con)
  close(con)
  writeLines(sprintf('%s\t%s\t%s', ind$ind, ind$sex, ind$pop), paste0(pref, '.ind'))
  writeLines(sprintf('%s\t%s\t%s\t%d\t%s\t%s', snp$SNP, snp$CHR, format(snp$cm, scientific = FALSE),
                     as.integer(snp$POS), snp$A1, snp$A2), paste0(pref, '.snp'))
  invisible(pref)
}

# Deterministic random dataset with a few populations, some missingness and
# one pseudohaploid population (codes only 0/2/NA, no hets). Returns a list
# with geno (nsnp x nind, 0/1/2/NA), ind, snp.
make_packedancestrymap_dataset = function(nsnp, nind, seed = 172, npop = 3, pseudohaploid_pop = 'popC') {
  withr::with_seed(seed, {
    pops = paste0('pop', LETTERS[seq_len(npop)])
    pop = rep(pops, length.out = nind)
    geno = matrix(sample(c(0L, 1L, 2L, NA), nsnp * nind, replace = TRUE,
                         prob = c(0.35, 0.25, 0.3, 0.1)), nsnp, nind)
    ph = which(pop == pseudohaploid_pop)
    if(length(ph)) {
      sub = geno[, ph, drop = FALSE]
      het = which(sub == 1L)
      sub[het] = sample(c(0L, 2L), length(het), replace = TRUE)
      geno[, ph] = sub
    }
    ind = data.frame(ind = sprintf('ind%03d', seq_len(nind)), sex = 'U', pop = pop,
                     stringsAsFactors = FALSE)
    snp = data.frame(SNP = sprintf('rs%d', seq_len(nsnp)), CHR = '1',
                     cm = seq(0, 1, length.out = nsnp), POS = seq_len(nsnp) * 1000L,
                     A1 = 'A', A2 = 'C', stringsAsFactors = FALSE)
    list(geno = geno, ind = ind, snp = snp)
  })
}

# Write the same dataset in both layouts under `dir`; returns both prefixes.
write_both_layouts = function(dir, dat, name = 'fixture') {
  geno_pref = file.path(dir, paste0(name, '_geno'))
  tgeno_pref = file.path(dir, paste0(name, '_tgeno'))
  write_packedancestrymap_fixture(geno_pref, dat$geno, dat$ind, dat$snp, transposed = FALSE)
  write_packedancestrymap_fixture(tgeno_pref, dat$geno, dat$ind, dat$snp, transposed = TRUE)
  list(geno = geno_pref, tgeno = tgeno_pref)
}

# Copy a triplet to a new prefix so a test can corrupt one file.
copy_triplet = function(from, to) {
  for(ext in c('.geno', '.ind', '.snp')) file.copy(paste0(from, ext), paste0(to, ext), overwrite = TRUE)
  invisible(to)
}
