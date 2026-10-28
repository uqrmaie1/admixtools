# Regression tests for uqrmaie1/admixtools#172.
#
# Recent AADR releases (v62.0.p1, v66 onwards) ship their .geno files in the
# transposed packedancestrymap layout ("TGENO" header, one row per
# individual) instead of the original layout ("GENO" header, one row per
# SNP). The C++ readers used to infer the row size as filesize / (nsnp + 1)
# and always read SNP-major rows, so a TGENO file produced chance-level
# genotypes with the right dimensions and no error, and (when the inferred
# size rounded down) a heap overrun. The readers now take the layout from
# the header, validate the file size and the .ind/.snp row counts against
# it, and decode TGENO rows per individual.
#
# Fixtures are written in pure R by helper-packedancestrymap.R. The shapes
# below deliberately include the reporter's shape (482 SNPs x 266
# individuals, where the old inference was one byte short), GENO headers
# that are padded to 48 bytes (nind < 192), exactly 48 bytes (nind = 192)
# and larger than 48 bytes, counts that are and are not multiples of 4, and
# single-row / single-column matrices.

.read_geno_matrix = function(pref, ...) {
  read_packedancestrymap(pref, verbose = FALSE, ...)$geno
}

.read_afs = function(pref, ...) {
  packedancestrymap_to_afs(pref, verbose = FALSE, ...)
}

# R reference for packedancestrymap_to_afs on a 0/1/2/NA matrix.
.reference_afs = function(geno, pop, ploidy, pops) {
  afs = counts = matrix(NA_real_, nrow(geno), length(pops), dimnames = list(NULL, pops))
  for(p in seq_along(pops)) {
    cols = which(pop == pops[p])
    alt = rowSums(sweep(geno[, cols, drop = FALSE], 2, 3 - ploidy[cols], '/'), na.rm = TRUE)
    obs = as.numeric((!is.na(geno[, cols, drop = FALSE])) %*% ploidy[cols])
    afs[, p] = ifelse(obs > 0, alt / obs, NA_real_)
    counts[, p] = obs
  }
  list(afs = afs, counts = counts)
}

.truth = function(dat) {
  truth = dat$geno
  dimnames(truth) = list(dat$snp$SNP, dat$ind$ind)
  truth
}

shapes = list(reporter = c(nsnp = 482L, nind = 266L),
              small_padded_header = c(nsnp = 20L, nind = 6L),
              wide = c(nsnp = 57L, nind = 226L),
              unaligned = c(nsnp = 23L, nind = 13L),
              geno_rows_exact = c(nsnp = 9L, nind = 200L),      # nind %% 4 == 0, header 50 bytes
              header_tie = c(nsnp = 8L, nind = 192L),           # GENO row size exactly 48
              big_header = c(nsnp = 3L, nind = 768L),           # GENO header 192 bytes
              one_by_one = c(nsnp = 1L, nind = 1L),
              one_ind = c(nsnp = 5L, nind = 1L),
              one_snp = c(nsnp = 1L, nind = 5L))

test_that("fixtures have the file sizes the format specification implies", {
  for(nm in names(shapes)) {
    sh = shapes[[nm]]
    withr::with_tempdir({
      dat = make_packedancestrymap_dataset(sh[['nsnp']], sh[['nind']], seed = sh[['nsnp']] + sh[['nind']])
      prefs = write_both_layouts(getwd(), dat)
      expect_equal(file.size(paste0(prefs$geno, '.geno')),
                   packedancestrymap_expected_size(sh[['nsnp']], sh[['nind']], transposed = FALSE), info = nm)
      expect_equal(file.size(paste0(prefs$tgeno, '.geno')),
                   packedancestrymap_expected_size(sh[['nsnp']], sh[['nind']], transposed = TRUE), info = nm)
    })
  }
})

test_that("GENO and TGENO layouts of the same matrix read identically and match the truth", {
  for(nm in names(shapes)) {
    sh = shapes[[nm]]
    withr::with_tempdir({
      dat = make_packedancestrymap_dataset(sh[['nsnp']], sh[['nind']], seed = sh[['nsnp']] + sh[['nind']])
      prefs = write_both_layouts(getwd(), dat)
      truth = .truth(dat)
      expect_identical(.read_geno_matrix(prefs$geno), truth, info = paste(nm, 'GENO'))
      expect_identical(.read_geno_matrix(prefs$tgeno), truth, info = paste(nm, 'TGENO'))
      # `transpose` only changes the output shape, for both layouts.
      expect_identical(.read_geno_matrix(prefs$geno, transpose = TRUE), t(truth), info = nm)
      expect_identical(.read_geno_matrix(prefs$tgeno, transpose = TRUE), t(truth), info = nm)
    })
  }
})

test_that("randomised shapes, windows and subsets decode identically in both layouts", {
  # Twenty-five shapes with nsnp and nind in 1..300 cover the combinations of
  # nind %% 4, nsnp %% 4 and GENO headers on both sides of 48 bytes that the
  # fixed shapes above may miss. Seeds are fixed so failures are reproducible.
  pick = function(x) x[sample.int(length(x), 1)]
  shp = withr::with_seed(2026, data.frame(nsnp = sample(1:300, 25, replace = TRUE),
                                          nind = sample(1:300, 25, replace = TRUE)))
  withr::with_tempdir({
    for(i in seq_len(nrow(shp))) {
      dat = make_packedancestrymap_dataset(shp$nsnp[i], shp$nind[i], seed = 1000 + i, npop = 2)
      prefs = write_both_layouts(getwd(), dat, name = paste0('rand', i))
      truth = .truth(dat)
      lab = paste0(shp$nsnp[i], 'x', shp$nind[i])
      expect_identical(.read_geno_matrix(prefs$geno), truth, info = lab)
      expect_identical(.read_geno_matrix(prefs$tgeno), truth, info = lab)
      withr::with_seed(2000 + i, {
        first = pick(seq_len(shp$nsnp[i]))
        last = pick(first:shp$nsnp[i])
        keep = sample(dat$ind$ind, max(1, shp$nind[i] %/% 3))
      })
      keep_cols = sort(match(keep, dat$ind$ind))
      for(lay in c('geno', 'tgeno')) {
        expect_identical(.read_geno_matrix(prefs[[lay]], first = first, last = last, inds = keep),
                         truth[first:last, keep_cols, drop = FALSE], info = paste(lab, lay, first, last))
      }
      a_geno = .read_afs(prefs$geno, first = first, last = last)
      a_tgeno = .read_afs(prefs$tgeno, first = first, last = last)
      expect_identical(a_tgeno$afs, a_geno$afs, info = lab)
      expect_identical(a_tgeno$counts, a_geno$counts, info = lab)
    }
  })
})

test_that("TGENO rows padded to 48 bytes (plink2 --export eigt) are accepted", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(10L, 7L)
    padded = write_packedancestrymap_fixture(file.path(getwd(), 'padded'), dat$geno, dat$ind, dat$snp,
                                             transposed = TRUE, min_rowbytes = 48L)
    expect_equal(file.size(paste0(padded, '.geno')), 48 + 7 * 48)
    expect_identical(.read_geno_matrix(padded), .truth(dat))
    unpadded = write_packedancestrymap_fixture(file.path(getwd(), 'unpadded'), dat$geno, dat$ind, dat$snp,
                                               transposed = TRUE)
    expect_identical(.read_afs(padded)$afs, .read_afs(unpadded)$afs)
  })
})

test_that("plink2 round trip: files written by plink2 --export eig / eigt read back identically", {
  plink2 = .plink2_bin()   # skips when plink2 is not on PATH
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(300L, 40L)
    prefs = write_both_layouts(getwd(), dat)
    # 'nohash' because the fixtures carry dummy hash fields.
    outs = list(eig = file.path(getwd(), 'plink2_eig'), eigt = file.path(getwd(), 'plink2_eigt'))
    for(fmt in names(outs)) {
      status = system2(plink2, c('--eigfile', prefs$geno, 'nohash', '--export', fmt, '--out', outs[[fmt]]),
                       stdout = FALSE, stderr = FALSE)
      testthat::skip_if(status != 0, paste('this plink2 cannot export', fmt))
    }
    expect_equal(admixtools:::read_packedancestrymap_header(paste0(outs$eig, '.geno'))$format, 'GENO')
    expect_equal(admixtools:::read_packedancestrymap_header(paste0(outs$eigt, '.geno'))$format, 'TGENO')
    g_eig = .read_geno_matrix(outs$eig)
    g_eigt = .read_geno_matrix(outs$eigt)
    expect_identical(g_eigt, g_eig)
    # plink2 may count the other allele; accept the matrix or its complement.
    truth = unname(.truth(dat))
    expect_true(identical(unname(g_eigt), truth) || identical(unname(2L - g_eigt), truth))
    a_eig = .read_afs(outs$eig); a_eigt = .read_afs(outs$eigt)
    expect_identical(a_eigt$afs, a_eig$afs)
    expect_identical(a_eigt$counts, a_eig$counts)
  })
})

test_that("the header reports the layout and both wrappers see the same counts", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(23L, 13L)
    prefs = write_both_layouts(getwd(), dat)
    hg = admixtools:::read_packedancestrymap_header(paste0(prefs$geno, '.geno'))
    ht = admixtools:::read_packedancestrymap_header(paste0(prefs$tgeno, '.geno'))
    expect_equal(hg[c('format', 'transposed', 'nind', 'nsnp')],
                 list(format = 'GENO', transposed = FALSE, nind = 13, nsnp = 23))
    expect_equal(ht[c('format', 'transposed', 'nind', 'nsnp')],
                 list(format = 'TGENO', transposed = TRUE, nind = 13, nsnp = 23))
    # A header written with single spaces (convertf style) parses the same way.
    con = file(paste0(prefs$geno, '.geno'), 'r+b')
    writeBin(c(charToRaw('GENO 13 23 1a 2b'), raw(48 - 16)), con)
    close(con)
    expect_equal(admixtools:::read_packedancestrymap_header(paste0(prefs$geno, '.geno'))$nsnp, 23)
    # The R helper and the C++ readers share one parser, so any whitespace
    # between header fields is accepted by both, and the same file reads
    # identically through read_packedancestrymap and packedancestrymap_to_afs.
    con = file(paste0(prefs$geno, '.geno'), 'r+b')
    writeBin(c(charToRaw('GENO\t13\t23\t1a\t2b'), raw(48 - 16)), con)
    close(con)
    expect_equal(admixtools:::read_packedancestrymap_header(paste0(prefs$geno, '.geno'))$nind, 13)
    expect_identical(.read_geno_matrix(prefs$geno), .truth(dat))
    expect_identical(.read_afs(prefs$geno)$afs, .read_afs(prefs$tgeno)$afs)
    # Counts that are not non-negative integers are rejected as a bad header.
    con = file(paste0(prefs$geno, '.geno'), 'r+b')
    writeBin(c(charToRaw('GENO abc 23'), raw(48 - 11)), con)
    close(con)
    expect_error(admixtools:::read_packedancestrymap_header(paste0(prefs$geno, '.geno')),
                 'valid individual and SNP counts')
    expect_error(read_packedancestrymap(prefs$geno, verbose = FALSE), 'valid individual and SNP counts')
    # Progress messages go through cat(), so they are output rather than conditions.
    expect_output(read_packedancestrymap(prefs$tgeno, verbose = TRUE), 'TGENO layout')
  })
})

test_that("individual and population subsets and unaligned SNP ranges, both layouts", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 266L)
    prefs = write_both_layouts(getwd(), dat)
    truth = .truth(dat)
    keep = dat$ind$ind[c(266, 3, 100, 4, 199, 1, 57)]   # out of file order on purpose
    keep_cols = sort(match(keep, dat$ind$ind))
    popB = dat$ind$ind[dat$ind$pop == 'popB']
    windows = list(c(3, 11), c(5, 5), c(482, 482), c(479, 482), c(1, 4), c(4, 9), c(400, Inf), c(1, Inf))

    for(lay in c('geno', 'tgeno')) {
      pref = prefs[[lay]]
      expect_identical(.read_geno_matrix(pref, inds = keep), truth[, keep_cols], info = lay)
      expect_identical(.read_geno_matrix(pref, pops = 'popB'), truth[, popB], info = lay)
      for(rng in windows) {
        first = rng[1]; last = min(rng[2], 482)
        lab = paste(lay, first, last)
        expect_identical(.read_geno_matrix(pref, first = first, last = last),
                         truth[first:last, , drop = FALSE], info = lab)
        expect_identical(.read_geno_matrix(pref, first = first, last = last, inds = keep, transpose = TRUE),
                         t(truth[first:last, keep_cols, drop = FALSE]), info = lab)
      }
    }
  })
})

test_that("packedancestrymap_to_afs: GENO and TGENO agree with each other and with an R reference", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 266L)
    prefs = write_both_layouts(getwd(), dat)
    pops = c('popA', 'popB', 'popC')

    for(adjust in c(TRUE, FALSE)) {
      a_geno = .read_afs(prefs$geno, adjust_pseudohaploid = adjust)
      a_tgeno = .read_afs(prefs$tgeno, adjust_pseudohaploid = adjust)
      expect_identical(a_tgeno$afs, a_geno$afs, info = paste('adjust', adjust))
      expect_identical(a_tgeno$counts, a_geno$counts, info = paste('adjust', adjust))
      # popC was made pseudohaploid (no hets), so with adjustment it
      # contributes one observed allele per non-missing genotype.
      ploidy = if(adjust) ifelse(dat$ind$pop == 'popC', 1, 2) else rep(2, nrow(dat$ind))
      ref = .reference_afs(dat$geno, dat$ind$pop, ploidy, pops)
      expect_equal(as.matrix(a_tgeno$afs[, pops]), ref$afs, ignore_attr = TRUE, info = paste('adjust', adjust))
      expect_equal(as.matrix(a_tgeno$counts[, pops]), ref$counts, ignore_attr = TRUE, info = paste('adjust', adjust))
    }

    # SNP window and population subset through the R wrapper.
    w_geno = .read_afs(prefs$geno, pops = c('popA', 'popC'), first = 7, last = 250)
    w_tgeno = .read_afs(prefs$tgeno, pops = c('popA', 'popC'), first = 7, last = 250)
    expect_identical(w_tgeno$afs, w_geno$afs)
    expect_identical(w_tgeno$counts, w_geno$counts)
    expect_equal(nrow(w_tgeno$afs), 250 - 7 + 1)

    # Direct C++ calls: transpose flag and boundary windows (0-based first,
    # exclusive last, as .make_block_reader and the wrappers use them).
    indvec = as.integer(factor(dat$ind$pop, levels = pops)) - 1L
    ploidy = rep(2L, 266L)
    for(win in list(c(0L, 482L), c(2L, 11L), c(4L, 5L), c(481L, 482L), c(0L, 1L), c(3L, 3L))) {
      for(tr in c(FALSE, TRUE)) {
        r_geno = admixtools:::cpp_packedancestrymap_to_afs(paste0(prefs$geno, '.geno'), 482L, 266L, indvec,
                                                            win[1], win[2], ploidy, tr, FALSE)
        r_tgeno = admixtools:::cpp_packedancestrymap_to_afs(paste0(prefs$tgeno, '.geno'), 482L, 266L, indvec,
                                                             win[1], win[2], ploidy, tr, FALSE)
        lab = paste('window', win[1], win[2], 'transpose', tr)
        expect_identical(r_tgeno, r_geno, info = lab)
        expect_equal(dim(r_tgeno$afs), if(tr) c(3L, win[2] - win[1]) else c(win[2] - win[1], 3L), info = lab)
        if(tr) {
          r_plain = admixtools:::cpp_packedancestrymap_to_afs(paste0(prefs$tgeno, '.geno'), 482L, 266L, indvec,
                                                               win[1], win[2], ploidy, FALSE, FALSE)
          expect_identical(t(r_tgeno$afs), r_plain$afs, info = lab)
          expect_identical(t(r_tgeno$counts), r_plain$counts, info = lab)
        }
      }
    }
  })
})

test_that("cpp_packedancestrymap_ploidy detects pseudohaploid samples in both layouts", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 266L)
    prefs = write_both_layouts(getwd(), dat)
    # indvec: population index (0-based) or -1 for unused; mark popA unused.
    indvec = as.integer(factor(dat$ind$pop, levels = c('popA', 'popB', 'popC'))) - 1L
    indvec[indvec == 0L] = -1L
    for(ntest in c(1000L, 482L, 17L, 4L, 1L, 0L)) {
      pg = admixtools:::cpp_packedancestrymap_ploidy(paste0(prefs$geno, '.geno'), 482L, 266L, indvec, ntest)
      pt = admixtools:::cpp_packedancestrymap_ploidy(paste0(prefs$tgeno, '.geno'), 482L, 266L, indvec, ntest)
      expect_identical(pt, pg, info = paste('ntest', ntest))
      # R reference: 2 if any het among the first ntest SNPs, else 1; unused stay 1.
      tested = seq_len(min(ntest, nrow(dat$geno)))
      ref = ifelse(indvec == -1L, 1,
                   ifelse(colSums(dat$geno[tested, , drop = FALSE] == 1L, na.rm = TRUE) > 0, 2, 1))
      expect_equal(pt, ref, ignore_attr = TRUE, info = paste('ntest', ntest))
    }
    expect_error(admixtools:::cpp_packedancestrymap_ploidy(paste0(prefs$tgeno, '.geno'), 482L, 266L, indvec, -1L),
                 'ntest')
  })
})

test_that("extract_f2 on a TGENO prefix reproduces the GENO result end to end", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 60L)
    prefs = write_both_layouts(getwd(), dat)
    out_geno = file.path(getwd(), 'f2_geno')
    out_tgeno = file.path(getwd(), 'f2_tgeno')
    args = list(pops = c('popA', 'popB', 'popC'), maxmiss = 1, blgsize = 0.1, verbose = FALSE)
    expect_no_warning(do.call(extract_f2, c(list(prefs$geno, out_geno), args)))
    expect_no_warning(do.call(extract_f2, c(list(prefs$tgeno, out_tgeno), args)))
    f2_geno = read_f2(out_geno, verbose = FALSE)
    f2_tgeno = read_f2(out_tgeno, verbose = FALSE)
    expect_identical(dim(f2_tgeno), dim(f2_geno))
    expect_identical(f2_tgeno, f2_geno)
    expect_true(all(is.finite(f2_geno)))

    # qpfstats = TRUE goes through f4blockdat_from_geno and .make_block_reader;
    # it needs at least four populations.
    dat4 = make_packedancestrymap_dataset(482L, 60L, npop = 4)
    prefs4 = write_both_layouts(getwd(), dat4, name = 'four_pops')
    out_geno_q = file.path(getwd(), 'f2q_geno')
    out_tgeno_q = file.path(getwd(), 'f2q_tgeno')
    args_q = list(pops = c('popA', 'popB', 'popC', 'popD'), maxmiss = 1, blgsize = 0.1,
                  qpfstats = TRUE, verbose = FALSE)
    # The qpfstats path prints progress through cat() regardless of verbose.
    invisible(utils::capture.output(suppressMessages(suppressWarnings(
      do.call(extract_f2, c(list(prefs4$geno, out_geno_q), args_q))))))
    invisible(utils::capture.output(suppressMessages(suppressWarnings(
      do.call(extract_f2, c(list(prefs4$tgeno, out_tgeno_q), args_q))))))
    expect_identical(read_f2(out_tgeno_q, verbose = FALSE), read_f2(out_geno_q, verbose = FALSE))
  })
})

test_that("read_anygeno / anygeno_to_aftable / f4blockdat_from_geno auto-detect a TGENO prefix", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 60L)
    prefs = write_both_layouts(getwd(), dat)
    expect_identical(suppressMessages(admixtools:::read_anygeno(prefs$tgeno, verbose = FALSE))$geno,
                     suppressMessages(admixtools:::read_anygeno(prefs$geno, verbose = FALSE))$geno)
    at_geno = suppressMessages(admixtools:::anygeno_to_aftable(prefs$geno, verbose = FALSE))
    at_tgeno = suppressMessages(admixtools:::anygeno_to_aftable(prefs$tgeno, verbose = FALSE))
    expect_identical(at_tgeno$afs, at_geno$afs)
    expect_identical(at_tgeno$counts, at_geno$counts)
    popcombs = tibble::tibble(pop1 = 'popA', pop2 = 'popB', pop3 = 'popC', pop4 = 'popA')
    fb_geno = suppressMessages(f4blockdat_from_geno(prefs$geno, popcombs = popcombs, blgsize = 0.1, verbose = FALSE))
    fb_tgeno = suppressMessages(f4blockdat_from_geno(prefs$tgeno, popcombs = popcombs, blgsize = 0.1, verbose = FALSE))
    expect_identical(fb_tgeno, fb_geno)
  })
})

test_that("malformed files fail loudly instead of returning garbage", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 266L)
    prefs = write_both_layouts(getwd(), dat)
    header_msg = "'GENO' \\(packedancestrymap\\) or 'TGENO'"

    # Unknown header tag.
    bad = copy_triplet(prefs$tgeno, file.path(getwd(), 'badtag'))
    con = file(paste0(bad, '.geno'), 'r+b'); writeBin(charToRaw('XGENO'), con); close(con)
    expect_error(read_packedancestrymap(bad, verbose = FALSE), header_msg)
    expect_error(packedancestrymap_to_afs(bad, verbose = FALSE), header_msg)

    # Arbitrary binary content (not valid UTF-8), a text EIGENSTRAT file, and
    # a file shorter than a header.
    binary = copy_triplet(prefs$tgeno, file.path(getwd(), 'binary'))
    writeBin(as.raw(withr::with_seed(1, sample(0:255, 4000, replace = TRUE))), paste0(binary, '.geno'))
    expect_error(read_packedancestrymap(binary, verbose = FALSE), header_msg)
    expect_error(packedancestrymap_to_afs(binary, verbose = FALSE), header_msg)
    text = copy_triplet(prefs$tgeno, file.path(getwd(), 'text'))
    writeLines(apply(ifelse(is.na(dat$geno), 9L, dat$geno), 1, paste, collapse = ''), paste0(text, '.geno'))
    expect_error(read_packedancestrymap(text, verbose = FALSE), header_msg)
    tiny = copy_triplet(prefs$tgeno, file.path(getwd(), 'tiny'))
    writeBin(charToRaw('TGENO 266'), paste0(tiny, '.geno'))
    expect_error(read_packedancestrymap(tiny, verbose = FALSE), header_msg)
    expect_error(packedancestrymap_to_afs(tiny, verbose = FALSE), 'too small')

    # File size inconsistent with the header (one byte missing / one extra).
    for(lay in c('geno', 'tgeno')) {
      genofile = paste0(prefs[[lay]], '.geno')
      bytes = readBin(genofile, 'raw', file.size(genofile))
      trunc = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('trunc_', lay)))
      writeBin(bytes[-length(bytes)], paste0(trunc, '.geno'))
      expect_error(read_packedancestrymap(trunc, verbose = FALSE), 'implies', info = lay)
      expect_error(packedancestrymap_to_afs(trunc, verbose = FALSE), 'implies', info = lay)
      extra = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('extra_', lay)))
      writeBin(c(bytes, as.raw(0)), paste0(extra, '.geno'))
      expect_error(read_packedancestrymap(extra, verbose = FALSE), 'implies', info = lay)
    }

    # .ind / .snp files that disagree with the header are rejected by both
    # wrappers (the readers would otherwise index past the sample vector).
    for(lay in c('geno', 'tgeno')) {
      short_ind = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('short_ind_', lay)))
      writeLines(readLines(paste0(short_ind, '.ind'))[1:100], paste0(short_ind, '.ind'))
      expect_error(read_packedancestrymap(short_ind, verbose = FALSE), 'declares', info = lay)
      expect_error(packedancestrymap_to_afs(short_ind, verbose = FALSE), 'declares', info = lay)
      long_ind = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('long_ind_', lay)))
      cat('extra\tU\tpopA\n', file = paste0(long_ind, '.ind'), append = TRUE)
      expect_error(read_packedancestrymap(long_ind, verbose = FALSE), 'declares', info = lay)
      expect_error(packedancestrymap_to_afs(long_ind, verbose = FALSE), 'declares', info = lay)
      long_snp = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('long_snp_', lay)))
      cat('rs_extra\t1\t1.5\t999999\tA\tC\n', file = paste0(long_snp, '.snp'), append = TRUE)
      expect_error(read_packedancestrymap(long_snp, verbose = FALSE), 'declares', info = lay)
      expect_error(packedancestrymap_to_afs(long_snp, verbose = FALSE), 'declares', info = lay)
      short_snp = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('short_snp_', lay)))
      writeLines(readLines(paste0(short_snp, '.snp'))[1:400], paste0(short_snp, '.snp'))
      expect_error(read_packedancestrymap(short_snp, verbose = FALSE), 'declares', info = lay)
      expect_error(packedancestrymap_to_afs(short_snp, verbose = FALSE), 'declares', info = lay)
    }

    # A window past the end of the file.
    expect_error(admixtools:::cpp_read_packedancestrymap(paste0(prefs$tgeno, '.geno'), 482L, 266L,
                                                         rep(1L, 266), 480L, 490L, FALSE, FALSE), 'outside')
    # A sample vector of the wrong length never reaches the genotype bytes.
    expect_error(admixtools:::cpp_read_packedancestrymap(paste0(prefs$tgeno, '.geno'), 482L, 266L,
                                                         rep(1L, 100), 0L, 482L, FALSE, FALSE), 'indvec')
    expect_error(admixtools:::cpp_packedancestrymap_ploidy(paste0(prefs$geno, '.geno'), 482L, 266L,
                                                           rep(0L, 10), 100L), 'indvec')
  })
})

test_that("invalid sample vectors and SNP windows are rejected before any genotype byte is read", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 266L)
    prefs = write_both_layouts(getwd(), dat)
    for(lay in c('geno', 'tgeno')) {
      genofile = paste0(prefs[[lay]], '.geno')
      # cpp_read_packedancestrymap sizes its output from the entries equal to
      # 1, so any other nonzero entry (or NA) would write past the matrix.
      for(bad in list(rep(2L, 266), c(NA_integer_, rep(1L, 265)), c(-1L, rep(1L, 265)))) {
        expect_error(admixtools:::cpp_read_packedancestrymap(genofile, 482L, 266L, bad, 0L, 482L, FALSE, FALSE),
                     'indvec', info = lay)
      }
      # The ploidy and allele-frequency readers take -1 (unused) or a
      # population index; NA (INT_MIN) is neither.
      expect_error(admixtools:::cpp_packedancestrymap_ploidy(genofile, 482L, 266L, c(NA_integer_, rep(0L, 265)), 100L),
                   'indvec', info = lay)
      expect_error(admixtools:::cpp_packedancestrymap_to_afs(genofile, 482L, 266L, c(-2L, rep(0L, 265)), 0L, 482L,
                                                             rep(2L, 266), FALSE, FALSE), 'indvec', info = lay)
      # An inverted or non-positive window is an error, not a corrupt .snp
      # file; a window that starts past the last SNP is empty (tested below).
      expect_error(read_packedancestrymap(prefs[[lay]], first = 10, last = 5, verbose = FALSE), 'is invalid', info = lay)
      expect_error(read_packedancestrymap(prefs[[lay]], first = 500, last = 10, verbose = FALSE), 'is invalid', info = lay)
      expect_error(read_packedancestrymap(prefs[[lay]], first = 0, verbose = FALSE), 'is invalid', info = lay)
      expect_error(read_packedancestrymap(prefs[[lay]], first = c(1, 2), verbose = FALSE), 'single whole', info = lay)
      expect_error(read_packedancestrymap(prefs[[lay]], first = 2.5, verbose = FALSE), 'single whole', info = lay)
      expect_error(read_packedancestrymap(prefs[[lay]], last = NA, verbose = FALSE), 'single whole', info = lay)
      # A missing .geno is named in the error.
      nogeno = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('nogeno_', lay)))
      file.remove(paste0(nogeno, '.geno'))
      expect_error(read_packedancestrymap(nogeno, verbose = FALSE), 'Could not open genotype file .*nogeno', info = lay)
      expect_error(packedancestrymap_to_afs(nogeno, verbose = FALSE), 'Could not open genotype file .*nogeno', info = lay)
      # A .snp file that disagrees with the header is caught even when the
      # requested window fits inside it.
      short_snp = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('short_snp_win_', lay)))
      writeLines(readLines(paste0(short_snp, '.snp'))[1:100], paste0(short_snp, '.snp'))
      expect_error(read_packedancestrymap(short_snp, first = 1, last = 50, verbose = FALSE), 'has 100 rows', info = lay)
      expect_error(read_packedancestrymap(short_snp, first = 90, verbose = FALSE), 'has 100 rows', info = lay)
      expect_error(read_packedancestrymap(short_snp, first = 200, verbose = FALSE), 'has 100 rows', info = lay)
      long_snp = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('long_snp_win_', lay)))
      cat('rs_extra\t1\t1.5\t999999\tA\tC\n', file = paste0(long_snp, '.snp'), append = TRUE)
      expect_error(read_packedancestrymap(long_snp, first = 1, last = 50, verbose = FALSE), 'has 483 rows', info = lay)
      expect_error(read_packedancestrymap(long_snp, first = 400, verbose = FALSE), 'has 483 rows', info = lay)
      # Blank lines in the .snp are not rows wherever they sit (readr's `skip`
      # counts raw lines, so these exercise the whole-file fallback), in
      # either read mode, and the window is taken from the right rows.
      truth = .truth(dat)
      lines = readLines(paste0(prefs[[lay]], '.snp'))
      for(at in c(0, 5, 250, 481, 482)) {
        blank = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('blank_snp_', at, '_', lay)))
        writeLines(append(lines, '', after = at), paste0(blank, '.snp'))
        lab = paste(lay, 'blank after line', at)
        expect_identical(.read_geno_matrix(blank), truth, info = lab)
        expect_identical(.read_geno_matrix(blank, first = 3, last = 11), truth[3:11, ], info = lab)
        expect_identical(.read_geno_matrix(blank, first = 300, last = 482), truth[300:482, ], info = lab)
        expect_identical(.read_geno_matrix(blank, first = 400), truth[400:482, ], info = lab)
      }
    }
    # A header that declares no SNPs still reads as an empty matrix. Both
    # layouts are then a bare 48-byte header, so the TGENO file is the GENO
    # file with its tag rewritten.
    dat5 = make_packedancestrymap_dataset(4L, 5L)
    geno0 = write_packedancestrymap_fixture(file.path(getwd(), 'nosnps_geno'), dat5$geno[0, , drop = FALSE],
                                            dat5$ind, dat5$snp[0, ], transposed = FALSE)
    expect_equal(file.size(paste0(geno0, '.geno')), 48)
    tgeno0 = copy_triplet(geno0, file.path(getwd(), 'nosnps_tgeno'))
    con = file(paste0(tgeno0, '.geno'), 'r+b')
    writeBin(c(charToRaw('TGENO 5 0 0 0'), raw(48 - 13)), con)
    close(con)
    for(pref0 in c(geno0, tgeno0)) {
      g0 = .read_geno_matrix(pref0)
      expect_identical(dim(g0), c(0L, 5L), info = basename(pref0))
      expect_identical(colnames(g0), dat5$ind$ind, info = basename(pref0))
    }
  })
})

test_that("windows past the end are empty, so a chunked loop can stop on an empty result", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 266L)
    prefs = write_both_layouts(getwd(), dat)
    truth = .truth(dat)
    windows = list(list(first = 483), list(first = 500), list(first = 483, last = 482),
                   list(first = 600, last = 700), list(first = 5, last = 4), list(first = 1, last = 0))
    for(lay in c('geno', 'tgeno')) {
      for(w in windows) {
        lab = paste(lay, paste(names(w), unlist(w), collapse = ' '))
        res = do.call(read_packedancestrymap, c(list(prefs[[lay]], verbose = FALSE), w))
        expect_identical(dim(res$geno), c(0L, 266L), info = lab)
        expect_identical(colnames(res$geno), dat$ind$ind, info = lab)
        expect_identical(nrow(res$snp), 0L, info = lab)
        expect_identical(names(res$snp), c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2'), info = lab)
        tr = do.call(read_packedancestrymap, c(list(prefs[[lay]], transpose = TRUE, verbose = FALSE), w))
        expect_identical(dim(tr$geno), c(266L, 0L), info = lab)
      }
      # Fixed-size chunks until the first empty one reassemble the matrix.
      chunks = list()
      first = 1
      repeat {
        g = .read_geno_matrix(prefs[[lay]], first = first, last = first + 99)
        if(nrow(g) == 0) break
        chunks[[length(chunks) + 1]] = g
        first = first + 100
      }
      expect_identical(do.call(rbind, chunks), truth, info = lay)
    }
  })
})

test_that("whitespace-only .snp lines are blank lines, and CRLF files take the fast path", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 266L)
    prefs = write_both_layouts(getwd(), dat)
    truth = .truth(dat)
    for(lay in c('geno', 'tgeno')) {
      lines = readLines(paste0(prefs[[lay]], '.snp'))
      spaces = c('   ', '\t', ' \t ')
      for(k in seq_along(spaces)) for(at in c(0, 5, 481, 482)) {
        f = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('ws', k, '_', at, '_', lay)))
        writeLines(append(lines, spaces[k], after = at), paste0(f, '.snp'))
        lab = paste(lay, 'whitespace line', k, 'after line', at)
        expect_identical(.read_geno_matrix(f), truth, info = lab)
        expect_identical(.read_geno_matrix(f, first = 3, last = 11), truth[3:11, ], info = lab)
        expect_identical(.read_geno_matrix(f, first = 400), truth[400:482, ], info = lab)
      }
      # A whitespace-only line does not hide a real extra row.
      f = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('ws_extra_', lay)))
      writeLines(c(lines, '   ', 'rs_extra\t1\t1.5\t999999\tA\tC'), paste0(f, '.snp'))
      expect_error(read_packedancestrymap(f, first = 1, last = 10, verbose = FALSE), 'has 483 rows', info = lay)
      # CRLF line endings, no blank lines: clean, so only the window is parsed.
      f = copy_triplet(prefs[[lay]], file.path(getwd(), paste0('crlf_', lay)))
      con = file(paste0(f, '.snp'), 'wb')
      writeLines(lines, con, sep = '\r\n')
      close(con)
      expect_true(admixtools:::cpp_count_text_rows(paste0(f, '.snp'))$clean, info = lay)
      expect_identical(.read_geno_matrix(f, first = 3, last = 11), truth[3:11, ], info = lay)
      expect_identical(.read_geno_matrix(f), truth, info = lay)
    }
  })
})

test_that("cpp_count_text_rows counts the lines that hold text", {
  withr::with_tempdir({
    count = function(text) {
      writeBin(charToRaw(text), 'f.txt')
      admixtools:::cpp_count_text_rows('f.txt')
    }
    expect_identical(count(''), list(rows = 0, clean = TRUE))
    expect_identical(count('a\nb\n'), list(rows = 2, clean = TRUE))
    expect_identical(count('a\nb'), list(rows = 2, clean = TRUE))
    expect_identical(count('a\r\nb\r\n'), list(rows = 2, clean = TRUE))
    expect_identical(count('  a\n\tb\n'), list(rows = 2, clean = TRUE))
    expect_identical(count('a\n\nb\n'), list(rows = 2, clean = FALSE))
    expect_identical(count('a\n  \nb\n'), list(rows = 2, clean = FALSE))
    expect_identical(count('a\nb\n \t'), list(rows = 2, clean = FALSE))
    expect_identical(count('\r\n'), list(rows = 0, clean = FALSE))
    # Lines that straddle the 1 MiB read buffer, in both states.
    expect_identical(count(paste0(strrep('x', 3e6), '\n', strrep(' ', 2^20 - 1), '\n', 'y')),
                     list(rows = 2, clean = FALSE))
    expect_identical(count(paste0(strrep(' ', 2^20 + 7), 'z\n', strrep('w', 2^20), '\n')),
                     list(rows = 2, clean = TRUE))
    expect_error(admixtools:::cpp_count_text_rows('missing.txt'), 'Could not open')
  })
})

test_that("counts in error messages are printed in full, not in scientific notation", {
  withr::with_tempdir({
    # The .geno files are bare headers declaring round counts; the R-level
    # checks fire before the C++ size check would.
    header = function(text) writeBin(c(charToRaw(text), raw(48 - nchar(text))), 'big.geno')
    writeLines('ind1\tU\tpopA', 'big.ind')
    writeLines(c('rs1\t1\t0\t1\tA\tC', 'rs2\t1\t0\t2\tA\tC'), 'big.snp')
    header('GENO 1 2000000 0 0')
    expect_error(read_packedancestrymap('big', verbose = FALSE), 'declares 2000000 SNPs, but big.snp has 2 rows')
    expect_error(read_packedancestrymap('big', first = 0, verbose = FALSE), 'invalid for the 2000000 SNPs')
    header('GENO 1000000 2 0 0')
    expect_error(read_packedancestrymap('big', verbose = FALSE), 'declares 1000000 individuals')
  })
})

test_that("the .snp row count is cached per file and refreshed when the file changes", {
  withr::with_tempdir({
    dat = make_packedancestrymap_dataset(482L, 266L)
    prefs = write_both_layouts(getwd(), dat)
    truth = .truth(dat)
    snp = normalizePath(paste0(prefs$geno, '.snp'))
    calls = 0
    real = admixtools:::cpp_count_text_rows
    testthat::local_mocked_bindings(cpp_count_text_rows = function(file) { calls <<- calls + 1; real(file) },
                                    .package = 'admixtools')
    for(first in c(1, 101, 201, 301, 401)) {
      expect_identical(.read_geno_matrix(prefs$geno, first = first, last = first + 99),
                       truth[first:min(first + 99, 482), ], info = first)
    }
    expect_identical(calls, 1)
    # Appending a row changes the file, so the next read counts again and
    # reports the mismatch instead of trusting the cached count.
    cat('rs_extra\t1\t1.5\t999999\tA\tC\n', file = snp, append = TRUE)
    expect_error(read_packedancestrymap(prefs$geno, first = 1, last = 10, verbose = FALSE), 'has 483 rows')
    expect_identical(calls, 2)
  })
})
