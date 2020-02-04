
#' Compute all pairwise f2 statistics
#'
#' This function takes a data frame with allele frequencies (see \code{\link{packedancestrymap_to_aftable}})
#' and computes blocked f2 statistics for all population pairs. \eqn{f2} for each SNP is computed as
#' \eqn{(p1 - p2)^2 - p1 (1 - p1)/(n1 - 1) - p2 (1 - p2)/(n2 - 1)}, where \eqn{p1} and \eqn{p2} are
#' allele frequencies in populations \eqn{1} and \eqn{2}, and \eqn{n1} and \eqn{n2} is the average number of
#' non-missing haplotypes in populations \eqn{1} and \eqn{2}. See \code{details}
#' @export
#' @param afmat matrix of allele frequencies for all populations (columns) and SNPs (rows).
#' @param countmat matrix of allele counts for all populations (columns) and SNPs (rows).
#' @param block_lengths vector with lengths of each block. \code{sum(block_lengths)}
#' has to match \code{nrow(afmat)}. See \code{\link{get_block_lengths}}
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param maxmem split up allele frequency data into blocks, if memory requirements exceed \code{maxmem} MB.
#' @param outdir directory into which to write f2 data (if \code{NULL}, data is returned instead)
#' @param overwrite should existing files be overwritten? only relevant if \code{outdir} is not \code{NULL}
#' @param verbose print progress updates
#' @details For each population pair, each of the \eqn{i = 1, \ldots, n} resutling values
#' (\eqn{n} is around 700 in practice) is the mean \eqn{f2} estimate across all SNPs except the ones in block \eqn{i}.
#'
#' \eqn{- p1 (1 - p1)/(2 n1 - 1) - p2 (1 - p2)/(2 n2 - 1)} is a correction term which makes the estimates
#' unbiased at low sample sizes.
#' @return a 3d array of dimensions npop x npop x nblocks with all pairwise leave-one-out f2 stats,
#' leaving each block out at a time
#' @examples
#' \dontrun{
#' f2_blocks = afs_to_f2_blocks(afmat, countmat, block_lengths)
#' }
afs_to_f2_blocks = function(afmat, countmat, block_lengths, f2_denom=1, maxmem = 8000,
                            outdir = NULL, overwrite = FALSE, verbose = TRUE) {

  #if('data.frame' %in% class(afs)) afs %<>% select(-seq_len(infocols)) %>% as.matrix

  nc = ncol(afmat)
  mem1 = lobstr::obj_size(afmat)
  mem2 = mem1*nc
  numsplits = ceiling(mem2/1e6/maxmem)
  width = ceiling(nc/numsplits)
  starts = seq(1, nc, width)
  numsplits2 = length(starts)
  ends = c(lead(starts)[-numsplits2]-1, nc)

  if(verbose) {
    alert_info(paste0('allele frequency matrix for ', nrow(afmat), ' SNPs and ',
                      nc, ' populations is ', round(mem1/1e6), ' MB\n'))
    alert_warning(paste0('computing pairwise f2 for all SNPs and population pairs requires ',
                         round(mem2/1e6), ' MB RAM without splitting\n'))
    if(numsplits2 > 1) alert_info(paste0('splitting into ', numsplits2, ' blocks of ',
                                         width, ' populations and up to ', maxmem, ' MB (',
                                         choose(numsplits2+1,2), ' block pairs)\n'))
  }

  f2_blocks = get_split_f2_blocks(afmat, countmat, block_lengths, starts=starts,
                                  ends=ends, outdir = outdir, overwrite = overwrite,
                                  f2_denom = f2_denom, verbose = verbose)
  f2_blocks
}


get_split_f2_blocks = function(afmat, countmat, block_lengths, starts, ends, outdir = NULL,
                               overwrite = FALSE, f2_denom = 1, verbose = TRUE) {
  # splits afmat into blocks by column, computes snp blocks on each pair of population blocks,
  #   and combines into 3d array
  numsplits2 = length(starts)
  cmb = combn(0:numsplits2, 2)+(1:0)
  arrlist = replicate(numsplits2, list())
  nsnp = nrow(afmat)
  numblocks = length(block_lengths)

  for(i in 1:ncol(cmb)) {
    if(numsplits2 > 1 & verbose) cat(paste0('\rpop pair block ', i, ' out of ', ncol(cmb)))
    c1 = cmb[1,i]
    c2 = cmb[2,i]
    s1 = starts[c1]:ends[c1]
    s2 = starts[c2]:ends[c2]
    b1 = afmat[, s1, drop=F]
    b2 = afmat[, s2, drop=F]
    nam1 = colnames(b1)
    nam2 = colnames(b2)
    #numer = mats_to_f2arr(afmat[,s1,drop=F], afmat[,s2,drop=F], countmat[,s1,drop=F], countmat[,s2,drop=F])
    #f2_subblock = bj_arr_lo_mean(numer, block_lengths, loo = loo) / f2_denom
    f2_subblock = mats_to_f2arr(b1, b2, countmat[,s1, drop=F], countmat[,s2, drop=F]) %>%
      block_arr_mean(block_lengths) %>%
      `/`(f2_denom) %>%
      `dimnames<-`(list(nam1, nam2, paste0('l', block_lengths)))

    if(c1 == c2) for(j in 1:dim(f2_subblock)[1]) f2_subblock[j,j,] = 0
    if(!is.null(outdir)) {
      write_f2(f2_subblock, outdir = outdir, overwrite = overwrite)
    } else {
      arrlist[[c1]][[c2]] = f2_subblock
      arrlist[[c2]][[c1]] = aperm(arrlist[[c1]][[c2]], c(2,1,3))
    }
  }
  if(numsplits2 > 1 & verbose) cat('\n')
  if(is.null(outdir)) {
    f2_blocks = do.call(abind, list(lapply(arrlist, function(x) do.call(abind, list(x, along=2))), along=1))
    dimnames(f2_blocks)[[1]] = dimnames(f2_blocks)[[2]] = colnames(afmat)
    dimnames(f2_blocks)[[3]] = paste0('l', block_lengths)
    return(f2_blocks)
  }
}


mats_to_f2arr = function(afmat1, afmat2, countmat1, countmat2) {
  # Compute f2 stats for all SNPs and all population pairs from two af matrices

  stopifnot(all.equal(nrow(afmat1), nrow(afmat2), nrow(countmat1), nrow(countmat2)))
  stopifnot(all.equal(ncol(afmat1), ncol(countmat1)))
  stopifnot(all.equal(ncol(afmat2), ncol(countmat2)))

  nsnp = nrow(afmat1)
  d1 = c(ncol(afmat1), 1, nsnp)
  d2 = c(1, ncol(afmat2), nsnp)
  afrr1 = rray(t(afmat1), dim = d1)
  afrr2 = rray(t(afmat2), dim = d2)
  denom1 = pmax(1, colMeans(countmat1, na.rm = T)-1)
  denom2 = t(pmax(1, colMeans(countmat2, na.rm = T)-1))
  pq1 = afrr1*(1-afrr1)
  pq2 = afrr2*(1-afrr2)
  pqarr = pq1/denom1 + pq2/denom2
  arr = as.array((afrr1 - afrr2)^2 - pqarr)
  dimnames(arr)[[1]] = colnames(afmat1)
  dimnames(arr)[[2]] = colnames(afmat2)
  arr
}

xmats_to_pairarrs = function(xmat1, xmat2) {
  # Compute aa and nn for all SNPs and all population pairs

  nr = nrow(xmat1)
  stopifnot(nr == nrow(xmat2))
  xmat1 %<>% fix_ploidy
  xmat2 %<>% fix_ploidy
  ploidy1 = attr(xmat1, 'ploidy')
  ploidy2 = attr(xmat2, 'ploidy')
  n1 = (!is.na(xmat1)) * rep(ploidy1, each = nr)
  n2 = (!is.na(xmat2)) * rep(ploidy2, each = nr)
  xmat1 %<>% replace_na(0)
  xmat2 %<>% replace_na(0)
  aa = prodarray(xmat1, xmat2)
  nn = prodarray(n1, n2)
  dimnames(aa)[1:2] = dimnames(nn)[1:2] = list(colnames(xmat1), colnames(xmat2))
  namedList(aa, nn)
}


indpairs_to_f2blocks = function(indivs, pairs, poplist, block_lengths) {
  # creates f2_blocks from per individual data
  # make fast version in Rcpp and use tibbles here for readability
  # indivs: data frame with columns ind, bl, a, n
  # pairs: data frame with columns ind1, ind2, bl, aa, nn
  # poplist: data frame with columns ind, pop
  # block_lengths
  # a is number of alt alleles, n number of ref + alt alleles.
  stopifnot('ind' %in% names(poplist) && 'pop' %in% names(poplist))

  indsums = indivs %>% left_join(poplist, by='ind') %>%
    group_by(pop, bl) %>% summarize(a = sum(a), n = sum(n)) %>% ungroup

  pairsums = pairs %>%
    left_join(poplist %>% transmute(ind1 = ind, pop1 = pop), by = 'ind1') %>%
    left_join(poplist %>% transmute(ind2 = ind, pop2 = pop), by = 'ind2') %>%
    group_by(pop1, pop2, bl) %>%
    summarize(aa = sum(aa), nn = sum(nn), pp = aa/nn) %>% ungroup
  # check if it should rather be pp = sum(aa/nn)

  pairsums_samepop = pairsums %>% filter(pop1 == pop2) %>% transmute(pop = pop1, bl, aa, nn, pp)

  main = pairsums %>%
    left_join(pairsums_samepop %>% transmute(pop1 = pop, bl, pp1 = pp), by = c('pop1', 'bl')) %>%
    left_join(pairsums_samepop %>% transmute(pop2 = pop, bl, pp2 = pp), by = c('pop2', 'bl')) %>%
    mutate(f2uncorr = pp1 + pp2 - 2*pp)

  corr = pairsums_samepop %>% left_join(indsums, by=c('bl', 'pop')) %>%
    mutate(den1 = pmax(nn - n, 1),
           n3unfix = n*nn, n3fix = nn/n^2, n3 = n3unfix * n3fix,
           den2 = pmax(n3 - nn, 1),
           corr = a/den1 - aa/den2)

  f2dat = main %>%
    left_join(corr %>% transmute(pop1 = pop, bl, corr1 = corr), by = c('bl', 'pop1')) %>%
    left_join(corr %>% transmute(pop2 = pop, bl, corr2 = corr), by = c('bl', 'pop2')) %>%
    mutate(f2 = f2uncorr - corr1 - corr2, f2 = ifelse(pop1 == pop2, 0, f2),
           pp = paste(pmin(pop1, pop2), pmax(pop1, pop2))) %>%
    group_by(pp, bl) %>% mutate(cnt = n()) %>% ungroup %>%
    bind_rows(filter(., pop1 != pop2 & cnt == 1) %>%
                rename(pop1 = pop2, pop2 = pop1, pp1 = pp2,
                       pp2 = pp1, corr1 = corr2, corr2 = corr1)) %>%
    arrange(bl, pop2, pop1)

  popnames = unique(poplist$pop)
  popnames2 = unique(pairsums_samepop$pop)
  npop = length(popnames)
  nblock = length(block_lengths)
  nsnp = sum(block_lengths)

  array(f2dat$f2, dim = c(npop, npop, nblock),
        dimnames = list(pop1 = popnames2,
                        pop2 = popnames2,
                        bl = paste0('l', block_lengths))) %>%
    `[`(popnames, popnames, ) %>%
    ifelse(is.nan(.), NA, .)
}


fix_ploidy = function(xmat) {
  # divides pseudohaploid columns by 2
  # returns list with updated xmat, with named 1/2 ploidy vector as attribute

  ploidy = apply(xmat, 2, function(x) length(unique(na.omit(x))))-1
  maxgt = apply(xmat, 2, max, na.rm = TRUE)
  fixcols = ploidy == 1 & maxgt == 2
  xmat[, fixcols] = xmat[, fixcols]/2
  xmat %<>% structure(ploidy = ploidy)
  xmat
}


# turns f2_data (f2 dir) into f2_blocks; divides by denom
# returns f2_blocks (a)rray with block_lengths in 3rd dimension names
get_f2 = function(f2_data, pops, f2_denom = 1, rray = FALSE) {

  stopifnot(!is.character(f2_data) || dir.exists(f2_data))
  if(is.character(f2_data)) {
    f2_blocks = f2_from_precomp(f2_data, pops = pops)
  } else {
    f2_blocks = f2_data
  }

  stopifnot(all(pops %in% dimnames(f2_blocks)[[1]]))
  f2_blocks = f2_blocks[pops, pops,] / f2_denom
  if(rray) f2_blocks = rray::rray(f2_blocks, dim_names = dimnames(f2_blocks))
  f2_blocks
}


