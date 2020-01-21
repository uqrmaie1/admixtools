
bj_arr_lo_mean = function(arr, block_lengths, loo = TRUE, verbose = TRUE) {
  # returns leave-one-block-out array means
  # arr is 3d
  # group over third dimension
  # used to return rray

  blockids = rep(seq_along(block_lengths), block_lengths)
  numblocks = length(block_lengths)
  tot = apply(arr, 1:2, sum, na.rm = TRUE)
  tot_nonmiss = apply(arr, 1:2, function(x) sum(!is.na(x)))
  parts = parts_nonmiss = array(NA, c(dim(arr)[1], dim(arr)[2], numblocks))
  for(i in 1:numblocks) {
    #if(verbose) cat(paste0('\r',  i, ' out of ', numblocks, ' SNP blocks processed...'))
    bid = which(blockids == i)
    parts[,,i] = apply(arr[,, bid, drop=F], 1:2, sum, na.rm = TRUE)
    parts_nonmiss[,,i] = apply(arr[,, bid, drop=F], 1:2, function(x) sum(!is.na(x)))
  }
  #if(verbose) cat('\n')
  if(loo) {
    sums = rray(tot) - rray(parts)
    nonmiss = rray(tot_nonmiss) - rray(parts_nonmiss)
  } else {
    sums = parts
    nonmiss = parts_nonmiss
  }
  out = as.array(sums/nonmiss)
  if(!is.null(dimnames(arr)[[1]])) dimnames(out)[[1]] = dimnames(arr)[[1]]
  if(!is.null(dimnames(arr)[[2]])) dimnames(out)[[2]] = dimnames(arr)[[2]]
  out
}

block_arr_mean = function(arr, block_lengths, verbose = TRUE) {
  # returns grouped array means
  # arr is 3d
  # group over third dimension

  blockids = rep(seq_along(block_lengths), block_lengths)
  numblocks = length(block_lengths)
  sums = nonmiss = arr[,,1:numblocks]*NA
  for(i in 1:numblocks) {
    bid = which(blockids == i)
    sums[,,i] = apply(arr[,, bid, drop=F], 1:2, sum, na.rm = TRUE)
    nonmiss[,,i] = apply(arr[,, bid, drop=F], 1:2, function(x) sum(!is.na(x)))
  }
  sums/nonmiss
}



bj_mat_stats = function(bj_lo_mat, block_lengths) {
  # input is matrix (one block per column)
  # output is list with vector of jackknife means and matrix of pairwise jackknife covariances
  # uses mean jackknife estimate instead of overall mean; probably makes very little difference
  # lengths normalization used to have -1, removed this because it resulted in negative values; where did this come from?
  numblocks = length(block_lengths)
  jest = rowMeans(bj_lo_mat)
  mnc = t(jest - bj_lo_mat) * sqrt((sum(block_lengths)/block_lengths-1)/numblocks)
  jvar = crossprod(mnc)
  namedList(jest, jvar)
}

bj_arr_stats = function(bj_lo_arr, block_lengths) {
  # input is 3d array (n x n x m)
  # output is list with jackknife means and jackknife variances
  # uses mean jackknife estimate instead of overall mean; probably makes very little difference

  numblocks = length(block_lengths)
  jest = apply(bj_lo_arr, 1:2, weighted.mean, block_lengths)
  # changed the previous line without testing. used to be:
  # jest = apply(bj_lo_arr, 1:2, mean)
  xtau = (rray(jest) - bj_lo_arr)^2 * rray(sum(block_lengths)/block_lengths-1, c(1, 1, numblocks))
  jvar = apply(xtau, 1:2, weighted.mean, block_lengths)
  # changed the previous line without testing. used to be:
  # jvar = apply(xtau, 1:2, mean)
  namedList(jest, jvar)
}

bj_pairarr_stats = function(bj_lo_arr, block_lengths) {
  # input is 3d array (m x n x p)
  # output is list with jackknife means and jackknife covariances
  # uses mean jackknife estimate instead of overall mean; probably makes very little difference
  numblocks = length(block_lengths)
  jest = c(t(apply(bj_lo_arr, 1:2, weighted.mean, block_lengths)))
  # changed the previous line without testing. used to be:
  # jest = c(t(apply(bj_lo_arr, 1:2, mean)))
  bj_lo_mat = bj_lo_arr %>% aperm(c(2,1,3)) %>% arr3d_to_mat
  mnc = t(jest - bj_lo_mat) * sqrt((sum(block_lengths)/block_lengths-1)/numblocks)
  jvar = crossprod(mnc)
  jest = t(matrix(jest, dim(bj_lo_arr)[2]))
  rownames(jest) = dimnames(bj_lo_arr)[[1]]
  colnames(jest) = dimnames(bj_lo_arr)[[2]]
  namedList(jest, jvar)
}

#' Compute all pairwise f2 statistics
#'
#' This function takes a data frame with allele frequencies (see \code{\link{packedancestrymap_to_aftable}}) and computes block jackknife f2 statistics for all population pairs. \eqn{f2} for each SNP is computed as \eqn{(P_1 - P_2)^2 - P_1 (1 - P_1)/(2 C_1 - 1) - P_2 (1 - P_2)/(2 C_2 - 1)}, where \eqn{P_1} and \eqn{P_2} are allele frequencies in populations \eqn{1} and \eqn{2}, and \eqn{C_1} and \eqn{C_2} is the average number of non-missing individuals in populations \eqn{1} and \eqn{2}. See \code{details}
#' @export
#' @param afs data frame of allele frequencies for each population. column 1 to \code{infocols} are SNP annotation columns, the other columns are allele frequencies.
#' @param popcounts named vector with number of samples for each population.
#' @param block_lengths vector with lengths of each jackknife block. \code{sum(block_lengths)} has to match \code{nrow(afs)}. See \code{\link{get_block_lengths}}
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param maxmem split up allele frequency data into blocks, if memory requirements exceed \code{maxmem} MB.
#' @param infocols number of initil columns with meta data
#' @param outdir directory into which to write f2 data (if \code{NULL}, data is returned instead)
#' @param overwrite should existing files be overwritten? only relevant if \code{outdir} is not \code{NULL}
#' @param loo if \code{TRUE} (the default), return the total estimates minus the estimates from each block. if \code{FALSE}, return the estimate from each block.
#' @param verbose print progress updates
#' @details For each population pair, each of the \eqn{i = 1, \ldots, n} resutling values (\eqn{n} is around 700 in practice) is the mean \eqn{f2} estimate across all SNPs except the ones in block \eqn{i}.
#'
#' \eqn{- P_1 (1 - P_1)/(2 C_1 - 1) - P_2 (1 - P_2)/(2 C_2 - 1)} is a correction term which makes the estimates unbiased at low sample sizes.
#' @return a 3d array of dimensions npop x npop x nblocks with all pairwise leave-one-out f2 stats, leaving each block out at a time
#' @examples
#' \dontrun{
#' f2_blocks = afs_to_f2_blocks(afs, popcounts, block_lengths)
#' }
afs_to_f2_blocks = function(afs, popcounts, block_lengths, f2_denom=1, maxmem=1e3,
                            infocols = 6, outdir = NULL, overwrite = FALSE, loo = TRUE, verbose = TRUE) {

  if('data.frame' %in% class(afs)) afs %<>% select(-seq_len(infocols)) %>% as.matrix
  #popcounts = as.vector(popcounts[colnames(afs)])

  mem1 = lobstr::obj_size(afs)
  mem2 = mem1*ncol(afs)
  numsplits = ceiling(mem2/1e6/maxmem)
  width = ceiling(ncol(afs)/numsplits)
  starts = seq(1, ncol(afs), width)
  numsplits2 = length(starts)
  ends = c(lead(starts)[-numsplits2]-1, ncol(afs))

  if(verbose) {
    alert_info(paste0('allele frequency matrix for ', nrow(afs), ' SNPs and ', ncol(afs), ' populations is ', round(mem1/1e6), ' MB\n'))
    alert_warning(paste0('matrix of pairwise f2 for all SNPs and population pairs requires ', round(mem2/1e6), ' MB\n'))
    if(numsplits2 > 1) alert_info(paste0('splitting into ', numsplits2, ' blocks of ', width, ' populations and up to ', maxmem, ' MB (', choose(numsplits2+1,2), ' block pairs)\n'))
  }

  f2_blocks = get_split_f2_blocks(afs, popcounts, block_lengths, starts=starts, ends=ends, outdir = outdir, overwrite = overwrite, f2_denom = f2_denom, loo = loo, verbose = verbose)
  f2_blocks
}



set_blocks = function(dat, dist = 0.05, distcol = 'cm') {

  sb = function(cumpos, CHR, POS) {
    o = cumpos
    cumpos[3] = POS
    cumpos[1] = o[2]
    cumpos[2] = pmax(0, POS-o[3]+o[2])
    if(o[2] %% dist < o[1] %% dist && o[2] > dist) cumpos[2] = POS-o[3]
    cumpos
  }

  newdat = do.call(rbind,
                   accumulate2(.x=dat$CHR, .y=dat[[distcol]], .f=sb,
                               .init=c(0, 0, dat[[distcol]][1]))) %>%
    as_tibble(.name_repair = NULL)
  dat %<>% bind_cols(newdat %>% slice(-1))
  dat %>% mutate(newblock = .data$V2 > lead(.data$V2, default=0) &
                   .data$CHR == lead(.data$CHR, default=0) |
                   .data$CHR > lag(.data$CHR, default=0),
                 block=cumsum(.data$newblock))

}

#' Find LD-independent blocks
#'
#' A new block begins at the SNP after the first SNP which is not within `dist` of the start of the last block. `afdat` needs to be ordered first by 'CHR', then by 'POS' or 'cm'
#' @export
#' @param afdat data frame with columns 'CHR' and either 'POS' or 'cm'
#' @param dist minimum distance between blocks
#' @param distcol which column should be used as distance column
#' @return a numeric vector where the ith element lists the number of SNPs in the ith block.
#' @examples
#' \dontrun{
#' prefix = 'path/to/packedancestrymap_prefix'
#' pops = c('Zerg', 'Protoss', 'Terran')
#' afdat = packedancestrymap_to_aftable(prefix, pops)
#' block_lengths = get_block_lengths(afdat)
#' }
get_block_lengths = function(afdat, dist = 0.05, distcol = 'cm') {
  afdat %>%
    set_blocks(dist = dist, distcol = 'cm') %$% block %>%
    rle %$% lengths
}

get_split_f2_blocks = function(afmat, popcounts, block_lengths, starts, ends, outdir = NULL,
                               overwrite = FALSE, f2_denom = 1, loo = TRUE, verbose = TRUE) {
  # splits afmat into blocks by column, computes lo jackknife blocks on each pair of blocks, and combines into 3d array
  numsplits2 = length(starts)
  cmb = combn(0:numsplits2, 2)+(1:0)
  arrlist = replicate(numsplits2, list())
  nsnp = nrow(afmat)

  for(i in 1:ncol(cmb)) {
    if(numsplits2 > 1 & verbose) cat(paste0('\rpop pair block ', i, ' out of ', ncol(cmb)))
    c1 = cmb[1,i]
    c2 = cmb[2,i]
    s1 = starts[c1]:ends[c1]
    s2 = starts[c2]:ends[c2]
    arrs = mats_to_f2arr(afmat[,s1,drop=F], afmat[,s2,drop=F], popcounts[,s1,drop=F], popcounts[,s2,drop=F])
    numer = arrs[[1]]
    denom = arrs[[2]]
    f2_subblock = bj_arr_lo_mean(numer, block_lengths, loo = loo)
    #f2_denom = bj_arr_lo_mean(denom, block_lengths, loo = loo)
    f2_subblock = f2_subblock / f2_denom
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
    f2_blocks = do.call(abind, list(lapply(arrlist, function(x) do.call(abind, list(x, along=2))), along=1)) %>%
      structure(block_lengths = block_lengths)
    dimnames(f2_blocks)[[1]] = dimnames(f2_blocks)[[2]] = colnames(afmat)
    return(f2_blocks)
  }
}

mats_to_f2arr = function(afmat1, afmat2, countmat1, countmat2) {

  stopifnot(all.equal(nrow(afmat1), nrow(afmat2), nrow(countmat1), nrow(countmat2)))
  stopifnot(all.equal(ncol(afmat1), ncol(countmat1)))
  stopifnot(all.equal(ncol(afmat2), ncol(countmat2)))

  nsnp = nrow(afmat1)
  d1 = c(ncol(afmat1), 1, nsnp)
  d2 = c(1, ncol(afmat2), nsnp)
  afrr1 = rray(t(afmat1), dim = d1)
  afrr2 = rray(t(afmat2), dim = d2)
  #denom1 = rray(t(pmax(-Inf, countmat1-1)), dim = d1)
  #denom2 = rray(t(pmax(-Inf, countmat2-1)), dim = d2)
  denom1 = pmax(1, colMeans(countmat1, na.rm = T)-1)
  denom2 = t(pmax(1, colMeans(countmat2, na.rm = T)-1))
  #denom1 = pmax(1, apply(countmat1, 2, max, na.rm = T)-1)
  #denom2 = t(pmax(1, apply(countmat2, 2, max, na.rm = T)-1))
  pq1 = afrr1*(1-afrr1)
  pq2 = afrr2*(1-afrr2)
  pqarr = pq1/denom1 + pq2/denom2
  arr = (afrr1 - afrr2)^2 - pqarr
  denom = NULL
  #denom = arr + pq1*countmat1/(countmat1-1) + pq2*countmat2/(countmat2-1)
  arr = as.array(arr)
  #denom = as.array(denom)
  # maybe set pop pairs with < 2 haplotypes to missing here
  #valid = rray(countmat1 > 1, dim = d1) & rray(countmat2 > 1, dim = d2)
  #arr[which(!valid)] = NA
  #denom[which(!valid)] = NA
  dimnames(arr)[[1]] = colnames(afmat1)
  dimnames(arr)[[2]] = colnames(afmat2)
  #dimnames(denom)[[1]] = colnames(afmat1)
  #dimnames(denom)[[2]] = colnames(afmat2)
  list(arr, denom)
}

#' Compute block jackknife f2 blocks and write them to disk
#'
#' This is intended for computing f2-statistics for a large number of populations, which cannot be stored in memory. It assumes that the allele frequencies have already been computed and are stored in .RData files, split into consecutive blocks for a set of populations. This function calls \code{\link{write_f2}}, which takes a (sub-)block of pairwise f2-statistics, and writes to disk one pair at a time.
#' @export
#' @param afmatprefix prefix of the allele frequency \code{.RData} files created by \code{\link{split_afmat}}
#' @param outdir directory where the f2 blocks will be stored
#' @param block1 block number of the first block of populations
#' @param block2 block number of the second block of populations
#' @param popcounts named vector with number of samples for each population.
#' @param block_lengths vector with lengths of each jackknife block. \code{sum(block_lengths)} has to match the number of SNPs.
#' @param f2_denom scaling factor applied to f2-statistics. If set to 0.278, this will be approximately equal to Fst.
#' @param verbose print progress updates
#' @return a numeric vector where the ith element lists the number of SNPs in the ith block.
#' @seealso \code{\link{split_afmat}} for creating split allele frequency data, \code{\link{write_f2}} for writing split f2 block jackknife estimates
#' @examples
#' \dontrun{
#' afmatall = packedancestrymap_to_aftable('path/to/packedancestrymap_prefix', allpopulations,
#'                                         na.action = 'none', return_matrix = TRUE)
#' split_afmat(afmatall, pops_per_block = 20, outprefix = 'afmatall_split_v42.1/afmatall_')
#' numblocks = 185 # this should be the number of split allele frequency files
#' for(j in 1:numblocks) {
#'   for(j in i:numblocks) {
#'     write_split_f2_block('afmatall_split_v42.1/afmatall_', 'f2blocks_v42.1/',
#'                          block1 = i, block2 = j, popcounts, block_lengths)
#'     }
#'   }
#' }
write_split_f2_block = function(afmatprefix, outdir, block1, block2, popcounts, block_lengths, f2_denom = 1, verbose = TRUE) {
  # reads two afmat blocks, computes f2 jackknife blocks, and writes output to outdir

  # load(paste0(afmatprefix, block1, '.RData'))
  # b1 = afs
  # load(paste0(afmatprefix, block2, '.RData'))
  # b2 = afs
  # rm(afs)
  b1 = readRDS(paste0(afmatprefix, block1, '.rds'))
  b2 = readRDS(paste0(afmatprefix, block2, '.rds'))
  # continue here: need to read count matrices in a similar way
  nam1 = colnames(b1)
  nam2 = colnames(b2)
  nsnp = nrow(b1)
  filenames = expand_grid(nam1, nam2) %>%
    transmute(nam = paste0(outdir, '/', pmin(nam1, nam2), '/', pmax(nam1, nam2), '.rds')) %$% nam
    #transmute(nam = paste0(outdir, '/', pmin(nam1, nam2), '/', pmax(nam1, nam2), '.RData')) %$% nam
  if(all(file.exists(filenames))) return()

  # replace below with this:
  # arr = mats_to_f2arr(b1, b2, popcounts[,s1,drop=F], popcounts[,s2,drop=F])
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  afrr1 = rray::rray(t(b1), dim=c(ncol(b1), 1, nsnp))
  afrr2 = rray::rray(t(b2), dim=c(1, ncol(b2), nsnp))
  pqarr = afrr1*(1-afrr1)/(2*c(popcounts[nam1])-1) + afrr2*(1-afrr2)/t(2*c(popcounts[nam2])-1)
  arr = as.array((afrr1 - afrr2)^2 - pqarr)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rm(afrr1, afrr2, pqarr)
  f2_subblock = bj_arr_lo_mean(arr / f2_denom, block_lengths, verbose = verbose)
  rm(arr)
  dimnames(f2_subblock)[[1]] = nam1
  dimnames(f2_subblock)[[2]] = nam2
  if(block1 == block2) for(i in 1:dim(f2_subblock)[1]) f2_subblock[i,i,] = 0
  write_f2(f2_subblock, outdir = outdir)
}


est_to_loo = function(arr) {
  # turns block estimates into leave-one-block-out estimates
  # assumes blocks are along 3rd dimension, and block_lengths is attribute

  block_lengths = attr(arr, 'block_lengths')
  tot = rray(apply(arr, 1:2, weighted.mean, block_lengths, na.rm=T))
  rel_bl = rray(block_lengths/sum(block_lengths), dim = c(1,1,dim(arr)[3]))
  out = (tot - arr*rel_bl) / (1-rel_bl)
  attributes(out) = attributes(arr)
  out
}


loo_to_est = function(arr) {
  # inverse of est_to_res

  block_lengths = attr(arr, 'block_lengths')
  rel_bl = rray(block_lengths/sum(block_lengths), dim = c(1,1,dim(arr)[3]))
  tot = rray(apply(arr, 1:2, weighted.mean, 1-rel_bl, na.rm=T))
  out = as.array((tot - arr * (1-rel_bl))/rel_bl)
  attributes(out) = attributes(arr)
  out
}
