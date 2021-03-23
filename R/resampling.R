block_arr_sum = function(arr, block_lengths, na.rm = TRUE) {
  # returns grouped array sums
  # arr is 3d, group over third dimension

  mat = arr3d_to_mat(arr)
  bm = block_mat_sum(mat, block_lengths, na.rm = na.rm)
  out = mat_to_arr3d(bm, dim(arr)[1])
  dimnames(out)[1:2] = dimnames(arr)[1:2]
  out
}


block_arr_mean = function(arr, block_lengths, na.rm = TRUE) {
  # returns grouped array means
  # arr is 3d, group over third dimension

  mat = arr3d_to_mat(arr)
  bm = block_mat_mean(mat, block_lengths, na.rm = na.rm)
  out = mat_to_arr3d(bm, dim(arr)[1])
  dimnames(out)[1:2] = dimnames(arr)[1:2]
  out
}

block_mat_sum = function(mat, block_lengths, na.rm = TRUE) {
  # returns grouped matrix sums
  # group over first dimension
  blockids = rep(seq_along(block_lengths), block_lengths)
  rowsum(mat, blockids, na.rm = na.rm)
}

block_mat_mean = function(mat, block_lengths, na.rm = TRUE) {
  # returns grouped matrix means
  # group over first dimension
  blockids = rep(seq_along(block_lengths), block_lengths)
  sums = rowsum(mat, blockids, na.rm = na.rm)
  nonmiss = rowsum((!is.na(mat))+0, blockids)
  replace_nan_with_na(sums/nonmiss)
}


jack_vec_stats = function(loo_vec, block_lengths, tot = NULL, na.rm = TRUE) {
  # input is a vector of leave-one-out estimates
  # output is list with jackknife mean and covariance
  # should give same results as 'jack_arr_stats' and 'jack_mat_stats'

  if(is.null(tot)) return(jack_vec_stats2(loo_vec, block_lengths, na.rm=na.rm))
  if(na.rm) {
    fin = is.finite(loo_vec)
    if(sum(fin) == 0) stop("Too many missing values!")
    block_lengths = block_lengths[fin]
    loo_vec = loo_vec[fin]
  }
  n = sum(block_lengths)
  numblocks = length(block_lengths)
  est = mean(tot - loo_vec, na.rm=na.rm)*numblocks + weighted.mean(loo_vec, block_lengths, na.rm=na.rm)
  h = n/block_lengths
  tau = h*tot - (h-1)*loo_vec
  var = mean((tau - est)^2/(h-1))

  namedList(est, var)
}

jack_vec_stats2 = function(loo_vec, block_lengths, na.rm = TRUE) {
  # input is a vector of leave-one-out estimates
  # output is list with jackknife mean and covariance

  # est == tot; only valid when estimates are additive

  if(na.rm) {
    fin = is.finite(loo_vec)
    if(sum(fin) == 0) stop("Too many missing values!")
    block_lengths = block_lengths[fin]
    loo_vec = loo_vec[fin]
  }
  n = sum(block_lengths)
  h = n/block_lengths
  est = weighted.mean(loo_vec, 1-1/h, na.rm=na.rm)
  var = mean((est - loo_vec)^2*(h-1))

  namedList(est, var)
}

jack_vec_stats3 = function(est_vec, block_lengths, na.rm = TRUE) {
  # input is a vector of blocked estimates
  # output is list with jackknife mean and covariance
  # should give same results as 'jack_arr_stats' and 'jack_mat_stats'

  # est == tot; only valid when estimates are additive

  if(na.rm) {
    fin = is.finite(est_vec)
    if(sum(fin) == 0) stop("Too many missing values!")
    block_lengths = block_lengths[fin]
    est_vec = est_vec[fin]
  }
  est = weighted.mean(est_vec, block_lengths, na.rm=na.rm) # only valid when estimates are additive
  var = mean((est - est_vec)^2/(sum(block_lengths)/block_lengths-1))

  namedList(est, var)
}


jack_mat_stats_proper_na = function(loo_mat, block_lengths, tot = NULL, na.rm = TRUE) {
  # input is matrix (one block per column)
  # output is list with vector of jackknife means and matrix of pairwise jackknife covariances
  # should give same results as 'jack_arr_stats'
  # handles NA correctly

  blmat = matrix(block_lengths, nrow(loo_mat), ncol(loo_mat), byrow = TRUE)
  blmat[is.na(loo_mat)] = NA
  hmat = rowSums(blmat,na.rm=TRUE)/blmat
  if(is.null(tot)) tot = weighted_row_means(loo_mat, 1-1/hmat, na.rm=na.rm) # only valid when estimates are additive
  numblocks = length(block_lengths)
  est = rowSums(tot - loo_mat, na.rm=na.rm) + weighted_row_means(loo_mat, block_lengths, na.rm=na.rm)
  estmat = matrix(est, nrow(loo_mat), numblocks)
  totmat = matrix(tot, nrow(loo_mat), numblocks)
  taumat = hmat*totmat - (hmat-1)*loo_mat
  xtau = (taumat - estmat) / sqrt(hmat-1)
  var = tcrossprod(replace_na(xtau, 0))/tcrossprod(!is.na(xtau))

  namedList(est, var)
}

jack_mat_stats = function(loo_mat, block_lengths, tot = NULL, na.rm = TRUE) {
  # input is matrix (one block per column)
  # output is list with vector of jackknife means and matrix of pairwise jackknife covariances
  # should give same results as 'jack_arr_stats'

  if(is.null(tot)) return(jack_mat_stats2(loo_mat, block_lengths, na.rm=na.rm))
  if(na.rm) {
    fin = apply(loo_mat, 2, function(x) sum(!is.finite(x))==0)
    if(sum(fin) == 0) stop("Too many missing values!")
    block_lengths = block_lengths[fin]
    loo_mat = loo_mat[,fin,drop=FALSE]
  }
  n = sum(block_lengths)
  numblocks = length(block_lengths)
  est = rowMeans(tot - loo_mat, na.rm=na.rm)*numblocks + weighted_row_means(loo_mat, block_lengths, na.rm=na.rm)
  estmat = matrix(est, nrow(loo_mat), numblocks)
  totmat = matrix(tot, nrow(loo_mat), numblocks)
  h = rep(n/block_lengths, each = nrow(loo_mat))
  taumat = h*totmat - (h-1)*loo_mat
  xtau = (taumat - estmat) / sqrt(h-1)
  var = tcrossprod(replace_na(xtau, 0))/tcrossprod(!is.na(xtau))

  namedList(est, var)
}

jack_mat_stats_old = function(loo_mat, block_lengths, na.rm = TRUE) {
  # input is matrix (one block per column)
  # output is list with vector of jackknife means and matrix of pairwise jackknife covariances

  # est == tot; only valid when estimates are additive

  if(na.rm) {
    fin = apply(loo_mat, 2, function(x) sum(!is.finite(x))==0)
    if(sum(fin) == 0) stop("Too many missing values!")
    block_lengths = block_lengths[fin]
    loo_mat = loo_mat[,fin,drop=FALSE]
  }
  n = sum(block_lengths)
  h = n/block_lengths
  est = weighted_row_means(loo_mat, 1-1/h, na.rm=na.rm)
  xtau = (est - loo_mat) * sqrt(rep(h, each = nrow(loo_mat))-1)
  var = tcrossprod(replace_na(xtau, 0)) / tcrossprod(!is.na(xtau))

  namedList(est, var)
}

jack_mat_stats2 = function(loo_mat, block_lengths, na.rm = TRUE) {
  # input is matrix (one block per column)
  # output is list with vector of jackknife means and matrix of pairwise jackknife covariances

  # est == tot; only valid when estimates are additive

  nr = nrow(loo_mat)
  bl_mat = matrix(rep(block_lengths, each = nr), nr)
  bl_mat[which(is.na(loo_mat))] = NA
  n = rowSums(bl_mat, na.rm=T)
  h = n/bl_mat
  est = sapply(1:nr, function(i) weighted.mean(loo_mat[i,], 1-1/h[i,], na.rm=na.rm))
  xtau = (est - loo_mat) * sqrt(h-1)
  #est = weighted_row_means(loo_mat, 1-1/h, na.rm=na.rm)
  #xtau = (est - loo_mat) * sqrt(rep(h, each = nrow(loo_mat))-1)
  var = tcrossprod(replace_na(xtau, 0)) / tcrossprod(!is.na(xtau))

  namedList(est, var)
}

jack_arr_stats_old = function(loo_arr, block_lengths, tot = NULL, na.rm = TRUE) {
  # input is 3d array (n x n x m) with leave-one-out statistics
  # output is list with jackknife means and jackknife variances
  # should give same results as 'jack_mat_stats'

  if(na.rm) {
    fin = apply(loo_arr, 3, function(x) sum(!is.finite(x))==0)
    if(sum(fin) == 0) stop("Too many missing values!")
    block_lengths = block_lengths[fin]
    loo_arr = loo_arr[,,fin,drop=FALSE]
  }
  n = sum(block_lengths)
  if(is.null(tot)) tot = apply(loo_arr, 1:2, weighted.mean, 1-block_lengths/n, na.rm=na.rm)
  # above line only valid when estimates are additive
  numblocks = length(block_lengths)
  totarr = replicate(numblocks, tot)
  d1 = dim(loo_arr)[1]
  d2 = dim(loo_arr)[2]
  est = apply(totarr - loo_arr, 1:2, mean, na.rm=na.rm)*numblocks +
    apply(loo_arr, 1:2, weighted.mean, block_lengths, na.rm=na.rm)
  estarr = replicate(numblocks, est)
  h = rep(n/block_lengths, each = d1*d2)
  tau = h*totarr - (h-1)*loo_arr
  xtau = (tau - estarr)^2 / (h-1)
  var = apply(xtau, 1:2, mean, na.rm = na.rm)

  namedList(est, var)
}

jack_arr_stats = function(loo_arr, block_lengths, tot = NULL, na.rm = TRUE) {
  # input is 3d array (n x n x m) with leave-one-out statistics
  # output is list with jackknife means and jackknife variances
  # should give same results as 'jack_mat_stats'

  d1 = dim(loo_arr)[1]
  d2 = dim(loo_arr)[2]
  d3 = dim(loo_arr)[3]
  bl_arr = (loo_arr*0+1) * rep(block_lengths, each = prod(dim(loo_arr)[1:2]))
  n = apply(bl_arr, 1:2, sum, na.rm=T)
  narr = replicate(d3, n)
  w = 1-bl_arr/narr
  if(is.null(tot)) tot = apply(loo_arr * w, 1:2, sum, na.rm=T)/apply(w, 1:2, sum, na.rm=T)
  #if(is.null(tot)) tot = apply(loo_arr, 1:2, weighted.mean, 1-block_lengths/n, na.rm=na.rm)
  # above line only valid when estimates are additive
  numblocks = apply(bl_arr, 1:2, function(x) sum(!is.na(x)))
  totarr = replicate(d3, tot)

  est = apply(totarr - loo_arr, 1:2, mean, na.rm=na.rm) * numblocks +
    apply(loo_arr*bl_arr, 1:2, sum, na.rm=T)/n
    #apply(loo_arr, 1:2, weighted.mean, block_lengths, na.rm=na.rm)

  estarr = replicate(dim(loo_arr)[3], est)
  h = narr/bl_arr
  tau = h*totarr - (h-1) * loo_arr
  xtau = (tau - estarr)^2 / (h-1)
  var = apply(xtau, 1:2, mean, na.rm = na.rm)

  namedList(est, var)
}


jack_dat_stats = function(dat, na.rm = TRUE) {
  # input is a grouped data frame with columns 'loo', 'n', and 'block'

  if(!'tot' %in% names(dat)) dat %<>% mutate(tot = weighted.mean(loo, 1-n/sum(n), na.rm=na.rm))
  dat %>%
    filter(is.finite(loo)) %>%
    mutate(h = sum(n)/n,
           est = mean(tot - loo, na.rm = na.rm)*n() + weighted.mean(loo, n, na.rm = na.rm),
           xtau = (h*tot - (h-1)*loo - est)^2/(h-1)) %>%
    summarize(est = est[1], var = mean(xtau, na.rm = na.rm), cnt = sum(n), n = sum(!is.na(xtau)))
}


boot_dat_stats = function(dat) {
  jack_dat_stats(dat) %>%
    mutate(var = var/n)
}


jack_pairarr_stats = function(loo_arr, block_lengths, tot = NULL, na.rm = TRUE) {
  # input is 3d array (m x n x p)
  # output is list with jackknife means and jackknife covariances

  if(na.rm) {
    fin = apply(loo_arr, 3, function(x) sum(!is.finite(x))==0)
    if(sum(fin) == 0) stop("Too many missing values!")
    block_lengths = block_lengths[fin]
    loo_arr = loo_arr[,,fin,drop=FALSE]
  }
  n = sum(block_lengths)
  numblocks = length(block_lengths)
  if(is.null(tot)) tot = c(t(apply(loo_arr, 1:2, weighted.mean, 1-block_lengths/n, na.rm=na.rm)))
  # above line only valid when estimates are additive
  totmat = replicate(numblocks, tot)
  loomat = loo_arr %>% aperm(c(2,1,3)) %>% arr3d_to_mat %>% t
  est = rowMeans(totmat - loomat, na.rm=na.rm)*numblocks +
    c(t(apply(loo_arr, 1:2, weighted.mean, block_lengths, na.rm=na.rm)))
  h = rep(n/block_lengths, each = nrow(loomat))
  xtau = (est * h - loomat * (h-1) - totmat) / sqrt(h-1)
  var = tcrossprod(xtau)/numblocks

  est = t(matrix(est, dim(loo_arr)[2]))
  rownames(est) = dimnames(loo_arr)[[1]]
  colnames(est) = dimnames(loo_arr)[[2]]
  namedList(est, var)
}


make_bootfun = function(jackfun) {
  function(x, y, ...) {
    out = jackfun(x, y, ...)
    out[[2]] = out[[2]]*length(y)/(length(y)-1)^2
    out
  }
}


boot_vec_stats = make_bootfun(jack_vec_stats)
boot_mat_stats = make_bootfun(jack_mat_stats)
boot_arr_stats = make_bootfun(jack_arr_stats)
boot_pairarr_stats = make_bootfun(jack_pairarr_stats)

cpp_boot_vec_stats = make_bootfun(cpp_jack_vec_stats)



#' Find LD-independent blocks
#'
#' A new block begins at the SNP after the first SNP which is not within `blgsize` of the start of the last block.
#' `dat` needs to be ordered first by 'CHR', then by 'POS' or 'cm'
#' @export
#' @param dat Data frame with columns 'CHR' and either 'POS' or 'cm'
#' @param blgsize SNP block size in Morgan. Default is 0.05 (50 cM). If `blgsize` is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
#' @param cpp Should the faster C++ version be used?
#' @param verbose Print progress updates
#' @return A numeric vector where the ith element lists the number of SNPs in the ith block.
#' @examples
#' \dontrun{
#' prefix = 'path/to/packedancestrymap_prefix'
#' pops = c('pop1', 'pop2', 'pop3')
#' afdat = packedancestrymap_to_afs(prefix, pops = pops)
#' block_lengths = get_block_lengths(afdat)
#' }
get_block_lengths = function(dat, blgsize = 0.05, cpp = TRUE, verbose = TRUE) {

  distcol = 'cm'
  if(blgsize >= 100) {
    distcol = 'POS'
    if(!distcol %in% names(dat)) stop(paste0(distcol, ' not found in column names!'))
    if(verbose) alert_warning("'blgsize' is >= 100 and interpreted as base pair distance!")
  } else if(length(unique(dat$cm)) < 2) {
    if(length(unique(dat$POS)) > 2) {
      blgsize = 2e6
      warning(paste0("No genetic linkage map found! Defining blocks by base pair distance of ", blgsize))
    } else {
      warning(paste0("No genetic linkage map or base positions found! ",
                     "Each chromosome will be its own block, which can make standard error estimates inaccurate."))
    }
  }

  numchr = parse_number(dat$CHR)
  if(any(is.na(numchr))) numchr = as.integer(as.factor(dat$CHR))
  if(cpp) return(cpp_get_block_lengths(numchr, dat[[distcol]], blgsize))

  fpos = -1e20
  lchrom = -1
  xsize = 0
  n = 0
  bsize = c()
  for(i in 1:nrow(dat)) {
    chrom = numchr[i]
    gpos = dat[[distcol]][i]
    dis = gpos - fpos
    if ((chrom != lchrom) || (dis >= blgsize)) {
      if (xsize > 0) {
        bsize[n+1] = xsize
        n = n+1
      }
      lchrom = chrom
      fpos = gpos
      xsize = 0
    }
    xsize = xsize + 1
  }
  if (xsize > 0) {
    bsize[n+1] = xsize
  }
  bsize
}



est_to_loo_old = function(arr, block_lengths = NULL) {
  # turns block estimates into leave-one-block-out estimates
  # assumes blocks are along 3rd dimension

  if(is.null(block_lengths)) block_lengths = parse_number(dimnames(arr)[[3]])
  nam = dimnames(arr)
  dm = length(dim(arr))
  if(dm == 0) arr = array(arr, c(1,1,length(arr)))
  else if(dm == 2) arr = array(c(arr), c(nrow(arr),1,ncol(arr)))
  tot = apply(arr, 1:2, weighted.mean, block_lengths, na.rm=T)
  rel_bl = rep(block_lengths/sum(block_lengths), each = length(tot))
  out = (replicate(dim(arr)[3], tot) - arr*rel_bl) / (1-rel_bl)
  if(dm == 0) out %<>% c
  else if(dm == 2) out = matrix(c(out), dim(arr)[1])
  dimnames(out) = nam
  out
}

#' Turn per-block estimates into leave-one-out estimates
#'
#' This works for any statistics which, when computed across `N` blocks, are equal
#' to the weighted mean of the statistics across the `N` blocks.
#' @export
#' @param arr 3d array with blocked estimates, with blocks in the 3rd dimension
#' @param block_lengths Optional block lengths. If `NULL`, will be parsed from 3rd dimnames in blocks
#' @return A 3d array with leave-one-out estimates for jackknife. Dimensions are equal to those of `arr`.
#' @seealso \code{\link{loo_to_est}} \code{\link{est_to_boo}}
est_to_loo = function(arr, block_lengths = NULL) {
  # turns block estimates into leave-one-block-out estimates
  # assumes blocks are along 3rd dimension

  if(is.null(block_lengths)) block_lengths = parse_number(dimnames(arr)[[3]])
  nam = dimnames(arr)
  dm = length(dim(arr))
  if(dm == 0) arr = array(arr, c(1,1,length(arr)))
  else if(dm == 2) arr = array(c(arr), c(nrow(arr),1,ncol(arr)))
  tot = apply(arr, 1:2, weighted.mean, block_lengths, na.rm=T)
  bl = rep(block_lengths, each = length(tot))
  if(any(is.na(arr))) rel_bl = bl / c(apply((arr*0+1)*bl, 1:2, sum, na.rm=T))
  else rel_bl = bl / sum(block_lengths)
  out = (replicate(dim(arr)[3], tot) - arr*rel_bl) / (1-rel_bl)
  if(dm == 0) out %<>% c
  else if(dm == 2) out = matrix(c(out), dim(arr)[1])
  dimnames(out) = nam
  out
}


#' Turn leave-one-out estimates to per-block estimates
#'
#' Inverse of \code{\link{est_to_loo}}
#' This works for any statistics which, when computed across `N` blocks, are equal
#' to the weighted mean of the statistics across the `N` blocks.
#' @export
#' @param arr 3d array with blocked estimates, with blocks in the 3rd dimension.
#' @param block_lengths Optional block lengths. If `NULL`, will be parsed from 3rd dimnames in blocks
#' @return A 3d array with leave-one-out estimates for jackknife. Dimensions are equal to those of `arr`.
#' @seealso \code{\link{est_to_loo}}
loo_to_est = function(arr, block_lengths = NULL) {
  # inverse of est_to_loo

  if(is.null(block_lengths)) block_lengths = parse_number(dimnames(arr)[[3]])
  nam = dimnames(arr)
  dm = length(dim(arr))
  if(dm == 0) arr = array(arr, c(1,1,length(arr)))
  else if(dm == 2) arr = array(c(arr), c(nrow(arr),1,ncol(arr)))
  bl = (arr*0+1) * rep(block_lengths, each = prod(dim(arr)[1:2]))
  if(any(is.na(arr))) rel_bl = bl / c(apply(bl, 1:2, sum, na.rm=T))
  else rel_bl = bl / sum(block_lengths)
  tot = apply(arr*(1-rel_bl), 1:2, sum, na.rm=T) / apply(array(1-rel_bl, dim(arr)), 1:2, sum, na.rm=T)
  out = (replicate(dim(arr)[3], tot) - arr * (1-rel_bl))/rel_bl
  if(dm == 0) out %<>% c
  else if(dm == 2) out = matrix(c(out), dim(arr)[1])
  dimnames(out) = nam
  out
}


est_to_loo_nafix = function(arr) {
  # turns block estimates into leave-one-block-out estimates
  # assumes blocks are along 3rd dimension

  block_lengths = parse_number(dimnames(arr)[[3]])
  tot = apply(arr, 1:2, weighted.mean, block_lengths, na.rm=T)
  rel_bl = rep(block_lengths/sum(block_lengths), each = length(tot))
  if(any(is.na(arr))) warning(paste0('Replacing ', sum(is.na(arr)), ' NAs with 0!'))
  arr %<>% replace_na(0)
  out = (replicate(dim(arr)[3], tot) - arr*rel_bl) / (1-rel_bl)
  dimnames(out) = dimnames(arr)
  out
}

#' Turn per-block estimates into bootstrap estimates
#'
#' This works for any statistics which, when computed across `N` blocks, are equal
#' to the weighted mean of the statistics across the `N` blocks.
#' @export
#' @param arr 3d array with blocked estimates, with blocks in the 3rd dimension.
#' @param nboot Number of bootstrap iterations
#' @return A 3d array with bootstrap estimates. The first two dimensions are equal to those of `arr`.
#'   The 3rd dimension is equal to `nboot`.
#' @seealso \code{\link{est_to_loo}}
est_to_boo = function(arr, nboot = dim(arr)[3], block_lengths = NULL) {
  # turns block estimates into bootstrap estimates
  # assumes blocks are along 3rd dimension

  dm = length(dim(arr))
  if(dm == 0) arr = array(arr, c(1,1,length(arr)))
  else if(dm == 2) arr = array(c(arr), c(nrow(arr),1,ncol(arr)))

  nboot = max(nboot, dim(arr)[3])
  if(is.null(block_lengths)) block_lengths = parse_number(dimnames(arr)[[3]])
  numblocks = length(block_lengths)
  sel = sample(seq_len(numblocks), numblocks*nboot, replace = TRUE)
  lengths = block_lengths[sel]
  grp = rep(seq_len(nboot), each = numblocks)
  arr %>%
    magrittr::multiply_by(rep(block_lengths, each = prod(dim(arr)[1:2]))) %>%
    matrix(numblocks, byrow=T) %>%
    `[`(sel,) %>%
    rowsum(grp, na.rm=T) %>%
    #`/`(numblocks) %>%
    magrittr::divide_by(c(tapply(lengths, grp, sum))) %>%
    t %>%
    array(c(dim(arr)[1:2], nboot),
          dimnames = list(dimnames(arr)[[1]], dimnames(arr)[[2]], rep('l1', nboot))) %>%
    `[`(,,)
}

#' Generate a list of leave-one-out arrays
#'
#' @export
#' @param arr 3d array with blocked estimates, with blocks in the 3rd dimension
#' @return A list of leave-one-out arrays, each of which is 1 element shorter than `arr` along the 3rd dimension
#' @seealso \code{\link{boo_list}} \code{\link{est_to_loo}}
loo_list = function(arr) map(1:dim(arr)[3], ~arr[,,-.])

#' Generate a list of bootstrap resampled arrays
#'
#' @export
#' @param arr 3d array with blocked estimates, with blocks in the 3rd dimension
#' @param nboot Number of bootstrap iterations
#' @return A list of bootstrap resampled arrays, with 3rd dimension equal to `nboot`
#' @seealso \code{\link{loo_list}} \code{\link{est_to_boo}} \code{\link{qpgraph_resample_snps2}}
boo_list = function(arr, nboot = dim(arr)[3]) {
  # returns list of bootstrap resampled arrays
  # arr is k x k x n; output is length nboot, output arrs are k x k x n
  sel = rerun(nboot, sample(1:dim(arr)[3], replace = TRUE))
  list(boo = map(sel, ~arr[,,.]), sel = sel, test = map(sel, ~arr[,,-.]))
}

boo_list_cons = function(arr, nboot = dim(arr)[3]) {
  # returns list of bootstrap resampled arrays
  # arr is k x k x n; output is length nboot, output arrs are k x k x n
  len = dim(arr)[3]
  sel = rerun(nboot, ((1:round(len/2))+sample(1:len, 1)) %% len + 1)
  list(boo = map(sel, ~arr[,,.]), sel = sel, test = map(sel, ~arr[,,-.]))
}


est_to_loo_dat = function(dat) {
  # like est_to_loo, but for a grouped data frame with columns 'est', 'block', and 'n'
  # adds column 'loo'
  dat %>%
    mutate(.rel_bl = n/sum(n),
           .tot = weighted.mean(est, n, na.rm=TRUE),
           loo = (.tot - est*.rel_bl) / (1-.rel_bl)) %>%
    select(-.rel_bl, -.tot)
}


loo_to_est_dat = function(dat) {
  # like loo_to_est, but for a grouped data frame with columns 'loo', 'block', and 'n'
  # adds column 'est'
  dat %>%
    mutate(.rel_bl = n/sum(n),
           .tot = weighted.mean(loo, 1-.rel_bl, na.rm=TRUE),
           est = (.tot - loo*(1-.rel_bl)) / .rel_bl) %>%
    select(-.rel_bl, -.tot)
}

est_to_boo_dat = function(dat, nboot = 1000) {
  # like est_to_loo, but for a grouped data frame with columns 'est', 'block', and 'length'
  # adds column 'boo'

  numblocks = length(unique(dat$block))
  boodat = map(1:nboot, ~tibble(.rep = ., block = sample(1:numblocks, numblocks, replace = T))) %>%
    bind_rows
  dat %>%
    right_join(boodat, by = 'block') %>%
    group_by(.rep, .add = T) %>%
    #summarize(loo = mean(est, na.rm = TRUE)) %>%
    summarize(loo = weighted.mean(est, length, na.rm = TRUE)) %>%
    select(-.rep) %>%
    mutate(length = 1)
}

# Takes a function `qpfun` which takes f2_blocks as input
# Returns a function which will repeadetly evaluate `qpfun` on
# jackknife or bootstrap resamplings of the f2_blocks, returning a nested data frame
make_resample_snps_fun = function(qpfun) {
  function(f2_blocks, boot = FALSE, ...) {
    verbose = TRUE
    if(boot) {
      if(boot == 1) boot = dim(f2_blocks)[3]
      boo = boo_list(f2_blocks, boot)
      f2dat = boo$boo
      sel = boo$sel
    } else {
      f2dat = loo_list(f2_blocks)
      sel = map(seq_len(length(f2dat)), ~setdiff(seq_len(length(f2dat)), .))
    }
    if(verbose) alert_info(paste0('Running models...\n'))
    fun = function(x) safely(qpfun)(x, verbose = FALSE, ...)
    tibble(id = seq_len(length(f2dat)), f2dat, sel) %>%
      mutate(out = furrr::future_map(f2dat, fun, .progress = verbose, .options = furrr::furrr_options(seed = TRUE)),
             result = map(out, 'result', .null = tibble()), error = map(out, 'error')) %>%
      select(-out) %>% unnest_wider(result)
  }
}

# Takes a function `qpfun` which takes f2_blocks as input
# Returns a function which will repeadetly evaluate `qpfun` on
# a subset of all samples, returning a nested data frame
make_resample_inds_fun = function(qpfun) {
  function(dir, inds, pops, verbose = TRUE, ...) {
    stopifnot(length(pops) == length(inds))
    poplist = tibble(ind = inds, pop = pops)
    lo_samples = poplist %>% group_by(pop) %>% mutate(cnt = n()) %>%
      ungroup %>% filter(cnt > 1) %>% select(-cnt)

    if(verbose) alert_info('Reading data...\n')
    f2dat = lo_samples$ind %>% rlang::set_names() %>%
      furrr::future_map(~f2_from_precomp_indivs(dir, filter(poplist, ind != .), verbose = FALSE)$f2_blocks,
                        .progress = verbose)

    if(verbose) alert_info(paste0('Running models...\n'))
    fun = function(x) safely(qpfun)(x, verbose = FALSE, ...)
    lo_samples %>% add_column(f2dat) %>%
      mutate(out = furrr::future_map(f2dat, fun, .progress = verbose),
             result = map(out, 'result', .null = tibble()), error = map(out, 'error')) %>%
      select(-out) %>% unnest_wider(result)
  }
}


#' Run models while leaving out SNP blocks
#'
#' These are wrapper functions around various `qp` functions, which will evaluate a model on multiple subdivisions of the data.
#' The subdivisions are formed by leaving out one or more SNP block at a time (see `boot` for details).
#' @name resample_snps
#' @param f2_blocks a 3d array of blocked f2 statistics
#' @param boot If `FALSE` (the default), block-jackknife resampling will be used to compute standard errors.
#' Otherwise, block-bootstrap resampling will be used to compute standard errors. If `boot` is an integer, that number
#' will specify the number of bootstrap resamplings. If `boot = TRUE`, the number of bootstrap resamplings will be
#' equal to the number of SNP blocks.
#' @param verbose print progress updates
#' @param ... named arguments which are passed to the `qp` function.
#' @return a nested data frame where each model is a row, and the columns are model parameters and model outputs
#' @examples
#' \dontrun{
#' res = qpadm_resample_snps(example_f2_blocks, target = target, left = left, right = right)
#' unnest(res, weights)
#' }
NULL


#' Run models while leaving out individuals
#'
#' These are wrapper functions around various `qp` functions, which will evaluate many models at once.
#' The models are formed by leaving out every individual which belongs to a population with more than one sample,
#' one at a time.
#' @name resample_inds
#' @param dir directory with precomputed data
#' @param inds vector of individual names
#' @param pops vector of population names. Should be the same length as `inds`
#' @param verbose print progress updates
#' @param ... named arguments passed to the `qp` function.
#' @return a nested data frame where each model is a row, and the columns are model parameters and model outputs
#' @examples
#' \dontrun{
#' res = qpadm_resample_inds
#' unnest(res, weights)
#' }
NULL


#' @rdname resample_snps
#' @export
qp3pop_resample_snps = make_resample_snps_fun(qp3pop)
#' @rdname resample_snps
#' @export
qpdstat_resample_snps = make_resample_snps_fun(qpdstat)
#' @rdname resample_snps
#' @export
qpwave_resample_snps = make_resample_snps_fun(qpwave)
#' @rdname resample_snps
#' @export
qpadm_resample_snps = make_resample_snps_fun(function(...) qpadm(..., getcov = FALSE))
#' @rdname resample_snps
#' @export
qpgraph_resample_snps = make_resample_snps_fun(qpgraph)


#' @rdname resample_inds
#' @export
qp3pop_resample_inds = make_resample_inds_fun(qp3pop)
#' @rdname resample_inds
#' @export
qpdstat_resample_inds = make_resample_inds_fun(qpdstat)
#' @rdname resample_inds
#' @export
qpwave_resample_inds = make_resample_inds_fun(qpwave)
#' @rdname resample_inds
#' @export
qpadm_resample_inds = make_resample_inds_fun(function(...) qpadm(..., getcov = FALSE))
#' @rdname resample_inds
#' @export
qpgraph_resample_inds = make_resample_inds_fun(qpgraph)


#' Summarize graph fits
#'
#' This function takes the output of \code{\link{qpgraph_resample_snps}} and creates a data frame with summaries of the estimated graph parameters. `weight` has the mean of all fits, while `low` and `high` have lower and upper quantiles.
#' @export
#' @param fits A nested data frame where each line is the output of `qpgraph()`
#' @return A data frame summarizing the fits
#' @seealso \code{\link{qpgraph_resample_snps}}
summarize_fits = function(fits, q_low = 0.05, q_high = 0.95) {
  fits$edges %>% bind_rows(.id = 'i') %>% group_by(from, to) %>%
    summarize(type = type[1], low = quantile(weight, q_low), high = quantile(weight, q_high), weight = mean(weight))
}


snpdat_to_jackest = function(dat) {
  # dat is a (grouped) data frame with one row per SNP, and a column `block` (or `CHR` and `cm`)
  # output has columns est and se or each other column in input

  if(! 'block' %in% names(dat)) {
    bl = get_block_lengths(dat)
    gr = group_vars(dat)
    dat %<>% ungroup %>% mutate(block = rep(1:length(bl), bl)) %>% group_by(gr)
  }

  dat %>% select(-any_of(c('CHR', 'cm', 'POS'))) %>%
    select(where(is.numeric), block, group_cols()) %>%
    pivot_longer(-c(block, group_cols()), names_to = '.col', values_to = 'v') %>%
    group_by(.col, block, .add = TRUE) %>% summarize(est = mean(v, na.rm=T), n = sum(!is.na(v))) %>%
    est_to_loo_dat() %>%
    group_by(.col, .add = TRUE) %>% jack_dat_stats() %>% ungroup %>% suppressMessages()
}


