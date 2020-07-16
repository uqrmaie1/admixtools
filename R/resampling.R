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
  sums/nonmiss
}


jack_vec_stats = function(loo_vec, block_lengths) {
  # input is a vector of leave-one-out estimates
  # output is list with jackknife mean and covariance
  # should give same results as 'jack_arr_stats' and 'jack_mat_stats'

  fin = is.finite(loo_vec)
  block_lengths = block_lengths[fin]
  loo_vec = loo_vec[fin]
  tot = weighted.mean(loo_vec, 1-block_lengths/sum(block_lengths), na.rm = TRUE)
  est = mean(loo_vec, na.rm = TRUE)
  y = sum(block_lengths)/block_lengths
  xtau = (tot * y - loo_vec * (y-1) - est) / sqrt(y-1)
  var = mean(xtau^2, na.rm = TRUE)

  namedList(est, var)
}


jack_mat_stats = function(loo_mat, block_lengths) {
  # input is matrix (one block per column)
  # output is list with vector of jackknife means and matrix of pairwise jackknife covariances
  # should give same results as 'jack_arr_stats'

  numblocks = length(block_lengths)
  tot = weighted_row_means(loo_mat, 1-block_lengths/sum(block_lengths))
  est = rowMeans(loo_mat, na.rm = TRUE)
  totmat = replicate(numblocks, tot)
  estmat = replicate(numblocks, est)
  y = rep(sum(block_lengths)/block_lengths, each = nrow(loo_mat))
  xtau = (totmat * y - loo_mat * (y-1) - estmat) / sqrt(y-1)
  #var = tcrossprod(xtau)/numblocks
  var = tcrossprod(replace_na(xtau, 0))/tcrossprod(!is.na(xtau))

  namedList(est, var)
}



jack_arr_stats = function(loo_arr, block_lengths) {
  # input is 3d array (n x n x m) with leave-one-out statistics
  # output is list with jackknife means and jackknife variances
  # should give same results as 'jack_mat_stats'

  numblocks = length(block_lengths)
  d1 = dim(loo_arr)[1]
  d2 = dim(loo_arr)[2]
  est = apply(loo_arr, 1:2, mean, na.rm = TRUE)
  tot = apply(loo_to_est(loo_arr), 1:2, weighted.mean, block_lengths, na.rm=T)
  estarr = replicate(numblocks, est)
  totarr = replicate(numblocks, tot)
  y = rep(sum(block_lengths)/block_lengths, each = d1*d2)
  xtau = (totarr * y - loo_arr * (y-1) - estarr) / sqrt(y-1)
  var = apply(xtau^2, 1:2, mean, na.rm = TRUE)

  namedList(est, var)
}


jack_dat_stats = function(dat) {
  # input is a grouped data frame with columns 'loo', 'length', and 'block'
  dat %>%
    mutate(est = mean(loo, na.rm = TRUE),
           y = sum(length)/length,
           tot = weighted.mean(loo, 1-1/y, na.rm = TRUE),
           xtau = (tot*y - loo*(y-1) - est)/sqrt(y-1)) %>%
    summarize(est = est[1], var = mean(xtau^2, na.rm = TRUE), n = sum(!is.na(xtau)))
}

boot_dat_stats = function(dat) {
  jack_dat_stats(dat) %>%
    mutate(var = var/n)
}


jack_pairarr_stats = function(loo_arr, block_lengths) {
  # input is 3d array (m x n x p)
  # output is list with jackknife means and jackknife covariances
  # todo: make equivalent to jack_mat_stats

  numblocks = length(block_lengths)
  est = c(t(apply(loo_arr, 1:2, mean)))
  tot = c(t(apply(loo_arr, 1:2, weighted.mean, 1-block_lengths/sum(block_lengths))))
  bj_lo_mat = loo_arr %>% aperm(c(2,1,3)) %>% arr3d_to_mat %>% t
  totmat = replicate(numblocks, tot)
  estmat = replicate(numblocks, est)
  y = rep(sum(block_lengths)/block_lengths, each = nrow(bj_lo_mat))
  xtau = (totmat * y - bj_lo_mat * (y-1) - estmat) / sqrt(y-1)
  var = tcrossprod(xtau)/numblocks

  est = t(matrix(est, dim(loo_arr)[2]))
  rownames(est) = dimnames(loo_arr)[[1]]
  colnames(est) = dimnames(loo_arr)[[2]]
  namedList(est, var)
}

make_bootfun = function(jackfun) {
  function(x, y) {
    out = jackfun(x, y)
    out[[2]] = out[[2]]/length(y)
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
#' A new block begins at the SNP after the first SNP which is not within `dist` of the start of the last block.
#' `afdat` needs to be ordered first by 'CHR', then by 'POS' or 'cm'
#' @export
#' @param afdat data frame with columns 'CHR' and either 'POS' or 'cm'
#' @param dist minimum distance between blocks
#' @param distcol which column should be used as distance column
#' @return a numeric vector where the ith element lists the number of SNPs in the ith block.
#' @examples
#' \dontrun{
#' prefix = 'path/to/packedancestrymap_prefix'
#' pops = c('pop1', 'pop2', 'pop3')
#' afdat = packedancestrymap_to_aftable(prefix, pops = pops)
#' block_lengths = get_block_lengths(afdat)
#' }
get_block_lengths = function(dat, dist = 0.05, distcol = 'cm') {

  cpp_get_block_lengths(as.integer(as.factor(dat$CHR)), dat[[distcol]], dist)

  # fpos = -1e20
  # lchrom = -1
  # xsize = 0
  # n = 0
  # bsize = c()
  # for(i in 1:nrow(dat)) {
  #   chrom = dat$CHR[i]
  #   gpos = dat[[distcol]][i]
  #   dis = gpos - fpos
  #   if ((chrom != lchrom) || (dis >= dist)) {
  #     if (xsize > 0) {
  #       bsize[n+1] = xsize
  #       n = n+1
  #     }
  #     lchrom = chrom
  #     fpos = gpos
  #     xsize = 0
  #   }
  #   xsize = xsize + 1
  # }
  # if (xsize > 0) {
  #   bsize[n+1] = xsize
  # }
  # bsize
}






#' @export
est_to_loo = function(arr) {
  # turns block estimates into leave-one-block-out estimates
  # assumes blocks are along 3rd dimension

  block_lengths = parse_number(dimnames(arr)[[3]])
  tot = apply(arr, 1:2, weighted.mean, block_lengths, na.rm=T)
  rel_bl = rep(block_lengths/sum(block_lengths), each = length(tot))
  out = (replicate(dim(arr)[3], tot) - arr*rel_bl) / (1-rel_bl)
  dimnames(out) = dimnames(arr)
  out
}



#' @export
loo_to_est = function(arr) {
  # inverse of est_to_res

  block_lengths = parse_number(dimnames(arr)[[3]])
  rel_bl = block_lengths/sum(block_lengths)
  tot = apply(arr, 1:2, weighted.mean, 1-rel_bl, na.rm=T)
  rel_bl = rep(rel_bl, each = length(tot))
  out = (replicate(dim(arr)[3], tot) - arr * (1-rel_bl))/rel_bl
  dimnames(out) = dimnames(arr)
  out
}


#' @export
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

#' @export
est_to_boo = function(arr, nboot = dim(arr)[3]) {
  # turns block estimates into bootstrap estimates
  # assumes blocks are along 3rd dimension

  nboot = max(nboot, dim(arr)[3])
  block_lengths = parse_number(dimnames(arr)[[3]])
  numblocks = length(block_lengths)
  sel = sample(seq_len(numblocks), numblocks*nboot, replace = TRUE)
  arr %>%
    matrix(numblocks, byrow=T) %>%
    `[`(sel,) %>%
    rowsum(rep(seq_len(nboot), each = numblocks), na.rm=T) %>%
    `/`(length(block_lengths)) %>%
    t %>%
    array(c(dim(arr)[1:2], nboot),
          dimnames = c(dimnames(arr)[1:2], list(rep('l1', nboot))))
}

# returns list of arrays, with each block left out at a time.
# arr is k x k x n; output is length n, output arrs are k x k x (n-1)
#' @export
loo_list = function(arr) map(1:dim(arr)[3], ~arr[,,-.])

#' @export
boo_list = function(arr, nboot = dim(arr)[3]) {
  # returns list of arrays, with each block left out at a time.
  # arr is k x k x n; output is length nboot, output arrs are k x k x n
  sel = rerun(nboot, sample(1:dim(arr)[3], replace = TRUE))
  list(boo = map(sel, ~arr[,,.]), sel = sel, test = map(sel, ~arr[,,-.]))
}


est_to_loo_dat = function(dat) {
  # like est_to_loo, but for a grouped data frame with columns 'est', 'block', and 'length'
  # adds column 'loo'
  dat %>%
    mutate(.rel_bl = length/sum(length),
           .tot = weighted.mean(est, length, na.rm=TRUE),
           loo = (.tot - est*.rel_bl) / (1-.rel_bl)) %>%
    select(-.rel_bl, -.tot)
}

loo_to_est_dat = function(dat) {
  # like loo_to_est, but for a grouped data frame with columns 'loo', 'block', and 'length'
  # adds column 'est'
  dat %>%
    mutate(.rel_bl = length/sum(length),
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
    summarize(loo = mean(est, na.rm = TRUE)) %>%
    select(-.rep) %>%
    mutate(length = 1)
}

#' Takes a function `qpfun` which takes f2_blocks as input
#' Returns a function which will repeadetly evaluate `qpfun` on
#' jackknife or bootstrap resamplings of the f2_blocks, returning a nested data frame
#' @export
make_resample_snps_fun = function(qpfun) {
  function(f2_blocks, boot = FALSE, verbose = TRUE, ...) {
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
      mutate(out = furrr::future_map(f2dat, fun, .progress = verbose),
             result = map(out, 'result', .null = tibble()), error = map(out, 'error')) %>%
      select(-out) %>% unnest_wider(result)
  }
}

#' Takes a function `qpfun` which takes f2_blocks as input
#' Returns a function which will repeadetly evaluate `qpfun` on
#' a subset of all samples, returning a nested data frame
#' @export
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
#' @param boot If `FALSE` (the default), each block will be left out at a time.
#' Otherwise bootstrap resampling is performed `n` times, where `n` is either equal to `boot` if it is an integer,
#' or equal to the number of blocks if `boot` is `TRUE`.
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
#' @param multiprocess If `TRUE` (the default), models will run on multiple cores in parallel.
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
qpadm_resample_inds = make_resample_inds_fun(function(...) qpadm(..., getcov = FALSE))
#' @rdname resample_inds
#' @export
qpgraph_resample_inds = make_resample_inds_fun(qpgraph)





