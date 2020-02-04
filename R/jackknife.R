
block_arr_mean = function(arr, block_lengths) {
  # returns grouped array means
  # arr is 3d, group over third dimension

  blockids = rep(seq_along(block_lengths), block_lengths)
  mat = arr3d_to_mat(arr)
  sums = rowsum(mat, blockids, na.rm = TRUE)
  nonmiss = rowsum((!is.na(mat))+0, blockids)
  out = mat_to_arr3d(sums/nonmiss, dim(arr)[1])
  dimnames(out)[1:2] = dimnames(arr)[1:2]
  out
}


jack_mat_stats = function(loo_mat, block_lengths) {
  # input is matrix (one block per column)
  # output is list with vector of jackknife means and matrix of pairwise jackknife covariances
  # uses mean jackknife estimate instead of overall mean; probably makes very little difference
  # lengths normalization used to have -1, removed this because it resulted in negative values; where did this come from?

  # todo: make weighted mean; test how big the difference is

  numblocks = length(block_lengths)
  est = rowMeans(loo_mat)
  mnc = t(est - loo_mat) * sqrt((sum(block_lengths)/block_lengths-1)/numblocks)
  var = crossprod(mnc)
  namedList(est, var)
}


jack_arr_stats = function(loo_arr, block_lengths) {
  # input is 3d array (n x n x m)
  # output is list with jackknife means and jackknife variances
  # uses mean jackknife estimate instead of overall mean; probably makes very little difference

  numblocks = length(block_lengths)
  est = apply(loo_arr, 1:2, weighted.mean, block_lengths)
  xtau = (rray(est) - loo_arr)^2 * rray(sum(block_lengths)/block_lengths-1, c(1, 1, numblocks))
  var = apply(xtau, 1:2, weighted.mean, block_lengths)
  namedList(est, var)
}


jack_pairarr_stats = function(loo_arr, block_lengths) {
  # input is 3d array (m x n x p)
  # output is list with jackknife means and jackknife covariances
  # uses mean jackknife estimate instead of overall mean; probably makes very little difference

  # todo: make weighted mean for var; test how big the difference is

  numblocks = length(block_lengths)
  est = c(t(apply(loo_arr, 1:2, weighted.mean, block_lengths)))
  bj_lo_mat = loo_arr %>% aperm(c(2,1,3)) %>% arr3d_to_mat %>% t
  mnc = t(est - bj_lo_mat) * sqrt((sum(block_lengths)/block_lengths-1)/numblocks)
  var = crossprod(mnc)
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

boot_mat_stats = make_bootfun(jack_mat_stats)
boot_arr_stats = make_bootfun(jack_arr_stats)
boot_pairarr_stats = make_bootfun(jack_pairarr_stats)


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
#' pops = c('Zerg', 'Protoss', 'Terran')
#' afdat = packedancestrymap_to_aftable(prefix, pops = pops)
#' block_lengths = get_block_lengths(afdat)
#' }
get_block_lengths = function(afdat, dist = 0.05, distcol = 'cm') {
  afdat %>%
    set_blocks(dist = dist, distcol = 'cm') %$% block %>%
    rle %$% lengths
}


#' @export
est_to_loo = function(arr) {
  # turns block estimates into leave-one-block-out estimates
  # assumes blocks are along 3rd dimension

  block_lengths = parse_number(dimnames(arr)[[3]])
  tot = rray(apply(arr, 1:2, weighted.mean, block_lengths, na.rm=T))
  rel_bl = rray(block_lengths/sum(block_lengths), dim = c(1,1,dim(arr)[3]))
  out = as.array((tot - arr*rel_bl) / (1-rel_bl))
  dimnames(out) = dimnames(arr)
  out
}


#' @export
loo_to_est = function(arr) {
  # inverse of est_to_res

  block_lengths = parse_number(dimnames(arr)[[3]])
  rel_bl = rray(block_lengths/sum(block_lengths), dim = c(1,1,dim(arr)[3]))
  tot = rray(apply(arr, 1:2, weighted.mean, 1-as.vector(rel_bl), na.rm=T))
  out = as.array((tot - arr * (1-rel_bl))/rel_bl)
  dimnames(out) = dimnames(arr)
  out
}

#' @export
est_to_loo_nafix = function(arr) {
  # turns block estimates into leave-one-block-out estimates
  # assumes blocks are along 3rd dimension

  block_lengths = parse_number(dimnames(arr)[[3]])
  tot = rray(apply(arr, 1:2, weighted.mean, block_lengths, na.rm=T))
  rel_bl = rray(block_lengths/sum(block_lengths), dim = c(1,1,dim(arr)[3]))
  if(any(is.na(arr))) warning(paste0('Replacing ', sum(is.na(arr)), ' NAs with 0!'))
  arr %<>% replace_na(0)
  out = as.array((tot - arr*rel_bl) / (1-rel_bl))
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
loo_list = function(arr) map(1:dim(arr)[3], ~arr[,,-.])

boo_list = function(arr, nboot = dim(arr)[3]) {
  # returns list of arrays, with each block left out at a time.
  # arr is k x k x n; output is length nboot, output arrs are k x k x n

  rerun(nboot, arr[,,sample(1:dim(arr)[3], replace = TRUE)])
}


make_resample_snps_fun = function(qpfun) {
  function(f2_blocks, boot = FALSE, multicore = TRUE, verbose = TRUE, ...) {
    if(boot) {
      f2dat = boo_list(f2_blocks, max(boot, dim(f2_blocks)[3]))
    } else {
      f2dat = loo_list(f2_blocks)
    }

    if(multicore) {
      oplan = future::plan()
      future::plan('multicore')
      on.exit(future::plan(oplan))
    }
    if(verbose) alert_info(paste0('Running models...\n'))
    fun = function(x) safely(qpfun)(x, verbose = FALSE, ...)
    tibble(id = seq_len(length(f2dat)), f2dat) %>%
      mutate(out = furrr::future_map(f2dat, fun, .progress = verbose),
             result = map(out, 'result', .null = tibble()), error = map(out, 'error')) %>%
      select(-out) %>% unnest_wider(result)
  }
}


make_resample_inds_fun = function(qpfun) {
  function(dir, inds, pops, multicore = TRUE, verbose = TRUE, ...) {
    stopifnot(length(pops) == length(inds))
    poplist = tibble(ind = inds, pop = pops)
    lo_samples = poplist %>% group_by(pop) %>% mutate(cnt = n()) %>%
      ungroup %>% filter(cnt > 1) %>% select(-cnt)

    if(multicore) {
      oplan = future::plan()
      future::plan('multicore')
      on.exit(future::plan(oplan))
    }
    if(verbose) alert_info('Reading data...\n')
    f2dat = lo_samples$ind %>% set_names %>%
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
#' These are wrapper functions around various `qp` functions, which will evaluate many models at once.
#' The models are formed by leaving out one or more SNP block at a time (see `boot` for details).
#' @name resample_snps
#' @param f2_blocks a 3d array of blocked f2 statistics
#' @param boot If `FALSE` (the default), each block will be left out at a time.
#' Otherwise bootstrap resampling is performed `n` times, where `n` is either equal to `boot` if it is an integer,
#' or equal to the number of blocks if `boot` is `TRUE`.
#' @param multicore If `TRUE` (the default), models will run on multiple cores in parallel.
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
#' @param multicore If `TRUE` (the default), models will run on multiple cores in parallel.
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






