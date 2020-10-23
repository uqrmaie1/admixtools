
# a, b: num left - 1, num right - 1 (dimensions of f4_est)
qpwave_dof = function(a, b, r) a*b - qpadm_dof(a, b, r)
qpwave_dof = function(a, b, r) r*(a+b-r)
qpwave_dof = function(a, b, r) (a+b)*r - r^2
qpwave_dof = function(a, b, r) a*b - (a-r)*(b-r)

qpadm_dof = function(a, b, r) a*b - qpwave_dof(a, b, r)
qpadm_dof = function(a, b, r) a*b - (a+b)*r + r^2
qpadm_dof = function(a, b, r) (a-r)*(b-r)



opt_A = function(B, xmat, qinv, fudge = 0.0001) {
  # return A which minimizes covariance-weighted t(c(E)) %*% qinv %*% c(E), where E = xmat - (A %*% B)
  # A: nr * rnk
  # B: rnk * nc
  # xmat: nr * nc
  # tdim: rnk * max(nr, nc)
  # coeffs: tdim * tdim
  # rhs: rnk * nc
  # qinv: nr*nc * nr*nc
  B2 = diag(nrow(xmat)) %x% B
  coeffs = B2 %*% qinv %*% t(B2)
  rhs = c(t(xmat)) %*% qinv %*% t(B2)
  diag(coeffs) = diag(coeffs) + fudge*sum(diag(coeffs))
  A2 = solve(coeffs, rhs[1,])
  matrix(A2, nrow(xmat), byrow = TRUE)
}


opt_B = function(A, xmat, qinv, fudge = 0.0001) {
  # return B which minimizes covariance-weighted t(c(E)) %*% qinv %*% c(E), where E = xmat - (A %*% B)
  # A: nr * rnk
  # B: rnk * nc
  # xmat: nr * nc
  # tdim: rnk * max(nr, nc)
  # coeffs: tdim * tdim
  # rhs: rnk * nc
  # qinv: nr*nc * nr*nc
  A2 = A %x% diag(ncol(xmat))
  coeffs = t(A2) %*% qinv %*% A2
  rhs = c(t(xmat)) %*% qinv %*% A2
  diag(coeffs) = diag(coeffs) + fudge*sum(diag(coeffs))
  B2 = solve(coeffs, rhs[1,])
  matrix(B2, ncol = ncol(xmat), byrow = TRUE)
}


get_weights_covariance = function(f4_lo, qinv, block_lengths, fudge = 0.0001, boot = FALSE,
                                  constrained = FALSE, qpsolve = NULL) {
  rnk = dim(f4_lo)[1]-1
  numreps = dim(f4_lo)[3]
  wmat = matrix(NA, numreps, dim(f4_lo)[1])
  for(i in 1:numreps) {
    wmat[i,] = qpadm_weights(as.matrix(f4_lo[,,i]), qinv, rnk, fudge = fudge,
                             constrained = constrained, qpsolve = qpsolve)$weights
  }
  if(!boot) wmat = wmat * (numreps-1) / sqrt(numreps)
  cov(wmat)
}


qpadm_weights = function(xmat, qinv, rnk, fudge = 0.0001, iterations = 20,
                         constrained = FALSE, qpsolve = NULL) {
  if(rnk == 0) return(list(weights = 1))
  nr = nrow(xmat)
  f4svd = svd(xmat)
  B = t(f4svd$v[, seq_len(rnk), drop=FALSE])
  A = xmat %*% t(B)
  for(i in seq_len(iterations)) {
    A = opt_A(B, xmat, qinv, fudge = fudge)
    B = opt_B(A, xmat, qinv, fudge = fudge)
  }
  x = t(cbind(A, 1))
  y = c(rep(0, rnk), 1)
  rhs = crossprod(x)
  lhs = crossprod(x, y)
  if(constrained) {
    Amat = cbind(1, -1, diag(nr), -diag(nr))
    bvec = c(1, -1, rep(0, nr), -rep(1, nr))
    w = qpsolve(rhs, lhs, Amat, bvec)
  } else w = as.matrix(solve(rhs, lhs))[,1]
  weights = w/sum(w)
  namedList(weights, A, B)
}


#' Estimate admixture weights
#'
#' `qpadm` models a target population as a mixture of left (source) populations, given a set of right (outgroup) populations.
#' It can be used to estimate whether the left populations explain all genetic variation in the target population, relative to the right populations, and to estimate admixture proportions of the left populations to the target population.
#' @export
#' @param data The input data in the form of:
#' \itemize{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}}
#' \item A directory with f2 statistics
#' \item The prefix of a genotype file
#' }
#' @param left Left populations (sources)
#' @param right Right populations (outgroups)
#' @param target Target population
#' @param f4blocks Instead of f2 blocks, f4 blocks can be supplied. This is used by \code{\link{qpadm_multi}}
#' @param fudge Value added to diagonal matrix elements before inverting
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix
#' of f4 statistics, as well as the weight standard errors, will be computed using block-jackknife.
#' Otherwise bootstrap resampling is performed `n` times, where `n` is either equal to `boot` if it is an integer,
#' or equal to the number of blocks if `boot` is `TRUE`. The the covariance matrix of f4 statistics,
#' as well as the weight standard errors, will be computed using bootstrapping.
#' @param getcov Compute weights covariance. Setting `getcov = FALSE` will speed up the computation.
#' @param constrained Constrain admixture weights to be non-negative
#' @param cpp Use C++ functions. Setting this to `FALSE` will be slower but can help with debugging.
#' @param verbose Print progress updates
#' @param ... If `data` is the prefix of genotype files, additional arguments will be passed to \code{\link{f4blockdat_from_geno}}
#' @return `qpadm` returns a list with up to four data frames describing the model fit:
#' \enumerate{
#' \item `weights`: A data frame with estimated admixture proportions where each row is a left population.
#' \item `f4`: A data frame with estimated and fitted f4-statistics
#' \item `rankdrop`: A data frame describing model fits with different ranks, including *p*-values for the overall fit
#' and for nested models (comparing two models with rank difference of one). A model with `L` left populations and `R` right populations has an f4-matrix of dimensions `(L-1)*(R-1)`. If no two left population form a clade with respect to all right populations, this model will have rank `(L-1)*(R-1)`.
#'   \itemize{
#'     \item `f4rank`: Tested rank
#'     \item `dof`: Degrees of freedom of the chi-squared null distribution: `(L-1-f4rank)*(R-1-f4rank)`
#'     \item `chisq`: Chi-sqaured statistic, obtained as `E'QE`, where `E` is the difference between estimated and fitted f4-statistics, and `Q` is the f4-statistic covariance matrix.
#'     \item `p`: p-value obtained from `chisq` as `pchisq(chisq, df = dof, lower.tail = FALSE)`
#'     \item `dofdiff`: Difference in degrees of freedom between this model and the model with one less rank
#'     \item `chisqdiff`: Difference in chi-squared statistics
#'     \item `p_nested`: *p*-value testing whether the difference between two models of rank difference 1 is significant
#'   }
#' \item `popdrop`: A data frame describing model fits with different populations. Note that all models with fewer populations use the same set of SNPs as the first model.
#'   \itemize{
#'     \item `pat`: A binary code indicating which populations are present in this model. A `1` represents dropped populations. The full model is all zeros.
#'     \item `wt`: Number of populations dropped
#'     \item `dof`: Degrees of freedom of the chi-squared null distribution: `(L-1-f4rank)*(R-1-f4rank)`
#'     \item `chisq`: Chi-sqaured statistic, obtained as `E'QE`, where `E` is the difference between estimated and fitted f4-statistics, and `Q` is the f4-statistic covariance matrix.
#'     \item `p`: *p*-value obtained from `chisq` as `pchisq(chisq, df = dof, lower.tail = FALSE)`
#'     \item `f4rank`: Tested rank
#'     \item `feasible`: A model is feasible if all weights fall between 0 and 1
#'     \item `<population name>`: The weights for each population in this model
#'   }
#' }
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 10)
#' @seealso \code{\link{qpwave}}, \code{\link{lazadm}}
#' @examples
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' target = 'Denisova.DG'
#' qpadm(example_f2_blocks, left, right, target)
#' \dontrun{
#' qpadm("/my/geno/prefix", left, right, target)
#' }
qpadm = function(data, left, right, target, f4blocks = NULL,
                 fudge = 0.0001, boot = FALSE, getcov = TRUE,
                 constrained = FALSE, cpp = TRUE, verbose = TRUE, ...) {

  #----------------- prepare f4 stats -----------------
  f2_blocks = NULL
  if(is.null(f4blocks)) {
    if(all(file.exists(left, right))) {
      left %<>% readLines
      right %<>% readLines
    }
    if(!is.null(target)) left = c(target, setdiff(left, target))
    if(is_geno_prefix(data)) {
      f4blockdat = f4blockdat_from_geno(data, left = left, right = right, verbose = verbose, ...)
      f4blocks = f4blockdat %>% f4blockdat_to_f4blocks()
    } else {
      if(verbose) alert_info('Computing f4 stats...\n')
      f2_blocks = get_f2(data, pops = c(left, right), afprod = TRUE, verbose = verbose)
      f4blocks = f2blocks_to_f4blocks(f2_blocks, left, right)
    }
    f4blocks = f4blocks[left[-1], right[-1],,drop=FALSE]
  } else {
    if(!all(map_lgl(list(data, left, right), is.null)))
      stop("Can't provide 'f4blocks' and data/left/right at the same time!")
    naml = dimnames(f4blocks)[1]
    namr = dimnames(f4blocks)[2]
    left = c(names(naml), naml[[1]])
    right = c(names(namr), namr[[1]])
    if(!is.null(target)) target = left[1]
  }
  pops = c(left, right)
  if(any(duplicated(pops)))
    stop(paste0('Duplicated pops: ', paste0(unique(pops[duplicated(pops)]), collapse=', ')))

  f4stats = f4blocks %>% f4blocks_to_f4stats(boot = boot)

  f4_est = f4stats$est
  f4_var = f4stats$var
  f4_lo = f4stats$f4_lo
  block_lengths = f4stats$block_lengths
  diag(f4_var) = diag(f4_var) + fudge*sum(diag(f4_var))
  qinv = solve(f4_var)
  out = list()

  #----------------- compute admixture weights -----------------
  if(!is.null(target)) {
    if(verbose) alert_info('Computing admixture weights...\n')
    if(cpp) {
      qpadm_weights = cpp_qpadm_weights
      get_weights_covariance = cpp_get_weights_covariance
    }
    weight = qpadm_weights(f4_est, qinv, rnk = length(left)-2, fudge = fudge, constrained = constrained,
                           qpsolve = qpsolve)$weights %>% c
    if(getcov) {
      if(verbose) alert_info('Computing standard errors...\n')
      se = sqrt(diag(get_weights_covariance(f4_lo, qinv, block_lengths, fudge = fudge, boot = boot,
                                            constrained = constrained, qpsolve = qpsolve)))
    } else se = NA
    out$weights = tibble(target, left = left[-1], weight, se) %>% mutate(z = weight/se)

    wvec = out$weights %>% select(left, weight) %>% deframe
    if(!is.null(f2_blocks)) out$f4 = fitted_f4(f2_blocks, wvec, target, left[-1], right)
    #else out$f4 = f4blockdat %>% f4blockdat_to_f4out(boot = boot)
    else out$f4 = NULL
  } else {
    if(!is.null(f2_blocks)) out$f4 = f4(f2_blocks, left[1], left[-1], right[1], right[-1], verbose = FALSE)
    #else out$f4 = f4blockdat %>% f4blockdat_to_f4out(boot = boot)
    else out$f4 = NULL
  }

  #----------------- compute number of admixture waves -----------------
  if(verbose) alert_info('Computing number of admixture waves...\n')
  out$rankdrop = drop_ranks(f4_est, qinv, fudge, constrained, cpp)
  if(!is.null(target)) out$popdrop = drop_pops(f4_est, qinv, fudge, constrained, cpp, left[-1])

  out
}

#' Estimate admixture waves
#'
#' `qpwave` compares two sets of populations (`left` and `right`) to each other. It estimates a lower bound on the number of admixtue waves that went from `left` into `right`, by comparing a matrix of f4-statistics to low-rank approximations. For a rank of 0 this is equivalent to testing whether `left` and `right` form clades relative to each other.
#' @export
#' @inheritParams qpadm
#' @return `qpwave` returns a list with up to two data frames describing the model fit:
#' \enumerate{
#' \item `f4` A data frame with estimated f4-statistics
#' \item `rankdrop`: A data frame describing model fits with different ranks, including *p*-values for the overall fit
#' and for nested models (comparing two models with rank difference of one). A model with `L` left populations and `R` right populations has an f4-matrix of dimensions `(L-1)*(R-1)`. If no two left population form a clade with respect to all right populations, this model will have rank `(L-1)*(R-1)`.
#'   \itemize{
#'     \item `f4rank`: Tested rank
#'     \item `dof`: Degrees of freedom of the chi-squared null distribution: `(L-1-f4rank)*(R-1-f4rank)`
#'     \item `chisq`: Chi-sqaured statistic, obtained as `E'QE`, where `E` is the difference between estimated and fitted f4-statistics, and `Q` is the f4-statistic covariance matrix.
#'     \item `p`: p-value obtained from `chisq` as `pchisq(chisq, df = dof, lower.tail = FALSE)`
#'     \item `dofdiff`: Difference in degrees of freedom between this model and the model with one less rank
#'     \item `chisqdiff`: Difference in chi-squared statistics
#'     \item `p_nested`: *p*-value testing whether the difference between two models of rank difference 1 is significant
#'   }
#' }
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 10)
#' @seealso \code{\link{qpadm}}
#' @examples
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' qpwave(example_f2_blocks, left, right)
qpwave = function(data, left, right,
                  fudge = 0.0001, boot = FALSE,
                  constrained = FALSE, cpp = TRUE, verbose = TRUE)
  qpadm(data = data, left = left, right = right, target = NULL,
        fudge = fudge, boot = boot,
        constrained = constrained, cpp = cpp, verbose = verbose)



f2blocks_to_f4stats = function(f2_blocks, left, right, boot = FALSE) {

  f2blocks_to_f4blocks(f2_blocks, left, right) %>%
    f4blocks_to_f4stats(boot = boot)
}

f2blocks_to_f4blocks = function(f2_blocks, left, right) {

  if(length(right) < 2) stop('Not enough right populations!')
  if(length(left) < 2) stop('Not enough left populations!')

  nr = length(right) - 1
  nl = length(left) - 1
  f4_blocks = (f2_blocks[left[-1], rep(right[1], nr), , drop = FALSE] +
                 f2_blocks[rep(left[1], nl), right[-1], , drop = FALSE] -
                 f2_blocks[rep(left[1], nl), rep(right[1], nr), , drop = FALSE] -
                 f2_blocks[left[-1], right[-1], , drop = FALSE])/2
  dimnames(f4_blocks)[[2]] = right[-1]
  names(dimnames(f4_blocks))[1:2] = c(left[1], right[1])
  f4_blocks
}


f4blocks_to_f4stats = function(f4_blocks, boot = FALSE) {

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo_nafix)
  statfun = ifelse(boot, boot_pairarr_stats, jack_pairarr_stats)

  f4_lo = f4_blocks %>% samplefun
  block_lengths = parse_number(dimnames(f4_lo)[[3]])

  out = f4_lo %>% statfun(block_lengths)
  #out = f4_lo %>% jack_pairarr_stats(block_lengths)
  out$f4_lo = f4_lo
  out$block_lengths = block_lengths
  out
}

#' Turn f4 block data to 3d array
#'
#' @param f4blockdat f4 block data frame generated by \code{\link{f4blockdat_from_geno}}
#' @param remove_na Remove blocks with missing values
#' @export
f4blockdat_to_f4blocks = function(f4blockdat, remove_na = TRUE) {
  # assumes pop1 = left[1], pop2 = left[-1], pop3 = right[1], pop4 = right[-1]
  stopifnot(length(unique(f4blockdat$pop1)) == 1 && length(unique(f4blockdat$pop3)) == 1)
  stopifnot(all(c('block', 'n', paste0('pop', 1:4)) %in% names(f4blockdat)))

  f4blockdat %<>% arrange(block, pop4, pop2)
  p1 = f4blockdat$pop1[1]
  p3 = f4blockdat$pop3[1]
  p2 = unique(f4blockdat$pop2)
  p4 = unique(f4blockdat$pop4)
  bl = f4blockdat %>% select(block, n) %>%
    group_by(block) %>% slice_max(n, with_ties = F) %$% n %>% paste0('l', .)

  stopifnot(nrow(f4blockdat) == length(p2)*length(p4)*length(bl))

  arr = array(f4blockdat$est, c(length(p2), length(p4), length(bl)), list(p2, p4, bl) %>% purrr::set_names(c(p1, p3, 'bl')))
  if(remove_na) {
    keep = apply(arr, 3, function(x) sum(is.na(x)) == 0)
    arr = arr[,,keep, drop = FALSE]
    if(sum(!keep) > 0) warning(paste0('Discarding ', sum(!keep), ' block(s) due to missing values!\n',
                                        'Discarded block(s): ', paste0(which(!keep), collapse = ', ')))
  }
  arr
}

f4blockdat_to_f4out = function(f4blockdat, boot) {

  samplefun = ifelse(boot, function(...) est_to_boo_dat(...), est_to_loo_dat)
  datstatfun = ifelse(boot, boot_dat_stats, jack_dat_stats)
  totn = f4blockdat %>%
    group_by(pop1, pop2, pop3, pop4) %>%
    summarize(n = sum(n, na.rm=T))

  if(!'den' %in% names(f4blockdat)) {
    f4out = f4blockdat %>%
      group_by(pop1, pop2, pop3, pop4) %>%
      samplefun()
  } else {
    f4out = f4blockdat %>%
      rename(num = est) %>%
      pivot_longer(c(num, den), names_to = 'type', values_to = 'est') %>%
      group_by(type, pop1, pop2, pop3, pop4) %>%
      samplefun() %>%
      pivot_wider(names_from = type, values_from = c(est, loo)) %>%
      mutate(loo = loo_num/loo_den)
  }
  f4out %>%
    datstatfun %>%
    ungroup %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    transmute(pop1, pop2, pop3, pop4, est, se, z, p) %>%
    left_join(totn, by = c('pop1', 'pop2', 'pop3', 'pop4'))
}



lazadm_old = function(f2_data, left, right, target = NULL,
                      boot = FALSE, constrained = TRUE) {

  #----------------- prepare f4 stats -----------------
  if(is.null(target)) {
    left %<>% readLines
    right %<>% readLines
    target = left[1]
    left = left[-1]
  }
  pops = c(target, left, right)
  stopifnot(!any(duplicated(pops)))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo_nafix)
  f2_blocks = get_f2(f2_data, pops, afprod = TRUE) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  f2_mat = apply(f2_blocks, 1:2, weighted.mean, block_lengths)

  r = 1:length(right)
  og_indices = expand.grid(r, r, r) %>%
    filter(Var1 != Var2, Var1 != Var3, Var2 < Var3)

  pos1 = target
  pos2 = right[og_indices[,1]]
  pos3 = right[og_indices[,2]]
  pos4 = right[og_indices[,3]]

  y = f2_mat[cbind(pos1, pos4)] +
      f2_mat[cbind(pos2, pos3)] -
      f2_mat[cbind(pos1, pos3)] -
      f2_mat[cbind(pos2, pos4)]

  x = f2_mat[pos4, left] +
      f2_mat[cbind(pos2, pos3)] -
      f2_mat[pos3, left] -
      f2_mat[cbind(pos2, pos4)]

  #----------------- compute admixture weights -----------------
  lhs = crossprod(x, y)
  rhs = crossprod(x)
  nc = length(left)

  if(constrained) weight = qpsolve(rhs, lhs, diag(nc), rep(0, nc))
  else weight = solve(rhs, lhs)[,1]
  weight = weight/sum(weight)

  tibble(target, left, weight)
}


#' Estimate admixture weights
#'
#' Models target as a mixture of left populations, and outgroup right populations. Uses Lazaridis method
#' based non-negative least squares of f4 matrix.
#' @export
#' @inheritParams qpadm
#' @return `lazadm` returns a data frame with weights and standard errors for each left population
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 9)
#' @seealso \code{\link{qpadm}}
#' @examples
#' target = 'Denisova.DG'
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' lazadm(example_f2_blocks, left, right, target)
#' lazadm(example_f2_blocks, left, right, target, constrained = FALSE)
lazadm = function(data, left, right, target,
                  boot = FALSE, constrained = TRUE) {

  #----------------- prepare f4 stats -----------------
  if(is.null(target)) {
    left %<>% readLines
    right %<>% readLines
    target = left[1]
    left = left[-1]
  }
  pops = c(target, left, right)
  stopifnot(!any(duplicated(pops)))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo_nafix)
  statfun = ifelse(boot, boot_mat_stats, jack_mat_stats)
  f2_blocks = get_f2(data, pops, afprod = TRUE)
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])
  loo_blocks = f2_blocks %>% samplefun
  loo_blocks %<>% abind::abind(apply(f2_blocks, 1:2, weighted.mean, block_lengths))

  numblocks = length(block_lengths) + 1

  r = match(right, pops)
  og_indices = expand.grid(r, r, r) %>%
    filter(Var1 != Var2, Var1 != Var3, Var2 < Var3)
  ncomb = nrow(og_indices)
  nleft = length(left)
  blocknums = rep(1:numblocks, each = ncomb)

  pos1 = match(target, pops)
  pos2 = og_indices[,1]
  pos3 = og_indices[,2]
  pos4 = og_indices[,3]

  ymat = matrix(loo_blocks[cbind(pos1, pos4, blocknums)] +
                loo_blocks[cbind(pos2, pos3, blocknums)] -
                loo_blocks[cbind(pos1, pos3, blocknums)] -
                loo_blocks[cbind(pos2, pos4, blocknums)], ncomb)

  xarr = loo_blocks[pos4, left,] +
    array(loo_blocks[cbind(pos2, pos3, rep(blocknums, each = nleft))], c(ncomb, nleft, numblocks)) -
    loo_blocks[pos3, left,] -
    array(loo_blocks[cbind(pos2, pos4, rep(blocknums, each = nleft))], c(ncomb, nleft, numblocks))

  #----------------- compute admixture weights -----------------
  lhs = matrix(NA, nleft, numblocks)
  rhs = array(NA, c(nleft, nleft, numblocks))
  for(i in 1:numblocks) {
    lhs[,i] = crossprod(xarr[,,i], ymat[,i])
    rhs[,,i] = crossprod(xarr[,,i])
  }

  if(constrained) {
    Amat = cbind(1, -1, diag(nleft), -diag(nleft))
    bvec = c(1, -1, rep(0, nleft), -rep(1, nleft))
    fun = function(rhs, lhs) qpsolve(rhs, lhs, Amat, bvec)
  } else fun = function(rhs, lhs) solve(rhs, lhs)[,1]

  wmat = sapply(1:numblocks, function(i) {
    w = fun(rhs[,,i], lhs[,i,drop = FALSE])
    if(sum(w) > 1e-10) w = w/sum(w)
    w
  })
  stats = wmat[,-numblocks] %>% statfun(block_lengths, tot = wmat[,numblocks])
  weight = stats$est
  se = sqrt(diag(stats$var))

  tibble(target, left, weight, se) %>% mutate(z = weight/se)
}



drop_ranks = function(f4_est, qinv, fudge, constrained, cpp) {
  # drops rank and fits qpadm model

  rnk = nrow(f4_est) - 1
  fitrank = function(x) qpadm_fit(f4_est, qinv, x, fudge = fudge,
                                  constrained = constrained, cpp = cpp, addweights = FALSE)
  #rankdrop = map_dfr(rev(seq_len(max(1, rnk))-(rnk==0)), fitrank) %>%
  rankdrop = map_dfr(rnk:0, fitrank) %>%
    mutate(dofdiff = lead(dof)-dof, chisqdiff = lead(chisq)-chisq,
           p_nested = pchisq(chisqdiff, dofdiff, lower.tail = FALSE))
  rankdrop
}


drop_pops = function(f4_est, qinv, fudge, constrained, cpp, left) {
  # drops each subset of left populations and fits qpadm model
  #if(nrow(f4_est) == 1) return()
  fitpop = function(x, y) qpadm_fit(x, y, nrow(x)-1, fudge = fudge, constrained = constrained,
                                    cpp = cpp, addweights = TRUE)
  nc = ncol(f4_est)
  popdrop = left %>% power_set %>% enframe(name = 'i', value = 'pop') %>%
    #filter(map_int(pop, length) > 1) %>%
    mutate(f4mat = map(pop, ~f4_est[., , drop=F]),
           ind = map(pop, ~((rep(match(., left), each=nc)-1)*nc+(1:nc))),
           qinvs = map(ind, ~qinv[., ., drop=F]),
           fd = map2(f4mat, qinvs, fitpop)) %$% fd %>% bind_rows

  wmat = popdrop %>% select(-1:-4) %>% as.matrix
  pat = apply(is.na(wmat)+0, 1, paste0, collapse='')
  feasible = apply(wmat, 1, function(x) all(between(na.omit(x), 0, 1)))
  popdrop %<>% add_column(pat, feasible)

  child_patterns = function(pat) {
    if(length(pat) == 0) return(list())
    ones = str_locate_all(pat, '0')[[1]][,1]
    map(ones, ~{pat2 = pat; str_sub(pat2, ., .) = 1; pat2})
  }

  rnk = nrow(f4_est) - 1
  nested = popdrop %>% filter(f4rank == rnk-1) %>% select(pat, dof, chisq, feasible)
  for(i in rev(seq_len(max(0,rnk-1)))) {
    children = child_patterns(nested$pat)
    thisrnk = popdrop %>% filter(f4rank == i, pat %in% children, feasible) %>%
      top_n(1, -jitter(chisq)) %>% select(pat, dof, chisq, feasible)
    nested %<>% bind_rows(thisrnk)
  }
  nested %<>% arrange(-dof) %>%
    mutate(dofdiff = dof-lead(dof), chisqdiff = chisq-lead(chisq),
           p_nested = pchisq(chisqdiff, dofdiff, lower.tail = FALSE)) %>%
    transmute(best = TRUE, pat, dofdiff, chisqdiff, p_nested)

  popdrop %<>% left_join(nested, by = 'pat') %>% mutate(wt = rnk-f4rank) %>%
    select(pat, wt, dof, chisq, p, everything(),
           feasible, best, dofdiff, chisqdiff, p_nested) %>%
    arrange(dof, pat)
  popdrop
}



qpadm_evaluate_fit = function(xmat, qinv, A, B, f4rank) {

  res = t(xmat - A %*% B)
  chisq = (t(c(res)) %*% qinv %*% c(res))[,1]
  dof = qpadm_dof(nrow(A), ncol(B), f4rank)
  p = if(dof == 0) NA else pchisq(chisq, df = dof, lower.tail = FALSE)
  tibble(f4rank, dof, chisq, p)
}


qpadm_fit = function(xmat, qinv, rnk, fudge = 0.0001, iterations = 20,
                     constrained = FALSE, cpp = TRUE, addweights = FALSE) {
  # returns one-row data frame with fit of one qpadm model

  if(cpp) qpadm_weights = cpp_qpadm_weights
  if(rnk == 0) {
    fit = list(A = matrix(0, nrow(xmat)), B = t(matrix(0, ncol(xmat))), weights = 1)
  } else {
    fit = qpadm_weights(xmat, qinv, rnk, fudge = fudge, constrained = constrained, qpsolve = qpsolve)
  }
  out = qpadm_evaluate_fit(xmat, qinv, fit$A, fit$B, rnk)
  if(addweights && length(fit$weights) == nrow(xmat)) {
    #out %<>% bind_cols(fit$weights %>% t %>% as_data_frame %>% set_names(rownames(xmat)))
    out %<>% bind_cols(fit$weights %>% t %>% as_tibble(.name_repair = ~rownames(xmat)))
  }
  out
}

add_weighted_f2 = function(f2_blocks, weights) {

  nam = dimnames(f2_blocks)[[1]]
  npop = dim(f2_blocks)[1]
  nblocks = dim(f2_blocks)[3]
  stopifnot(all(names(weights) %in% nam))
  matchedweights = rep(0, npop)
  matchedweights[match(names(weights), nam)] = weights
  add = apply(f2_blocks * matchedweights, 2:3, sum, na.rm=T)
  f2_blocks = abind::abind(f2_blocks, add, along = 1)
  f2_blocks = abind::abind(f2_blocks, rbind(add, 0), along = 2)
  dimnames(f2_blocks)[[1]] = dimnames(f2_blocks)[[2]] = c(nam, 'fit')
  f2_blocks
}


fitted_f4 = function(f2_blocks, weights, target, left, right) {

  weights = weights/sum(weights)
  pops = names(weights)
  f2_blocks_plus = add_weighted_f2(f2_blocks, weights)
  fitf4 = f4(f2_blocks_plus, target, c(left, 'fit'), right, right, verbose = FALSE) %>% filter(pop3 != pop4)
  fitf4 %>% left_join(enframe(weights, name = 'pop2', value = 'weight'), by = 'pop2') %>%
    arrange(pop1, pop3, pop4, pop2)
}

#' Faster version of \code{\link{qpadm}} with reduced output
#'
#' Models target as a mixture of left populations, given a set of outgroup right populations.
#' Can be used to estimate admixture proportions, and to estimate the number of independent
#' admixture events.
#' @export
#' @inheritParams qpadm
#' @param rnk Rank of f4-matrix. Defaults to one less than full rank.
#' @return Data frame with `f4rank`, `dof`, `chisq`, `p`, `feasible`
#' @seealso \code{\link{qpadm}}
#' @examples
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' target = 'Denisova.DG'
#' qpadm_p(example_f2_blocks, left, right, target)
qpadm_p = function(f2_data, left, right, target = NULL, fudge = 0.0001, boot = FALSE,
                   constrained = FALSE, rnk = length(setdiff(left, target)) - 1, cpp = TRUE) {

  #force(rnk)
  if(!is.null(target)) left = c(target, setdiff(left, target))
  f2_blocks = get_f2(f2_data, pops = left, pops2 = right, afprod = TRUE)
  f4dat = f2blocks_to_f4stats(f2_blocks, left, right, boot = boot)
  f4_est = f4dat$est
  f4_var = f4dat$var
  diag(f4_var) = diag(f4_var) + fudge*sum(diag(f4_var))
  qinv = solve(f4_var)
  out = qpadm_fit(f4_est, qinv, rnk, fudge = fudge,
                  constrained = constrained, cpp = cpp, addweights = TRUE)
  w = out %>% select(-1:-4) %>% as.matrix
  out %>% select(1:4) %>% mutate(feasible = all(between(w, 0, 1)))
}


#' Test if two sets of populations form two clades
#'
#' A thin wrapper around \code{\link{qpadm_p}} with `rnk` set to zero
#' @export
#' @inheritParams qpadm
test_cladality = function(f2_data, left, right, fudge = 0.0001, boot = FALSE, cpp = TRUE) {
  qpadm_p(f2_data, left, right, fudge = fudge, boot = boot, rnk = 0, cpp = cpp)
}


find_right = function(f2_blocks, target, pops) {
  f4p = f4(f2_blocks, verbose = FALSE) %>%
    select(pop1:pop4, p)
}


all_comb = function(pops, target = NULL, left = NULL, right = NULL) {
  stopifnot(length(union(pops, target)) > 5)
  stopifnot(is.null(left) || is.null(right) || is.null(target))
  if(is.null(target)) target = setdiff(pops, union(left, right))
  if(!is.null(left)) return(bind_rows(map_dfr(target, ~all_comb_right(., left, right = setdiff(pops, c(., left))))))
  if(!is.null(right)) return(bind_rows(map_dfr(target, ~all_comb_left(., left = setdiff(pops, c(., right)), right))))
  lr = map(target, ~setdiff(pops, .))
  lrcomb = map(lr, all_lr)
  len = map_dbl(lrcomb, ~length(.$left))
  t = rep(target, len)
  l = flatten(map(lrcomb, 1))
  r = flatten(map(lrcomb, 2))
  tibble(target = t, left = l, right = r)
}

all_comb_right = function(target, left, right) {
  nleft = length(left)
  r = power_set(right) %>% discard(~length(.) <= nleft)
  tibble(target, left = list(left), right = r)
}
all_comb_left = function(target, left, right) {
  nright = length(right)
  l = power_set(left) %>% keep(~between(length(.), 2, nright-1))
  tibble(target, left = l, right = list(right))
}

all_lr = function(pops) {
  # splits pops into all possible combinations of left and right, and unused
  npops = length(pops)
  left = power_set(pops) %>% discard(~length(.) >= npops/2 | length(.) < 2)
  right = map(left, ~power_set(setdiff(pops, .)))
  right2 = map2(left, right, ~discard(.y, function(x) length(x) <= length(.x)))
  l = rep(left, map_dbl(right2, length))
  r = flatten(right2)
  list(left = l, right = r)
}

all_lr2 = function(pops, rightfix = 0) {
  # splits pops into all possible combinations of left and right
  maxleft = min(length(pops), floor((length(pops) + rightfix - 1)/2))
  left = map(1:maxleft, ~combn(pops, ., simplify = F)) %>% flatten
  right = map(left, ~setdiff(pops, .))
  namedList(left, right)
}


#' Compute all pairwise qpwave p-values
#'
#' For all pairs of left populations qpwave rank 0 Chi-squared statistics and p-values will be computed. This tests for each pair of left populations whether it forms a clade with respect to the right populations.
#' @export
#' @param f2_data A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}. Alternatively, a directory with f2 statistics.
#' @param left Left populations
#' @param right Right populations
#' @examples
#' \dontrun{
#' left = pops[5:7]
#' right = pops[1:4]
#' f2_blocks = f2_from_precomp('/my/f2/dir/', pops = left, pops2 = right, afprod = TRUE)
#' qpwave_pairs(f2_blocks, left, right)
#' # If f2-stats are too big to load them into memory,
#' # the following will read the data for each pair from disk:
#' qpwave_pairs('/my/f2/dir/', left, right)
#' }
qpwave_pairs = function(f2_data, left, right) {
  expand_grid(pop1 = left, pop2 = left) %>%
    filter(pop1 < pop2) %>%
    mutate(out = furrr::future_map2(pop1, pop2, ~qpadm_p(f2_data, .y, right, .x, rnk = 0))) %>%
    unnest_wider(out) %>%
    select(pop1, pop2, chisq, p) %>%
    bind_rows(rename(., pop1 = pop2, pop2 = pop1)) %>%
    arrange(pop1, pop2)
}



#' Compute all pairwise qpadm p-values
#'
#' For all pairs of left populations qpadm Chi-squared statistics and p-values will be computed
#' @export
#' @param f2_blocks 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' @param left Left populations
#' @param right Right populations
#' @param target Target population
#' @param verbose Print progress updates
#' @return A data frame with Chi-squared statistics and p-values for each population combination
qpadm_rotate = function(f2_blocks, left, right, target, verbose = TRUE) {

  lr = all_lr2(left, length(right))
  if(verbose) alert_info(paste0('Evaluating ', length(lr[[1]]), ' models...\n'))
  qpadm_eval_rotate(f2_blocks, target, lr, right, verbose = verbose)

}

qpadm_eval_rotate = function(f2_blocks, target, leftright_dat, rightfix, verbose = TRUE) {
  leftright_dat %>%
    as_tibble %>%
    mutate(res = map2(left, right, ~qpadm_p(f2_blocks, .x, c(.y, rightfix), target),
                                    .progress = verbose)) %>%
    unnest_wider(res) #%>%
  #mutate(chisq = map(rankdrop, 'chisq') %>% map_dbl(1)) %>%
  #arrange(chisq)
}



qpadmmodels_to_popcombs = function(models) {
  # takes data frame of left, right, target populations
  # returns all f4 popcombs numbered by model

  if('target' %in% names(models)) models$left = map2(models$left, models$target, ~union(.y, .x))
  models %>% as_tibble %>% mutate(model = 1:n()) %>% rowwise %>%
    transmute(model, pop1 = left[1], pop2 = list(left[-1]), pop3 = right[1], pop4 = list(right[-1])) %>%
    unnest_longer(pop2) %>% unnest_longer(pop4)
}


#' Run multiple qpadm models
#'
#' This function runs multiple qpadm models, re-using f4-statistics where possible. Supports parallel evaluation of models, which can be turned on with `future::plan('multisession')` or similar, and turned off with `future::plan('sequential')`.
#' @export
#' @param data The input data in the form of:
#' \itemize{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}}
#' \item A directory with f2 statistics
#' \item The prefix of a genotype file
#' }
#' @param models A nested list (or data frame) with qpadm models. It should consist of two or three other named lists (or columns) containing `left`, `right`, (and `target`) populations.
#' @param allsnps Use all SNPs with allele frequency estimates in every population of any given population quadruple. If `FALSE` (the default) only SNPs which are present in all populations in `popcombs` (or any given model in it) will be used. When there are populations with lots of (non-randomly) missing data, `allsnps = TRUE` can lead to false positive results. This option only has an effect if `data` is the prefix of a genotype file. If `data` are f2-statistics, the behavior will be determined by the options that were used in computing the f2-statistics.
#' @param verbose Print progress updates
#' @param ... Further arguments passed to \code{\link{qpadm}}
#' @return A list where each element is the output of one qpadm model.
#' @examples
#' \dontrun{
#' pops = paste0('pop', 1:10)
#' # the following specifies two models: one with 2/3/1 and one with 1/2/1 left/right/target populations
#' models = tibble(
#'            left = list(c('pop1', 'pop2'), c('pop3')),
#'            right = list(c('pop5', 'pop6', 'pop7'), c('pop7', 'pop8')),
#'            target = c('pop10', 'pop10'))
#' results = qpadm_multi('/my/geno/prefix', models)
#' }
qpadm_multi = function(data, models, allsnps = FALSE, verbose = TRUE, ...) {

  if(!'left' %in% names(models) || !'right' %in% names(models))
    stop("'models' should have elements 'left' and 'right'!")
  if(length(unique(map_dbl(models, length))) != 1)
    stop("'left', 'right' (and 'target') should be of equal length!")

  if('target' %in% names(models)) models$left = map2(models$left, models$target, ~union(.y, .x))
  if(!'tibble' %in% class(models)) models %<>% as_tibble

  if(is_geno_prefix(data)) {
    popcombs = qpadmmodels_to_popcombs(models)
    f4blockdat = data %>% f4blockdat_from_geno(popcombs, allsnps = allsnps, verbose = verbose)
    f4blocks = f4blockdat %>% split(.$model) %>% furrr::future_map(quietly(f4blockdat_to_f4blocks)) %>% map('result')
  } else {
    if(verbose && !allsnps) alert_warning('allsnps = FALSE is not effective when using precomputed f2 statistics\n')
    pops = models %>% select(any_of(c('target', 'left'))) %>% unlist %>% unique
    pops2 = models %>% select(right) %>% unlist %>% unique
    f2blocks = data %>% get_f2(pops, pops2, afprod = TRUE)
    f4blocks = models %>% rowwise %>% mutate(f4blocks = list(f2blocks_to_f4blocks(f2blocks, left, right))) %>% pull
  }

  if(verbose) alert_info('Running models...\n')
  f4blocks %>%
    furrr::future_map(function(.x, ...)
      qpadm(NULL, NULL, NULL, target = TRUE, f4blocks = .x, verbose = FALSE, ...),
      .progress = verbose, ...)
}


