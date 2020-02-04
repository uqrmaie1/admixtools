
# a, b: num left - 1, num right - 1 (dimensions of f4_est)
qpwave_dof = function(a, b, r) r*(a+b-r)
qpwave_dof = function(a, b, r) (a+b)*r - r^2
qpwave_dof = function(a, b, r) a*b - (a-r)*(b-r)

qpadm_dof = function(a, b, r) a*b - (a+b)*r + r^2
qpadm_dof = function(a, b, r) (a-r)*(b-r)
qpadm_dof = function(a, b, r) a*b - qpwave_dof(a, b, r)


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
  if(!boot) wmat = wmat * sqrt(numreps-1)
  cov(wmat)
}


qpadm_weights = function(xmat, qinv, rnk, fudge = 0.0001, iterations = 20,
                         constrained = FALSE, qpsolve = NULL) {
  nc = nrow(xmat)
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
  #w = solve(t(cbind(A, 1)), c(rep(0, rnk), 1))
  if(constrained) w = qpsolve(rhs, lhs, diag(nc), rep(0, nc))
  else w = as.matrix(solve(rhs, lhs))[,1]
  weights = w/sum(w)
  namedList(weights, A, B)
}


#' Estimate admixture weights
#'
#' Models target as a mixture of left populations, given a set of outgroup right populations.
#' Can be used to estimate admixture proportions, and to estimate the number of independent
#' admixture events (see `wave`).
#' @export
#' @param f2_data a 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' alternatively, a directory with precomputed data. see \code{\link{extract_f2}} and \code{\link{extract_indpairs}}.
#' @param target target population
#' @param left source populations
#' @param right outgroup populations
#' @param wave If `TRUE` (the default), lower rank models and models with fewer source populations will also be tested.
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param fudge value added to diagonal matrix elements before inverting
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix
#' of f4 statistics, as well as the weight standard errors, will be computed using block-jackknife.
#' Otherwise bootstrap resampling is performed `n` times, where `n` is either equal to `boot` if it is an integer,
#' or equal to the number of blocks if `boot` is `TRUE`. The the covariance matrix of f4 statistics,
#' as well as the weight standard errors, will be computed using bootstrapping.
#' @param getcov should standard errors be returned? Setting this to `FALSE` makes this function much faster.
#' @param constrained if `FALSE` (default), admixture weights can be negative.
#' if `TRUE`, they will all be non-negative, as in \code{\link{lazadm}}
#' @param cpp should optimization be done using C++ or R function? `cpp = TRUE` is much faster.
#' @param verbose print progress updates
#' @return a list with output describing the model fit:
#' \enumerate{
#' \item `weights` a data frame with estimated admixture proportions where each row is left population.
#' estimated edge length, and for admixture edges, it is the estimated admixture weight.
#' \item `rankdrop` a data frame describing model fits with different ranks, including p-values for the overall fit
#' and for nested models (comparing two models with rank difference of one). returned only if `wave = TRUE`
#' \item `popdrop` a data frame describing model fits with different populations. returned only if `wave = TRUE`
#' }
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 10)
#' @aliases qpwave
#' @section Alias:
#' `qpwave`
#' @seealso \code{\link{lazadm}}
#' @examples
#' target = 'Denisova.DG'
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' qpadm(example_f2_blocks, target, left, right)
qpadm = function(f2_data, target = NULL, left = NULL, right = NULL,
                 wave = TRUE, f2_denom = 1, fudge = 0.0001, boot = FALSE,
                 getcov = TRUE, constrained = FALSE, cpp = TRUE, verbose = TRUE) {

  stopifnot(is.null(left) && is.null(right) || length(right) > length(left))

  if(is.null(target)) {
    left %<>% readLines
    target = left[1]
    left = left[-1]
    right %<>% readLines
  }
  pops = c(target, left, right)
  stopifnot(!any(duplicated(pops)))

  #----------------- prepare f4 stats -----------------
  if(verbose) alert_info('Computing f4 stats...\n')
  f2_blocks = get_f2(f2_data, pops, f2_denom)
  f2_blocks %<>% rray(dim_names = dimnames(f2_blocks))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, boot_pairarr_stats, jack_pairarr_stats)

  f4_blocks = (f2_blocks[target, right[-1], ] +
               f2_blocks[left, right[1], ] -
               f2_blocks[target, right[1], ] -
               f2_blocks[left, right[-1], ])/2

  f4_lo = f4_blocks %>% samplefun
  block_lengths = parse_number(dimnames(f4_lo)[[3]])
  f4dat = f4_lo %>% statfun(block_lengths)
  f4_est = f4dat$est
  f4_var = f4dat$var

  #----------------- compute admixture weights -----------------
  if(verbose) alert_info('Computing admixture weights...\n')
  diag(f4_var) = diag(f4_var) + fudge*sum(diag(f4_var))
  qinv = solve(f4_var)
  rnk = length(left)-1

  if(cpp) {
    qpadm_weights = cpp_qpadm_weights
    get_weights_covariance = cpp_get_weights_covariance
  }
  weight = qpadm_weights(f4_est, qinv, rnk, fudge = fudge, constrained = constrained,
                         qpsolve = qpsolve)$weights %>% c
  if(getcov) {
    if(verbose) alert_info('Computing standard errors...\n')
    se = sqrt(diag(get_weights_covariance(f4_lo, qinv, block_lengths, fudge = fudge, boot = boot,
                                          constrained = constrained, qpsolve = qpsolve)))
  } else se = rep(NA, length(weight))

  out = list(weights = tibble(target, left, weight, se) %>% mutate(z = weight/se))

  #----------------- compute number of admixture waves -----------------
  if(wave) {
    if(verbose) alert_info('Computing number of admixture waves...\n')

    rankdrop = drop_ranks(f4_est, qinv, fudge, constrained, cpp)
    popdrop = drop_pops(f4_est, qinv, fudge, constrained, cpp)
    out = c(out, namedList(rankdrop, popdrop))
  }
  out
}

#' @export
qpwave = qpadm


#' Estimate admixture weights
#'
#' Models target as a mixture of left populations, and outgroup right populations. Uses Lazaridis method
#' based non-negative least squares of f4 matrix.
#' @export
#' @param f2_data a 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' alternatively, a directory with precomputed data. see \code{\link{extract_f2}} and \code{\link{extract_indpairs}}.
#' @param target target population
#' @param left source populations (or `leftlist` file)
#' @param right outgroup populations (or `rightlist` file)
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param getcov should standard errors be returned? Currently not implemented.
#' @param constrained if `TRUE` (default), admixture weights will all be non-negative.
#' if `FALSE`, they can be negative, as in \code{\link{qpadm}}
#' @return a data frame with weights and standard errors for each left population
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 9)
#' @seealso \code{\link{qpadm}}
#' @examples
#' target = 'Denisova.DG'
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' lazadm(example_f2_blocks, target, left, right)
#' lazadm(example_f2_blocks, target, left, right, constrained = FALSE)
lazadm = function(f2_data, target = NULL, left = NULL, right = NULL,
                  f2_denom = 1, boot = FALSE, getcov = FALSE, constrained = TRUE) {

  #----------------- prepare f4 stats -----------------
  if(is.null(target)) {
    left %<>% readLines
    target = left[1]
    left = left[-1]
    right %<>% readLines
  }
  pops = c(target, left, right)
  stopifnot(!any(duplicated(pops)))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  f2_blocks = get_f2(f2_data, pops, f2_denom) %>% samplefun
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

  if(constrained) weight = -qpsolve(rhs, lhs, -diag(nc), rep(0, nc))
  else weight = solve(rhs, lhs)[,1]
  weight = weight/sum(weight)

  # todo: implement this
  if(getcov) se = rep(NA, length(weight))
  else se = rep(NA, length(weight))

  tibble(target, left, weight, se) %>% mutate(z = weight/se)
}


drop_ranks = function(f4_est, qinv, fudge, constrained, cpp) {
  # drops rank and fits qpadm model

  rnk = nrow(f4_est) - 1
  fitrank = function(x) qpadm_fit(f4_est, qinv, x, fudge = fudge,
                                  constrained = constrained, cpp = cpp, addweights = FALSE)
  rankdrop = map_dfr(rev(seq_len(rnk)), fitrank) %>%
    mutate(dofdiff = lead(dof)-dof, chisqdiff = lead(chisq)-chisq,
           p_nested = pchisq(chisqdiff, dofdiff, lower.tail = FALSE))
  rankdrop
}


drop_pops = function(f4_est, qinv, fudge, constrained, cpp) {
  # drops each subset of left populations and fits qpadm model

  fitpop = function(x, y) qpadm_fit(x, y, nrow(x)-1, fudge = fudge, constrained = constrained,
                                    cpp = cpp, addweights = TRUE)
  nc = ncol(f4_est)
  popdrop = left %>% power_set %>% enframe(name = 'i', value = 'pop') %>%
    filter(map_int(pop, length) > 1) %>%
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
  for(i in rev(seq_len(rnk-1))) {
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



qpadm_evaluate_fit = function(xmat, qinv, A, B) {

  f4rank = ncol(A)
  res = t(xmat - A %*% B)
  chisq = (t(c(res)) %*% qinv %*% c(res))[,1]
  dof = qpadm_dof(nrow(A), ncol(B), f4rank)
  p = pchisq(chisq, df = dof, lower.tail = FALSE)
  tibble(f4rank, dof, chisq, p)
}


qpadm_fit = function(xmat, qinv, rnk, fudge = 0.0001, iterations = 20,
                     constrained = FALSE, cpp = TRUE, addweights = FALSE) {
  # returns one-row data frame with fit of one qpadm model

  if(cpp) qpadm_weights = cpp_qpadm_weights
  fit = qpadm_weights(xmat, qinv, rnk, fudge = fudge, constrained = constrained, qpsolve = qpsolve)
  out = qpadm_evaluate_fit(xmat, qinv, fit$A, fit$B)
  if(addweights) {
    out %<>% bind_cols(fit$weights %>% t %>% as_data_frame %>% set_names(rownames(xmat)))
  }
  out
}

