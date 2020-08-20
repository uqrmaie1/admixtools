
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
  if(constrained) w = qpsolve(rhs, lhs, diag(nc), rep(0, nc))
  else w = as.matrix(solve(rhs, lhs))[,1]
  weights = w/sum(w)
  namedList(weights, A, B)
}


#' Estimate admixture weights
#'
#' Models target as a mixture of left populations, given a set of outgroup right populations.
#' Can be used to estimate admixture proportions, and to estimate the number of independent
#' admixture events.
#' @export
#' @param f2_data A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}}. Alternatively, a directory with f2 statistics.
#' @param left Left populations (sources)
#' @param right Right populations (outgroups)
#' @param target Target population
#' @param fudge value added to diagonal matrix elements before inverting
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix
#' of f4 statistics, as well as the weight standard errors, will be computed using block-jackknife.
#' Otherwise bootstrap resampling is performed `n` times, where `n` is either equal to `boot` if it is an integer,
#' or equal to the number of blocks if `boot` is `TRUE`. The the covariance matrix of f4 statistics,
#' as well as the weight standard errors, will be computed using bootstrapping.
#' @param getcov Compute weights covariance. If not needed, turning this off will speed things up.
#' @param constrained Constrain admixture weights to be non-negative
#' @param cpp Use C++ functions. Setting this to `FALSE` will be slower but can help with debugging.
#' @param verbose Print progress updates
#' @return A list with output describing the model fit:
#' \enumerate{
#' \item `weights` A data frame with estimated admixture proportions where each row is a left population.
#' \item `f4` A data frame with estimated and fitted f4-statistics
#' \item `rankdrop` A data frame describing model fits with different ranks, including p-values for the overall fit
#' and for nested models (comparing two models with rank difference of one).
#' \item `popdrop` A data frame describing model fits with different populations.
#' }
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 10)
#' @seealso \code{\link{qpwave}}, \code{\link{lazadm}}
#' @examples
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' target = 'Denisova.DG'
#' qpadm(example_f2_blocks, left, right, target)
qpadm = function(f2_data, left, right, target,
                 fudge = 0.0001, boot = FALSE, getcov = TRUE,
                 constrained = FALSE, cpp = TRUE, verbose = TRUE) {

  if(all(file.exists(left, right))) {
    left %<>% readLines
    right %<>% readLines
  }
  if(!is.null(target)) left = c(target, setdiff(left, target))
  pops = c(left, right)
  if(any(duplicated(pops))) stop(paste0('Duplicated pops: ', paste0(unique(pops[duplicated(pops)]),collapse=', ')))

  #----------------- prepare f4 stats -----------------
  if(verbose) alert_info('Computing f4 stats...\n')
  f2_blocks = get_f2(f2_data, pops = pops, afprod = TRUE)
  f4dat = f2_to_f4(f2_blocks, left, right, boot = boot)

  f4_est = f4dat$est
  f4_var = f4dat$var
  f4_lo = f4dat$f4_lo
  block_lengths = f4dat$block_lengths
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
    #if(isTRUE(all.equal(dimnames(f2_blocks)[[1]], dimnames(f2_blocks)[[2]]))) {
    out$f4 = fitted_f4(f2_blocks, wvec, target, left[-1], right)
    #}
  } else {
    out$f4 = f4(f2_blocks, left[1], left[-1], right[1], right[-1], verbose = FALSE)
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
#' @return A list with output describing the model fit:
#' \enumerate{
#' \item `f4` A data frame with estimated f4-statistics
#' \item `rankdrop` A data frame describing model fits with different ranks, including p-values for the overall fit
#' and for nested models (comparing two models with rank difference of one).
#' }
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 10)
#' @seealso \code{\link{qpadm}}
#' @examples
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' qpwave(example_f2_blocks, left, right)
qpwave = function(f2_data, left, right,
                  fudge = 0.0001, boot = FALSE,
                  constrained = FALSE, cpp = TRUE, verbose = TRUE)
  qpadm(f2_data = f2_data, left = left, right = right, target = NULL,
        fudge = fudge, boot = boot,
        constrained = constrained, cpp = cpp, verbose = verbose)


f2_to_f4 = function(f2_blocks, left, right, boot = FALSE) {

  if(length(right) < 2) stop('Not enough right populations!')
  if(length(left) < 2) stop('Not enough left populations!')
  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo_nafix)
  statfun = ifelse(boot, boot_pairarr_stats, jack_pairarr_stats)

  # f4_blocks = (f2_blocks[left, right[1], ] +
  #              f2_blocks[target, right[-1], ] -
  #              f2_blocks[target, right[1], ] -
  #              f2_blocks[left, right[-1], ])/2

  nr = length(right) - 1
  nl = length(left) - 1
  f4_blocks = (f2_blocks[left[-1], rep(right[1], nr), , drop = FALSE] +
               f2_blocks[rep(left[1], nl), right[-1], , drop = FALSE] -
               f2_blocks[rep(left[1], nl), rep(right[1], nr), , drop = FALSE] -
               f2_blocks[left[-1], right[-1], , drop = FALSE])/2

  f4_lo = f4_blocks %>% samplefun
  block_lengths = parse_number(dimnames(f4_lo)[[3]])

  out = f4_lo %>% statfun(block_lengths)
  #out = f4_lo %>% jack_pairarr_stats(block_lengths)
  out$f4_lo = f4_lo
  out$block_lengths = block_lengths
  out
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
#' @param f2_data A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}}.
#' @param left Left populations (sources)
#' @param right Right populations (outgroups)
#' @param target Target population
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix
#' of f4 statistics, as well as the weight standard errors, will be computed using block-jackknife.
#' Otherwise bootstrap resampling is performed `n` times, where `n` is either equal to `boot` if it is an integer,
#' or equal to the number of blocks if `boot` is `TRUE`. The the covariance matrix of f4 statistics,
#' as well as the weight standard errors, will be computed using bootstrapping.
#' @param constrained If `TRUE` (default), admixture weights will all be non-negative.
#' If `FALSE`, they can be negative, as in \code{\link{qpadm}}
#' @return A data frame with weights and standard errors for each left population
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
lazadm = function(f2_data, left, right, target,
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
  f2_blocks[f2_blocks < 0] = 0
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])
  numblocks = length(block_lengths)

  r = 1:length(right)
  og_indices = expand.grid(r, r, r) %>%
    filter(Var1 != Var2, Var1 != Var3, Var2 < Var3)
  ncomb = nrow(og_indices)
  nleft = length(left)
  blocknums = rep(1:numblocks, each = ncomb)

  pos1 = 1
  pos2 = og_indices[,1]
  pos3 = og_indices[,2]
  pos4 = og_indices[,3]

  ymat = matrix(f2_blocks[cbind(pos1, pos4, blocknums)] +
                f2_blocks[cbind(pos2, pos3, blocknums)] -
                f2_blocks[cbind(pos1, pos3, blocknums)] -
                f2_blocks[cbind(pos2, pos4, blocknums)], ncomb)

  xarr = f2_blocks[pos4, left,] +
    array(f2_blocks[cbind(pos2, pos3, rep(blocknums, each = nleft))], c(ncomb, nleft, numblocks)) -
    f2_blocks[pos3, left,] -
    array(f2_blocks[cbind(pos2, pos4, rep(blocknums, each = nleft))], c(ncomb, nleft, numblocks))

  #----------------- compute admixture weights -----------------
  lhs = matrix(NA, nleft, numblocks)
  rhs = array(NA, c(nleft, nleft, numblocks))
  for(i in 1:numblocks) {
    lhs[,i] = crossprod(xarr[,,i], ymat[,i])
    rhs[,,i] = crossprod(xarr[,,i])
  }

  if(constrained) fun = function(rhs, lhs) qpsolve(rhs, lhs, diag(nleft), rep(0, nleft))
  else fun = function(rhs, lhs) solve(rhs, lhs)[,1]

  wmat = sapply(1:numblocks, function(i) {w = fun(rhs[,,i], lhs[,i,drop = FALSE]); w/sum(w)})
  weight = rowMeans(wmat)
  se = sqrt(diag(cov(t(wmat))))
  if(!boot) se = se * (numblocks-1) / sqrt(numblocks)

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
  f4dat = f2_to_f4(f2_blocks, left, right, boot = boot)
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


