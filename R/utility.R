

catfun = function(colfun, ...) {
  if(shiny::isRunning()) {
    l = list(..., appendLF = FALSE)
    do.call(message, l[-1])
  } else {
    cat(colfun(...))
  }
}

alert_success = function(msg) catfun(crayon::green, cli::symbol$tick, msg)
alert_info = function(msg) catfun(crayon::cyan, cli::symbol$info, msg)
alert_warning = function(msg) catfun(crayon::yellow, '!', msg)
alert_danger = function(msg) catfun(crayon::red, cli::symbol$cross, msg)

gg_color_hue = function(n, l=65, c=100) {
  hues = seq(15, 375, length=n+1)
  grDevices::hcl(h=hues, l=l, c=c)[1:n]
}

qpsolve = function(...) quadprog::solve.QP(...)$solution
#qpsolve = function(cc, q1, a3, a4) -quadprogpp::QP.Solve(cc, q1, -a3, -a4)
#qpsolve = function(cc, q1, a3, a4) -rcppeigen_quadratic_solve(cc, q1, matrix(0, nrow(cc), 0), vector(), -a3, -a4)[,1]


# qpsolve = function(...) tryCatch({quadprog::solve.QP(...)$solution},
#                           error = function(e) {
#                             if(str_detect(e$message, 'constraints are inconsistent')) {
#                               warning(e$message)
#                               #browser()
#                               # set lower bounds to -1e3 for solve.QP, then truncate at 0
#                               ell = list(...)
#                               ilow = (1:(length(ell[[4]])/2))
#                               low = ell[[4]][ilow]
#                               ell[[4]][ilow] = -1e3
#                               cc = solve((ell[[1]]+t(ell[[1]]))/2)
#                               #diag(cc) = diag(cc) + 0.0001 # this will result in extremely high scores for some graphs
#                               #diag(cc) = diag(cc) + 0.0001*mean(diag(cc)) # this will sometimes lead to 'constraints inconsistent
#                               diag(cc) = diag(cc) + 0.01*mean(diag(cc))
#                               cc = (cc+t(cc))/2
#                               ell[[1]] = solve(cc)
#                               return(do.call(quadprog::solve.QP, ell)$solution)
#                               #return(pmax(low, do.call(quadprog::solve.QP, ell)$solution))
#                             } else {
#                               stop(e)
#                             }
#                           })

# recurisvely increases fudge if constraints are inconsistent
# qpsolve = function(...) tryCatch({quadprog::solve.QP(...)$solution},
#                           error = function(e) {
#                             if(str_detect(e$message, 'constraints are inconsistent')) {
#                               ell = list(...)
#                               cc = solve(ell[[1]]); diag(cc) = diag(cc)+sum(diag(cc))*0.1; ell[[1]] = solve(cc)
#                               return(do.call(qpsolve, ell))
#                               #return(pmax(low, do.call(quadprog::solve.QP, ell)$solution))
#                             } else {
#                               stop(e)
#                             }
#                           })

make_pseudohaploid = function(xmat) round(jitter(xmat/2))*2

#' @export
namedList = function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}


mean_impute = function(mat, by=2) {
  out = apply(mat, by, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
  if(by == 1) out = t(out)
  out
}


arr3d_to_pairmat = function(arr, diag = TRUE) {
  # input is arr of dimension n x n x m
  # output is choose(n+1, 2) x m
  indx1 = which(lower.tri(arr[,,1], diag = diag))
  d1 = dim(arr)[1]
  d3 = dim(arr)[3]
  npair = choose(d1+diag, 2)
  indx = rep(indx1, d3) + rep((0:(d3-1))*d1^2, each=npair)
  matrix(as.array(arr)[indx], npair, d3)
}


prodarray = function(m1, m2) {
  (array(m1, c(nrow(m1), ncol(m1), 1))[,,rep(1, ncol(m2))] *
    array(m2, c(nrow(m1), 1, ncol(m2)))[,rep(1, ncol(m1)),]) %>%
    aperm(c(2,3,1))
}

outer_array = function(m1, m2, FUN = `*`) {
  nr = nrow(m1)
  nc1 = ncol(m1)
  nc2 = ncol(m2)
  FUN(array(m1, c(nr, nc1, 1))[,,rep(1, nc2)],
      array(m2, c(nr, 1, nc2))[,rep(1, nc1),]) %>%
    aperm(c(2,3,1))
}

arr3d_to_mat = function(arr) {
  # input is arr of dimension m x n x p
  # output is p x (m x n) (1..m, 1..m, ... times n)
  matrix(arr, nrow = dim(arr)[3], byrow = TRUE)
}

mat_to_arr3d = function(mat, dim1) {
  # input is p x (m x n) (1..m, 1..m, ... times n)
  # output is arr of dimension m x n x p
  stopifnot(ncol(mat) %% dim1 == 0)
  array(t(mat), c(dim1, ncol(mat)/dim1, nrow(mat)))
}

mat_to_rep21 = function(mat, rep) {
  # input is m x n matrix
  # output is rep x n x m array
  array(rep(t(mat), each = rep), c(rep, ncol(mat), nrow(mat)))
}
mat_to_2rep1 = function(mat, rep) {
  # input is m x n matrix
  # output is n x rep x m array
  array(rep(t(mat), each = rep), c(rep, ncol(mat), nrow(mat))) %>%
    aperm(c(2,1,3))
}

# for a m*m*n 3d array, return the m*n diagonal matrix
diag_3d = function(arr) {
  arr %<>% as.array
  d1 = dim(arr)[1]
  d3 = dim(arr)[3]
  matrix(arr[as.matrix(expand.grid(1:d1, 1:d3)[,c(1,1,2)])], d1)
}

weighted_row_means = function(mat, weights, na.rm = TRUE) {

  if(!is.matrix(weights)) weights = matrix(weights, nrow(mat), ncol(mat), byrow = TRUE)
  rowSums(mat * weights, na.rm = na.rm) / rowSums((!is.na(mat)) * weights, na.rm = na.rm)
}


power_set = function(l, nmax=length(l)) flatten(map(seq_along(l[seq_len(nmax)]), ~combn(l, ., simplify=F)))

ztop = function(z) 2*pnorm(-abs(z))


interleave = function(v1,v2) {
  # zips two vectors
  ord1 = 2*(1:length(v1))-1
  ord2 = 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}


multistart = function (parmat, fn, args, gr = NULL, lower = -Inf, upper = Inf, method = NULL,
                       hessian = FALSE, control = list(), verbose=TRUE, ...) {
  # same as from optimx library, but can be quiet
  nset = nrow(parmat)
  npar = ncol(parmat)
  if(nset < 1) stop("multistart: no starting parameters!")
  ansret = matrix(NA, nset, npar + 4)
  for (i in 1:nset) {
    # ans = cpp_lbfgsb(parmat[i,], args[[1]], args[[2]], args[[3]], args[[4]], args[[5]], args[[6]], args[[8]])
    # cpp_lbfgsb is ~ 10% faster for large graphs (20 pops, 12 admixnodes), but very similar in speed for smaller ones.
    # only returns 'par' and 'value' at the moment
    tryCatch({
      ans = optim(par = parmat[i, ], fn = fn, gr = gr, lower = lower,
                  upper = upper, method = method, hessian = hessian,
                  control = control, args = args, ...)
      ansret[i, ] = c(ans$par, ans$value, ans$counts[1], ans$counts[2], ans$convergence)
    }, error = function(e) {if(!str_detect(e$message, 'constraints are inconsistent')) stop(e)})
    if(verbose) cat(paste0('\r', i, ' out of ', nset))
  }
  if(verbose) cat('\n')
  if(all(is.na(ansret[,1]))) stop("Optimization unsuccessful (constraints inconsistent). Try increasing 'diag'.")
  as_tibble(ansret, .name_repair = ~c(paste0('p', seq_len(npar)), 'value', 'fevals', 'gevals', 'convergence'))
}

# on OSX, multiprocess multistart is only slightly faster than sequential in qpGraph. Need to test on Linux.
multistart2 = function (parmat, fn, args, gr = NULL, lower = -Inf, upper = Inf, method = NULL,
                       hessian = FALSE, control = list(), verbose = TRUE, ...) {

  nset = nrow(parmat)
  npar = ncol(parmat)
  if(nset < 1) stop("multistart: no starting parameters!")
  ansret = matrix(NA, nset, npar + 4)
  furrr::future_map(seq_len(nset), function(.x, ...) {
    ans = optim(par = parmat[.x, ], fn = fn, gr = gr, lower = lower,
                upper = upper, method = method, hessian = hessian,
                control = control, args = args, ...)
    c(ans$par, ans$value, ans$counts[1], ans$counts[2], ans$convergence)
  }, .progress = verbose, ...) %>%
    do.call(rbind, .) %>%
    as_tibble %>%
    set_colnames(c(paste0('p', seq_len(npar)), 'value', 'fevals', 'gevals', 'convergence'))
}

#' Return shortest unique prefixes
#'
#' @export
#' @param strings A character vector
#' @param min_length Minimum length of prefixes
#' @return Character vector with shortest unique prefixes
#' @examples
#' strings = c("Abra", "Abracadabra", "Simsalabim")
#' shortest_unique_prefixes(strings)
shortest_unique_prefixes = function(strings, min_length = 1) {
  # given a vector of strings, return a vector of the same length with the shortest unique prefixes
  if(length(strings) == 0) return(strings)
  len = rep(0, length(strings))
  for(i in 1:max(nchar(strings))) {
    pref = str_sub(strings, 1, i)
    tab = which(pref %in% names(which(table(pref) == 1)))
    tab = intersect(tab, which(len == 0))
    len[tab] = i
  }
  str_sub(strings, 1, pmax(len, min_length))
}


dflist_to_confints = function(dflist, boot = FALSE) {
  # takes a list of data frames with identical columns
  # groups by all non-numeric variables
  if(boot) {
    out = dflist %>%
      bind_rows %>%
      group_by_if(function(x) !is.numeric(x)) %>%
      summarize_all(list(
        cnt = ~n(),
        mean = ~mean(., na.rm=T),
        se2 = ~sum((mean(., na.rm=T) - .)^2, na.rm=T),
        lo = ~quantile(., 0.005),
        hi = ~quantile(., 0.995))) %>%
        ungroup %>% mutate(se2 = se2/cnt, se = sqrt(se2))
  } else {
    out = dflist %>%
      bind_rows %>%
      group_by_if(function(x) !is.numeric(x)) %>%
      summarize_all(list(
        cnt = ~n(),
        mean = ~mean(., na.rm=T),
        se2 = ~sum((mean(., na.rm=T) - .)^2, na.rm=T))) %>%
      ungroup %>% mutate(se = sqrt(se2)) %>%
      mutate(lo = mean - 1.96*se, hi = mean + 1.96*se)
  }
  out %>%
    pivot_longer(c(mean, lo, hi), names_to = 'k', values_to = 'v') %>%
    mutate(lab = ifelse(type == 'admix', paste0(round(v*100), '%'), round(v*1000))) %>%
    pivot_wider(names_from = k, values_from = c(v, lab)) %>%
    rename(mean = v_mean, lo = v_lo, hi = v_hi) %>%
    mutate(label = paste0(lab_mean, '\n[', lab_lo, '-', lab_hi, ']')) %>%
    select(-se2, -lab_mean, -lab_lo, -lab_hi)
}

replace_nan_with_na = function(arr) {
  arr[is.nan(arr)] = NA
  arr
}


