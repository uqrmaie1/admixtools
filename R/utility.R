

alert_success = function(msg) cat(crayon::green(cli::symbol$tick, msg))
alert_info = function(msg) cat(crayon::cyan(cli::symbol$info, msg))
alert_warning = function(msg) cat(crayon::yellow('!', msg))
alert_danger = function(msg) cat(crayon::red(cli::symbol$cross, msg))


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


arr3d_to_pairmat = function(arr, diag=TRUE) {
  # input is arr of dimension n x n x m
  # output is choose(n+1, 2) x m
  indx1 = which(lower.tri(arr[,,1], diag = diag))
  d1 = dim(arr)[1]
  d3 = dim(arr)[3]
  npair = choose(d1+diag, 2)
  indx = rep(indx1, d3) + rep((0:(d3-1))*d1^2, each=npair)
  matrix(as.array(arr)[indx], npair, d3)
}


arr3d_to_mat = function(arr) {
  # input is arr of dimension m x n x p
  # output is p x (m x n) (1..m, 1..m, ... times n)
  matrix(arr, ncol=dim(arr)[3])
}


power_set = function(l, nmax=length(l)) purrr::flatten(purrr::map(seq_along(l[seq_len(nmax)]), ~combn(l, ., simplify=F)))



gg_color_hue = function(n, l=65, c=100) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=c)[1:n]
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
    ans = optim(par = parmat[i, ], fn = fn, gr = gr, lower = lower,
                upper = upper, method = method, hessian = hessian,
                control = control, args=args, ...)
    ansret[i, ] = c(ans$par, ans$value, ans$counts[1], ans$counts[2], ans$convergence)
    if(verbose) cat(paste0('\r', i, ' out of ', nset))
  }
  if(verbose) cat('\n')
  as_tibble(ansret) %>% set_colnames(c(paste0('p', seq_len(npar)), 'value', 'fevals', 'gevals', 'convergence'))
}
