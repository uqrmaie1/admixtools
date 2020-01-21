

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
  matrix(arr, ncol = dim(arr)[3])
}

#' @export
indpairs_to_f2blocks = function(indivs, pairs, poplist, block_lengths) {
  # creates f2_blocks from per individual data
  # make fast version in Rcpp and use tibbles here for readability
  # indivs: data frame with columns ind, bl, a, n
  # pairs: data frame with columns ind1, ind2, bl, aa, nn
  # poplist: data frame with columns ind, pop
  # block_lengths
  # a is number of alt alleles, n number of ref + alt alleles.

  pairsums = pairs %>%
    left_join(poplist %>% transmute(ind1 = ind, pop1 = pop), by = 'ind1') %>%
    left_join(poplist %>% transmute(ind2 = ind, pop2 = pop), by = 'ind2') %>%
    left_join(indivs %>% transmute(ind1 = ind, bl), by = c('ind1', 'bl')) %>%
    left_join(indivs %>% transmute(ind2 = ind, bl), by = c('ind2', 'bl')) %>%
    group_by(pop1, pop2, bl) %>% summarize(aa = sum(aa), nn = sum(nn), pp = aa/nn) %>% ungroup
  # check if it should rather be pp = sum(aa/nn)

  pairsums_samepop = pairsums %>% filter(pop1 == pop2) %>% transmute(pop = pop1, bl, aa, nn, pp)

  main = pairsums %>%
    left_join(pairsums_samepop %>% transmute(pop1 = pop, bl, pp1 = pp), by = c('pop1', 'bl')) %>%
    left_join(pairsums_samepop %>% transmute(pop2 = pop, bl, pp2 = pp), by = c('pop2', 'bl')) %>%
    mutate(f2uncorr = pp1 + pp2 - 2*pp)

  indsums = indivs %>% left_join(poplist, by='ind') %>%
    group_by(pop, bl) %>% summarize(a = sum(a), n = sum(n)) %>% ungroup
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
    bind_rows(filter(., pop1 != pop2 & cnt == 1) %>% rename(pop1 = pop2, pop2 = pop1, pp1 = pp2, pp2 = pp1, corr1 = corr2, corr2 = corr1)) %>% arrange(bl, pop2, pop1)

  popnames = unique(poplist$pop)
  popnames2 = unique(pairsums_samepop$pop)
  npop = length(popnames)
  nblock = length(block_lengths)
  nsnp = sum(block_lengths)

  array(f2dat$f2, dim = c(npop, npop, nblock),
        dimnames = list(pop1 = popnames2,
                        pop2 = popnames2,
                        bl = paste0('bl', 1:nblock))) %>%
    `[`(popnames, popnames, ) %>%
    ifelse(is.nan(.), NA, .) %>%
    structure(block_lengths = block_lengths) %>%
    est_to_loo
}


update_f2blocks = function(indivs, pairs, poplist, f2_blocks) {
  # creates f2_blocks from per individual data
  # make fast version in Rcpp and use tibbles here for readability
  # indivs: data frame with columns ind, bl, a, n
  # pairs: data frame with columns ind1, ind2, bl, aa, nn
  # poplist: data frame with columns ind, pop
  # block_lengths
  # a is number of alt alleles, n number of ref + alt alleles.


}




fix_ploidy = function(xmat) {
  # divides pseudohaploid columns by 2
  # returns list with updated xmat, with named 1/2 ploidy vector as attribute

  ploidy = apply(xmat, 2, function(x) length(unique(na.omit(x))))-1
  maxgt = apply(xmat, 2, max, na.rm = TRUE)
  xmat[, ploidy == 1 && maxgt == 2] = xmat[, ploidy == 1 && maxgt == 2]/2
  xmat %<>% structure(ploidy = ploidy)
  xmat
}


treat_missing = function(afmatrix, countmatrix = NULL, snpfile = NULL,
                         na.action = 'none', verbose = TRUE) {
  discard = FALSE
  if(na.action == 'impute') {
    afmatrix = mean_impute(afmatrix, by=1)
    isnan = is.nan(afmatrix)
    afmatrix[isnan] = NA
    discard = apply(afmatrix, 1, function(x) all(is.na(x)))
    if(verbose) alert_danger(paste0(sumna-sum(isnan), ' allele frequencies were imputed, (on average ', round((sumna-sum(isnan))/numpop),' per population) ', sum(discard), ' SNPs were removed\n'))
  }
  if(na.action == 'remove') {
    discard = apply(afmatrix, 1, function(x) any(is.na(x)))
    if(verbose) alert_danger(paste0(sum(discard), ' SNPs were removed, ', sum(!discard), ' SNPs remain\n'))
  }
  afmatrix = afmatrix[!discard,]
  countmatrix = countmatrix[!discard,]
  snpfile = snpfile[!discard,]

  namedList(afmatrix, countmatrix, snpfile)
}



power_set = function(l, nmax=length(l)) purrr::flatten(purrr::map(seq_along(l[seq_len(nmax)]), ~combn(l, ., simplify=F)))

ztop = function(z) 2*pnorm(-abs(z))

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
