

graph_to_pwts = function(graph, leaves) {
  # input igraph object
  # output: numedges*(numpops-1) matrix 'pwts' which indicates all paths from pop to outpop; admixture edges are mapped onto parent edges and weighted

  root = get_root(graph)
  pwts = matrix(0, length(E(graph)), length(leaves))
  colnames(pwts) = leaves
  rownames(pwts) = attr(E(graph), 'vnames')

  admixnodes = which(degree(graph, mode='in') == 2)
  admixedges = unlist(incident_edges(graph, admixnodes, mode='in'))

  allpaths = all_simple_paths(graph, root, leaves, mode='out')
  pathcounts = table(names(sapply(allpaths, function(x) tail(x, 1))))
  for(i in seq_len(length(allpaths))) {
    target = names(tail(allpaths[[i]],1))
    ln = length(allpaths[[i]])
    pth2 = allpaths[[i]][c(1, 1+rep(seq_len(ln-2), each=2), ln)]
    rowind = as.vector(E(graph)[get.edge.ids(graph, pth2)])
    pwts[rowind,target] = pwts[rowind,target] + 1/pathcounts[target]
  }

  if(!is.null(admixedges)) pwts = pwts[-admixedges,]
  #pwts[,-1] - pwts[,1]
  pwts
}


tree_to_pwts = function(graph, leaves) {
  # should give same output as graph_to_pwts for trees, but faster
  # returns numedges*(numpops-1) matrix
  # in tree, each edge maps to one node

  pwts = matrix(0, length(E(graph))+1, length(leaves))
  colnames(pwts) = leaves
  rownames(pwts) = names(V(graph))
  enames = attr(E(graph), 'vnames') %>% str_replace('.+\\|', '')

  paths = shortest_paths(graph, 'R', leaves)$vpath
  indmat = paths %>% map(as.numeric) %>% map(tail, -1) %>% imap(~cbind(.x, .y)) %>% do.call(rbind, .)
  pwts[indmat] = 1
  pwts = pwts[,-1] - pwts[,1]
  pwts[enames,]
}



expand_path = function(path) {
  # igraph path (sequence of vertices)
  # duplicates inner vertices, so that this works with igraph functions that process vertex sequences pairwise
  ln = length(path)
  path[c(1, 1+rep(seq_len(ln-2), each=2), ln)]
}


graph_to_weightind = function(graph) {
  # input igraph object
  # output:
  # map (leaf, edge) -> paths
  # map path -> weights
  # ultimately: indices for weights into paths, indices for paths into pwts

  # room for improvement here...

  leaves = get_leaves(graph)
  root = get_root(graph)
  admixnodes = which(degree(graph, mode='in') == 2)
  admixedges = unlist(incident_edges(graph, admixnodes, mode='in'))
  normedges = setdiff(1:length(E(graph)), admixedges)
  paths = all_simple_paths(graph, root, leaves, mode='out')
  ends = sapply(paths, tail, 1)
  edge_per_path = paths %>% map(expand_path) %>% map(~get.edge.ids(graph, .))
  weight_per_path = edge_per_path %>% map(~(which(admixedges %in% .)))

  path_edge_table = do.call(rbind, lapply(seq_len(length(weight_per_path)),
                                          function(i) tibble(path=i, edge=c(edge_per_path[[i]])))) %>%
    mutate(edge2 = match(edge, normedges)) %>% filter(!is.na(edge2)) %>%
    mutate(leaf = as.vector(ends[path]), leaf2 = match(leaf, leaves)) %>%
    left_join(enframe(c(table(ends))) %>% transmute(leaf=as.numeric(name), numpaths=value), by='leaf') %>%
    group_by(leaf2, edge2) %>% mutate(cnt = n(), keep = cnt < numpaths) %>% filter(keep) %>% as.matrix

  path_admixedge_table = do.call(rbind, lapply(seq_len(length(weight_per_path)),
                                               function(i) tibble(path=i, admixedge=c(weight_per_path[[i]])))) %>%
    as.matrix
  list(path_edge_table, path_admixedge_table, length(paths))
}



fill_pwts = function(pwts, weights, path_edge_table, path_admixedge_table, numpaths = NULL) {
  # puts weights onto pwts, using index matrix and vectors

  if(length(weights)==0) return(pwts)
  wts2 = rep(weights, each=2)*c(1,-1) + (0:1)
  path_weights = path_admixedge_table %>% as_tibble %>% mutate(w = wts2[admixedge]) %>%
    group_by(path) %>% summarize(w = prod(w))
  pwts_weights = path_edge_table %>% as_tibble %>% left_join(path_weights, by='path') %>%
    mutate(w = ifelse(is.na(w), 1, w)) %>% group_by(edge2, leaf2) %>% summarize(w = sum(w)) %>% ungroup %>% as.matrix
  pwts[pwts_weights[,1:2]] = pwts_weights[,3]
  pwts
}


optimweightsfun = function(weights, args) {
  # likelihood function used in optimizing admixture weights
  # weights is vector of admixture weights to be optmized; only values for first incoming edge; 2nd is 1 - first

  pwts = args[[1]]
  ppinv = args[[2]]
  f3_est = args[[3]]
  path_edge_table = args[[4]] # indices into pwts with weight positions
  path_admixedge_table = args[[5]]
  #index3 = args[[6]]
  cmb = args[[7]]
  qpsolve = args[[8]]
  lower = args[[9]]
  upper = args[[10]]
  fudge = args[[11]]

  pwts = fill_pwts(pwts, weights, path_edge_table, path_admixedge_table)
  pwts = pwts[,-1] - pwts[,1]
  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])
  branch_lengths = opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve, lower, upper, fudge = fudge)
  f3_fit = ppwts_2d %*% branch_lengths
  get_score(f3_fit, f3_est, ppinv)
}



opt_edge_lengths = function(ppwts_2d, ppinv, f3_est, qpsolve, lower, upper, fudge = 1e-4) {
  # finds optimal edge lengths
  # pwts2d: npair x nedge design matrix with paths to outpop
  # ppinv: inverse of npair x npair matrix of varianc-covariance matrix of jackknife f3 stats
  # f3_est: estimated f3 stats

  pppp = t(ppwts_2d) %*% ppinv
  cc = pppp %*% ppwts_2d
  #diag(cc) = diag(cc) + fudge
  diag(cc) = diag(cc) + fudge*mean(diag(cc))
  cc = (cc+t(cc))/2
  q1 = (pppp %*% f3_est)[,1]
  nc = ncol(cc)
  #tryCatch({
  qpsolve(cc, q1, cbind(diag(nc), -diag(nc)), c(lower, -upper))
  #}, error = function(e) browser())
}

get_score = function(f3_fit, f3_est, ppinv) {
  res = f3_fit - f3_est
  lik = t(res) %*% ppinv %*% res
  lik[1,1]
}


#' Compute the fit of an admixture graph
#'
#' Computes the fit of an admixturegraph for a given graph topology and empirical f2-block-jackknife statistics.
#' @export
#' @param f2_blocks A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' @param graph An admixture graph represented as a matrix of edges, an \code{\link{igraph}} object, or the path to a qpGraph graph file. Edges can be constrained by providing a matrix or data frame of edges with columns titled `lower` and `upper` with lower and upper bounds, respectively. By default, admixture edges are constrained to be between zero and one (with paired edges summing to one), and drift edges have a lower bound at zero.
#' @param f2_denom Scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix of
#' f3 statistics will be computed using block-jackknife. Otherwise bootstrap resampling is performed `n` times,
#' where `n` is either equal to `boot` if it is an integer, or equal to the number of blocks if `boot` is `TRUE`.
#' The covariance matrix of f3 statistics will be computed using bootstrapping.
#' @param fudge Regularization term added to the diagonal elements of the covariance matrix of fitted branch lengths (after scaling by the matrix trace).
#' @param fudge_cov Regularization term added to the diagonal elements of the covariance matrix of estimated f3 statistics (after scaling by the matrix trace).
#' @param fnscale Optimization parameter passed to `control` in \code{\link{optim}}
#' @param lsqmode Least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @param numstart Number of random initializations. Defaults to 10.
#' @param seed Seed for generating starting weights.
#' @param cpp Use C++ functions. Setting this to `FALSE` will be slower but can help with debugging.
#' @param return_f4 Return all f4-statistics. Defaults to `FALSE` because this can be slow.
#' @param f3precomp Precomputed f3-statistics. This should be the output of `qpgraph_precompute_f3` and can be provided instead of `f2_blocks`. This can speed things up if many graphs are evaluated using the same set of f3-statistics.
#' @param f2_blocks_test An optional 3d array of f2-statistics used for computing an out-of-sample score. Ideally this contains SNP blocks which are not part of `f2_blocks`. This allows to estimate the fit of a graph without overfitting and will not be used during the optimization step
#' @param low_q Reported lower quantile of fitted admixture edge weights across random initializations
#' @param high_q Reported upper quantile of fitted admixture edge weights across random initializations
#' @param verbose Print progress updates
#' @return A list with data describing the model fit:
#' \enumerate{
#' \item `edges` a data frame where each row is an edge in the graph. For regular edges,
#' the column `weight` is the estimated edge length, and for admixture edges, it is the estimated admixture weight.
#' \item `score` the likelihood score of the fitted graph. lower values correspond to better fits.
#' the score is calculated as the inner product of the residuals (difference between estimated and
#' fitted f3 statistics), weighted by the inverse of the f3 covariance matrix.
#' \item `f2` estimated and fitted f2 statistics
#' \item `f3` estimated and fitted f3 statistics
#' \item `f4` estimated and fitted f4 statistics (if `return_f4 = TRUE`)
#' \item `opt` a data frame with details of the weight-fitting step, including the randomly sampled starting weights.
#' }
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @seealso \code{\link{qpgraph_wrapper}} for a wrapper functions which call the original qpGraph program.
#' @examples
#' out = qpgraph(example_f2_blocks, example_graph)
#' plot_graph(out$edges)
qpgraph = function(f2_blocks, graph, f2_denom = 1, boot = FALSE, fudge = 1e-4, fudge_cov = 1e-5, fnscale = 1e-6, lsqmode = FALSE,
                   numstart = 10, seed = NULL, cpp = TRUE, return_f4 = FALSE, f3precomp = NULL, f2_blocks_test = NULL,
                   low_q = 0, high_q = 1, verbose = FALSE) {

  #----------------- process graph -----------------
  if('matrix' %in% class(graph)) {
    edges = as.data.frame(graph, stringsAsFactors = FALSE)
  } else if('character' %in% class(graph)) {
    edges = parse_qpgraph_graphfile(graph)
  } else if('igraph' %in% class(graph)) {
    edges = igraph::as_edgelist(graph) %>% as.data.frame(stringsAsFactors = FALSE)
  } else if('data.frame' %in% class(graph)) {
    edges = graph
  } else stop(paste0('Cannot parse graph of class ', class(graph),'!'))

  if(cpp) {
    optimweightsfun = cpp_optimweightsfun
    opt_edge_lengths = cpp_opt_edge_lengths
    fill_pwts = cpp_fill_pwts
    get_pairindex = cpp_get_pairindex
  }

  if(class(graph)[1] != 'igraph') graph = graph_from_edgelist(as.matrix(edges[,1:2]))
  nedges = length(E(graph))
  admixnodes = which(degree(graph, mode='in') == 2)
  nadmix = length(admixnodes)
  admixedgesfull = sapply(seq_len(nadmix), function(i) incident_edges(graph, admixnodes, mode='in')[[i]][1:2])
  normedges = setdiff(1:nedges, admixedgesfull)

  pops = get_leafnames(graph)
  npop = length(pops)
  cmb = combn(0:(npop-1), 2)+(1:0)

  if(!is.null(f3precomp)) {
    precomp = f3precomp
    f3pops = attr(precomp$f3_est, 'pops')
    pairmatch = get_pairindex(match(pops, f3pops))
    precomp$f3_est = precomp$f3_est[pairmatch]
    precomp$ppinv = precomp$ppinv[pairmatch, pairmatch]
    stopifnot(all(!is.na(precomp$ppinv)))
  } else {
    precomp = qpgraph_precompute_f3(f2_blocks, pops, f2_denom = f2_denom, boot = boot,
                                    seed = seed, fudge_cov = fudge_cov, lsqmode = lsqmode)
  }

  f3_est = precomp$f3_est
  f3_est = pmax(0, f3_est)
  ppinv = precomp$ppinv

  weight = low = high = rep(NA, nedges)
  pwts = graph_to_pwts(graph, pops)
  opt = NULL

  mim = .Machine$integer.max
  if('lower' %in% names(edges)) {
    elower = replace_na(edges$lower[normedges], 0)
    eupper = replace_na(edges$upper[normedges], mim)
  } else {
    elower = rep(0, length(normedges))
    eupper = rep(mim, length(normedges))
  }

  if(nadmix > 0) {
    #if(is.null(numstart)) numstart = 10*nadmix
    if(!is.null(seed)) set.seed(seed)
    if('lower' %in% names(edges)) {
      alower = replace_na(pmax(edges$lower[admixedgesfull[1,]], 1-edges$upper[admixedgesfull[2,]]), 0)
      aupper = replace_na(pmin(edges$upper[admixedgesfull[1,]], 1-edges$lower[admixedgesfull[2,]]), 1)
      aupper = pmin(1, aupper) + 1e-9
    } else {
      alower = rep(0, nadmix)
      aupper = rep(1, nadmix)
    }
    parmat = matrix(runif(numstart*nadmix), numstart)
    if(verbose) alert_info(paste0('testing ', nrow(parmat), ' combinations of admixture weight starting values\n'))
    weightind = graph_to_weightind(graph)
    arglist = list(pwts, ppinv, f3_est, weightind[[1]], weightind[[2]], weightind[[3]], cmb, qpsolve, elower, eupper, fudge)
    oo = multistart(parmat, optimweightsfun, args = arglist, method = 'L-BFGS-B',
                    lower = alower, upper = aupper, control=list(maxit = 1e4, fnscale = fnscale),
                    verbose = verbose)
    best = oo %>% top_n(1, -value)
    opt = data.frame(parmat, oo, stringsAsFactors = F) #%>% filter(!is.na(convergence))

    #if(nrow(opt) == 0) stop('Optimization not successful! Increase fudge or numstart!')
    #else if(verbose) alert_info(paste0('Optimiztion successful for ', nrow(opt), ' combinations\n'))

    admnames = names(V(graph))[admixnodes]
    colnames(opt)[1:(nadmix*2)] = paste0(rep(c('i.', 'e.'), each = nadmix), rep(admnames, 2))
    hilo = apply(as.matrix(oo[,1:nadmix]), 2, function(x) quantile(x, c(low_q, high_q), na.rm = TRUE))

    wts = as.matrix(best[,1:nadmix])[1,]
    weight[admixedgesfull[1,]] = wts
    weight[admixedgesfull[2,]] = 1-wts
    low[admixedgesfull[1,]] = pmin(hilo[1,], hilo[2,])
    low[admixedgesfull[2,]] = pmin(1-hilo[1,], 1-hilo[2,])
    high[admixedgesfull[1,]] = pmax(hilo[1,], hilo[2,])
    high[admixedgesfull[2,]] = pmax(1-hilo[1,], 1-hilo[2,])
    pwts = fill_pwts(pwts, wts, weightind[[1]], weightind[[2]], weightind[[3]])
  }
  pwts = pwts[,-1] - pwts[,1]
  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])
  branch_lengths = opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve, elower, eupper, fudge = fudge)
  f3_fit = ppwts_2d %*% branch_lengths
  score = get_score(f3_fit, f3_est, ppinv)
  if(!is.null(f2_blocks_test)) {
    precomp_test = qpgraph_precompute_f3(f2_blocks_test, pops, f2_denom = f2_denom, boot = boot,
                                         seed = seed, fudge_cov = fudge_cov, lsqmode = lsqmode)
    score_test = get_score(f3_fit, precomp_test$f3_est, ppinv)
  } else {
    score_test = NULL
  }

  #weight[normedges] = branch_lengths
  weight[normedges] = pmax(0, branch_lengths)
  edges %<>% select(1:2) %>% set_colnames(c('from', 'to')) %>%  as_tibble %>%
    mutate(type = ifelse(1:n() %in% normedges, 'edge', 'admix'), weight = weight, low = low, high = high)
  f2 = precomp$f2out
  f3 = precomp$f3out %>% mutate(fit = c(f3_fit), diff = fit - est, z = diff/se, p = ztop(z))

  out = namedList(edges, score, score_test, f2, f3, opt, ppinv)
  #if(return_f4) out$f4 = f4(f2_blocks)
  if(return_f4) {
    if(verbose) alert_info(paste0('Computing f4\n'))
    out$f4 = fitf4(f2_blocks[pops, pops, ], f2, f3, cmb)
  }
  out
}



#' Compute f3-statistics from f2-statistics.
#'
#' Takes a 3d array of f2 block jackknife estimates and computes f3-statistics between the
#' first population \eqn{p1} and all population pairs \eqn{i, j}: \eqn{f3(p1; p_i, p_j)}
#' @export
#' @param f2_data A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' alternatively, a directory with precomputed data. see \code{\link{extract_f2}} and \code{\link{extract_indpairs}}.
#' @param pops Populations for which to compute f3-statistics
#' @param outpop Outgroup population. used as the basis of the f3-statistics. If `NULL` (the default),
#' the first population in `pops` will be used as the basis.
#' @param f2_denom Scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix of f3 statistics
#' will be computed using block-jackknife. Otherwise bootstrap resampling is performed `n` times, where `n` is either
#' equal to `boot` if it is an integer, or equal to the number of blocks if `boot` is `TRUE`. The the covariance matrix
#' of f3 statistics will be computed using bootstrapping.
#' @param fudge_cov Regularization term added to the diagonal elements of the covariance matrix of estimated f3 statistics (after scaling by the matrix trace).
#' @param lsqmode Least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @return A list with four items
#' \enumerate{
#' \item `f3_est` a matrix with f3-statistics for all population pairs with the output
#' \item `ppinv` a matrix with the inverse of the f3-statistic covariance matrix
#' \item `f2out` a data frame with f2 estimates
#' \item `f3out` a data frame with f3 estimates
#' }
#' @examples
#' pops = get_leafnames(example_igraph)
#' qpgraph_precompute_f3(example_f2_blocks, pops)$f3_est
#' \dontrun{
#' qpgraph_precompute_f3(f2_dir, pops, f2_denom = 0.278)
#' }
qpgraph_precompute_f3 = function(f2_data, pops, outpop = NULL, f2_denom = 1, boot = FALSE,
                                 seed = NULL, fudge_cov = 1e-5, lsqmode = FALSE) {
  # returns list of f3_est and ppinv for subset of populations.
  # f3_est and ppinv are required for qpgraph_slim; f2out and f3out are extra output
  # f2_blocks may contain more populations than the ones used in qpgraph
  # f2_blocks input here should be subset which is used by qpgraph function

  #----------------- read f-stats -----------------
  if(!is.null(outpop)) pops = c(outpop, setdiff(pops, outpop))

  if(!is.null(seed)) set.seed(seed)
  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  matstatfun = ifelse(boot, boot_mat_stats, jack_mat_stats)
  arrstatfun = ifelse(boot, boot_arr_stats, jack_arr_stats)
  f2_blocks = get_f2(f2_data, pops, f2_denom) %>% samplefun
  #f2_blocks = array(pmax(0, f2_blocks), dim(f2_blocks), dimnames(f2_blocks))
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  npop = length(pops)
  npair = choose(npop, 2)
  cmb = combn(0:(npop-1), 2)+(1:0)

  f2 = arrstatfun(f2_blocks, block_lengths)
  f2out = tibble(pop1 = combn(pops, 2)[1,],
                 pop2 = combn(pops, 2)[2,],
                 est = f2[[1]][lower.tri(f2[[1]])],
                 se = sqrt(f2[[2]][lower.tri(f2[[2]])]))

  f3_blocks = (f2_blocks[,rep(1, npop),] + f2_blocks[rep(1, npop),,] - f2_blocks)/2
  #f3_blocks = array(pmax(0, f3_blocks), dim(f3_blocks), dimnames(f3_blocks))
  f3_blocks_2d = arr3d_to_pairmat(f3_blocks[-1,-1,])
  f3dat = matstatfun(f3_blocks_2d, block_lengths)
  #f3dat = jack_mat_stats(f3_blocks_2d, block_lengths)
  f3_est = f3dat$est
  f3_var = f3dat$var
  f3out = tibble(pop1 = pops[1],
                 pop2 = pops[cmb[1,]+1],
                 pop3 = pops[cmb[2,]+1],
                 est = f3_est, se = sqrt(diag(f3_var)))
  diag(f3_var) = diag(f3_var) + sum(diag(f3_var))*fudge_cov
  # in qpGraph fudge_cov is 1e-5; sometimes quadprog doesn't converge unless this is larger
  #   has large effect on magnitude of likelihood score
  if(lsqmode) ppinv = diag(1/diag(f3_var))
  else ppinv = solve(f3_var)

  f3_est %<>% structure(pops = pops)
  ppinv %<>% structure(pops = pops)
  namedList(f3_est, ppinv, f2out, f3out, f3_blocks_2d)
}

get_pairindex = function(perm) {
  # returns index vector that matches population pairs in
  # f3_est and ppinv (which were computed using pops) to graph populations
  cmb = combn(0:(length(perm)-1), 2)+(1:0)
  popind = setdiff(perm, 1)
  orig_order = apply(cmb+1, 2, paste0, collapse='')
  new_order = apply(matrix(popind[c(cmb)], 2), 2, function(x) paste0(sort(x), collapse=''))
  match(new_order, orig_order)
}

#' @export
qpgraph_anorexic = function(f3precomp, graph, fudge = 1e-4, fnscale = 1e-6,
                            numstart = 10, seed = NULL, verbose = FALSE, cpp = TRUE) {

  # only works for trees at the moment, because weightind order is coupled to pwts order
  admixnodes = which(degree(graph, mode='in') == 2)
  stopifnot(length(admixnodes) == 0)

  if(cpp) {
    opt_edge_lengths = cpp_opt_edge_lengths
  }

  graphpops = get_leafnames(graph)
  f3_est = f3precomp$f3_est
  ppinv = f3precomp$ppinv
  f3pops = attr(f3_est, 'pops')
  pwts = tree_to_pwts(graph, graphpops)
  pwts = pwts[, match(f3pops[-1], setdiff(graphpops, f3pops[1]))]

  cmb = combn(0:(length(f3pops)-1), 2)+(1:0)
  ppwts_2d = t(pwts[,cmb[1,]] * pwts[,cmb[2,]])

  branch_lengths = opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve,
                                    lower = rep(0, nrow(pwts)), upper = rep(.Machine$integer.max, nrow(pwts)),
                                    fudge = fudge)
  f3_fit = ppwts_2d %*% branch_lengths
  score = get_score(f3_fit, f3_est, ppinv)

  edges = as_tibble(as_edgelist(graph), .name_repair = ~c('from', 'to')) %>%
    mutate(type = 'edge', weight = c(branch_lengths))

  namedList(edges, score, opt = NULL)
}

# not used
f3out_to_fittedf2out = function(f2out, f3out) {
  # computes fitted f2 statistics data frame from f2 and f3 data frames
  # will not include f2(outgroup, X)

  f2out %>%
    right_join(f3out %>% filter(pop2 != pop3) %>% transmute(pop1=pop2, pop2=pop3, f3 = fit), by = c('pop1', 'pop2')) %>%
    left_join(f3out %>% filter(pop2 == pop3) %>% transmute(pop1=pop2, f2_1 = fit), by = c('pop1')) %>%
    left_join(f3out %>% filter(pop2 == pop3) %>% transmute(pop2=pop3, f2_2 = fit), by = c('pop2')) %>%
    transmute(pop1, pop2, est, se, fit = (f2_1 + f2_2 - f3*2), diff = fit - est, z = diff/se, p = ztop(z))
}




fitf4 = function(f2_blocks, f2, f3, cmb) {
  # returns a tibble with estimated and fitted f4-statistics

  f2_out = f3 %>% filter(pop2 == pop3) %$% fit
  f2_fit = f3 %>% mutate(f21 = f2_out[cmb[1,]], f22 = f2_out[cmb[2,]], f2fit = (f21 + f22 - fit*2))
  f2_fit2 = f2 %>%
    left_join(f2_fit, by = c('pop1'='pop2', 'pop2'='pop3')) %>%
    filter(!is.na(f2fit)) %>%
    select(pop1, pop2, f2fit) %>%
    bind_rows(f2_fit %>% filter(pop2 == pop3) %>% transmute(pop1, pop2, f2fit = fit)) %>%
    bind_rows(rename(., pop1 = pop2, pop2 = pop1)) %>%
    bind_rows(tibble(pop1 = unique(.$pop1), pop2 = pop1, f2fit = 0))
  x = f4(f2_blocks, unique_only = F) %>% select(-z, -p)
  x %>%
    left_join(f2_fit2 %>% rename(c1 = f2fit), by = c('pop1' = 'pop1', 'pop4' = 'pop2')) %>%
    left_join(f2_fit2 %>% rename(c2 = f2fit), by = c('pop2' = 'pop1', 'pop3' = 'pop2')) %>%
    left_join(f2_fit2 %>% rename(c3 = f2fit), by = c('pop1' = 'pop1', 'pop3' = 'pop2')) %>%
    left_join(f2_fit2 %>% rename(c4 = f2fit), by = c('pop2' = 'pop1', 'pop4' = 'pop2')) %>%
    mutate(fit = (c1 + c2 - c3 - c4)/2, diff = fit - est, z = diff/se, p = ztop(z)) %>%
    select(-c1:-c4)
}



#' Compare the fit of two qpgraph models
#'
#' Takes two data frames with model fits computed on two graphs for on the same populations and tests whether the scores of one graph are significantly better than the scores of the other
#' @export
#' @param fits1 The fits of the first graph
#' @param fits2 The fits of the second graph
#' @param boot should match the `boot` parameter in `qpgraph_resample_snps` (`FALSE` by default)
#' @examples
#' \dontrun{
#' nblocks = dim(example_f2_blocks)[3]
#' train = sample(1:nblocks, round(nblocks/2))
#' fits1 = qpgraph_resample_snps(example_f2_blocks[,,train], graph = graph1, f2_blocks_test = example_f2_blocks[,,-train])
#' fits2 = qpgraph_resample_snps(example_f2_blocks[,,train], graph = graph2, f2_blocks_test = example_f2_blocks[,,-train])
#' compare_fits2(fit1, fit2)
#' }
compare_fits2 = function(fits1, fits2, boot = FALSE) {

  matstatfun = ifelse(boot, boot_mat_stats, jack_mat_stats)
  stats = matstatfun(t(fits1$score_test - fits2$score_test), rep(1, length(fits1$score_test)))

  diff = stats$est
  se = sqrt(stats$var)
  z = diff/se
  p = ztop(z)
  namedList(diff, se, z, p, scores1 = fits1$score_test, scores2 = fits2$score_test, boot)

}

#' Compare the fit of two qpgraph models
#'
#' Takes the bootstrap score distribution of two fits on the same populations and tests whether the scores of one graph are significantly better than the scores of the other.
#' @export
#' @param fits1 Fitted models for the first graph
#' @param fits2 Fitted models for the second graph
#' @examples
#' \dontrun{
#' boo = boo_list(f2_blocks, nboot = 100)
#' fits1 = qpgraph_resample_snps2(boo$boo, graph1, boo$test)
#' fits2 = qpgraph_resample_snps2(boo$boo, graph2, boo$test)
#' compare_fits3(fits1$score_test, fits2$score_test)
#' }
#' \dontrun{
#' boo = boo_list(f2_blocks, nboot = 100)
#' f3precomp1 = qpgraph_precompute_f3(f2_blocks, get_leafnames(graph1))
#' f3precomp2 = qpgraph_precompute_f3(f2_blocks, get_leafnames(graph2))
#' fits1 = qpgraph_resample_snps2(boo$boo, graph1, boo$test, f3precomp = f3precomp1)
#' fits2 = qpgraph_resample_snps2(boo$boo, graph2, boo$test, f3precomp = f3precomp2)
#' compare_fits3(fits1$score_test, fits2$score_test)
#' }
compare_fits3 = function(scores1, scores2) {

  scorediff = scores1 - scores2
  ci_low = unname(quantile(scorediff, 0.025, na.rm = T))
  ci_high = unname(quantile(scorediff, 0.975, na.rm = T))

  scorediff = na.omit(scores1 - scores2)
  stats = boot_mat_stats(t(scorediff), rep(1, length(scorediff)))

  diff = stats$est
  se = sqrt(stats$var)
  z = diff/se
  p = ztop(z)
  frac = mean(scorediff < 0)
  p_emp = min(frac, 1-frac)*2
  p_emp = max(p_emp, 1/length(scorediff))
  c(diff=diff, se=se, z=z, p=p, p_emp=p_emp, ci_low=ci_low, ci_high=ci_high)
}

#' @export
qpgraph_resample_snps2 = function(f2_blocks, graph, f2_blocks_test, verbose = TRUE, ...) {

  ell = list(...)
  fun = function(f2dat, f2dat_test, g) function() safely(qpgraph)(f2_blocks = f2dat, graph = g, f2_blocks_test = f2dat_test, verbose = FALSE, ...)

  tibble(id = seq_len(length(f2_blocks)), graph = list(graph), f2_blocks, f2_blocks_test) %>%
    mutate(fun2 = pmap(list(f2_blocks, f2_blocks_test, graph), fun)) %>%
    mutate(out = furrr::future_invoke_map(fun2, .progress = verbose),
           result = map(out, 'result', .null = tibble()), error = map(out, 'error')) %>%
    select(-out, -fun2) %>% unnest_wider(result)
}



#' Compare the fit of two qpgraph models
#'
#' Takes two data frames with model fits computed on two graphs for on the same populations and tests whether the scores of one graph are significantly better than the scores of the other.
#' @export
#' @param fit1 The fit of the first graph
#' @param fit2 The fit of the second graph
#' @param f2_blocks f2 blocks used for fitting `fit1` and `fit2`. Used in combination with `f2_blocks_test` to compute f-statistics covariance matrix.
#' @param f2_blocks_test f2 blocks which were not used for fitting `fit1` and `fit2`
#' @param boot if `TRUE`, bootstrap resampling will be used on `f2_blocks_test` with the number of resamplings equal to the number of blocks. If `FALSE` jackknife will be used. If set to a number, bootstrap resampling will be used on `f2_blocks_test` with the number of resamplings equal to `boot`. If bootstrap resampling is enabled, empirical p-values (`p_emp`) and 95\% confidence intervals (`ci_low` and `ci_high`) will be reported.
#' @param seed random seed used if `boot` is `TRUE`. does not need to match a seed used in fitting the models
#' @examples
#' \dontrun{
#' nblocks = dim(example_f2_blocks)[3]
#' train = sample(1:nblocks, round(nblocks/2))
#' fit1 = qpgraph(example_f2_blocks[,,train], graph1)
#' fit2 = qpgraph(example_f2_blocks[,,train], graph2)
#' compare_fits4(fit1, fit2, example_f2_blocks[,,train], example_f2_blocks[,,-train])
#' }
compare_fits4 = function(fit1, fit2, f2_blocks, f2_blocks_test, boot = FALSE, seed = NULL) {

  stopifnot(all.equal(sort(attr(fit1$ppinv, 'pops')), sort(attr(fit2$ppinv, 'pops'))))
  matstatfun = ifelse(boot, boot_mat_stats, jack_mat_stats)

  pops = attr(fit1$ppinv, 'pops')
  ppinv = qpgraph_precompute_f3(abind::abind(f2_blocks, f2_blocks_test), pops, boot = boot, seed = seed)$ppinv
  f3_test = qpgraph_precompute_f3(f2_blocks_test, pops, boot = boot, seed = seed)$f3_blocks_2d
  f3_fit = fit1$f3 %>%
    left_join(fit2$f3 %>% bind_rows(rename(., pop2=pop3, pop3=pop2) %>% filter(pop2 != pop3)),
              by = c('pop1', 'pop2', 'pop3'))
  scores1 = map_dbl(1:dim(f2_blocks_test)[3], ~get_score(f3_fit$fit.x, f3_test[,.], ppinv))
  scores2 = map_dbl(1:dim(f2_blocks_test)[3], ~get_score(f3_fit$fit.y, f3_test[,.], ppinv))

  scorediff = na.omit(scores1 - scores2)
  stats = matstatfun(t(scorediff), rep(1, length(scorediff)))

  diff = stats$est
  se = sqrt(stats$var)
  z = diff/se
  p = ztop(z)
  frac = mean(scorediff < 0)
  p_emp = ci_low = ci_high = NA
  if(boot) {
    p_emp = min(frac, 1-frac)*2
    ci_low = unname(quantile(scorediff, 0.025, na.rm = T))
    ci_high = unname(quantile(scorediff, 0.975, na.rm = T))
  }
  namedList(diff, se, z, p, p_emp, ci_low, ci_high, scores1, scores2)
}



