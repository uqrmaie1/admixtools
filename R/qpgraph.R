

graph_to_pwts = function(graph) {
  # input igraph object
  # output: numedges*(numpops-1) matrix 'pwts' which indicates all paths from pop to outpop; admixture edges are mapped onto parent edges and weighted
  # assumes first pop is outpop, and 'R' is root

  leaves = get_leafnames(graph)
  #leaves = V(graph)$name[degree(graph, v = V(graph), mode = c('out')) == 0]
  pwts = matrix(0, length(E(graph)), length(leaves))
  colnames(pwts) = leaves
  rownames(pwts) = attr(E(graph), 'vnames')

  admixnodes = which(degree(graph, mode='in')==2)
  admixedges = unlist(incident_edges(graph, admixnodes, mode='in'))

  allpaths = all_simple_paths(graph, V(graph)[1], leaves, mode='out')
  pathcounts = table(names(sapply(allpaths, function(x) tail(x,1))))
  for(i in seq_len(length(allpaths))) {
    target = names(tail(allpaths[[i]],1))
    ln = length(allpaths[[i]])
    pth2 = allpaths[[i]][c(1, 1+rep(seq_len(ln-2), each=2), ln)]
    rowind = as.vector(E(graph)[get.edge.ids(graph, pth2)])
    pwts[rowind,target] = pwts[rowind,target] + 1/pathcounts[target]
  }

  if(!is.null(admixedges)) pwts = pwts[-admixedges,]
  pwts[,-1] - pwts[,1]
}

tree_to_pwts = function(graph) {
  # should give same output as graph_to_pwts for trees, but faster
  # returns numedges*(numpops-1) matrix
  # in tree, each edge maps to one node

  leaves = get_leafnames(graph)
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
  admixnodes = which(degree(graph, mode='in')==2)
  admixedges = unlist(incident_edges(graph, admixnodes, mode='in'))
  normedges = setdiff(1:length(E(graph)), admixedges)
  paths = all_simple_paths(graph, V(graph)[1], leaves[-1], mode='out')
  ends = sapply(paths, tail, 1)
  edge_per_path = paths %>% map(expand_path) %>% map(~get.edge.ids(graph, .))
  weight_per_path = edge_per_path %>% map(~(which(admixedges %in% .)))

  path_edge_table = do.call(rbind, lapply(seq_len(length(weight_per_path)),
                                          function(i) tibble(path=i, edge=c(edge_per_path[[i]])))) %>%
    mutate(edge2 = match(edge, normedges)) %>% filter(!is.na(edge2)) %>%
    mutate(leaf = as.vector(ends[path]), leaf2 = match(leaf, leaves[-1])) %>%
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

  pwts = fill_pwts(pwts, weights, path_edge_table, path_admixedge_table)
  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])
  branch_lengths = opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve, lower, upper)
  f3_fit = ppwts_2d %*% branch_lengths
  get_score(f3_fit, f3_est, ppinv)
}



opt_edge_lengths = function(ppwts_2d, ppinv, f3_est, qpsolve, lower, upper) {
  # finds optimal edge lengths
  # pwts2d: npair x nedge design matrix with paths to outpop
  # ppinv: inverse of npair x npair matrix of varianc-covariance matrix of jackknife f3 stats
  # f3_est: estimated f3 stats

  pppp = t(ppwts_2d) %*% ppinv
  cc = pppp %*% ppwts_2d
  diag(cc) = diag(cc) + mean(diag(cc))*0.0001
  cc = (cc+t(cc))/2
  q1 = (pppp %*% f3_est)[,1]
  nc = ncol(cc)
  qpsolve(cc, q1, cbind(diag(nc), -diag(nc)), c(lower, -upper))
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
#' @param f2_data a 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' alternatively, a directory with precomputed data. see \code{\link{extract_f2}} and \code{\link{extract_indpairs}}.
#' @param graph an admixture graph represented as a matrix of edges, an \code{\link{igraph}} object, or the path to a qpGraph graph file.
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix of
#' f3 statistics will be computed using block-jackknife. Otherwise bootstrap resampling is performed `n` times,
#' where `n` is either equal to `boot` if it is an integer, or equal to the number of blocks if `boot` is `TRUE`.
#' The covariance matrix of f3 statistics will be computed using bootstrapping.
#' @param fudge try increasing this, if you get the error message \code{constraints are inconsistent, no solution!}.
#' @param fnscale optimization parameter passed to `control` in \code{\link{optim}}
#' @param lsqmode least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @param numstart number of random initializations. defaults to 10 times the number of admixture nodes.
#' @param seed seed for generating starting weights.
#' @param cpp should optimization be done using C++ or R function? \code{cpp = TRUE} is much faster.
#' @param return_f4 return all f4 statistics? Can take a while.
#' @param verbose print progress updates
#' @return a list with output describing the model fit:
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
#' @seealso \code{\link{qpgraph_wrapper}} for a wrapper functions which call the original qpGraph program;
#' \code{\link{qpgraph_slim}} for a faster function function which requires f3 estimates
#' and an inverted covariance matrix as input instead.
#' @examples
#' out = qpgraph(example_f2_blocks, example_graph)
#' plot_graph(out$edges)
qpgraph = function(f2_data, graph, f2_denom = 1, boot = FALSE, fudge = 1e-3, fnscale = 1e-6, lsqmode = FALSE,
                   numstart = NULL, seed = NULL, cpp = TRUE, return_f4 = FALSE, f3precomp = NULL,
                   low_q = 0, high_q = 1, verbose = FALSE) {
  # modelled after AdmixTools qpGraph

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

  if(class(graph)[1] != 'igraph') graph = graph_from_edgelist(as.matrix(edges[,1:2]))
  nedges = length(E(graph))
  admixnodes = which(degree(graph, mode='in')==2)
  nadmix = length(admixnodes)
  admixedgesfull = sapply(seq_len(nadmix), function(i) incident_edges(graph, admixnodes, mode='in')[[i]][1:2])
  normedges = setdiff(1:nedges, admixedgesfull)

  pops = get_leafnames(graph)
  npop = length(pops)
  cmb = combn(0:(npop-1), 2)+(1:0)

  if(!is.null(f3precomp)) precomp = f3precomp
  else precomp = qpgraph_precompute_f3(f2_data, pops, f2_denom = f2_denom, boot = boot, fudge = fudge, lsqmode = lsqmode)
  f3_est = precomp$f3_est
  ppinv = precomp$ppinv

  weightind = graph_to_weightind(graph)
  weight = low = high = rep(NA, nedges)
  pwts = graph_to_pwts(graph)
  opt = NULL

  if(cpp) {
    optimweightsfun = cpp_optimweightsfun
    opt_edge_lengths = cpp_opt_edge_lengths
    fill_pwts = cpp_fill_pwts
  }

  mim = .Machine$integer.max
  if('lower' %in% names(edges)) {
    elower = replace_na(edges$lower[normedges], 0)
    alower = replace_na(pmax(edges$lower[admixedgesfull[1,]], 1-edges$upper[admixedgesfull[2,]]), 0)
  } else {
    elower = rep(0, length(normedges))
    alower = rep(0, nadmix)
  }
  if('upper' %in% names(edges)) {
    eupper = replace_na(edges$upper[normedges], mim)
    aupper = replace_na(pmin(edges$upper[admixedgesfull[1,]], 1-edges$lower[admixedgesfull[2,]]), 1)
    aupper = pmin(1, aupper)
  } else {
    eupper = rep(mim, length(normedges))
    aupper = rep(1, nadmix)
  }

  if(nadmix > 0) {
    if(is.null(numstart)) numstart = 10*nadmix
    set.seed(seed)
    parmat = matrix(runif(numstart*nadmix), numstart)
    if(verbose) alert_info(paste0('testing ', nrow(parmat), ' combinations of admixture weight starting values\n'))
    arglist = list(pwts, ppinv, f3_est, weightind[[1]], weightind[[2]], weightind[[3]], cmb, qpsolve, elower, eupper)
    oo = multistart(parmat, optimweightsfun, args=arglist, method='L-BFGS-B',
                    lower=alower, upper=aupper, control=list(maxit=1e4, fnscale=fnscale), verbose=verbose)

    best = oo %>% top_n(1, -.data$value)
    admnames = names(V(graph))[admixnodes]
    opt = data.frame(parmat, oo, stringsAsFactors = F)
    colnames(opt)[1:(nadmix*2)] = paste0(rep(c('i.', 'e.'), each = nadmix), rep(admnames, 2))
    hilo = apply(as.matrix(oo[,1:nadmix]), 2, function(x) quantile(x, c(low_q, high_q)))

    wts = as.matrix(best[,1:nadmix])[1,]
    weight[admixedgesfull[1,]] = wts
    weight[admixedgesfull[2,]] = 1-wts
    low[admixedgesfull[1,]] = pmin(hilo[1,], hilo[2,])
    low[admixedgesfull[2,]] = pmin(1-hilo[1,], 1-hilo[2,])
    high[admixedgesfull[1,]] = pmax(hilo[1,], hilo[2,])
    high[admixedgesfull[2,]] = pmax(1-hilo[1,], 1-hilo[2,])
    pwts = fill_pwts(pwts, wts, weightind[[1]], weightind[[2]], weightind[[3]])
  }

  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])
  branch_lengths = opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve, elower, eupper)
  f3_fit = ppwts_2d %*% branch_lengths
  score = get_score(f3_fit, f3_est, ppinv)

  weight[normedges] = branch_lengths
  edges %<>% select(1:2) %>% set_colnames(c('from', 'to')) %>%  as_tibble %>%
    mutate(type = ifelse(1:n() %in% normedges, 'edge', 'admix'), weight = weight, low = low, high = high)
  f2 = precomp$f2out
  f3 = precomp$f3out %>% mutate(fit = c(f3_fit), diff = fit - est, z = diff/se, p.value = ztop(z))

  out = namedList(edges, score, f2, f3, opt)
  if(return_f4) out$f4 = f4(f2_data)
  out
}

#' Compute the fit of an admixturegraph.
#'
#' Computes the fit of an admixture graph for a given graph topology given f3 estimates and an inverted covariance matrix.
#' @export
#' @param graph an admixture graph represented as an \code{\link{igraph}} object
#' @param f3_est matrix of f3 estimates. output of \code{\link{qpgraph_precompute_f3}}
#' @param ppinv inverted covariance matrix of f3-statistics
#' @param fnscale try increasing or decreasing this, if optimization does not converge.
#' @param numstart number of random initializations. defaults to 10 times the number of admixture nodes.
#' @param seed seed for generating starting weights.
#' @param cpp should optimization be run in C++ or in R? C++ is faster.
#' @param verbose print progress updates
#' @return a list with output describing the model fit:
#' \enumerate{
#' \item `edges` a data frame where each row is an edge in the graph. For regular edges, the column `weight` is the
#' estimated edge length, and for admixture edges, it is the estimated admixture weight.
#' \item `score` the likelihood score of the fitted graph. lower values correspond to better fits.
#' the score is calculated as the inner product of the residuals (difference between estimated and fitted f3 statistics),
#' weighted by the inverse of the f3 covariance matrix.
#' \item `opt` a data frame with details of the weight-fitting step, including the randomly sampled starting weights.
#' }
#' @seealso \code{\link{qpgraph_wrapper}}for a wrapper functions which call the original qpGraph program;
#' \code{\link{qpgraph}} for a slower function which requires f2 block jackknife statistics as input instead.
#' \code{\link{qpgraph_precompute_f3}} computes the required input (`f3_est`, `ppinv`) from a 3d array of f2 statistics.
#' @examples
#' pops = get_leafnames(example_igraph)
#' precomp = qpgraph_precompute_f3(example_f2_blocks, pops)
#' f3_est = precomp$f3_est
#' ppinv = precomp$ppinv
#' out = qpgraph_slim(example_igraph, f3_est, ppinv)
#' plot_graph(out$edges)
qpgraph_slim = function(graph, f3_est, ppinv, fnscale = 1e-6, numstart = 10,
                        seed = NULL, cpp = TRUE, verbose = FALSE) {
  # modelled after AdmixTools qpGraph
  # optimised for testing many topologies for a given set of populations
  # graph is igraph

  if(cpp) {
    get_pairindex = cpp_get_pairindex
    optimweightsfun = cpp_optimweightsfun
    opt_edge_lengths = cpp_opt_edge_lengths
    fill_pwts = cpp_fill_pwts
  }

  nedges = length(E(graph))
  admixnodes = which(degree(graph, mode='in')==2)
  nadmix = length(admixnodes)
  admixedgesfull = sapply(seq_len(nadmix), function(i) incident_edges(graph, admixnodes, mode='in')[[i]][1:2])
  normedges = setdiff(1:nedges, admixedgesfull)
  graphpops = get_leafnames(graph)
  f3pops = attr(f3_est, 'pops')

  pairmatch = get_pairindex(match(graphpops, f3pops))
  f3_est = f3_est[pairmatch]
  ppinv = ppinv[pairmatch, pairmatch]
  stopifnot(all(!is.na(ppinv)))

  weightind = graph_to_weightind(graph)
  weight = rep(NA, nedges)
  pwts = graph_to_pwts(graph)
  opt = NULL
  cmb = combn(0:(length(f3pops)-1), 2)+(1:0)
  elower = rep(0, length(normedges))
  eupper = rep(.Machine$integer.max, length(normedges))

  if(nadmix > 0) {
    set.seed(seed)
    parmat = matrix(runif(numstart*nadmix), numstart)

    arglist = list(pwts, ppinv, f3_est, weightind[[1]], weightind[[2]], weightind[[3]], cmb, qpsolve, elower, eupper)
    opt = multistart(parmat, optimweightsfun, args=arglist, method='L-BFGS-B',
                     lower=0, upper=1, control=list(maxit=1e4, fnscale=fnscale), verbose=verbose)

    best = opt %>% top_n(1, -.data$value)
    opt = data.frame(parmat, opt, stringsAsFactors = F)

    wts = as.matrix(best[,1:nadmix])[1,]
    weight[admixedgesfull[1,]] = wts
    weight[admixedgesfull[2,]] = 1-wts
    pwts = fill_pwts(pwts, wts, weightind[[1]], weightind[[2]], weightind[[3]])
  }

  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])

  branch_lengths = opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve, elower, eupper)
  f3_fit = ppwts_2d %*% branch_lengths
  score = get_score(f3_fit, f3_est, ppinv)

  weight[normedges] = branch_lengths
  edges = as_tibble(as_edgelist(graph), .name_repair = ~c('from', 'to')) %>%
    mutate(type = ifelse(1:n() %in% normedges, 'edge', 'admix'), weight = weight)

  #namedList(edges, score, opt, f3_fit)
  namedList(edges, score, opt)
}


#' Compute f3-statistics from f2-statistics.
#'
#' Takes a 3d array of f2 block jackknife estimates and computes f3-statistics between the
#' first population \eqn{p1} and all population pairs \eqn{i, j}: \eqn{f3(p1; p_i, p_j)}
#' @export
#' @param f2_data a 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' alternatively, a directory with precomputed data. see \code{\link{extract_f2}} and \code{\link{extract_indpairs}}.
#' @param pops populations for which to compute f3-statistics
#' @param outpop outgroup population. used as the basis of the f3-statistics. If `NULL` (the default),
#' the first population in `pops` will be used as the basis.
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix of f3 statistics
#' will be computed using block-jackknife. Otherwise bootstrap resampling is performed `n` times, where `n` is either
#' equal to `boot` if it is an integer, or equal to the number of blocks if `boot` is `TRUE`. The the covariance matrix
#' of f3 statistics will be computed using bootstrapping.
#' @param fudge try increasing this, if you get the error message \code{constraints are inconsistent, no solution!}.
#' @param lsqmode least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @return a list with four items
#' \enumerate{
#' \item `f3_est` a matrix with f3-statistics for all population pairs with the output
#' \item `ppinv` a matrix with the inverse of the f3-statistic covariance matrix
#' \item `2out` a data frame with f2 estimates
#' \item `f3out` a data frame with f3 estimates
#' }
#' @examples
#' pops = get_leafnames(example_igraph)
#' qpgraph_precompute_f3(example_f2_blocks, pops)$f3_est
#' \dontrun{
#' qpgraph_precompute_f3(f2_dir, pops, f2_denom = 0.278)
#' }
qpgraph_precompute_f3 = function(f2_data, pops, outpop = NULL, f2_denom = 1, boot = FALSE,
                                 fudge = 1e-3, lsqmode = FALSE) {
  # returns list of f3_est and ppinv for subset of populations.
  # f3_est and ppinv are required for qpgraph_slim; f2out and f3out are extra output
  # f2_blocks may contain more populations than the ones used in qpgraph
  # f2_blocks input here should be subset which is used by qpgraph function

  #----------------- read f-stats -----------------
  if(!is.null(outpop)) pops = c(outpop, setdiff(pops, outpop))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo_nafix)
  matstatfun = ifelse(boot, boot_mat_stats, jack_mat_stats)
  arrstatfun = ifelse(boot, boot_arr_stats, jack_arr_stats)
  f2_blocks = get_f2(f2_data, pops, f2_denom) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])
  f2_blocks %<>% rray(dim_names = dimnames(f2_blocks))

  npop = length(pops)
  npair = choose(npop, 2)
  cmb = combn(0:(npop-1), 2)+(1:0)

  f2 = arrstatfun(f2_blocks, block_lengths)
  f2out = tibble(pop1 = combn(pops, 2)[1,],
                 pop2 = combn(pops, 2)[2,],
                 est = f2[[1]][lower.tri(f2[[1]])],
                 se = sqrt(f2[[2]][lower.tri(f2[[2]])]))

  f3_blocks = (f2_blocks[,1,] + f2_blocks[1,,] - f2_blocks)/2
  f3_blocks_2d = arr3d_to_pairmat(f3_blocks[-1,-1,])
  f3dat = matstatfun(f3_blocks_2d, block_lengths)
  f3_est = f3dat$est
  f3_var = f3dat$var
  f3out = tibble(pop1 = pops[1],
                 pop2 = pops[cmb[1,]+1],
                 pop3 = pops[cmb[2,]+1],
                 est = f3_est, se = sqrt(diag(f3_var)))
  diag(f3_var) = diag(f3_var) + sum(diag(f3_var))*fudge
  # in qpGraph fudge is 1e-5; sometimes quadprog doesn't converge unless this is larger
  #   has large effect on magnitude of likelihood score
  if(lsqmode) ppinv = diag(1/diag(f3_var))
  else ppinv = solve(f3_var)

  f3_est %<>% structure(pops = pops)
  ppinv %<>% structure(pops = pops)
  namedList(f3_est, ppinv, f2out, f3out)
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
qpgraph_anorexic = function(graph, f3_est, ppinv, fnscale = 1e-6,
                            numstart = 10, seed = NULL, verbose = FALSE, cpp = TRUE) {

  # only works for trees at the moment, because weightind order is coupled to pwts order
  admixnodes = which(degree(graph, mode='in') == 2)
  stopifnot(length(admixnodes) == 0)

  if(cpp) {
    opt_edge_lengths = cpp_opt_edge_lengths
  }

  graphpops = get_leafnames(graph)
  f3pops = attr(f3_est, 'pops')
  pwts = tree_to_pwts(graph)
  pwts = pwts[, match(f3pops[-1], setdiff(graphpops, f3pops[1]))]

  cmb = combn(0:(length(f3pops)-1), 2)+(1:0)
  ppwts_2d = t(pwts[,cmb[1,]] * pwts[,cmb[2,]])

  branch_lengths = opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve,
                                    lower = rep(0, nrow(pwts)), upper = rep(.Machine$integer.max, nrow(pwts)))
  f3_fit = ppwts_2d %*% branch_lengths
  score = get_score(f3_fit, f3_est, ppinv)

  edges = as_tibble(as_edgelist(graph), .name_repair = ~c('from', 'to')) %>%
    mutate(type = 'edge', weight = c(branch_lengths))

  namedList(edges, score, opt = NULL)
}

f3out_to_fittedf2out = function(f2out, f3out) {
  # computes fitted f2 statistics data frame from f2 and f3 data frames
  # will not include f2(outgroup, X)

  f2out %>%
    right_join(f3out %>% filter(pop2 != pop3) %>% transmute(pop1=pop2, pop2=pop3, f3 = fit), by = c('pop1', 'pop2')) %>%
    left_join(f3out %>% filter(pop2 == pop3) %>% transmute(pop1=pop2, f2_1 = fit), by = c('pop1')) %>%
    left_join(f3out %>% filter(pop2 == pop3) %>% transmute(pop2=pop3, f2_2 = fit), by = c('pop2')) %>%
    transmute(pop1, pop2, est, se, fit = (f2_1 + f2_2 - f3*2), diff = fit - est, z = diff/se, p.value = ztop(z))
}



