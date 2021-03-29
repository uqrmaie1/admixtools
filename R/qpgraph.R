

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
  baseind = args[[12]]
  constrained = args[[13]]

  pwts = fill_pwts(pwts, weights, path_edge_table, path_admixedge_table)
  pwts = pwts[,-baseind] - pwts[,baseind]
  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])
  branch_lengths = opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve, lower, upper, fudge = fudge, constrained = constrained)
  f3_fit = ppwts_2d %*% branch_lengths
  qpgraph_score(f3_fit, f3_est, ppinv)
}



opt_edge_lengths = function(ppwts_2d, ppinv, f3_est, qpsolve, lower, upper, fudge = 1e-4, constrained = TRUE) {
  # finds optimal edge lengths
  # pwts2d: npair x nedge design matrix with paths to outpop
  # ppinv: inverse of npair x npair covariance matrix of jackknife f3 stats
  # f3_est: estimated f3 stats

  pppp = t(ppwts_2d) %*% ppinv
  cc = pppp %*% ppwts_2d
  nc = ncol(cc)
  diag(cc) = diag(cc) + fudge*mean(diag(cc))
  sc = sqrt(diag(cc))
  q1 = (pppp %*% f3_est)[,1]/sc
  cc = cc/tcrossprod(sc)
  if(constrained) qpsolve(cc, q1, cbind(diag(nc), -diag(nc)), c(lower*sc, -upper*sc))/sc
  else solve(cc, q1)/sc
}

qpgraph_score = function(f3_fit, f3_est, ppinv = diag(length(f3_fit))) {
  res = f3_est - f3_fit
  lik = t(res) %*% ppinv %*% res
  lik[1,1]
}

treemix_score = function(f3_fit, f3_est, ppinv) {
  res = f3_est - f3_fit
  se = sqrt(diag(solve(ppinv)))
  sum(res^2/(2*se^2) + log(se * sqrt(2*pi)))
}


#' Compute the fit of an admixture graph
#'
#' Computes the fit of a given admixturegraph from f2-statistics. Drift edge weights and admixture edges weights are optimized until the (negative) likelihood score is minimized. The likelihood score is based on the squared difference between estimated and fitted f3-statistics.
#' @export
#' @param data Input data in one of three forms:
#' \enumerate{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}} (fastest option)
#' \item A directory which contains pre-computed f2-statistics
#' \item The prefix of genotype files (slowest option)
#' }
#' @param graph An admixture graph represented as a matrix of edges, an \code{\link{igraph}} object, or the path to a *qpGraph* graph file. Edges can be constrained by providing a matrix or data frame of edges with columns titled `lower` and `upper` with lower and upper bounds, respectively. By default, admixture edges are constrained to be between zero and one (with paired edges summing to one), and drift edges have a lower bound at zero.
#' @param lambdascale Scales f2-statistics. This has no effect on the fit, but is used in the original *qpGraph* program to display branch weights on a scale that corresponds to FST distances.
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix of
#' f3 statistics will be computed using block-jackknife. Otherwise bootstrap resampling is performed `n` times,
#' where `n` is either equal to `boot` if it is an integer, or equal to the number of blocks if `boot` is `TRUE`.
#' The covariance matrix of f3 statistics will be computed using bootstrap resampling.
#' @param diag Regularization term added to the diagonal elements of the covariance matrix of fitted branch lengths (after scaling by the matrix trace). Default is 0.0001.
#' @param diag_f3 Regularization term added to the diagonal elements of the covariance matrix of estimated f3 statistics (after scaling by the matrix trace). In the original *qpGraph* program, this is fixed at 0.00001.
#' @param lsqmode Least-squares mode. If `TRUE`, the likelihood score will be computed using a diagonal matrix with `1/(sum(diag(f3_var)) * diag_f3)`, in place of the inverse f3-statistic covariance matrix.
#'
#' `lsqmode = 2` will use the identity matrix instead, which is equivalent to computing the score as the sum of squared residuals (`sum((f3_est-f3_fit)^2)`).
#'
#' Both of these options do not take the covariance of f3-statistics into account. This can lead to bias, but is more stable in cases where the inverse f3-statistics covariance matrix can not be estimated precisely (for example because the number of populations is large). An alternative to `lsqmode = TRUE` that doesn't completely ignore the covariance of f3-statistics is to increase `diag_f3`.
#' @param numstart Number of random initializations of starting weights. Defaults to 10. Increasing this number will make the optimization slower, but reduce the risk of not finding the optimal weights. Check the `opt` output to see how much the optimization depends on the starting weights.
#' @param seed Random seed for generating starting weights.
#' @param cpp Use C++ functions. Setting this to `FALSE` will be slower but can help with debugging.
#' @param return_f4 Return all f4-statistics, as well as the z-score of the worst f4-statistic residual. Defaults to `FALSE` because this can be slow.
#' @param f3precomp Optional precomputed f3-statistics. This should be the output of \code{\link{qpgraph_precompute_f3}} and can be provided instead of `data`. This can speed things up if many graphs are evaluated using the same set of f3-statistics.
#' @param f3basepop Optional f3-statistics base population. Inference will be based on f3-statistics of the form `f3(f3basepop; i, j)` for all population pairs `(i, j)`. Defaults to the outgroup population if the graph has one. This option is ignored if `f3precomp` is provided. Changing `f3basepop` should make very little difference.
#' @param ppinv Optional inverse f3-statistics covariance matrix
#' @param f2_blocks_test An optional 3d array of f2-statistics used for computing an out-of-sample score. This should contain only SNP blocks which are not part of `f2_blocks`. This allows to estimate the fit of a graph without overfitting and will not be used during the optimization step
#' @param verbose Print progress updates
#' @return `qpgraph` returns a list with data describing the model fit:
#' \itemize{
#' \item `edges`: A data frame where each row is an edge in the graph. For regular edges,
#' the column `weight` is the estimated edge length, and for admixture edges, it is the estimated admixture weight.
#' \item `score`: The likelihood score of the fitted graph. Lower values correspond to better fits.
#' The score is calculated as the inner product of the residuals (difference between estimated and
#' fitted f3 statistics), weighted by the inverse of the f3 covariance matrix. See \code{\link{qpgraph_score}}
#' \item `f2`: Estimated and fitted f2 statistics. p-values and z-scores test the significance of the difference.
#' \item `f3`: Estimated and fitted f3 statistics. p-values and z-scores test the significance of the difference.
#' \item `f4`: Estimated and fitted f4 statistics (if `return_f4 = TRUE`). p-values and z-scores test the significance of the difference.
#' \item `opt`: A data frame with details of the weight-fitting step, including the randomly sampled starting weights. The column `value` contains the score for each set of starting weights. Columns starting with `x` denote initial weights, and columns starting with `y` denote fitted weights.
#' \item `worst_residual`: The highest absolute z-score of f4-statistics residuals (fitted - estimated f4); (returned if `return_f4 = TRUE`)
#' }
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @seealso \code{\link{qpgraph_wrapper}} for a wrapper functions which calls the original *qpGraph* program.
#' @examples
#' out = qpgraph(example_f2_blocks, example_graph)
#' plot_graph(out$edges)
qpgraph = function(data, graph, lambdascale = 1, boot = FALSE, diag = 1e-4, diag_f3 = 1e-5,
                   lsqmode = FALSE, numstart = 10, seed = NULL, cpp = TRUE, return_f4 = FALSE, f3precomp = NULL,
                   f3basepop = NULL, constrained = TRUE, ppinv = NULL, f2_blocks_test = NULL, verbose = FALSE) {

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
  if(lambdascale == -1) lambdascale = 1
  if(!lambdascale > 0) stop("'lambdascale' has to be > 0!")
  if(!isFALSE(return_f4) && is.null(data)) stop("Can't compute f4 without f2 data!")

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

  if(!is.null(data) && !is.null(f3precomp)) stop("'f2_blocks' and 'f3precomp' can't both be provided!")

  if(!is.null(f3precomp)) {
    precomp = f3precomp
    f3pops = attr(precomp$f3_est, 'pops')
    pairmatch = get_pairindex(match(pops, f3pops))
    precomp$f3_est = precomp$f3_est[pairmatch]
    precomp$ppinv = precomp$ppinv[pairmatch, pairmatch]
    precomp$f3out %<>% slice(pairmatch)
    baseind = which(pops == f3pops[1])
  } else {
    if(is.data.frame(data) || is.matrix(data)) {
      # sets f3 covariance matrix to identity matrix
      if(is.data.frame(data)) {
        data %<>% select(1:3) %>% set_colnames(c('pop1', 'pop2', 'f2')) %>% filter(pop1 < pop2) %>%
          bind_rows(rename(., pop1=pop2,pop2=pop1)) %>% bind_rows(tibble(pop1=pops,pop2=pops,f2=0)) %>%
          arrange(pop1,pop2)
        f2mat = data %>% pivot_wider(names_from=pop2, values_from=f2) %>% column_to_rownames('pop1') %>%
          as.matrix
      } else {
        f2mat = data
        data = f2mat %>% as_tibble(rownames = 'pop1') %>%
          pivot_longer(-pop1, names_to = 'pop2', values_to = 'f2')
      }
      f2mat = f2mat[pops,pops]
      precomp = list()
      f3mat = (t(t(-f2mat + f2mat[,1])+f2mat[,1])/2)[-1,-1]
      precomp$f3_est = c(f3mat[!upper.tri(f3mat)])
      precomp$ppinv = diag(choose(length(pops), 2))
      precomp$f3out = data %>% transmute(pop1,pop2,est=f2,se=1) %>% filter(pop1 < pop2)
    } else {
    f2_blocks = get_f2(data, pops, afprod = FALSE, verbose = verbose)
    precomp = qpgraph_precompute_f3(f2_blocks, pops, f3basepop = f3basepop, lambdascale = lambdascale, boot = boot,
                                    seed = seed, diag_f3 = diag_f3, lsqmode = lsqmode)
    }
    baseind = if(is.null(f3basepop)) 1 else which(pops == f3basepop)
  }
  stopifnot(all(!is.na(precomp$ppinv)))
  if(!is.null(ppinv)) {
    if(!is.null(f3precomp)) stop("'f3precomp' and 'ppinv' can't both be provided!")
    f3pops = attr(ppinv, 'pops')
    pairmatch = get_pairindex(match(pops, f3pops))
    precomp$ppinv = ppinv[pairmatch, pairmatch]
  }

  f3_est = precomp$f3_est
  #f3_est = pmax(0, f3_est)
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
    if(!is.null(seed)) set.seed(seed)
    if('lower' %in% names(edges)) {
      alower = replace_na(pmax(edges$lower[admixedgesfull[1,]], 1-edges$upper[admixedgesfull[2,]], na.rm=T), 0)
      aupper = replace_na(pmin(edges$upper[admixedgesfull[1,]], 1-edges$lower[admixedgesfull[2,]], na.rm=T), 1)
      aupper = pmin(1, aupper) + 1e-9
    } else if(constrained) {
      alower = rep(0, nadmix)
      aupper = rep(1, nadmix)
    } else {
      alower = rep(-Inf, nadmix)
      aupper = rep(Inf, nadmix)
    }
    parmat = matrix(runif(numstart*nadmix), numstart)
    if(verbose) alert_info(paste0('Testing ', nrow(parmat), ' combinations of admixture weight starting values...\n'))
    weightind = graph_to_weightind(graph)
    arglist = list(pwts, ppinv, f3_est, weightind[[1]], weightind[[2]], weightind[[3]],
                   cmb, qpsolve, elower, eupper, diag, baseind, constrained)
    oo = multistart(parmat, optimweightsfun, args = arglist, method = 'L-BFGS-B',
                    lower = alower, upper = aupper, control=list(maxit = 1e4, fnscale = 1e-6),
                    verbose = verbose)
    best = oo %>% top_n(1, -value)
    opt = data.frame(parmat, oo, stringsAsFactors = F)

    admnames = names(V(graph))[admixnodes]
    colnames(opt)[1:(nadmix*2)] = paste0(rep(c('i.', 'e.'), each = nadmix), rep(admnames, 2))
    hilo = apply(as.matrix(oo[,1:nadmix]), 2, function(x) quantile(x, c(0, 1), na.rm = TRUE))

    wts = as.matrix(best[,1:nadmix])[1,]
    ta = c(t(admixedgesfull))
    weight[ta] = c(wts, 1-wts)
    low[ta] = c(pmin(hilo[1,], hilo[2,]), pmin(1-hilo[1,], 1-hilo[2,]))
    high[ta] = c(pmax(hilo[1,], hilo[2,]), pmax(1-hilo[1,], 1-hilo[2,]))
    pwts = fill_pwts(pwts, wts, weightind[[1]], weightind[[2]], weightind[[3]])
  }
  pwts = pwts[,-baseind] - pwts[,baseind]
  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])
  branch_lengths = opt_edge_lengths(ppwts_2d, ppinv, f3_est, qpsolve, elower, eupper, fudge = diag, constrained = constrained)
  f3_fit = ppwts_2d %*% branch_lengths
  score = qpgraph_score(f3_fit, f3_est, ppinv)
  if(!is.null(f2_blocks_test)) {
    precomp_test = qpgraph_precompute_f3(f2_blocks_test, pops, lambdascale = lambdascale, boot = boot,
                                         seed = seed, diag_f3 = diag_f3, lsqmode = lsqmode)
    score_test = qpgraph_score(f3_fit, precomp_test$f3_est, ppinv)
  }

  # if(constrained) weight[normedges] = pmax(0, branch_lengths)
  # else weight[normedges] = branch_lengths
  weight[normedges] = branch_lengths
  edges %<>% select(1:2) %>% set_colnames(c('from', 'to')) %>%  as_tibble %>%
    mutate(type = ifelse(1:n() %in% normedges, 'edge', 'admix'), weight = weight, low = low, high = high)
  f2 = precomp$f2out
  f3 = precomp$f3out %>% mutate(fit = c(f3_fit), diff = est - fit, z = diff/se, p = ztop(z))

  out = namedList(edges, score, f2, f3, opt, ppinv)
  if(!is.null(f2_blocks_test)) {
    out[['score_test']] = score_test
    out[['f3_test']] = precomp_test$f3out
  }
  if(isTRUE(return_f4) || return_f4 == 'f4') {
    if(verbose) alert_info(paste0('Computing f4\n'))
    out$f4 = fitf4(f2_blocks[pops, pops, ], f2, f3)
    out$worst_residual = max(abs(out$f4$z))
  } else if(return_f4 == 'f2') {
    out$worst_residual = fitf2(f2_blocks[pops, pops, ], f2, f3)$z %>% abs %>% max
  } else if(return_f4 == 'f3') {
    out$worst_residual = max(abs(f3$z))
  }
  out
}





#' Compute f3-statistics from f2-statistics.
#'
#' Takes a 3d array of f2 block jackknife estimates and computes f3-statistics between the
#' first population \eqn{p1} and all population pairs \eqn{i, j}: \eqn{f3(p1; p_i, p_j)}
#' @export
#' @param data Input data in one of three forms:
#' \enumerate{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}} (fastest option)
#' \item A directory which contains pre-computed f2-statistics
#' \item The prefix of genotype files (slowest option)
#' }
#' @param pops Populations for which to compute f3-statistics
#' @param f3basepop f3-statistics base population. If `NULL` (the default),
#' the first population in `pops` will be used as the basis.
#' @param lambdascale Scales f2-statistics. This has no effect on the fit, but is used in the original *qpGraph* program to display branch weights on a scale that corresponds to FST distances.
#' @param boot If `FALSE` (the default), block-jackknife resampling will be used to compute standard errors.
#' Otherwise, block-bootstrap resampling will be used to compute standard errors. If `boot` is an integer, that number
#' will specify the number of bootstrap resamplings. If `boot = TRUE`, the number of bootstrap resamplings will be
#' equal to the number of SNP blocks.
#' @param seed Random seed used if `boot` is `TRUE`.
#' @param diag_f3 Regularization term added to the diagonal elements of the covariance matrix of estimated f3 statistics (after scaling by the matrix trace). In the original *qpGraph* program, this is fixed at 0.00001.
#' @param lsqmode Least-squares mode. If `TRUE`, the likelihood score will be computed using a diagonal matrix with `1/(sum(diag(f3_var)) * diag_f3)`, in place of the inverse f3-statistic covariance matrix. `lsqmode = 2` will use the identity matrix instead, which is equivalent to computing the score as the sum of squared residuals (`sum((f3_est-f3_fit)^2)`). Both of these options do not take the covariance of f3-statistics into account. This can lead to bias, but is more stable in cases where the inverse f3-statistics covariance matrix can not be estimated precisely (for example because the number of populations is large). An alternative to using `lsqmode = TRUE` which doesn't completely ignore the covariance of f3-statistics is to increase `diag_f3`.
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
#' qpgraph_precompute_f3(f2_dir, pops)
#' }
qpgraph_precompute_f3 = function(data, pops, f3basepop = NULL, lambdascale = 1, boot = FALSE,
                                 seed = NULL, diag_f3 = 1e-5, lsqmode = FALSE) {
  # returns list of f3_est and ppinv for subset of populations.
  # f3_est and ppinv are required for qpgraph_slim; f2out and f3out are extra output
  # f2_blocks may contain more populations than the ones used in qpgraph
  # f2_blocks input here should be subset which is used by qpgraph function

  #----------------- read f-stats -----------------
  if(!is.null(f3basepop)) pops = c(f3basepop, setdiff(pops, f3basepop))

  if(!is.null(seed)) set.seed(seed)
  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  matstatfun = ifelse(boot, boot_mat_stats, jack_mat_stats)
  arrstatfun = ifelse(boot, boot_arr_stats, jack_arr_stats)
  f2_blocks = (get_f2(data, pops) * lambdascale) %>% samplefun
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

  f3_blocks = (f2_blocks[,rep(1, npop),,drop=F] + f2_blocks[rep(1, npop),,,drop=F] - f2_blocks)/2
  #f3_blocks = array(pmax(0, f3_blocks), dim(f3_blocks), dimnames(f3_blocks))
  f3_blocks_2d = arr3d_to_pairmat(f3_blocks[-1,-1,,drop=F])
  f3dat = matstatfun(f3_blocks_2d, block_lengths)
  #f3dat = jack_mat_stats(f3_blocks_2d, block_lengths)
  f3_est = f3dat$est
  f3_var = f3dat$var
  f3out = tibble(pop1 = pops[1],
                 pop2 = pops[cmb[1,]+1],
                 pop3 = pops[cmb[2,]+1],
                 est = f3_est, se = sqrt(diag(f3_var)))
  add_diag = sum(diag(f3_var)) * diag_f3
  diag(f3_var) = diag(f3_var) + add_diag
  # in qpGraph diag_f3 is 1e-5; has large effect on magnitude of likelihood score
  if(lsqmode == 2) ppinv = diag(nrow(f3_var))
  else if(lsqmode) ppinv = diag(nrow(f3_var)) / add_diag
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

qpgraph_anorexic = function(f3precomp, graph, diag = 1e-4,
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
                                    fudge = diag, constrained = TRUE)
  f3_fit = ppwts_2d %*% branch_lengths
  score = qpgraph_score(f3_fit, f3_est, ppinv)

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
    transmute(pop1, pop2, est, se, fit = (f2_1 + f2_2 - f3*2), diff = est - fit, z = diff/se, p = ztop(z))
}

fitf2 = function(f2_blocks, f2, f3) {
  # returns a tibble with estimated and fitted f2-statistics
  cmb = combn(0:(dim(f2_blocks)[1]-1), 2)+(1:0)
  f2_out = f3 %>% filter(pop2 == pop3) %$% fit
  f2_fit = f3 %>% mutate(f21 = f2_out[cmb[1,]], f22 = f2_out[cmb[2,]], f2fit = (f21 + f22 - fit*2))
  f2_fit2 = f2 %>%
    left_join(f2_fit, by = c('pop1'='pop2', 'pop2'='pop3')) %>%
    filter(!is.na(f2fit)) %>%
    transmute(pop1, pop2, fit = f2fit) %>%
    bind_rows(f2_fit %>% filter(pop2 == pop3) %>% select(pop1, pop2, fit))
  f2(f2_blocks, verbose = FALSE) %>%
    inner_join(f2_fit2, by = c('pop1', 'pop2')) %>%
    mutate(diff = est - fit, z = diff/se, p = ztop(z))
}


fitf4 = function(f2_blocks, f2, f3) {
  # returns a tibble with estimated and fitted f4-statistics

  cmb = combn(0:(dim(f2_blocks)[1]-1), 2)+(1:0)
  f2_out = f3 %>% filter(pop2 == pop3) %$% fit
  f2_fit = f3 %>% mutate(f21 = f2_out[cmb[1,]], f22 = f2_out[cmb[2,]], f2fit = (f21 + f22 - fit*2))
  f2_fit2 = f2 %>%
    left_join(f2_fit, by = c('pop1'='pop2', 'pop2'='pop3')) %>%
    filter(!is.na(f2fit)) %>%
    select(pop1, pop2, f2fit) %>%
    bind_rows(f2_fit %>% filter(pop2 == pop3) %>% transmute(pop1, pop2, f2fit = fit)) %>%
    bind_rows(rename(., pop1 = pop2, pop2 = pop1)) %>%
    bind_rows(tibble(pop1 = unique(.$pop1), pop2 = pop1, f2fit = 0))
  x = f4(f2_blocks, unique_only = FALSE, verbose = FALSE) %>% select(-z, -p)
  x %>%
    left_join(f2_fit2 %>% rename(c1 = f2fit), by = c('pop1' = 'pop1', 'pop4' = 'pop2')) %>%
    left_join(f2_fit2 %>% rename(c2 = f2fit), by = c('pop2' = 'pop1', 'pop3' = 'pop2')) %>%
    left_join(f2_fit2 %>% rename(c3 = f2fit), by = c('pop1' = 'pop1', 'pop3' = 'pop2')) %>%
    left_join(f2_fit2 %>% rename(c4 = f2fit), by = c('pop2' = 'pop1', 'pop4' = 'pop2')) %>%
    mutate(fit = (c1 + c2 - c3 - c4)/2, diff = est - fit, z = diff/se, p = ztop(z)) %>%
    select(-c1:-c4)
}



#' Compare the fit of two qpgraph models
#'
#' Takes two data frames with model fits computed on two graphs for on the same populations and tests whether the scores of one graph are significantly better than the scores of the other
#' @param fits1 The fits of the first graph
#' @param fits2 The fits of the second graph
#' @param boot should match the `boot` parameter in `qpgraph_resample_snps` (`FALSE` by default)
#' @examples
#' \dontrun{
#' nblocks = dim(example_f2_blocks)[3]
#' train = sample(1:nblocks, round(nblocks/2))
#' fits1 = qpgraph_resample_snps(example_f2_blocks[,,train], graph = graph1,
#'                               f2_blocks_test = example_f2_blocks[,,-train])
#' fits2 = qpgraph_resample_snps(example_f2_blocks[,,train], graph = graph2,
#'                               f2_blocks_test = example_f2_blocks[,,-train])
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
#' @param scores1 Scores for the first graph
#' @param scores2 Scores for the second graph
#' @seealso \code{\link{qpgraph_resample_multi}}
#' @examples
#' \dontrun{
#' fits = qpgraph_resample_multi(f2_blocks, list(graph1, graph2), nboot = 100)
#' compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)
#' }
compare_fits = function(scores1, scores2) {

  scorediff = scores1 - scores2
  ci_low = unname(quantile(scorediff, 0.025, na.rm = T))
  ci_high = unname(quantile(scorediff, 0.975, na.rm = T))

  scorediff = na.omit(scores1 - scores2)
  stats = boot_mat_stats(t(scorediff), rep(1, length(scorediff)))

  diff = stats$est
  se = sqrt(stats$var)
  z = diff/se
  p = ztop(z)[1]
  frac = mean(scorediff < 0)
  p_emp_nocorr = min(frac, 1-frac)*2
  p_emp = max(p_emp_nocorr, 1/length(scorediff))
  namedList(diff, se, z, p, p_emp, p_emp_nocorr, ci_low, ci_high)
}




#' Evaluate a qpgraph model many times
#'
#' This function is used in combination with \code{\link{compare_fits}} in order to test whether one graph has a significantly better fit than another. using a list of bootstrap resampled f2-statistics and corresponding test-set f2-statistics
#' @export
#' @param f2_blocks A list of bootstrap resampled f2-statistics
#' @param graph An admixture graph
#' @param f2_blocks_test A list of f2-statistics from all blocks which were not used in the corresponding f2-array in `f2_blocks`
#' @param verbose Print progress updates
#' @param ... Arguments passed to \code{\link{qpgraph}}
#' @return A data frame with \code{\link{qpgraph}} results for each iteration of bootstrap resampled f2-statistics
#' @seealso \code{\link{compare_fits}} \code{\link{boo_list}}
qpgraph_resample_snps2 = function(f2_blocks, graph, f2_blocks_test, verbose = TRUE, ...) {

  #progressr::with_progress({
  #pb = progressr::progressor(steps = length(f2_blocks))

  ell = list(...)
  fun = function(f2dat, f2dat_test, g) function() {
    #pb$tick()
    #pb()
    safely(qpgraph)(data = f2dat, graph = g, f2_blocks_test = f2dat_test, verbose = FALSE, ...)
  }

  #pb = progress::progress_bar$new(total = length(f2_blocks))
  tibble(id = seq_len(length(f2_blocks)), graph = list(graph), f2_blocks, f2_blocks_test) %>%
    mutate(fun2 = pmap(list(f2_blocks, f2_blocks_test, graph), fun)) %>%
    mutate(out = furrr::future_invoke_map(fun2, .progress = verbose, .options = furrr::furrr_options(seed = TRUE)),
           result = map(out, 'result', .null = tibble()), error = map(out, 'error')) %>%
    select(-out, -fun2) %>% unnest_wider(result)
  #})
}

# qpgraph_resample_snps2 = function(f2_blocks, graph, f2_blocks_test, verbose = TRUE, ...) {
#
#   fun = function(f2dat, f2dat_test, g) function() safely(qpgraph)(data = f2dat, graph = g, f2_blocks_test = f2dat_test, verbose = FALSE, ...)
#
#   tibble(id = seq_len(length(f2_blocks)), graph = list(graph), f2_blocks, f2_blocks_test) %>%
#     mutate(fun2 = pmap(list(f2_blocks, f2_blocks_test, graph), fun)) %>%
#     mutate(out = furrr::future_invoke_map(fun2, .progress = verbose),
#            result = map(out, 'result', .null = tibble()), error = map(out, 'error')) %>%
#     select(-out, -fun2) %>% unnest_wider(result)
# }



#' Evaluate a qpgraph models many times
#'
#' This function is used in combination with \code{\link{compare_fits}} in order to test whether one graph has a significantly better fit than another. It creates bootstrap resampled SNP block training and test sets, and uses them to evaluate multiple graphs.
#' @export
#' @param f2_blocks 3d array of f2-statistics
#' @param graphlist A list of admixture graphs
#' @param nboot Number of bootstrap iterations
#' @param verbose Print progress updates
#' @param ... Arguments passed to \code{\link{qpgraph}}
#' @return A list of same length as `graphlist` with data frames with \code{\link{qpgraph}} results for each iteration of bootstrap resampled f2-statistics
#' @examples
#' \dontrun{
#' fits = qpgraph_resample_multi(f2_blocks, list(graph1, graph2), nboot = 100)
#' compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)
#' }
#' @seealso \code{\link{compare_fits}}
qpgraph_resample_multi = function(f2_blocks, graphlist, nboot, verbose = TRUE, ...) {

  boo = boo_list(f2_blocks, nboot = nboot)
  #f3pre = map(graphlist, ~qpgraph_precompute_f3(f2_blocks, get_leafnames(.))$ppinv)
  #map2(graphlist, f3pre, function(.x, .y, ...) qpgraph_resample_snps2(
  #  boo$boo, .x, boo$test, ppinv = .y, verbose = verbose, ...), ...)
  map(graphlist, function(.x, ...) qpgraph_resample_snps2(
    boo$boo, .x, boo$test, verbose = verbose, ...), ...)
}


#' Compare the fit of two qpgraph models
#'
#' Takes two data frames with model fits computed on two graphs for on the same populations and tests whether the scores of one graph are significantly better than the scores of the other.
#' @param fit1 The fit of the first graph
#' @param fit2 The fit of the second graph
#' @param f2_blocks f2 blocks used for fitting `fit1` and `fit2`. Used in combination with `f2_blocks_test` to compute f-statistics covariance matrix.
#' @param f2_blocks_test f2 blocks which were not used for fitting `fit1` and `fit2`
#' @param boot If `FALSE` (the default), block-jackknife resampling will be used to compute standard errors.
#' Otherwise, block-bootstrap resampling will be used to compute standard errors. If `boot` is an integer, that number
#' will specify the number of bootstrap resamplings. If `boot = TRUE`, the number of bootstrap resamplings will be
#' equal to the number of SNP blocks. If bootstrap resampling is enabled, empirical p-values (`p_emp`) and 95 confidence intervals (`ci_low` and `ci_high`) will be reported.
#' @param seed Random seed used if `boot` is `TRUE`. Does not need to match a seed used in fitting the models
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
  scores1 = map_dbl(1:dim(f2_blocks_test)[3], ~qpgraph_score(f3_fit$fit.x, f3_test[,.], ppinv))
  scores2 = map_dbl(1:dim(f2_blocks_test)[3], ~qpgraph_score(f3_fit$fit.y, f3_test[,.], ppinv))

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


fit_to_worstf3z = function(fit) {
  fit$f3 %>% left_join(fit$f3_test, by = c('pop1', 'pop2', 'pop3')) %>%
    mutate(z = (est.y-fit)/sqrt(se.x^2+se.y^2)) %>% slice_max(abs(z)) %>% pluck('z', 1)
}



