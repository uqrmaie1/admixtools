
f2_qpgraph = function(p1, p2, c1, c2) {
  # c1, c2: counts of individuals
  (p1-p2)^2 - p1*(1-p1)/(2*c1-1) - p2*(1-p2)/(2*c2-1)
}

f2_qpgraph1 = function(p1, p2) {
  # c1, c2: counts of individuals
  (p1-p2)^2 - p1*(1-p1) - p2*(1-p2)
}


fst_qpgraph = function(p1, p2, c1, c2) {
  # c1, c2: counts of individuals
  # is the same as fst_reich
  num = (p1-p2)^2 - p1*(1-p1)/(2*c1-1) - p2*(1-p2)/(2*c2-1)
  denom = p1 + p2 - 2*p1*p2
  num/denom
}


qpgraph_wrapper2 = function(bin='./qpGraph', parfile='./parfile', graphfile='./graphfile', outfile='./out', printonly=FALSE) {
  # wrapper around AdmixTools qpGraph
  # input is locations of parfile and graphfile
  # output is parsed output

  cmd = paste0(bin, ' -p ', parfile, ' -g ', graphfile, ' > ', outfile)
  if(printonly) {
    print(cmd)
  } else{
    system(cmd)
    return(parse_qpgraph_output(outfile))
  }
}



#' Wrapper function around the original qpGraph program
#' @export
#' @param graph an admixture graph or qpGraph graph file
#' @param bin location of the qpGraph binary
#' @param pref prefix of the packedancestrymap format genotype files.
#' @param parfile qpGraph parameter file
#' @param outdir output directory
#' @param printonly should output be executed or the command just be printed?
#' @param lambdascale lambdascale
#' @param diag diag
#' @param outpop outgroup population
#' @param lsqmode least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @param hires hires
#' @param forcezmode forcezmode
#' @param allsnps allsnps
#' @param bigiter bigiter
#' @return a list with parsed qpGraph output
#' \enumerate{
#' \item \code{edges}: data frame
#' \item \code{score}: scalar
#' \item \code{f2}: data frame
#' }
#' @examples
#' \dontrun{
#' qpgraph_wrapper(example_graph,
#'                  bin = 'path/to/qpGraph',
#'                  pref = 'path/to/packedancestrymap_prefix')
#' qpgraph_wrapper('path/to/graphfile',
#'                 bin = 'path/to/qpGraph',
#'                 parfile = 'path/to/parfile')
#' }
qpgraph_wrapper = function(graph, bin, pref = NULL, parfile = NULL, outdir='.',
                           printonly=FALSE, lambdascale=-1, diag=0.0001, outpop='NULL',
                           lsqmode='NO', hires='NO', forcezmode='NO', allsnps='NO', bigiter=100, env='') {
  # wrapper around AdmixTools qpGraph
  # makes parfile and graphfile
  stopifnot(!is.null(parfile) | !is.null(pref))

  if(is.null(parfile)) {
    pref = normalizePath(pref, mustWork = FALSE)
    parfile = paste0('indivname:       ', pref, '.ind\n',
                     'snpname:         ', pref, '.snp\n',
                     'genotypename:    ', pref, '.geno\n',
                     'outpop:         ', outpop, '\n',
                     'blgsize: 0.05\n',
                     'details: YES\n',
                     'fstdetails: YES\n',
                     'diag: ', diag, '\n',
                     'lsqmode: ', lsqmode, '\n',
                     'hires: ', hires, '\n',
                     'forcezmode: ', forcezmode, '\n',
                     'allsnps: ', allsnps, '\n',
                     'lambdascale: ', lambdascale, '\n',
                     'bigiter: ', bigiter, '\n')
    pf = paste0(outdir, '/parfile')
    write(parfile, pf)
  } else {
    pf = parfile
  }

  if(class(graph) != 'character') {

    if(class(graph)[1] == 'igraph') graph = igraph::as_edgelist(graph)
    edg = as_tibble(graph) %>% set_colnames(c('from', 'to'))
    edg %<>% group_by(.data$to) %>% mutate(type = ifelse(n()==1, 'edge', 'admix')) %>% ungroup
    e1 = (edg %>% filter(.data$type == 'edge'))$from
    e2 = (edg %>% filter(.data$type == 'edge'))$to
    a1 = (edg %>% filter(.data$type == 'admix'))$from
    a2 = (edg %>% filter(.data$type == 'admix'))$to
    leaves = setdiff(edg$to, edg$from)
    admix = tibble()
    for(m in unique(a2)) {
      admix %<>% bind_rows(tibble(v1='admix', v2=m, v3=(edg %>% filter(.data$to == m))$from[1], v4=(edg %>% filter(.data$to == m))$from[2]))
    }

    simfile = tibble(v1 = c('root'), v2 = c('R'), v3='', v4='') %>%
      bind_rows(tibble(v1 = 'label', v2=leaves, v3=leaves, v4='')) %>%
      bind_rows(tibble(v1='edge', v2=paste0('e', 1:length(e1)), v3=e1, v4=e2)) %>%
      bind_rows(admix)

    gf = paste0(outdir, '/graphfile')
    simfile %>% write_tsv(gf, col_names=F)
  } else {
    stopifnot(file.exists(graph))
    gf = graph
  }

  qpgraph_wrapper2(bin = paste0(env,' ', bin), parfile = pf, graphfile = gf,
                   outfile = paste0(outdir, '/out'), printonly = printonly)
}



#' @export
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


#' @export
expand_path = function(path) {
  # igraph path (sequence of vertices)
  # duplicates inner vertices, so that this works with igraph functions that process vertex sequences pairwise
  ln = length(path)
  path[c(1, 1+rep(seq_len(ln-2), each=2), ln)]
}

#' @export
graph_to_weightind = function(graph) {
  # input igraph object
  # output:
  # map (leaf, edge) -> paths
  # map path -> weights
  # ultimately: indices for weights into paths, indices for paths into pwts

  leaves = get_leaves(graph)
  admixnodes = which(degree(graph, mode='in')==2)
  admixedges = unlist(incident_edges(graph, admixnodes, mode='in'))
  normedges = setdiff(1:length(E(graph)), admixedges)
  paths = all_simple_paths(graph, V(graph)[1], leaves[-1], mode='out')
  ends = sapply(paths, tail, 1)
  edge_per_path = paths %>% map(expand_path) %>% map(~get.edge.ids(graph, .))
  weight_per_path = edge_per_path %>% map(~(which(admixedges %in% .)))

  path_edge_table = do.call(rbind, lapply(seq_len(length(weight_per_path)), function(i) tibble(path=i, edge=c(edge_per_path[[i]])))) %>% mutate(edge2 = match(edge, normedges)) %>% filter(!is.na(edge2)) %>% mutate(leaf = as.vector(ends[path]), leaf2 = match(leaf, leaves[-1])) %>% left_join(enframe(c(table(ends))) %>% transmute(leaf=as.numeric(name), numpaths=value), by='leaf') %>% group_by(leaf2, edge2) %>% mutate(cnt = n(), keep = cnt < numpaths) %>% filter(keep) %>% as.matrix

  path_admixedge_table = do.call(rbind, lapply(seq_len(length(weight_per_path)), function(i) tibble(path=i, admixedge=c(weight_per_path[[i]])))) %>% as.matrix
  list(path_edge_table, path_admixedge_table, length(paths))
}



#' @export
fill_pwts = function(pwts, weights, path_edge_table, path_admixedge_table, numpaths = NULL) {
  # puts weights onto pwts, using index matrix and vectors

  if(length(weights)==0) return(pwts)
  wts2 = rep(weights, each=2)*c(1,-1) + (0:1)
  path_weights = path_admixedge_table %>% as_tibble %>% mutate(w = wts2[admixedge]) %>% group_by(path) %>% summarize(w = prod(w))
  pwts_weights = path_edge_table %>% as_tibble %>% left_join(path_weights, by='path') %>% mutate(w = ifelse(is.na(w), 1, w)) %>% group_by(edge2, leaf2) %>% summarize(w = sum(w)) %>% ungroup %>% as.matrix
  pwts[pwts_weights[,1:2]] = pwts_weights[,3]
  pwts
}

#' @export
optimweightsfun = function(weights, args) {
  # likelihood function used in optimizing admixture weights
  # weights is vector of admixture weights to be optmized; only values for first incoming edge; 2nd is 1 - first

  pwts = args[[1]]
  ppinv = args[[2]]
  f3_jest = args[[3]]
  path_edge_table = args[[4]] # indices into pwts with weight positions
  path_admixedge_table = args[[5]]
  #index3 = args[[6]]
  cmb = args[[7]]
  pwts = fill_pwts(pwts, weights, path_edge_table, path_admixedge_table)
  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])
  q2 = opt_edge_lengths(ppwts_2d, ppinv, f3_jest, args[[8]])
  w2 = (ppwts_2d %*% q2) - f3_jest
  lik = t(w2) %*% ppinv %*% w2
  lik[1,1]
}

#' @export
opt_edge_lengths = function(ppwts_2d, ppinv, f3_jest, qpsolve) {
  # finds optimal edge lengths
  # pwts2d: npair x nedge design matrix with paths to outpop
  # ppinv: inverse of npair x npair matrix of varianc-covariance matrix of jackknife f3 stats
  # f3_jest: estimated f3 stats
  pppp = t(ppwts_2d) %*% ppinv
  cc = pppp %*% ppwts_2d
  diag(cc) = diag(cc) + mean(diag(cc))*0.0001
  cc = (cc+t(cc))/2
  q1 = -(pppp %*% f3_jest)[,1]
  nc = ncol(cc)
  -qpsolve(cc, q1, -diag(nc), rep(0, nc))

  # earlier version of qpsolve:
  # this is very slow, don't use it (takes 10 times as long):
  # pracma::quadprog(cc, q1, A=-diag(nc), b=rep(0, nc))$xmin
  # this is as fast as quadprogpp::QP.Solve, and on CRAN:
  # -quadprog::solve.QP(cc, q1, -diag(nc), rep(0, nc))$solution
  # this was used while most testing was done:
  #quadprogpp::QP.Solve(cc, q1, CI = diag(nc), ci0 = rep(0, nc))
}

#' @export
get_score = function(ppwts_2d, ppinv, f3_jest, q2) {
  w2 = (ppwts_2d %*% q2) - f3_jest
  lik = t(w2) %*% ppinv %*% w2
  lik[1,1]
}




#' Compute the fit of an admixturegraph
#'
#' Computes the fit of an admixturegraph for a given graph topology and empirical f2-block-jackknife statistics. f2-block-jackknife statistics can be provided via the arguments \code{f2_blocks} and \code{block_lengths}, or via \code{f2_dir}.
#' @export
#' @param graph an admixture graph represented as a matrix of edges, an \code{\link{igraph}} object, or the path to a qpGraph graph file.
#' @param f2_blocks 3d array of block-jackknife leave-one-block-out estimates of f2 statistics. output of \code{\link{afs_to_f2_blocks}}.
#' @param block_lengths the jackknife block lengths used in computing the f2 statistics. see \code{\link{get_block_lengths}}.
#' @param f2_dir a directory with f2 statistics for each population pair in the graph. must contain 'block_lengths.RData'.
#' @param lsqmode least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param fnscale try increasing or decreasing this, if optimization does not converge.
#' @param fudge try increasing this, if you get the error message \code{constraints are inconsistent, no solution!}.
#' @param numstart number of random initializations. defaults to 10 times the number of admixture nodes.
#' @param seed seed for generating starting weights.
#' @param verbose print optimization iterations
#' @param cpp should optimization be done using C++ or R function? \code{cpp = TRUE} is much faster.
#' @return a list of qpGraph output data
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @seealso \code{\link{qpgraph_wrapper}} for a wrapper functions which call the original qpGraph program; \code{\link{qpgraph_slim}} for a faster function function which requires f3 estimates and an inverted covariance matrix as input instead.
#' @examples
#' out = qpgraph(example_graph, example_f2_blocks, example_block_lengths)
#' plot_graph(out$edges)
qpgraph = function(graph, f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL, lsqmode=FALSE, f2_denom = 1, fnscale = 1e-6, fudge = 1e-3, numstart = NULL, seed = NULL, verbose = FALSE, cpp = TRUE) {
  # modelled after AdmixTools qpGraph

  #----------------- process graph -----------------
  if('data.frame' %in% class(graph) | 'matrix' %in% class(graph)) {
    edges = as.matrix(graph)
  } else if(class(graph) == 'character') {
    edges = parse_qpgraph_graphfile(graph)
  } else if(class(graph)[1] == 'igraph') {
    edges = igraph::as_edgelist(graph)
  } else {
    stop(paste0('Cannot parse graph of class ', class(graph),'!'))
  }
  grph = graph_from_edgelist(edges)
  if(class(graph)[1] == 'igraph') grph = graph
  nedges = length(E(grph))
  admixnodes = which(degree(grph, mode='in')==2)
  nadmix = length(admixnodes)
  admixedgesfull = sapply(seq_len(nadmix), function(i) incident_edges(grph, admixnodes, mode='in')[[i]][1:2])
  normedges = setdiff(1:nedges, admixedgesfull)

  pops = get_leafnames(grph)
  npop = length(pops)
  cmb = combn(0:(npop-1), 2)+(1:0)

  popind = match(pops, dimnames(f2_blocks)[[1]])

  precomp = qpgraph_precompute_f3(pops, f2_blocks = f2_blocks, block_lengths = block_lengths, f2_dir = f2_dir, f2_denom = f2_denom, fudge = fudge, lsqmode = lsqmode)

  f3_jest = precomp$f3_jest
  ppinv = precomp$ppinv
  f2out = precomp$f2out
  f3out = precomp$f3out

  if(verbose) alert_info('preprocessing done')
  weightind = graph_to_weightind(grph)
  weight = rep(NA, nedges)
  pwts = graph_to_pwts(grph)
  opt = NULL
  qpsolve = function(...) quadprog::solve.QP(...)$solution

  if(cpp) {
    optimweightsfun = cpp_optimweightsfun
    opt_edge_lengths = cpp_opt_edge_lengths
    fill_pwts = cpp_fill_pwts
  }

  if(nadmix > 0) {
    if(is.null(numstart)) numstart = 10*nadmix
    set.seed(seed)
    parmat = matrix(runif(numstart*nadmix), numstart)
    if(verbose) alert_info(paste0('testing ', nrow(parmat), ' combinations of admixture weight starting values\n'))
    arglist = list(pwts, ppinv, f3_jest, weightind[[1]], weightind[[2]], weightind[[3]], cmb, qpsolve)
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

  q2 = opt_edge_lengths(ppwts_2d, ppinv, f3_jest, qpsolve)
  score = get_score(ppwts_2d, ppinv, f3_jest, q2)
  weight[normedges] = q2
  edges = as_tibble(edges, .name_repair = ~c('from', 'to')) %>%
    mutate(type = ifelse(1:n() %in% normedges, 'edge', 'admix'), weight=weight, label=round(weight, 2))

  namedList(edges, score, f2=f2out, f3=f3out, opt)
}

#' Compute the fit of an admixturegraph.
#'
#' Computes the fit of an admixture graph for a given graph topology given f3 estimates and an inverted covariance matrix.
#' @export
#' @param graph an admixture graph represented as an \code{\link{igraph}} object
#' @param f3_jest matrix of f3 estimates. output of \code{\link{qpgraph_precompute_f3}}
#' @param ppinv inverted covariance matrix of f3-statistics
#' @param fnscale try increasing or decreasing this, if optimization does not converge.
#' @param fudge try increasing this, if you get the error message \code{constraints are inconsistent, no solution!}.
#' @param numstart number of random initializations. defaults to 10 times the number of admixture nodes.
#' @param seed seed for generating starting weights.
#' @param verbose print optimization iterations
#' @param cpp should optimization be run in C++ or in R? C++ is faster.
#' @return a list of qpGraph output data
#' @seealso \code{\link{qpgraph_wrapper}}for a wrapper functions which call the original qpGraph program; \code{\link{qpgraph}} for a slower function which requires f2 block jackknife statistics as input instead. \code{\link{qpgraph_precompute_f3}} computes the required input from a 3d array of \code{f2_statistics}
#' @examples
#' pops = get_leafnames(example_igraph)
#' precomp = qpgraph_precompute_f3(pops, example_f2_blocks, example_block_lengths)
#' f3_jest = precomp$f3_jest
#' ppinv = precomp$ppinv
#' out = qpgraph_slim(example_igraph, f3_jest, ppinv)
#' plot_graph(out$edges)
qpgraph_slim = function(graph, f3_jest, ppinv, fnscale = 1e-6, numstart = 10,
                        seed = NULL, verbose = FALSE, cpp = TRUE) {
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
  f3pops = attr(f3_jest, 'pops')

  pairmatch = get_pairindex(match(graphpops, f3pops))
  f3_jest = f3_jest[pairmatch]
  ppinv = ppinv[pairmatch, pairmatch]
  stopifnot(all(!is.na(ppinv)))

  weightind = graph_to_weightind(graph)
  weight = rep(NA, nedges)
  pwts = graph_to_pwts(graph)
  opt = NULL
  cmb = combn(0:(length(f3pops)-1), 2)+(1:0)
  qpsolve = function(...) quadprog::solve.QP(...)$solution

  if(nadmix > 0) {
    set.seed(seed)
    parmat = matrix(runif(numstart*nadmix), numstart)

    arglist = list(pwts, ppinv, f3_jest, weightind[[1]], weightind[[2]], weightind[[3]], cmb, qpsolve)
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

  q2 = opt_edge_lengths(ppwts_2d, ppinv, f3_jest, qpsolve)
  score = get_score(ppwts_2d, ppinv, f3_jest, q2)
  weight[normedges] = q2
  edges = as_tibble(as_edgelist(graph), .name_repair = ~c('from', 'to')) %>%
    mutate(type = ifelse(1:n() %in% normedges, 'edge', 'admix'), weight=weight, label=round(weight, 2))

  namedList(edges, score, opt)
}


#' Compute f3-statistics from f2-statistics.
#'
#' Takes a 3d array of f2 block jackknife estimates and computes f3-statistics between the first population \code{p1} and all population pairs \code{i, j}: \code{f3(p1; p_i, p_j)}
#' @export
#' @param pops populations for which to compute f3-statistics
#' @param f2_blocks 3d array of block-jackknife leave-one-block-out estimates of f2 statistics. output of \code{\link{afs_to_f2_blocks}}. they are weighted by inverse of outgroup heterozygosity, if outgroup was specified.
#' @param block_lengths the jackknife block lengths used in computing the f2 statistics. see \code{\link{get_block_lengths}}.
#' @param f2_dir a directory with f2 statistics for each population pair in the graph. must contain 'block_lengths.RData'.
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param fudge try increasing this, if you get the error message \code{constraints are inconsistent, no solution!}.
#' @param lsqmode least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @return a list with four items
#' \enumerate{
#' \item \code{f3_jest} a matrix with f3-statistics for all population pairs with the output
#' \item \code{ppinv} a matrix with the inverse of the f3-statistic covariance matrix
#' \item \code{f2out} a data frame with f2 estimates
#' \item \code{f3out} a data frame with f3 estimates
#' }
#' @examples
#' pops = get_leafnames(example_igraph)
#' qpgraph_precompute_f3(pops, example_f2_blocks, example_block_lengths)
#' \dontrun{
#' qpgraph_precompute_f3(pops, f2_dir = f2_dir, f2_denom = 0.278)
#' }
qpgraph_precompute_f3 = function(pops, f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL,
                                 f2_denom = 1, fudge=1e-3, lsqmode=FALSE) {
  # returns list of f3_jest and ppinv for subset of populations.
  # f3_jest and ppinv are required for qpgraph_slim; f2out and f3out are extra output
  # f2_blocks may contain more populations than the ones used in qpgraph
  # f2_blocks input here should be subset which is used by qpgraph function

  if(is.null(f2_dir) & (is.null(f2_blocks) | is.null(block_lengths))) stop('You have to provide an f2_dir argument, or f2_blocks and block_lengths!')
  if(!is.null(f2_dir) & (is.null(f2_blocks) | is.null(block_lengths))) {
    f2_blocks = read_f2(f2_dir, pops = pops)
    load(paste0(f2_dir, '/block_lengths.RData'))
  }
  stopifnot(all(pops %in% dimnames(f2_blocks)[[1]]))

  f2_blocks = f2_blocks / f2_denom
  f2_blocks = rray(f2_blocks[pops, pops, ])

  npop = length(pops)
  npair = choose(npop, 2)
  cmb = combn(0:(npop-1), 2)+(1:0)

  f2 = bj_arr_stats(f2_blocks, block_lengths)
  f2out = tibble(pop1=combn(pops, 2)[1,],
                 pop2=combn(pops, 2)[2,],
                 f2est = f2[[1]][lower.tri(f2[[1]])],
                 se = sqrt(f2[[2]][lower.tri(f2[[2]])]))

  f3_blocks = (f2_blocks[,1,] + f2_blocks[1,,] - f2_blocks)/2
  f3_blocks_2d = arr3d_to_pairmat(f3_blocks[-1,-1,])
  f3dat = bj_mat_stats(f3_blocks_2d, block_lengths)
  f3_jest = f3dat$jest
  f3_jvar = f3dat$jvar
  f3out = tibble(pop1=pops[cmb[1,]],
                 pop2=pops[cmb[2,]],
                 f3est = f3_jest, se = sqrt(diag(f3_jvar)))
  diag(f3_jvar) = diag(f3_jvar) + sum(diag(f3_jvar))*fudge
  # in qpGraph fudge is 1e-5; sometimes quadprog doesn't converge unless this is larger; has large effect on magnitude of likelihood score
  if(lsqmode) ppinv = diag(1/diag(f3_jvar))
  else ppinv = solve(f3_jvar)

  f3_jest = structure(f3_jest, pops = pops)
  ppinv = structure(ppinv, pops = pops)
  namedList(f3_jest, ppinv, f2out, f3out)
}

get_pairindex = function(perm) {
  # returns index vector that matches population pairs in
  # f3_jest and ppinv (which were computed using pops) to graph populations
  cmb = combn(0:(length(perm)-1), 2)+(1:0)
  popind = setdiff(perm, 1)
  orig_order = apply(cmb+1, 2, paste0, collapse='')
  new_order = apply(matrix(popind[c(cmb)], 2), 2, function(x) paste0(sort(x), collapse=''))
  match(new_order, orig_order)
}

#' @export
qpgraph_anorexic = function(graph, f3_jest, ppinv, fnscale = 1e-6,
                            numstart = 10, seed = NULL, verbose = FALSE, cpp = TRUE) {

  # only works for trees at the moment, because weightind order is coupled to pwts order
  admixnodes = which(degree(graph, mode='in')==2)
  stopifnot(length(admixnodes) == 0)

  if(cpp) {
    opt_edge_lengths = cpp_opt_edge_lengths
  }

  graphpops = get_leafnames(graph)
  f3pops = attr(f3_jest, 'pops')
  pwts = tree_to_pwts(graph)
  pwts = pwts[, match(f3pops[-1], setdiff(graphpops, f3pops[1]))]

  cmb = combn(0:(length(f3pops)-1), 2)+(1:0)
  ppwts_2d = t(pwts[,cmb[1,]] * pwts[,cmb[2,]])
  qpsolve = function(...) quadprog::solve.QP(...)$solution

  q2 = opt_edge_lengths(ppwts_2d, ppinv, f3_jest, qpsolve)
  score = get_score(ppwts_2d, ppinv, f3_jest, q2)
  edges = as_tibble(as_edgelist(graph), .name_repair = ~c('from', 'to')) %>%
    mutate(type = 'edge', weight = c(q2))

  namedList(edges, score, opt = NULL)
}




