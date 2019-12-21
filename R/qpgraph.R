
#' A wrapper around qpGraph which requires an existing 'parfile' and 'graphfile'
#' @export
#' @param bin location of the qpGraph binary file.
#' @param parfile qpGraph parameter file.
#' @param graphfile qpGraph graph file.
#' @param outfile output file.
#' @param printonly should output be executed or the command just be printed?
#' @return a list of qpGraph output data.
#' @examples
#' \dontrun{
#' qpgraph_wrapper(bin = 'path/to/qpGraph',
#'                 pref = 'path/to/packedancestrymap_prefix',
#'                 parfile = 'path/to/parfile',
#'                 graphfile = 'path/to/graphfile')
#' }
qpgraph_wrapper = function(bin='./qpGraph', parfile='./parfile', graphfile='./graphfile', outfile='./out', printonly=FALSE) {
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



#' A wrapper around qpGraph which creates 'parfile' and 'graphfile'
#' @export
#' @param graph a graph as two column edge matrix .
#' @param outpop an outgroup population.
#' @param bin location of the qpGraph binary file.
#' @param pref prefix of the packedancestrymap format genotype files.
#' @param outdir output directory
#' @param printonly should output be executed or the command just be printed?
#' @param lambdascale lambdascale
#' @param lsqmode least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @param diag diag
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
#' qpgraph_wrapper2(graph1,
#'                  bin = 'path/to/qpGraph',
#'                  pref = 'path/to/packedancestrymap_prefix')
#' }
qpgraph_wrapper2 = function(graph, bin, pref, outpop='NULL', outdir='.', printonly=FALSE, lambdascale=-1, lsqmode='NO', diag=0.0001, hires='NO', forcezmode='NO', allsnps='NO', bigiter=100) {
  # wrapper around AdmixTools qpGraph
  # makes parfile and graphfile

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

  pops = union(edg[[1]], edg[[2]])
  simfile = tibble(v1 = c('root'), v2 = c('R'), v3='', v4='') %>%
    bind_rows(tibble(v1 = 'label', v2=leaves, v3=leaves, v4='')) %>%
    bind_rows(tibble(v1='edge', v2=paste0('e', 1:length(e1)), v3=e1, v4=e2)) %>%
    bind_rows(admix)

  pf = paste0(outdir, '/parfile')
  gf = paste0(outdir, '/graphfile')
  of = paste0(outdir, '/out')

  write(parfile, pf)
  simfile %>% write_tsv(gf, col_names=F)

  qpgraph_wrapper(bin=bin, parfile=pf, graphfile=gf, outfile=of, printonly=printonly)
}



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



#' @export
graph_to_pwts = function(grph) {
  # input igraph object
  # output: (numpops-1)*numedges matrix 'pwts' which indicates all paths from pop to outpop; admixture edges are mapped onto parent edges and weighted
  # assumes first pop is outpop, and 'R' is root

  leaves = V(grph)$name[degree(grph, v = V(grph), mode = c('out')) == 0]
  pwts = matrix(0, length(E(grph)), length(leaves))
  colnames(pwts) = leaves
  rownames(pwts) = attr(E(grph), 'vnames')

  admixnodes = which(degree(grph, mode='in')==2)
  admixedges = unlist(incident_edges(grph, admixnodes, mode='in'))

  allpaths = all_simple_paths(grph, V(grph)[1], leaves, mode='out')
  pathcounts = table(names(sapply(allpaths, function(x) tail(x,1))))
  for(i in seq_len(length(allpaths))) {
    target = names(tail(allpaths[[i]],1))
    ln = length(allpaths[[i]])
    pth2 = allpaths[[i]][c(1, 1+rep(seq_len(ln-2), each=2), ln)]
    rowind = as.vector(E(grph)[get.edge.ids(grph, pth2)])
    pwts[rowind,target] = pwts[rowind,target] + 1/pathcounts[target]
  }

  if(!is.null(admixedges)) pwts = pwts[-admixedges,]
  pwts[,-1] - pwts[,1]
}


#' @export
expand_path = function(path) {
  # igraph path (sequence of vertices)
  # duplicates inner vertices, so that this works with igraph functions that process vertex sequences pairwise
  ln = length(path)
  path[c(1, 1+rep(seq_len(ln-2), each=2), ln)]
}

#' @export
graph_to_weightind = function(grph) {
  # input igraph object
  # output:
  # map (leaf, edge) -> paths
  # map path -> weights
  # ultimately: indices for weights into paths, indices for paths into pwts

  leaves = V(grph)[degree(grph, v = V(grph), mode = c('out')) == 0]
  admixnodes = which(degree(grph, mode='in')==2)
  admixedges = unlist(incident_edges(grph, admixnodes, mode='in'))
  normedges = setdiff(1:length(E(grph)), admixedges)
  paths = all_simple_paths(grph, V(grph)[1], leaves[-1], mode='out')
  ends = sapply(paths, tail, 1)
  edge_per_path = paths %>% map(expand_path) %>% map(~get.edge.ids(grph, .))
  weight_per_path = edge_per_path %>% map(~(which(admixedges %in% .)))

  path_edge_table = do.call(rbind, lapply(seq_len(length(weight_per_path)), function(i) tibble(path=i, edge=c(edge_per_path[[i]])))) %>% mutate(edge2 = match(edge, normedges)) %>% filter(!is.na(edge2)) %>% mutate(leaf = as.vector(ends[path]), leaf2 = match(leaf, leaves[-1])) %>% left_join(enframe(c(table(ends))) %>% transmute(leaf=as.numeric(name), numpaths=value), by='leaf') %>% group_by(leaf2, edge2) %>% mutate(cnt = n(), keep = cnt < numpaths) %>% filter(keep) %>% as.matrix

  path_admixedge_table = do.call(rbind, lapply(seq_len(length(weight_per_path)), function(i) tibble(path=i, admixedge=c(weight_per_path[[i]])))) %>% as.matrix
  list(path_edge_table, path_admixedge_table, length(paths))
}



#' @export
fill_pwts = function(pwts, weights, path_edge_table, path_admixedge_table) {
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
  q2 = opt_edge_lengths(ppwts_2d, ppinv, f3_jest)
  w2 = (ppwts_2d %*% q2) - f3_jest
  lik = t(w2) %*% ppinv %*% w2
  lik[1,1]
}

#' @export
opt_edge_lengths = function(ppwts_2d, ppinv, f3_jest) {
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

  # very slow, don't use (takes 10 times as long):
  # pracma::quadprog(cc, q1, A=-diag(nc), b=rep(0, nc))$xmin
  # as fast as quadprogpp::QP.Solve, and on CRAN. switch to this maybe:
  -quadprog::solve.QP(cc, q1, -diag(nc), rep(0, nc))$solution
  #quadprogpp::QP.Solve(cc, q1, CI = diag(nc), ci0 = rep(0, nc))
}

#' @export
get_score = function(ppwts_2d, ppinv, f3_jest, q2) {
  w2 = (ppwts_2d %*% q2) - f3_jest
  lik = t(w2) %*% ppinv %*% w2
  lik[1,1]
}




#' Compute the fit of an admixturegraph.
#'
#' Computes the fit of an admixturegraph for a given graph topology and empirical f2-block-jackknife statistics. f2-block-jackknife statistics can be provided via the arguments \code{f2_blocks} and \code{block_lengths}, or via \code{f2_dir}.
#' @export
#' @param graph an admixture graph represented as a matrix of edges, an \code{\link{igraph}} object, or the path to a qpGraph graph file.
#' @param f2_blocks 3d array of block-jackknife leave-one-block-out estimates of f2 statistics. output of \code{\link{afs_to_f2_blocks}}. they are weighted by inverse of outgroup heterozygosity, if outgroup was specified.
#' @param block_lengths the jackknife block lengths used in computing the f2 statistics. see \code{\link{get_block_lengths}}.
#' @param f2_dir a directory with f2 statistics for each population pair in the graph. must contain 'block_lengths.RData'.
#' @param lsqmode least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @param fnscale try increasing or decreasing this, if optimization does not converge.
#' @param fudge try increasing this, if you get the error message \code{constraints are inconsistent, no solution!}.
#' @param numstart number of random initializations. defaults to 10 times the number of admixture nodes.
#' @param seed seed for generating starting weights.
#' @param verbose print optimization iterations
#' @return a list of qpGraph output data
#' @examples
#' qpgraph(graph, f2_blocks, block_lengths)
#' plot_graph(qpgraph(graph, f2_blocks, block_lengths)$edges)
qpgraph = function(graph, f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL, lsqmode=FALSE, fnscale=1e-6, fudge=1e-3, numstart=NULL, seed=NULL, verbose=FALSE, cpp=TRUE) {
  # modelled after AdmixTools qpGraph

  #----------------- process graph -----------------
  if('data.frame' %in% class(graph) | 'matrix' %in% class(graph)) {
    edges = as.matrix(graph)
  } else if(class(graph) == 'character') {
    edges = parse_qpgraph_graphfile(graph)
  } else if(class(graph)[1] == 'igraph') {
    edges = igraph::as_edgelist(graph)
  } else {
    stop('Cannot parse graph!')
  }
  grph = graph_from_edgelist(edges)
  nedges = length(E(grph))
  admixnodes = which(degree(grph, mode='in')==2)
  nadmix = length(admixnodes)
  admixedgesfull = sapply(seq_len(nadmix), function(i) incident_edges(grph, admixnodes, mode='in')[[i]][1:2])
  normedges = setdiff(1:nedges, admixedgesfull)

  pops = V(grph)$name[degree(grph, v = V(grph), mode = c('out')) == 0]

  #----------------- read f-stats -----------------
  if(is.null(f2_dir) & (is.null(f2_blocks) | is.null(block_lengths))) stop('You have to provide an f2_dir argument, or f2_blocks and block_lengths!')
  if(!is.null(f2_dir) & (is.null(f2_blocks) | is.null(block_lengths))) {
    f2_blocks = read_f2(f2_dir, pops = pops)
    load(paste0(f2_dir, '/block_lengths.RData'))
  }
  stopifnot(all(pops %in% dimnames(f2_blocks)[[1]]))

  popind = match(pops, dimnames(f2_blocks)[[1]])
  f2_blocks = rray(f2_blocks)
  npop = length(pops)
  npair = choose(npop, 2)
  cmb = combn(0:(npop-1), 2)+(1:0)

  #----------------- process f-stats -----------------
  f2 = bj_arr_stats(f2_blocks[popind,popind,], block_lengths)
  f2out = tibble(pop1=combn(pops, 2)[1,],
                 pop2=combn(pops, 2)[2,],
                 f2est = f2[[1]][lower.tri(f2[[1]])],
                 se = sqrt(f2[[2]][lower.tri(f2[[2]])]))
  # todo: write function to compute all pairwise fitted f2 stats from edge weights (if that's an interesting output)

  f3_blocks = (f2_blocks[, popind[1],] + f2_blocks[popind[1],,] - f2_blocks)/2
  f3_blocks_2d = arr3d_to_pairmat(f3_blocks[popind[-1], popind[-1],])
  sts = bj_mat_stats(f3_blocks_2d, block_lengths)
  f3_jest = sts[[1]]
  f3_jvar = sts[[2]]
  f3out = tibble(pop1=pops[cmb[1,]],
                 pop2=pops[cmb[2,]],
                 f3est = f3_jest, se = sqrt(diag(f3_jvar)))
  diag(f3_jvar) = diag(f3_jvar) + sum(diag(f3_jvar))*fudge

  #----------------- compute fit -----------------
  # in qpGraph fudge is 1e-5; sometimes quadprog doesn't converge unless this is larger; has large effect on magnitude of likelihood score
  if(lsqmode) ppinv = diag(1/diag(f3_jvar))
  else ppinv = solve(f3_jvar)

  weightind = graph_to_weightind(grph)
  weight = rep(NA, nedges)
  pwts = graph_to_pwts(grph)
  opt = NULL

  if(nadmix > 0) {
    if(is.null(numstart)) numstart = 10*nadmix
    set.seed(seed)
    parmat = matrix(runif(numstart*nadmix), numstart)
    if(verbose) alert_info(paste0('testing ', nrow(parmat), ' combinations of admixture weight starting values\n'))
    #arglist = list(pwts, ppinv, f3_jest, weightind[[1]], weightind[[2]], weightind[[3]], cmb, quadprogpp::QP.Solve); # need to change sign of third argument in qpsolve, if switching back to this
    arglist = list(pwts, ppinv, f3_jest, weightind[[1]], weightind[[2]], weightind[[3]], cmb, function(...) quadprog::solve.QP(...)$solution)
    if(cpp) optimweightsfun = cpp_optimweightsfun
    opt = multistart(parmat, optimweightsfun, args=arglist, method='L-BFGS-B',
                     lower=0, upper=1, control=list(maxit=1e4, fnscale=fnscale), verbose=verbose)

    best = opt %>% top_n(1, -.data$value)
    opt = data.frame(parmat, opt, stringsAsFactors = F)

    wts = as.matrix(best[,1:nadmix])[1,]
    weight[admixedgesfull[1,]] = wts
    weight[admixedgesfull[2,]] = 1-wts
    pwts = fill_pwts(pwts, wts, weightind[[1]], weightind[[2]])
  }

  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])
  q2 = opt_edge_lengths(ppwts_2d, ppinv, f3_jest)
  score = get_score(ppwts_2d, ppinv, f3_jest, q2)
  weight[normedges] = q2
  edges = as_tibble(edges) %>% set_colnames(c('from', 'to')) %>%
    mutate(type = ifelse(1:n() %in% normedges, 'edge', 'admix'), weight=weight, label=round(weight, 2))

  namedList(edges, score, f2=f2out, f3=f3out, opt)
}

#' @export
qpgraph_slim = function(grph, f3_jest, ppinv, pops, fnscale=1e-6, numstart=10, seed=NULL, verbose=FALSE, cpp=TRUE) {
  # modelled after AdmixTools qpGraph
  # optimised for testing many topologies for a given set of populations
  # graph is igraph

  nedges = length(E(grph))
  admixnodes = which(degree(grph, mode='in')==2)
  nadmix = length(admixnodes)
  admixedgesfull = sapply(seq_len(nadmix), function(i) incident_edges(grph, admixnodes, mode='in')[[i]][1:2])
  normedges = setdiff(1:nedges, admixedgesfull)

  npop = length(pops)
  cmb = combn(0:(npop-1), 2)+(1:0)
  graphpops = V(grph)$name[degree(grph, v = V(grph), mode = c('out')) == 0]
  popind = setdiff(match(graphpops, pops), 1)
  orig_order = apply(cmb+1, 2, paste0, collapse='')
  new_order = apply(matrix(popind[c(cmb)], 2), 2, function(x) paste0(sort(x), collapse=''))
  pairmatch = match(new_order, orig_order)
  f3_jest = f3_jest[pairmatch]
  ppinv = ppinv[pairmatch, pairmatch]
  stopifnot(all(!is.na(ppinv)))

  weightind = graph_to_weightind(grph)
  weight = rep(NA, nedges)
  pwts = graph_to_pwts(grph)
  opt = NULL

  if(nadmix > 0) {
    set.seed(seed)
    parmat = matrix(runif(numstart*nadmix), numstart)

    arglist = list(pwts, ppinv, f3_jest, weightind[[1]], weightind[[2]], weightind[[3]], cmb, quadprogpp::QP.Solve)
    if(cpp) optimweightsfun = cpp_optimweightsfun
    opt = multistart(parmat, optimweightsfun, args=arglist, method='L-BFGS-B',
                          lower=0, upper=1, control=list(maxit=1e4, fnscale=fnscale), verbose=verbose)


    best = opt %>% top_n(1, -.data$value)
    opt = data.frame(parmat, opt, stringsAsFactors = F)

    wts = as.matrix(best[,1:nadmix])[1,]
    weight[admixedgesfull[1,]] = wts
    weight[admixedgesfull[2,]] = 1-wts
    pwts = fill_pwts(pwts, wts, weightind[[1]], weightind[[2]])
  }

  ppwts_2d = t(pwts[,cmb[1,]]*pwts[,cmb[2,]])
  #q2 = cpp_opt_edge_lengths(ppwts_2d, ppinv, f3_jest, quadprogpp::QP.Solve)
  q2 = opt_edge_lengths(ppwts_2d, ppinv, f3_jest)
  score = get_score(ppwts_2d, ppinv, f3_jest, q2)
  weight[normedges] = q2
  edges = as_tibble(igraph::as_edgelist(grph)) %>% set_colnames(c('from', 'to')) %>%
    mutate(type = ifelse(1:n() %in% normedges, 'edge', 'admix'), weight=weight, label=round(weight, 2))

  namedList(edges, score, opt)
}


#' @export
qpgraph_precompute = function(pops, f2_blocks, block_lengths, lsqmode=FALSE, fudge=1e-3) {
  # returns list of f3_jest and ppinv for subset of populations.
  # f3_jest and ppinv are required for qpgraph_slim; f2out and f3out are extra output
  # f2_blocks may contain more populations than the ones used in qpgraph
  # f2_blocks input here should be subset which is used by qpgraph function

  stopifnot(all(pops %in% dimnames(f2_blocks)[[1]]))
  f2_blocks = rray(f2_blocks[pops, pops, ])

  npop = length(pops)
  npair = choose(npop, 2)
  cmb = combn(0:(npop-1), 2)+(1:0)

  # process f-stats
  f2 = bj_arr_stats(f2_blocks, block_lengths)
  f2out = tibble(pop1=combn(pops, 2)[1,],
                 pop2=combn(pops, 2)[2,],
                 f2est = f2[[1]][lower.tri(f2[[1]])],
                 se = sqrt(f2[[2]][lower.tri(f2[[2]])]))

  f3_blocks = (f2_blocks[,1,] + f2_blocks[1,,] - f2_blocks)/2
  f3_blocks_2d = arr3d_to_pairmat(f3_blocks[-1,-1,])
  sts = bj_mat_stats(f3_blocks_2d, block_lengths)
  f3_jest = sts[[1]]
  f3_jvar = sts[[2]]
  f3out = tibble(pop1=pops[cmb[1,]],
                 pop2=pops[cmb[2,]],
                 f3est = f3_jest, se = sqrt(diag(f3_jvar)))
  diag(f3_jvar) = diag(f3_jvar) + sum(diag(f3_jvar))*fudge
  # in qpGraph fudge is 1e-5; sometimes quadprog doesn't converge unless this is larger; has large effect on magnitude of likelihood score
  if(lsqmode) ppinv = diag(1/diag(f3_jvar))
  else ppinv = solve(f3_jvar)
  #pairnam = apply(matrix(pops[cmb+1], 2), 2, paste0, collapse=' ')
  #names(f3_jest) = colnames(ppinv) = rownames(ppinv) = pairnam

  namedList(f3_jest, ppinv, f2out, f3out)
}

