
# number of all binary trees with a given number of leaf nodes
numtrees = function(n) factorial(2*n-2)/(2^(n-1)*factorial(n-1))

# number of all binary tree topologies
numtreestop = function(n) factorial(2*n)/factorial(n+1)/factorial(n)

# number of possible DAGs
numdags = function(n) {
  if(n <= 1) return(1)
  sum(sapply(1:n, function(k) (-1)^(k+1) * choose(n, k) * 2^(k*(n-k)) * numdags(n-k)))
}

# number of trees with x admiture events
numtreesadmix = function(n, nadmix) numtrees(n) * numadmixplacements(2*n-2, nadmix)

# number of unique ways to add 'nadmix' undirectd edges; each added edge increases the number of edges by 3
numadmixplacements = function(numedges, nadmix) {
  if(nadmix == 0) return(1)
  choose(numedges, 2) * numadmixplacements(numedges+3, nadmix-1)
}


numadmix = function(grph) {
  sum(degree(grph, mode='in') == 2)
}

#' @export
get_leafnames = function(grph) {
  grph %>% V %>% {names(which(degree(grph, ., mode='out') == 0))}
}

get_leaves = function(grph) {
  grph %>% V %>% {.[degree(grph, ., mode='out') == 0]}
}

get_leaves2 = function(grph) {
  # uses 'subcomponent' to return leaves in a consistent order
  grph %>% subcomponent(V(.)[1], mode='out') %>% igraph::intersection(get_leaves(grph))
}

unify_vertex_names_rec = function(grph, node, vnamemap) {
  # recursive function which returns a vector of new, unique node names for all nodes which are reachable from 'node'
  children = neighbors(grph, node, mode='out')
  oldnam = names(node)
  if(length(children) == 0) {
    vnamemap[oldnam] = oldnam
    return(vnamemap)
  }

  newnam = rep(NA, length(children))
  for(i in seq_along(children)) {
    vnamemap = unify_vertex_names_rec(grph, children[i], vnamemap)
    newnam[i] = vnamemap[names(children[i])]
  }
  stopifnot(all(!is.na(newnam)))
  vnamemap[oldnam] = paste0(sort(newnam), '.', collapse='_')
  vnamemap
}

shortest_unique_prefixes = function(strings) {
# given a vector of strings, return a vector of the same length with the shortest unique prefixes
  len = rep(0, length(strings))
  for(i in 1:max(nchar(strings))) {
    pref = str_sub(strings, 1, i)
    tab = which(pref %in% names(which(table(pref) == 1)))
    tab = intersect(tab, which(len == 0))
    len[tab] = i
  }
  str_sub(strings, 1, len)
}

unify_vertex_names = function(grph, keep_unique=TRUE) {
  # this is an aweful function I wrote when I was very tired which changes the vertex names of inner vertices, so that all isomorphic graphs with the same leaf nodes have equally labelled inner nodes
  leaves = get_leaves(grph)
  lv = shortest_unique_prefixes(names(leaves))
  g = set_vertex_attr(grph, 'name', leaves, lv)
  ormap = names(V(g))
  names(ormap) = ormap

  nammap = unify_vertex_names_rec(g, V(g)[1], ormap)
  nammap[lv] = names(leaves)
  names(nammap)[match(lv, names(nammap))] = names(leaves)
  nammap['R'] = 'R'
  if(!keep_unique) {
    changed = setdiff(names(nammap), c('R', names(leaves)))
    nammap[changed] = paste0('n', as.numeric(as.factor(nammap[changed])))
  }
  set_vertex_attr(grph, 'name', names(nammap), nammap)
}


#' Generate a random admixture graph
#'
#' This function randomly generates an admixture graph for a given set of leaf nodes
#' @export
#' @param leaves the names of the leaf nodes, or a number specifying how many leaf nodes there should be
#' @param numadmix the number of admixture events
#' @param simple should edges leading to admixture nodes consist of separate admix edges and normal edges
#' @param outpop the outgroup population
#' @examples
#' rand_graph = random_admixturegraph(10, numadmix = 5)
#' plot_graph(rand_graph)
random_admixturegraph = function(leaves, numadmix = 0, simple = FALSE, outpop = NA) {
  # makes a random admixture graph
  # returns an 'igraph' graph object
  # 'leaves' can be a number of leaf nodes, or a character vector of leaf names
  stopifnot(class(leaves)[1] %in% c('numeric', 'character'))
  if(length(leaves) == 1) leaves = paste0('l', 1:leaves)
  if(is.na(outpop)) outpop = sample(leaves, 1)
  newick = random_newick(sample(setdiff(leaves, outpop)), outpop = outpop)
  grph = graph_from_edgelist(newick_to_edges(newick))
  grph = insert_admix_igraph_random(grph, numadmix)
  if(simple) return(grph)
  desimplify_graph(grph)
}


#' Remove redundant edges
#'
#' Nodes with in and out degree of one will be removed.
#' @param grph an igraph object
#' @export
#' @examples
#' simple = simplify_graph(igraph1)
#' plot_graph(igraph1)
#' plot_graph(simple)
simplify_graph = function(grph) {
  # removes redundant nodes
  convmat = FALSE
  if(class(grph) == 'matrix') {
    convmat = TRUE
    grph = graph_from_edgelist(grph)
  }
  repeat({
    redundant = which(degree(grph, mode='in') == 1 & degree(grph, mode='out') == 1)
    if(length(redundant) == 0) break
    newfrom = names((neighbors(grph, redundant[1], mode='in')))
    newto = names((neighbors(grph, redundant[1], mode='out')))
    grph = igraph::delete_vertices(grph, redundant[1])
    grph = igraph::add_edges(grph, c(newfrom, newto))
  })
  if(convmat) grph = igraph::as_edgelist(grph)
  grph
}

#' Add two nodes before each admixture node
#'
#' This is used to revert simplify_graph.
#' @export
#' @param grph an admixture graph
#' @examples
#' simple = simplify_graph(igraph1)
#' desimple = desimplify_graph(simple)
#' plot_graph(simple)
#' plot_graph(desimple)
desimplify_graph = function(grph) {
  # adds two nodes before each admixture node
  convmat = FALSE
  if(class(grph)[1] == 'matrix') {
    convmat = TRUE
    grph = graph_from_edgelist(grph)
  }
  admix = which(degree(grph, mode='in') == 2)
  if(length(admix) == 0) return(grph)
  parents = names(unlist(unname(lapply(admix, function(v) neighbors(grph, v, mode='in')))))
  admixnamesrep = rep(names(admix), each=2)
  newnam = paste0(admixnamesrep, c('a', 'b'))
  grph = igraph::add_vertices(grph, length(admix)*2, name=newnam)
  newedges = c(rbind(parents, newnam, newnam, admixnamesrep))
  deledges = paste(parents, admixnamesrep, sep='|')
  grph = igraph::add_edges(grph, newedges)
  grph = igraph::delete_edges(grph, deledges)
  if(convmat) grph = igraph::as_edgelist(grph)
  grph
}

split_graph = function(grph) {
  # removes admixture nodes from an admixturegraph
  # returns a tree and a list of edges that can be used to reconstruct the original admixturegraph
  admixedges = list()
  repeat({
    admix = names(which(degree(grph, mode='in') == 2))
    if(length(admix) == 0) break
    parents = names(neighbors(grph, admix[1], mode='in'))
    parents = sample(parents)
    # stopifnot(!any(c(parents[1], admix[1]) %in% unlist(admixedges)))
    # too strict; probably won't be able to always reintroduce admixedges at appropriate locations after SPR
    fromedge1 = names(neighbors(grph, parents[1], mode = 'in'))
    fromedge2 = setdiff(names(neighbors(grph, parents[1], mode = 'out')), admix[1])
    toedge1 = setdiff(names(neighbors(grph, admix[1], mode = 'in')), parents[1])
    toedge2 = names(neighbors(grph, admix[1], mode = 'out'))
    admixedges[[length(admixedges)+1]] = c(fromedge1, fromedge2, toedge1, toedge2)
    grph = igraph::delete_vertices(grph, c(admix[1], parents[1]))
    grph = igraph::add_edges(grph, c(fromedge1, fromedge2, toedge1, toedge2))
  })
  list(grph, admixedges)
}




#' Generate a random, binary graph with n terminal nodes
#' @export
#' @param n the number of terminal nodes, or a vector of population labels.
#' @param start prefix.
#' @param end postfix.
#' @param outpop outgroup (if \code{n} is a vector of labels).
#' @return tree in newick format.
#' @examples
#' random_newick(5)
#' random_newick(c('a', 'b', 'c', 'd')) # toplogy random, but pop order fixed
#' random_newick(sample(c('a', 'b', 'c', 'd'))) # toplogy and pop order random
random_newick = function(n, start='', end='', outpop=NULL) {
  # recursive function which returns topology of a random, binary tree in newick format
  # redirects to 'random_newick_named' when length(n) > 1
  if(length(n) > 1) {
    if(!is.null(outpop)) {
      start = paste0(start, '(', outpop, ',')
      end = paste0(')', end)
      n = setdiff(n, outpop)
    }
    return(random_newick_named(n, start, end))
  }
  if(n == 1) return('')
  n1 = sample(n-1,1)
  n2 = n-n1
  return(paste0(start, '(',random_newick(n1, ''),',',random_newick(n2, ''),')', end))
}

random_newick_named = function(names, start='', end='') {
  # recursive function which returns a labelled random, binary tree in newick format with named leaves in order of input
  n = length(names)
  if(n == 1) return(names)
  n1 = sample(n-1,1)
  return(paste0(start, '(',random_newick_named(names[1:n1]),',',random_newick_named(names[(n1+1):n]),')', end))
}

#' Turn a newick format tree to a matrix of edges
#' @export
#' @param newick tree in newick format.
#' @param node root label of the tree.
#' @param edgemat argument used for recursive function calls.
#' @return tree as two column matrix of edges (adjacency list)
#' @examples
#' newick = random_newick(c('a', 'b', 'c', 'd'))
#' newick
#' newick_to_edges(newick)
newick_to_edges = function(newick, node='R', edgemat=matrix(NA,0,2)) {
  # turns binary tree in newick format into matrix of edges (adjacency list)

  newick = gsub('^\\(', '', gsub('\\)$', '', gsub(';$', '', newick)))
  opencount = 0
  for(i in 1:nchar(newick)) {
    char = substr(newick, i,i)
    if(char == '(') opencount = opencount+1
    if(char == ')') opencount = opencount-1
    if(char == ',' && opencount == 0) break
  }
  left = substr(newick, 1, i-1)
  right = substr(newick, i+1, nchar(newick))
  stopifnot(str_count(newick, ',') >= 1)
  if(str_count(left, ',') == 0) {
    nodel = left
    edgesleft = matrix(NA, 0, 2)
  } else {
    nodel = paste0(node, 'l')
    edgesleft = newick_to_edges(left, nodel)
  }
  if(str_count(right, ',') == 0) {
    noder = right
    edgesright = matrix(NA, 0, 2)
  } else {
    noder = paste0(node, 'r')
    edgesright = newick_to_edges(right, noder)
  }
  rbind(c(node, nodel), edgesleft, c(node, noder), edgesright, edgemat)
}

insert_admix_igraph = function(grph, admixedges, substitute_missing=FALSE) {
  # inverse of 'split_graph'
  # only uses 'fromedge2' and 'toedge2' (2nd and 4th elements of admixedges)
  # stopifnot(all(c(sapply(admixedges, `[`, c(2, 4))) %in% names(V(grph))))
  # some nodes will get lost when splitting; insert random amix edges instead
  miss = 0
  for(i in rev(seq_len(length(admixedges)))) {
    fromedge2 = admixedges[[i]][2]
    toedge2 = admixedges[[i]][4]
    if(!fromedge2 %in% names(V(grph)) ||
       !toedge2 %in% names(V(grph))) {
      miss = miss+1
      next
    }
    toedge2parent = neighbors(grph, toedge2, mode='in')
    toedge2sib = setdiff(names(neighbors(grph, toedge2parent, mode='out')), toedge2)
    if(degree(grph, fromedge2, mode='in') > 1 ||
       degree(grph, toedge2sib, mode='in') > 1) {
      miss = miss+1
      next
    }
    fromedge1 = names(neighbors(grph, fromedge2, mode='in'))
    toedge1 = names(neighbors(grph, toedge2, mode='in'))
    grph = igraph::delete_edges(grph, c(paste(fromedge1, fromedge2, sep='|'), paste(toedge1, toedge2, sep='|')))
    new1 = paste0(fromedge1, 'x')
    while(new1 %in% names(V(grph))) new1 = paste0(new1, sample(letters, 1))
    j = i
    new2 = paste0('admix', j)
    while(new2 %in% names(V(grph))) {j=j+1; new2 = paste0('admix', j)}
    grph = igraph::add_vertices(grph, 2, name=c(new1, new2))
    newedges = c(fromedge1, new1, new1, fromedge2, toedge1, new2, new2, toedge2, new1, new2)
    grph = igraph::add_edges(grph, newedges)
    stopifnot(max(degree(grph)) <= 3)
  }
  if(substitute_missing) grph = insert_admix_igraph_random(grph, miss)
  else if (miss > 0) warning('did not insert or substitute admixedge')
  grph
}


insert_admix_igraph_random = function(grph, nadmix) {
  # grph should be simplified
  for(i in seq_len(nadmix)) {
    firstgen = neighbors(grph, V(grph)[1], mode='out')
    admix = V(grph)[(degree(grph, mode='in') > 1)]
    admixchildren = do.call(c, igraph::ego(grph, 1, admix, mode='out'))
    admixparents = do.call(c, igraph::ego(grph, 1, admix, mode='in'))
    admixsiblings = do.call(c, igraph::ego(grph, 1, admixparents, mode='out'))
    exclude = c(admix, firstgen)
    if(!is.null(admixchildren)) exclude = c(exclude, admixchildren)
    cand_from = sample(igraph::difference(V(grph)[-1], exclude))
    for(fromnode in names(cand_from)) {
      upstream = subcomponent(grph, fromnode, mode='in')
      downstream1 = neighbors(grph, fromnode, mode='out')
      siblings = neighbors(grph, upstream[2], mode='out') # upstream[2] should be parent of fromnode
      uncles = neighbors(grph, neighbors(grph, upstream[2], mode='in'), mode='out')
      exclude = c(upstream, downstream1, siblings, uncles)
      if(!is.null(admixsiblings)) exclude = c(exclude, admixsiblings)
      cand_to = igraph::difference(cand_from, exclude)
      if(length(cand_to) > 0) break
    }
    if(length(cand_from) == 0 | length(cand_to) == 0) {
      warning(paste0('could only insert ', i-1, ' admixture nodes'))
      break
    }
    tonode = sample(names(cand_to), 1)
    stopifnot(all(c(fromnode, tonode) %in% names(V(grph))))
    grphnew = insert_admix_igraph(grph, list(c(NA, fromnode, NA, tonode)), substitute_missing=FALSE)
    stopifnot(numadmix(grphnew) > numadmix(grph))
    stopifnot(igraph::is_simple(grphnew))
    grph = grphnew
  }
  grph
}


subtree_prune_and_regraft = function(grph, only_leaves=FALSE) {
  # cuts of a parts of the tree and attaches it at a random location
  # root -> outgroup stays fixed
  stopifnot(degree(grph, V(grph)[1], mode='in') == 0)
  stopifnot(max(degree(grph, V(grph), mode='in')) == 1)

  firstgen = neighbors(grph, V(grph)[1], mode='out')
  repeat({
    excl = firstgen
    if(only_leaves) excl = c(excl, V(grph)[degree(grph, mode='out') > 0])
    cutnode = sample(igraph::difference(V(grph)[-1], excl), 1)
    cutnodes = subcomponent(grph, cutnode, mode='out')
    cutparent = neighbors(grph, cutnode, mode='in')
    cutgrandparent = neighbors(grph, cutparent, mode='in')
    cutsibling = igraph::difference(neighbors(grph, cutparent, mode='out'), cutnode)
    hostnodes = igraph::difference(V(grph)[-1], c(cutnodes, cutparent, cutsibling, firstgen))
    if(length(hostnodes) > 0) break
  })

  hostnode = sample(names(hostnodes), 1)
  hostparent = names(neighbors(grph, hostnode, mode='in'))

  newnam = paste(hostnode, sample(letters, 1), sep='_')
  while(newnam %in% names(V(grph))) newnam = paste0(newnam, sample(letters, 1))

  grph %<>% igraph::add_vertices(1, name=newnam) %>%
    igraph::add_edges(c(names(cutgrandparent), names(cutsibling),
                hostparent, newnam,
                newnam, hostnode,
                newnam, names(cutnode))) %>%
    igraph::delete_vertices(names(cutparent)) %>%
    igraph::delete_edges(paste(hostparent, hostnode, sep='|'))
  stopifnot(igraph::is_simple(grph))
  grph
}


admixturegraph_prune_and_regraft = function(grph, desimplify=TRUE, only_leaves = FALSE) {
  # 1. remove admixture edges (one randomly selected from each admixture node)
  # 2. runs subtree_prune_and_regraft on resulting tree
  # 3. add admixture edges back on
  grph = simplify_graph(grph)
  spl = split_graph(grph)
  grph = subtree_prune_and_regraft(spl[[1]], only_leaves = only_leaves)
  grph = insert_admix_igraph(grph, spl[[2]], substitute_missing = TRUE)
  stopifnot(igraph::is_simple(grph))
  if(desimplify) grph = desimplify_graph(grph)
  grph
}


move_admixedge_once = function(grph, desimplify=TRUE) {
  # selects random admixture edge, and moves it to next closest possible spot
  # if not possible, select other admix node
  # if none possible, print warning
  grph = simplify_graph(grph)
  admix = sample(names(V(grph)[(degree(grph, mode='in') > 1)]))
  nadmix = length(admix)
  firstgen = names(neighbors(grph, V(grph)[1], mode='out'))
  for(i in seq_len(nadmix)) {
    parents = sample(names(neighbors(grph, admix[i], mode='in')))
    for(j in 1:2) {
      sibling = setdiff(names(neighbors(grph, parents[j], mode='out')), admix[i])
      grandparent = names(neighbors(grph, parents[j], mode='in'))
      for(n in names(subcomponent(grph, parents[j]))) {
        # subcomponent will return nodes ordered by distance to n
        if(!n %in% names(V(grph)[1]) &
           !n %in% firstgen &
           !n %in% parents &
           !n %in% admix &
           !n %in% c(names(neighbors(grph, parents[1], mode='out')),
                     names(neighbors(grph, parents[2], mode='out'))) &
           !n %in% names(subcomponent(grph, admix[i], mode='out'))) {
          newgrandparent = names(neighbors(grph, n, mode='in'))
          grph = igraph::add_edges(grph, c(grandparent, sibling, newgrandparent, parents[j], parents[j], n))
          grph = igraph::delete_edges(grph, c(paste(newgrandparent, n, sep='|'), paste(grandparent, parents[j], sep='|'), paste(parents[j], sibling, sep='|')))
          stopifnot(igraph::is_simple(grph))
          if(desimplify) grph = desimplify_graph(grph)
          return(grph)
        }
      }
    }
  }
  if(nadmix > 0) warning('no suitable attachment points found for admixture edge')
  if(desimplify) grph = desimplify_graph(grph)
  return(grph)
}

permute_leaves = function(grph, fix_outgroup=TRUE) {

  leaves = V(grph)[degree(grph, v = V(grph), mode = c('out')) == 0]
  if(fix_outgroup) leaves = leaves[-1]
  nam = names(leaves)
  igraph::set_vertex_attr(grph, 'name', leaves, sample(nam))
}

swap_leaves = function(grph, fix_outgroup=TRUE) {

  leaves = V(grph)[degree(grph, v = V(grph), mode = c('out')) == 0]
  if(fix_outgroup) leaves = leaves[-1]
  leaves = sample(leaves)[1:2]
  nam = names(leaves)
  igraph::set_vertex_attr(grph, 'name', leaves, sample(nam))
}


add_generation = function(models, numsel, qpgfun) {

  lastgen = max(models$generation)
  numpergen = models %>% filter(generation == lastgen) %>% nrow
  numeach = round(numpergen/numsel)
  stopifnot(numpergen > numsel)
  oldmodels = models %>% filter(generation == lastgen)
  winners = oldmodels %>% mutate(prob = score_to_prob(score)) %>% sample_n(numsel-1, weight=prob)
  winners %<>% bind_rows(oldmodels %>% top_n(1, -jitter(score)))
  newmodels = tibble(generation = lastgen+1, index = 1:numpergen, igraph = rep(winners$igraph, numeach), oldscore = rep(winners$score, numeach))
  mutations = c(replicate(numsel, namedList(identity)), sample_mutations(numpergen-numsel))
  newmodels %<>% mutate(mutation = names(mutations), igraph = imap(igraph, ~exec(mutations[[.y]], .x)))
  newmodels %<>% mutate(out = furrr::future_map(igraph, qpgfun)) %>% unnest_wider(out)
  #newmodels %<>% mutate(out = map(igraph, qpgfun)) %>% unnest_wider(out)
  models %>% bind_rows(newmodels)
}

sample_mutations = function(n) {
  # returns a random mutation function
  mutations = list(sprleaves = function(x) (admixturegraph_prune_and_regraft(x, only_leaves = TRUE)),
                   sprall = function(x) (admixturegraph_prune_and_regraft(x, only_leaves = FALSE)),
                   #permleaves = permute_leaves,
                   swapleaves = swap_leaves,
                   moveadmix = move_admixedge_once)
  sample(mutations, n, replace = TRUE)
}


score_to_prob = function(score, fix_best=TRUE) {
  logscore = log10(score)
  sc = 1-logscore/max(logscore)
  sc = sc^20
  sc = sc/sum(sc)
}


evolve_topology = function(init, numgen, numsel, qpgfun, verbose=TRUE) {

  out = init
  for(i in seq_len(numgen)) {
    out = add_generation(out, numsel, qpgfun)

    if(verbose) {
      thisgen = out %>% filter(generation == i)
      cat(paste0('Generation ', i, '\n'))
      cat(paste0('  Selected scores: ', paste(round(thisgen %>% filter(index <= numsel) %>% top_n(numsel, -jitter(score)) %$% score), collapse=', '), '\n'))
      cat(paste0('  Best new scores: ', paste(round(thisgen %>% filter(index > numsel) %>% top_n(numsel, -jitter(score)) %$% score), collapse=', '), '\n'))
      cat(paste0('  Median score: ', round(median(thisgen$score)), '\n'))
    }
  }
  out
}


#' Find well fitting admixturegraphs
#'
#' This function generates and evaluates admixturegraphs in \code{numgen} iterations to find well fitting admixturegraphs.
#' @export
#' @param pops population names for which to fit admixturegraphs
#' @param f3_jest 3d array with all pairwise f3 estimates with the outgroup
#' @param ppinv inverse of the f3 block jackknife covariance matrix
#' @param numgraphs number of graphs in each generation
#' @param numgen number of generations
#' @param numsel number of graphs which are selected in each generation
#' @param numadmix number of admixture events within each graph
#' @param outpop outgroup population
#' @return a nested data frame with one model per line
#' @examples
#' pops = get_leafnames(igraph1)
#' optimize_admixturegraph_single(igraph1, f3_jest, ppinv, numgraphs=10, numgen=2, numsel=2, numadmix=1)
optimize_admixturegraph_single = function(pops, f3_jest, ppinv, numgraphs=50, numgen=5, numsel=5, numadmix=0, outpop=NA, verbose=TRUE, cpp=TRUE, debug=FALSE, ...) {

  qpgfun = function(grph) qpgraph_slim(grph, f3_jest, ppinv, pops, fnscale = 1e-4, cpp=cpp, verbose=FALSE, ...)
  if(!debug) qpgfun = possibly(qpgfun, otherwise=NULL)

  initigraphs = replicate(numgraphs, random_admixturegraph(pops, numadmix, outpop=outpop), simplify = F)
  init = tibble(generation=0, index = seq_len(numgraphs), igraph = initigraphs, mutation = 'random_admixturegraph') %>%
    mutate(out = map(igraph, qpgfun)) %>% unnest_wider(out) %>% mutate(oldscore = score)

  evolve_topology(init, numgen, numsel, qpgfun, verbose=verbose)
}

#' Find a well fitting admixturegraph
#'
#' This function generates and evaluates admixturegraphs in \code{numgen} iterations across \code{numrep} independent repeats to find well fitting admixturegraphs. It uses the function \code{\link{future_map}} from the \code{\link{furrr}} package to parallelize across the independent repeats. The function \code{\link{future::plan}} can be called to specify the details of the parallelization. This can be used to parallelize across cores or across nodes on a compute cluster.
#' @export
#' @param pops population names for which to fit admixturegraphs
#' @param f3_jest 3d array with all pairwise f3 estimates with the outgroup
#' @param ppinv inverse of the f3 block jackknife covariance matrix
#' @param numgraphs number of graphs in each generation
#' @param numgen number of generations
#' @param numsel number of graphs which are selected in each generation
#' @param numadmix number of admixture events within each graph
#' @param outpop outgroup population
#' @return a nested data frame with one model per line
#' @examples
#' \dontrun{
#' optimize_admixturegraph(pops, f3_jest, ppinv, numrep = 200, numgraphs = 100,
#'                         numgen = 20, numsel = 5, numadmix = 3)
#' }
optimize_admixturegraph = function(pops, f3_jest, ppinv, numrep=10, numgraphs=50, numgen=5, numsel=5, numadmix=0, outpop=NA, verbose=TRUE, cpp=TRUE, debug=FALSE, ...) {

  if(is.na(outpop)) outpop = pops[1]
  else if(!outpop %in% pops) pops = c(outpop, pops)

  oa = function(x) optimize_admixturegraph_single(pops, outpop=outpop, f3_jest, ppinv,
                                                  numgen = numgen, numsel = numsel,
                                                  numgraphs = numgraphs, numadmix = numadmix,
                                                  verbose = verbose, cpp = cpp, debug = debug, ...)
  if(!debug) oa = possibly(oa, otherwise=NULL)
  res = furrr::future_map(as.list(seq_len(numrep)), oa)
  bind_rows(res, .id='run') %>% mutate(run = as.numeric(run))
}


summarize_graph = function(grph, exclude_outgroup = TRUE) {

  leaves = V(grph)[degree(grph, v = V(grph), mode = c('out')) == 0]
  if(exclude_outgroup) leaves = igraph::difference(leaves, V(grph)[2])

  paths = all_simple_paths(grph, V(grph)[1], leaves, mode = 'out') %>% enframe(name = 'pathnum', value='path') %>% mutate(name1 = map_chr(path, ~(tail(names(.), 1))))

  get_iap = function(x, y) suppressWarnings(names(tail(which(cumsum(x != y) == 0), 1)))

  paths2 = paths %>% rename(path2 = path, pathnum2 = pathnum, name2 = name1)
  paths3 = paths %>% rename(path3 = path, pathnum3 = pathnum, name3 = name1)

  tripletopo = expand_grid(paths, paths2, paths3) %>%
    filter(pathnum != pathnum2, pathnum != pathnum3, pathnum2 != pathnum3,
           name1 != name2, name1 != name3, name2 != name3) %>%
    mutate(iap12 = map2_chr(path, path2, get_iap),
           iap13 = map2_chr(path, path3, get_iap),
           iap23 = map2_chr(path2, path3, get_iap)) %>%
    mutate(topo = case_when(iap13 == iap23 ~ 0, iap12 == iap23 ~ 1, iap12 == iap13 ~ 2),
           top = case_when(topo == 0 ~ iap13, topo == 1 ~ iap12, topo == 2 ~ iap12)) %>%
    group_by(name1, name2, name3, top) %>%
    summarize(x13 = any(topo == 1) & !any(topo == 0),
              x23 = any(topo == 2) & !any(topo == 0),
              x31 = any(topo == 1) & any(topo == 0),
              x32 = any(topo == 2) & any(topo == 0),
              x12 = any(topo == 0) & !any(topo == 1),
              x21 = any(topo == 0) & any(topo == 1)) %>%
    group_by(name1, name2, name3) %>%
    summarize(x13 = any(x13), x23 = any(x23), x31 = any(x31), x32 = any(x32), x12 = any(x12), x21 = any(x21),
              toposet = paste0(x13+0, x23+0, x31+0, x32+0, x12+0, x21+0, collapse='')) %>%
    ungroup %>%
    mutate(tlr = !x13 & !x31 & (x23 | x32) & (x12 | x21))

  tripletopo
}

#' Summarize triples across graphs
#'
#' This summarizes topologies of population triples across graphs
#'
#' @export
#' @param results the output of \code{\link{optimize_admixturegraph}}
#' @param maxscore restrict summary to graphs with score not larger than \code{maxscore}
#' @examples
#' \dontrun{
#' summarize_triples(opt_results)
#' }
summarize_triples = function(results, maxscore=NA) {
  # results is output from 'optimize_admixturegraph'
  # takes at most one graph from each independent run

  sel = results %>% group_by(run) %>% top_n(1, -jitter(score)) %>% ungroup
  if(!is.na(maxscore)) sel %<>% filter(score <= maxscore)

  sel %>% mutate(topo = map(igraph, summarize_graph)) %>%
    select(run, generation, index, igraph, score, topo) %>%
    unnest(topo) %>%
    group_by(name1, name2, name3, toposet) %>%
    mutate(cnt = n()) %>%
    group_by(name1, name2, name3) %>%
    summarize(numgraphs = n(),
              x13 = mean(x13),
              x23 = mean(x23),
              x31 = mean(x31),
              x32 = mean(x32),
              clade = mean(substr(toposet, 1, 4) == '0000'),
              x12 = mean(x12),
              x21 = mean(x21),
              toptopo = toposet[cnt = max(cnt)][1],
              toptopocnt = max(cnt),
              topos = list(setNames(cnt, toposet))) %>% ungroup
}

#' @export
isomorphism_classes = function(igraphlist) {
  # retuns integer vector with the same length as 'igraphlist', which assigns each graph to a class
  # only considers topology

  numgraph = length(igraphlist)
  if(numgraph == 0) return(numeric())
  if(numgraph == 0) return(1)

  leaflist = map(igraphlist, ~names(get_leaves(.)))

  cmb = combn(numgraph, 2)
  iso = map2_lgl(cmb[1,], cmb[2,], ~igraph::isomorphic(igraphlist[[.x]], igraphlist[[.y]]))

  sets = rep(NA, numgraph)
  sets[1] = 1
  for(i in seq_len(length(iso))) {
    if(iso[i]) {
      sets[cmb[2, i]] = sets[cmb[1, i]]
    } else if(is.na(sets[cmb[2, i]])) {
      sets[cmb[2, i]] = max(sets, na.rm=T) + 1
    }
  }
  sets
}

#' @export
isomorphism_classes2 = function(igraphlist) {

  # considers topology and leaf labels
  # in order to make sure a combination of leaves and topolgy is recognized as unique, this function uses 'unify_vertex_names'
  # 'unify_vertex_names' often leads to stack overflow; fix this
  # 'unify_vertex_names' should make consistent and unique node names for a given combination of topology and leaves, as long as no node has both more than one incoming edge and more than one outgoing edge

  numgraph = length(igraphlist)
  if(numgraph == 0) return(numeric())
  if(numgraph == 0) return(1)

  igraphlist %<>% map(unify_vertex_names) %>% map(as_edgelist) %>% map(as_tibble) %>% map(~arrange(., V1, V2))

  cmb = combn(numgraph, 2)
  iso = map2_lgl(cmb[1,], cmb[2,], ~identical(igraphlist[[.x]], igraphlist[[.y]]))

  sets = rep(NA, numgraph)
  sets[1] = 1
  for(i in seq_len(length(iso))) {
    if(iso[i]) {
      sets[cmb[2, i]] = sets[cmb[1, i]]
    } else if(is.na(sets[cmb[2, i]])) {
      sets[cmb[2, i]] = max(sets, na.rm=T) + 1
    }
  }
  sets
}

#' Return all valid qpAdm models for an admixturegraph
#'
#' For large admixturegraph, there may be a large number of valid qpAdm models!
#'
#' @param grph an admixture graph as igraph object
#' @param add_outgroup should the graph outgroup be added to the qpAdm right populations? Technically this shouldn't be an informative outgroup for qpAdm.
#' @param nested should a nested data frame be returned (\code{TRUE}), or should populations be concatenated into strings (\code{FALSE})?
#' @param abbr maximum number of characters to print for each population. The default (-1) doesn't abbreviate the names.
#' @export
#' @examples
#' \dontrun{
#' qpadm_models(igraph2, add_outgroup = TRUE)
#' qpadm_models(igraph2, add_outgroup = TRUE) %>% slice(1) %$% list(target, left, right)
#' }
qpadm_models = function(grph, add_outgroup=FALSE, nested = TRUE, abbr = -1) {
  # don't do this for large models

  outgroup = names(V(grph)[2])
  nleaves = grph %>% get_leaves %>% length

  # all triples which could form target - left pop - right pop
  tlr = summarize_graph(grph) %>% filter(tlr)
  if(nrow(tlr) == 0) return(tibble())
  # table of all possible right pops for each combination of target and right pop
  rightpops = tlr %>% group_by(name1, name2) %>% summarize(right = list(name3)) %>% ungroup %>% unnest_longer(right)
  # table of all possible combinations of target and power set of left pops (restricted to nleaves/2 because we need more rightpops than leftpops)
  leftpops = tlr %>% group_by(name1) %>% summarize(left = list(power_set(unique(name2), nmax=min(length(unique(name2)), floor(nleaves/2))))) %>% ungroup %>% unnest_longer(left) %>% group_by(name1) %>% mutate(powerset = seq_len(n())) %>% ungroup %>% unnest_longer(left)
  # table of all possible combinations of target, power set of left pops, and right, where right is right for every left pop of power set
  out = leftpops %>% group_by(name1, powerset) %>% mutate(numleft = n()) %>% filter(numleft > 1) %>% left_join(rightpops, by=c('name1', 'left'='name2')) %>% group_by(name1, powerset, right) %>% mutate(numoginset = n()) %>% group_by(name1, powerset) %>% filter(numoginset == numleft) %>% select(-numoginset) %>% ungroup %>% rename(target = name1)

  # group so that each qpadm model is one line
  out %<>% group_by(target, powerset) %>% summarize(left = list(sort(unique(left))), right = list(sort(unique(right)))) %>% ungroup %>% select(-powerset)

  # add outpop to right pops (even though it shouldn't be a right-group-outgroup in the qpAdm sense)
  if(add_outgroup) out %<>% mutate(right = modify(right, ~sort(c(., outgroup))))

  # only keep models where nright > nleft
  out %<>% filter(right %>% map_int(length) > left %>% map_int(length))

  # add hash to be able to join on qpadm models
  out %>%
    mutate(hash_tl = map2_chr(target, left, ~digest::digest(list(.x, .y)))) %>%
    mutate(hash_tlr = pmap_chr(list(target, left, right), ~digest::digest(list(..1, ..2, ..3))))

  if(abbr != -1) out %<>% mutate(target = str_sub(target, 1, abbr),
                                 left = map(left, ~str_sub(., 1, abbr)),
                                 right = map(right, ~str_sub(., 1, abbr)))
  if(!nested) out %<>% select(target, left, right) %>% mutate(left = map_chr(left, ~paste(., collapse=',')), right = map_chr(right, ~paste(., collapse=',')))

  out

}




split_graph = function(graph) {
  # splits an admixture graph into trees



}



