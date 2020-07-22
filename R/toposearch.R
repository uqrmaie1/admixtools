
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

#' @export
numadmix = function(graph) {
  sum(degree(graph, mode='in') == 2)
}

#' @export
get_leafnames = function(graph) {
  graph %>% V %>% {names(which(degree(graph, ., mode='out') == 0))}
}

get_leaves = function(graph) {
  graph %>% V %>% {.[degree(graph, ., mode='out') == 0]}
}

get_leaves2 = function(graph) {
  # uses 'subcomponent' to return leaves in a consistent order
  graph %>% subcomponent(V(.)[1], mode='out') %>% igraph::intersection(get_leaves(graph))
}

get_root = function(graph) {
  root = names(which(igraph::degree(graph, mode = 'in') == 0))
  if(length(root) != 1) stop(paste0('Root problem ', root))
  root
}

#' @export
get_outpop = function(graph) {
  # returns outpop, if there is an edge from the root to a leave, NULL otherwise
  #graph %>% V %>% names %>% pluck(2)
  stopifnot(class(graph) == 'igraph')
  root = get_root(graph)
  outpop = intersect(names(igraph::neighbors(graph, root)), get_leafnames(graph))
  if(length(outpop) == 1) return(outpop)
}

unify_vertex_names_rec = function(graph, node, vnamemap, sep1 = '.', sep2 = '_') {
  # recursive function which returns a vector of new, unique node names for all nodes which are reachable from 'node'
  children = neighbors(graph, node, mode='out')
  oldnam = names(node)
  if(length(children) == 0) {
    vnamemap[oldnam] = oldnam
    return(vnamemap)
  }

  newnam = rep(NA, length(children))
  for(i in seq_along(children)) {
    vnamemap = unify_vertex_names_rec(graph, children[i], vnamemap, sep1 = sep1, sep2 = sep2)
    newnam[i] = vnamemap[names(children[i])]
  }
  stopifnot(all(!is.na(newnam)))
  vnamemap[oldnam] = paste0(sort(newnam), sep1, collapse = sep2)
  vnamemap
}


unify_vertex_names = function(graph, keep_unique = TRUE, sep1 = '.', sep2 = '_') {
  # this is an aweful function I wrote when I was very tired which changes the vertex names of inner vertices,
  # so that all isomorphic graphs with the same leaf nodes have equally labelled inner nodes
  leaves = get_leaves(graph)
  lv = shortest_unique_prefixes(names(leaves))
  g = set_vertex_attr(graph, 'name', leaves, lv)
  ormap = names(V(g))
  names(ormap) = ormap

  nammap = unify_vertex_names_rec(g, V(g)[1], ormap, sep1 = sep1, sep2 = sep2)
  nammap[lv] = names(leaves)
  names(nammap)[match(lv, names(nammap))] = names(leaves)
  nammap['R'] = 'R'
  if(!keep_unique) {
    changed = setdiff(names(nammap), c('R', names(leaves)))
    nammap[changed] = paste0('n', as.numeric(as.factor(nammap[changed])))
  }
  nammap %<>% map_chr(digest::digest) %>% paste0('x', .)
  set_vertex_attr(graph, 'name', names(nammap), nammap)
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
random_admixturegraph = function(leaves, numadmix = 0, simple = FALSE, outpop = NULL) {
  # makes a random admixture graph
  # returns an 'igraph' graph object
  # 'leaves' can be a number of leaf nodes, or a character vector of leaf names
  stopifnot(class(leaves)[1] %in% c('numeric', 'character'))
  if(length(leaves) == 1) leaves = paste0('l', 1:leaves)
  if(is.null(outpop)) outpop = sample(leaves, 1)
  newick = random_newick(sample(setdiff(leaves, outpop)), outpop = outpop)
  graph = graph_from_edgelist(newick_to_edges(newick))
  graph = insert_admix_igraph_random(graph, numadmix)
  if(simple) return(graph)
  desimplify_graph(graph)
}


#' Remove redundant edges
#'
#' Nodes with in and out degree of one will be removed.
#' @param graph an igraph object
#' @export
#' @examples
#' simple = simplify_graph(example_igraph)
#' plot_graph(example_igraph)
#' plot_graph(simple)
simplify_graph = function(graph) {
  # removes redundant nodes
  convmat = FALSE
  if(class(graph) == 'matrix') {
    convmat = TRUE
    graph = graph_from_edgelist(graph)
  }
  repeat({
    redundant = which(degree(graph, mode='in') == 1 & degree(graph, mode='out') == 1)
    if(length(redundant) == 0) break
    newfrom = names((neighbors(graph, redundant[1], mode='in')))
    newto = names((neighbors(graph, redundant[1], mode='out')))
    graph %<>% igraph::delete_vertices(redundant[1])
    graph %<>% igraph::add_edges(c(newfrom, newto))
    graph %<>% igraph::simplify()
  })
  if(convmat) graph = igraph::as_edgelist(graph)
  graph
}

#' Add two nodes before each admixture node
#'
#' This is used to revert simplify_graph.
#' @export
#' @param graph an admixture graph
#' @examples
#' simple = simplify_graph(igraph1)
#' desimple = desimplify_graph(simple)
#' plot_graph(simple)
#' plot_graph(desimple)
desimplify_graph = function(graph) {
  # adds two nodes before each admixture node
  convmat = FALSE
  if(class(graph)[1] == 'matrix') {
    convmat = TRUE
    graph = graph_from_edgelist(graph)
  }
  admix = which(degree(graph, mode='in') == 2)
  if(length(admix) == 0) return(graph)
  parents = names(unlist(unname(lapply(admix, function(v) neighbors(graph, v, mode='in')))))
  admixnamesrep = rep(names(admix), each=2)
  newnam = paste0(admixnamesrep, c('a', 'b'))
  graph = igraph::add_vertices(graph, length(admix)*2, name=newnam)
  newedges = c(rbind(parents, newnam, newnam, admixnamesrep))
  deledges = paste(parents, admixnamesrep, sep='|')
  graph = igraph::add_edges(graph, newedges)
  graph = igraph::simplify(graph)
  graph = igraph::delete_edges(graph, deledges)
  if(max(degree(graph, mode='in')) > 2 || max(degree(graph, mode='out')) > 2) stop('Desimplify failed')
  if(convmat) graph = igraph::as_edgelist(graph)
  graph
}

split_graph = function(graph) {
  # removes admixture nodes from an admixturegraph
  # returns a tree and a list of edges that can be used to reconstruct the original admixturegraph
  fromnodes = tonodes = c()
  repeat({
    admix = names(which(degree(graph, mode='in') == 2))
    if(length(admix) == 0) break
    parents = names(neighbors(graph, admix[1], mode='in'))
    parents = sample(parents)
    # stopifnot(!any(c(parents[1], admix[1]) %in% unlist(admixedges)))
    # too strict; probably won't be able to always reintroduce admixedges at appropriate locations after SPR
    fromnode0 = names(neighbors(graph, parents[1], mode = 'in'))
    fromnode = setdiff(names(neighbors(graph, parents[1], mode = 'out')), admix[1])
    tonode0 = setdiff(names(neighbors(graph, admix[1], mode = 'in')), parents[1])
    tonode = names(neighbors(graph, admix[1], mode = 'out'))
    fromnodes = c(fromnodes, fromnode)
    tonodes = c(tonodes, tonode)
    graph = igraph::delete_vertices(graph, c(admix[1], parents[1]))
    graph = igraph::add_edges(graph, c(fromnode0, fromnode, tonode0, tonode))
  })
  tree = graph
  namedList(tree, fromnodes, tonodes)
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


insert_edge = function(graph, from, to) {
  # inserts edge from edge 'from' and to edge 'to'
  # more atomic than 'insert_admix_igraph', for better reusability
  n1 = as_ids(from)
  n2 = as_ids(to)
  v1 = str_split(n1, '\\|')[[1]]
  v2 = str_split(n2, '\\|')[[1]]
  o1 = str_replace(n1, '\\|', '_')
  o2 = str_replace(n2, '\\|', '_')
  graph %>%
    add_vertices(2, name = c(o1, o2)) %>%
    add_edges(c(v1[1], o1, o1, v1[2], v2[1], o2, o2, v2[2], o1, o2)) %>%
    delete_edges(c(n1, n2))
}

insert_edges = function(graph, from, to) {
  # from and to are vectors of node names
  # edges will be inserted above each node pair
  stopifnot(length(from) == length(to))
  reduce2(from, to, .init = graph, .f = function(g, x, y) {
    e = incident_edges(g, c(x, y), mode = 'in')
    insert_edge(g, e[[1]], e[[2]])
  })
}

#' @export
insert_admix_igraph = function(graph, fromnodes, tonodes, substitute_missing = FALSE,
                               allow_below_admix = FALSE, desimplify = FALSE) {
  # inserts edges fromnodes -> tonodes into graph
  # inverse of 'split_graph'
  # stopifnot(all(c(sapply(admixedges, `[`, c(2, 4))) %in% names(V(graph))))
  # some nodes will get lost when splitting; insert random amix edges instead

  stopifnot(length(fromnodes) == length(tonodes))
  miss = 0
  maxprev = max(degree(graph))
  for(i in rev(seq_len(length(fromnodes)))) {
    fromnode = fromnodes[i]
    tonode = tonodes[i]
    if(!fromnode %in% names(V(graph)) ||
       !tonode %in% names(V(graph))) {
      miss = miss+1
      next
    }
    toedge2parent = neighbors(graph, tonode, mode='in')
    toedge2sib = setdiff(names(neighbors(graph, toedge2parent, mode='out')), tonode)
    if(!allow_below_admix) {
      if(degree(graph, fromnode, mode='in') > 1 ||
         degree(graph, toedge2sib, mode='in') > 1) {
        miss = miss+1
        next
      }
    }
    fromnode_parent = names(neighbors(graph, fromnode, mode='in'))
    tonode_parent = names(neighbors(graph, tonode, mode='in'))
    graph = igraph::delete_edges(graph, c(paste(fromnode_parent, fromnode, sep='|'),
                                          paste(tonode_parent, tonode, sep='|')))
    new1 = paste0(fromnode_parent, 'x')
    while(new1 %in% names(V(graph))) new1 = paste0(new1, sample(letters, 1))
    j = i
    new2 = paste0('admix', j)
    while(new2 %in% names(V(graph))) {j=j+1; new2 = paste0('admix', j)}
    graph = igraph::add_vertices(graph, 2, name=c(new1, new2))
    newedges = c(fromnode_parent, new1, new1, fromnode, tonode_parent, new2, new2, tonode, new1, new2)
    graph = igraph::add_edges(graph, newedges)
    stopifnot(max(degree(graph)) <= maxprev)
  }
  if (miss > 0) {
    if(substitute_missing) graph = insert_admix_igraph_random(graph, miss)
    else warning('did not insert or substitute admixedge')
  }
  if(desimplify) graph %<>% simplify_graph %>% desimplify_graph
  graph
}


insert_admix_igraph_random = function(graph, nadmix) {
  # graph should be simplified
  for(i in seq_len(nadmix)) {
    firstgen = neighbors(graph, V(graph)[1], mode='out')
    admix = V(graph)[(degree(graph, mode='in') > 1)]
    admixchildren = do.call(c, igraph::ego(graph, 1, admix, mode='out'))
    admixparents = do.call(c, igraph::ego(graph, 1, admix, mode='in'))
    admixsiblings = do.call(c, igraph::ego(graph, 1, admixparents, mode='out'))
    exclude = c(admix, firstgen)
    if(!is.null(admixchildren)) exclude = c(exclude, admixchildren)
    cand_from = sample(igraph::difference(V(graph)[-1], exclude))
    for(fromnode in names(cand_from)) {
      upstream = subcomponent(graph, fromnode, mode='in')
      downstream1 = neighbors(graph, fromnode, mode='out')
      siblings = neighbors(graph, upstream[2], mode='out') # upstream[2] should be parent of fromnode
      uncles = neighbors(graph, neighbors(graph, upstream[2], mode='in'), mode='out')
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
    stopifnot(all(c(fromnode, tonode) %in% names(V(graph))))
    graphnew = insert_admix_igraph(graph, fromnode, tonode, substitute_missing = FALSE)
    stopifnot(numadmix(graphnew) > numadmix(graph))
    stopifnot(igraph::is_simple(graphnew))
    graph = graphnew
  }
  graph
}


insert_admix_igraph_random2 = function(graph, nadmix) {
  leaves = get_leafnames(graph)[-1]
  for(i in 1:nadmix) {
    sel = sample(leaves)[1:2]
    graph %<>% insert_edges(sel[1], sel[2])
  }
  graph
}


subtree_prune_and_regraft = function(graph, only_leaves = FALSE) {
  # cuts of a parts of the tree and attaches it at a random location
  # root -> outgroup stays fixed
  stopifnot(degree(graph, V(graph)[1], mode='in') == 0)
  stopifnot(max(degree(graph, V(graph), mode='in')) == 1)

  firstgen = neighbors(graph, V(graph)[1], mode='out')
  repeat({
    excl = firstgen
    if(only_leaves) excl = c(excl, V(graph)[degree(graph, mode='out') > 0])
    cutnode = sample(igraph::difference(V(graph)[-1], excl), 1)
    cutnodes = subcomponent(graph, cutnode, mode='out')
    cutparent = neighbors(graph, cutnode, mode='in')
    cutgrandparent = neighbors(graph, cutparent, mode='in')
    cutsibling = igraph::difference(neighbors(graph, cutparent, mode='out'), cutnode)
    hostnodes = igraph::difference(V(graph)[-1], c(cutnodes, cutparent, cutsibling, firstgen))
    if(length(hostnodes) > 0) break
  })

  hostnode = sample(names(hostnodes), 1)
  hostparent = names(neighbors(graph, hostnode, mode='in'))

  newnam = paste(hostnode, sample(letters, 1), sep='_')
  while(newnam %in% names(V(graph))) newnam = paste0(newnam, sample(letters, 1))

  graph %<>% igraph::add_vertices(1, name=newnam) %>%
    igraph::add_edges(c(names(cutgrandparent), names(cutsibling),
                hostparent, newnam,
                newnam, hostnode,
                newnam, names(cutnode))) %>%
    igraph::delete_vertices(names(cutparent)) %>%
    igraph::delete_edges(paste(hostparent, hostnode, sep='|'))
  stopifnot(igraph::is_simple(graph))
  graph
}


admixturegraph_prune_and_regraft = function(graph, desimplify = TRUE, only_leaves = FALSE) {
  # 1. remove admixture edges (one randomly selected from each admixture node)
  # 2. runs subtree_prune_and_regraft on resulting tree
  # 3. add admixture edges back on
  if(numadmix(graph) == 0) return(subtree_prune_and_regraft(graph, only_leaves = only_leaves))
  graph = simplify_graph(graph)
  spl = split_graph(graph)
  graph = subtree_prune_and_regraft(spl$tree, only_leaves = only_leaves)
  graph = insert_admix_igraph(graph, spl$fromnodes, spl$tonodes, substitute_missing = TRUE)
  stopifnot(igraph::is_simple(graph))
  if(desimplify) graph = desimplify_graph(graph)
  graph
}
#' @export
spr_leaves = function(graph) admixturegraph_prune_and_regraft(graph, only_leaves = TRUE)
#' @export
spr_all = function(graph) admixturegraph_prune_and_regraft(graph, only_leaves = FALSE)

#' @export
move_admixedge_once = function(graph, desimplify = TRUE) {
  # selects random admixture edge, and moves it to next closest possible spot
  # if not possible, select other admix node
  # if none possible, print warning
  graph = simplify_graph(graph)
  admix = sample(names(V(graph)[(degree(graph, mode='in') > 1)]))
  nadmix = length(admix)
  firstgen = names(neighbors(graph, V(graph)[1], mode='out'))
  for(i in seq_len(nadmix)) {
    parents = sample(names(neighbors(graph, admix[i], mode='in')))
    for(j in 1:2) {
      sibling = setdiff(names(neighbors(graph, parents[j], mode='out')), admix[i])
      grandparent = names(neighbors(graph, parents[j], mode='in'))
      for(n in names(subcomponent(graph, parents[j]))) {
        # subcomponent will return nodes ordered by distance to n
        if(!n %in% names(V(graph)[1]) &
           !n %in% firstgen &
           !n %in% parents &
           !n %in% admix &
           !n %in% c(names(neighbors(graph, parents[1], mode='out')),
                     names(neighbors(graph, parents[2], mode='out'))) &
           !n %in% names(subcomponent(graph, admix[i], mode='out'))) {
          newgrandparent = names(neighbors(graph, n, mode='in'))
          graph = igraph::add_edges(graph, c(grandparent, sibling, newgrandparent, parents[j], parents[j], n))
          graph = igraph::delete_edges(graph, c(paste(newgrandparent, n, sep='|'),
                                              paste(grandparent, parents[j], sep='|'),
                                              paste(parents[j], sibling, sep='|')))
          stopifnot(igraph::is_simple(graph))
          if(desimplify) graph = desimplify_graph(graph)
          return(graph)
        }
      }
    }
  }
  if(nadmix > 0) warning('no suitable attachment points found for admixture edge')
  if(desimplify) graph = desimplify_graph(graph)
  return(graph)
}

#' @export
permute_leaves = function(graph, fix_outgroup=TRUE) {

  leaves = V(graph)[degree(graph, v = V(graph), mode = c('out')) == 0]
  if(fix_outgroup) leaves = leaves[-1]
  nam = names(leaves)
  igraph::set_vertex_attr(graph, 'name', leaves, sample(nam))
}

#' @export
swap_leaves = function(graph, fix_outgroup=TRUE) {

  leaves = V(graph)[degree(graph, v = V(graph), mode = c('out')) == 0]
  if(fix_outgroup) leaves = leaves[-1]
  leaves = sample(leaves)[1:2]
  nam = names(leaves)
  igraph::set_vertex_attr(graph, 'name', leaves, sample(nam))
}

#' @export
flipadmix_random = function(graph, desimplify = TRUE) {
  graph %<>% simplify_graph
  admixedges = graph %>% find_admixedges %>% sample_frac(1)
  el = igraph::as_edgelist(graph)
  out = graph
  for(i in seq_len(nrow(admixedges))) {
    newg = admixedges %>% slice(i) %>% mutate(graph = map2(from, to, ~flipadmix(el, .x, .y))) %$% graph[[1]]
    if(!is.null(newg)) {
      out = newg
      break
    }
  }
  if(desimplify) out %<>% desimplify_graph
  out
}


add_generation = function(models, numgraphs, numsel, qpgfun, mutfuns, parallel = TRUE, verbose = TRUE) {

  space = paste0(paste(rep(' ', 50), collapse=''), '\r')
  if(verbose) alert_info(paste0('Selecting winners...', space))
  lastgen = max(models$generation)
  oldmodels = models %>% filter(generation == lastgen, !is.na(score))
  numsel = min(nrow(oldmodels), numsel)
  numeach = ceiling(numgraphs/numsel)
  winners = oldmodels %>% mutate(prob = score_to_prob(score)) %>% sample_n(numsel-1, weight=prob) %>% select(-prob)
  winners %<>% bind_rows(oldmodels %>% top_n(1, -jitter(score, amount = 1e-9)) %>% slice(1))
  sq = (numsel+1):numgraphs
  newmodels = tibble(generation = lastgen+1, index = sq,
                     igraph = rep(winners$igraph, numeach)[sq],
                     oldscore = rep(winners$score, numeach)[sq])
  mutations = sample(mutfuns, numgraphs-numsel, replace = TRUE)
  if(parallel) {
    map = furrr::future_map
    imap = furrr::future_imap
  }
  if(verbose) alert_info(paste0('Generating new graphs...', space))
  newmodels %<>% mutate(mutation = names(mutations),
                        igraph = imap(igraph, ~tryCatch(exec(mutations[[.y]], .x), error = function(e) .x)))
  if(verbose) alert_info(paste0('Evaluating graphs...', space))
  newmodels %<>% mutate(out = map(igraph, qpgfun))
  if(verbose) alert_info(paste0('Attaching to previous generations...', space))
  winners %>%
    mutate(generation = lastgen+1, index = 1:n()) %>%
    bind_rows(newmodels %>% unnest_wider(out))
}



score_to_prob = function(score) {
  logscore = log10(score)
  sc = 1-logscore/max(logscore)+1e-10
  sc = sc^20
  sc = replace_na(sc/sum(sc), 0)
  sc
}


evolve_topology = function(init, numgraphs, numgen, numsel, qpgfun, mutfuns, repnum, parallel = TRUE,
                           keep = 'all', store_intermediate = NULL, stop_at = NULL, verbose = TRUE) {
  out = init
  for(i in seq_len(numgen)) {
    if(!is.null(stop_at) && stop_at < Sys.time()) break
    newmodels = add_generation(out, numgraphs, numsel, qpgfun, mutfuns, parallel = parallel, verbose = verbose)
    #if(min(newmodels$score) > min(out$score)) browser()
    if(!is.null(store_intermediate)) {
      fl = paste0(store_intermediate, '_rep', repnum, '_gen', i, '.rds')
      saveRDS(newmodels, fl)
    }
    if(keep == 'all') out %<>% bind_rows(newmodels)
    else if(keep == 'best') out %<>% group_by(generation) %>% top_n(1, -jitter(score, amount = 1e-9)) %>%
      ungroup %>% bind_rows(newmodels)
    else out = newmodels
    if(verbose) {
      #best = newmodels %>% filter(index > numsel) %>% top_n(numsel, -jitter(score, amount = 1e-9)) %$% score %>% sort
      best = newmodels %>% top_n(numsel, -jitter(score, amount = 1e-9)) %$% score %>% sort
      alert_success(paste0('Generation ', i, '  Best scores: ', paste(round(best, 3), collapse=', '),
                    paste(rep(' ', 30), collapse=''), '\n'))
    }
  }
  if(verbose) cat('\n')
  if(keep == 'best') out %<>% group_by(generation) %>% top_n(1, -jitter(score, amount = 1e-9)) %>% ungroup
  out
}


optimize_admixturegraph_single = function(pops, precomp, repnum, numgraphs = 50, numgen = 5,
                                          numsel = 5, numadmix = 0, outpop = NA, initgraph = NULL,
                                          mutfuns = NULL, parallel = TRUE, stop_at = NULL,
                                          store_intermediate = NULL,
                                          cpp = TRUE, debug = FALSE, keep = 'all', verbose = TRUE) {
  if(numadmix == 0 && is.null(initgraph)) {
    qpgfun0 = qpgraph_anorexic
    mutfuns = mutfuns[setdiff(names(mutfuns), c('move_admixedge_once'))]
  } else {
    qpgfun0 = function(...) qpgraph(f2_blocks = NULL, ...)
  }

  qpgfun = function(graph) qpgfun0(graph = graph, f3precomp = precomp, cpp = cpp, numstart = 1, verbose = FALSE)
  if(!debug) qpgfun = possibly(qpgfun, otherwise = NULL)
  space = paste0(paste(rep(' ', 50), collapse=''), '\r')
  if(verbose) alert_info(paste0('Generate new graphs...', space))
  if(is.null(initgraph)) initgraphs = replicate(numgraphs, random_admixturegraph(pops, numadmix, outpop=outpop), simplify = FALSE)
  else {
    missing = numadmix - numadmix(initgraph)
    initgraphs = replicate(numgraphs, initgraph, simplify = FALSE)
    if(missing > 0) initgraphs %<>% map(~insert_admix_igraph_random(., missing))
  }
  if(verbose) alert_info(paste0('Evaluate graphs...', space))
  if(parallel) map = furrr::future_map
  init = tibble(generation=0, index = seq_len(numgraphs),
                igraph = initgraphs, mutation = 'random_admixturegraph') %>%
    mutate(out = map(igraph, qpgfun), isn = map_lgl(out, is.null))
  if(all(init$isn)) stop('All NULL!')

  init %<>% select(-isn) %>% unnest_wider(out) %>% mutate(oldscore = score)
  if(verbose) {
    best = init %>% filter(!is.na(score)) %>% top_n(min(numsel, n()), -jitter(score, amount = 1e-9)) %$% score
    alert_success(paste0('Generation 0  Best new scores: ', paste(round(best), collapse=', '), space))
  }

  evolve_topology(init, numgraphs, numgen, numsel, qpgfun, mutfuns, repnum, parallel = parallel,
                  keep = keep, store_intermediate = store_intermediate, stop_at = stop_at, verbose = verbose)
}



#' Find well fitting admixture graphs
#'
#' This function generates and evaluates admixture graphs in `numgen` iterations across `numrep` independent repeats
#' to find well fitting admixturegraphs. It uses the function \code{\link{future_map}} from the \code{\link{furrr}}
#' package to parallelize across the independent repeats. The function \code{\link{future::plan}} can be called
#' to specify the details of the parallelization. This can be used to parallelize across cores or across nodes on
#' a compute cluster. Setting `numadmix` to 0 will search for well fitting trees, which is much faster than searching
#' for admixture graphs with many admixture nodes.
#' @export
#' @param f2_data A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' alternatively, a directory with precomputed data. see \code{\link{extract_f2}} and \code{\link{extract_indpairs}}.
#' @param pops Populations for which to fit admixture graphs
#' @param outpop Outgroup population
#' @param numrep Number of independent repetitions (each repetition can be run in parallel)
#' @param numgraphs Number of graphs in each generation
#' @param numgen Number of generations
#' @param numsel Number of graphs which are selected in each generation. Should be less than `numgraphs`.
#' @param numadmix Number of admixture events within each graph
#' @param keep By default (`all`), all evaluated graphs will be returned. `best` will only return the best fitting
#' graph from each repeat and each generation. `last` will return all graphs from the last generation.
#' @param f2_denom Scales f2-statistics. A value of around 0.278 converts f2 to Fst.
#' @param cpp Use C++ functions. Setting this to `FALSE` will be slower but can help with debugging.
#' @param initgraph Optional graph to start with. If `NULL`, optimization will start with random graphs.
#' @param mutfuns The names of funcations used to modify graphs.
#' \itemize{
#' \item `spr_leaves` Subtree prune and regraft leaves. Cuts a leaf node and attaches it to a random other edge in the graph.
#' \item `spr_all` Subtree prune and regraft. Cuts any edge and attaches the new orphan node to a random other edge in the graph, keeping the number of admixture nodes constant.
#' \item `swap_leaves` Swaps two leaf nodes.
#' \item `move_admixedge_once` Moves an admixture edge to a nearby location.
#' \item `flipadmix_random` Flips the direction of an admixture edge (if possible).
#' }
#' @param store_intermediate Path and prefix of files for intermediate results to `.rds`. Can be useful if `find_graphs` doesn't finish sucessfully.
#' @param parallel Parallelize over repeats (if `numrep > 1`) or graphs (if `numrep == 1`) by replacing `purrr::map` with `furrr::future_map`. Will only be effective if `future::plan()` has been set.
#' @param stop_at Stop execution after finishing the generation running at `stop_at` seconds. Currently not working.
#' @param debug If `TRUE` each repeat is run sequentially in a loop and not via \code{\link{furrr::map}}).
#' Errors will interrupt execution. This is the default if `numrep = 1`
#' @param fudge_cov Regularization term added to the covariance matrix of estimated f3 statistics (after scaling by the matrix trace).
#' @param verbose Print progress updates
#' @param ... Additional arguments passed to `qpgraph`
#' @return a nested data frame with one model per line
#' @examples
#' \dontrun{
#' find_graphs(example_f2_blocks, numrep = 200, numgraphs = 100,
#'             numgen = 20, numsel = 5, numadmix = 3)
#' }
find_graphs = function(f2_data, pops = NULL, outpop = NULL, numrep = 10, numgraphs = 50,
                       numgen = 5, numsel = 5, numadmix = 0, keep = c('all', 'best', 'last'),
                       f2_denom = 1, cpp = TRUE, initgraph = NULL,
                       mutfuns = c('spr_leaves', 'spr_all', 'swap_leaves', 'move_admixedge_once', 'flipadmix_random'),
                       store_intermediate = NULL,
                       parallel = TRUE, stop_at = NULL,
                       debug = FALSE, fudge_cov = 1e-5, verbose = TRUE, ...) {

  keep = rlang::arg_match(keep)
  if(numsel >= numgraphs || numsel < 1) stop("'numsel' has to be smaller than 'numgraphs' and greater than 0!")
  if(!is.null(pops) && !is.null(initgraph)) stop("You can't provide 'pops' and 'initgraph' at the same time!")

  if(!is.null(initgraph)) {
    if('data.frame' %in% class(initgraph) || 'matrix' %in% class(initgraph)) {
      initgraph = graph_from_edgelist(as.matrix(initgraph[,1:2]))
    } else if('character' %in% class(initgraph)) {
      initgraph = graph_from_edgelist(as.matrix(parse_qpgraph_graphfile(initgraph)[,1:2]))
    } else if(!'igraph' %in% class(initgraph)) stop("'initgraph' format not recognized!")
  }

  if(is.character(f2_data) && file.exists(f2_data) && !isTRUE(file.info(f2_data)$isdir)) {
    # parse Nick's fstats
    precomp = parse_fstats(f2_data)
    precomp$f3_est = precomp$f3
    precomp$ppinv = solve(precomp$f3var)
    if(!is.null(initgraph)) {
      pops = get_leafnames(initgraph)
    } else pops = precomp$pops
  } else {
    if(is.null(pops)) {
      if(!is.null(initgraph)) {
        pops = get_leafnames(initgraph)
      } else if(is.character(f2_data)) {
        pops = list.dirs(f2_data, full.names=FALSE, recursive=FALSE)
      } else pops = dimnames(f2_data)[[1]]
    }
    #pops = get_leafnames(initgraph)
    precomp = qpgraph_precompute_f3(f2_data, pops, outpop = outpop, f2_denom = f2_denom, fudge_cov = fudge_cov)
  }

  if(is.null(outpop)) outpop = pops[1]
  else if(!outpop %in% pops) pops = c(outpop, pops)
  mutfuns = sapply(mutfuns, get)

  oa = function(i) optimize_admixturegraph_single(pops, precomp, repnum = i,
                                                 numgen = numgen, numsel = numsel,
                                                 numgraphs = numgraphs, numadmix = numadmix,
                                                 outpop = outpop, cpp = cpp, initgraph = initgraph,
                                                 mutfuns = mutfuns, parallel = parallel && numrep == 1,
                                                 stop_at = stop_at, store_intermediate = store_intermediate,
                                                 debug = debug, keep = keep, verbose = verbose && numrep == 1)
  if(!debug && numrep > 1) oa = possibly(oa, otherwise = NULL)
  if(parallel && numrep > 1) {
    res = furrr::future_map(as.list(seq_len(numrep)), oa)
  } else {
    res = list()
    for(i in seq_len(numrep)) {
      if(verbose && numrep > 1) alert_info(paste0('Repeat ', i, ' out of ', numrep, '...\n'))
      res[[i]] = oa(i)
      #if(is.null(res[[i]])) browser()
    }
  }
  res = bind_rows(res, .id='run')
  if(nrow(res) > 0) res %<>% mutate(run = as.numeric(run))
  res
}



summarize_graph = function(graph, exclude_outgroup = TRUE) {

  leaves = V(graph)[degree(graph, v = V(graph), mode = c('out')) == 0]
  if(exclude_outgroup) leaves = igraph::difference(leaves, V(graph)[2])

  paths = all_simple_paths(graph, V(graph)[1], leaves, mode = 'out') %>%
    enframe(name = 'pathnum', value='path') %>% mutate(name1 = map_chr(path, ~(tail(names(.), 1))))

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
#' @param maxscore restrict summary to graphs with score not larger than `maxscore`
#' @examples
#' \dontrun{
#' summarize_triples(opt_results)
#' }
summarize_triples = function(results, maxscore = NA) {
  # results is output from 'optimize_admixturegraph'
  # takes at most one graph from each independent run

  sel = results %>% group_by(run) %>% top_n(1, -jitter(score, amount = 1e-9)) %>% ungroup
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
  # in order to make sure a combination of leaves and topolgy is recognized as unique, this function uses
  # 'unify_vertex_names'
  # 'unify_vertex_names' often leads to stack overflow; fix this
  # 'unify_vertex_names' should make consistent and unique node names for a given combination of topology and leaves,
  # as long as no node has both more than one incoming edge and more than one outgoing edge

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
#' @export
#' @param graph an admixture graph as igraph object
#' @param add_outgroup should the graph outgroup be added to the qpAdm right populations?
#' Technically this shouldn't be an informative outgroup for qpAdm.
#' @param nested should a nested data frame be returned (`TRUE`), or should populations be concatenated
#' into strings (`FALSE`)?
#' @param abbr maximum number of characters to print for each population. The default (-1) doesn't abbreviate the names.
#' @examples
#' \dontrun{
#' qpadm_models(igraph2, add_outgroup = TRUE)
#' qpadm_models(igraph2, add_outgroup = TRUE) %>% slice(1) %$% list(target, left, right)
#' }
qpadm_models = function(graph, add_outgroup=FALSE, nested = TRUE, abbr = -1) {
  # don't do this for large models

  outgroup = names(V(graph)[2])
  nleaves = graph %>% get_leaves %>% length

  # all triples which could form target - left pop - right pop
  tlr = summarize_graph(graph) %>% filter(tlr)
  if(nrow(tlr) == 0) return(tibble())
  # table of all possible right pops for each combination of target and right pop
  rightpops = tlr %>% group_by(name1, name2) %>% summarize(right = list(name3)) %>% ungroup %>% unnest_longer(right)
  # table of all possible combinations of target and power set of left pops (restricted to nleaves/2
  # because we need more rightpops than leftpops)
  leftpops = tlr %>% group_by(name1) %>%
    summarize(left = list(power_set(unique(name2), nmax=min(length(unique(name2)), floor(nleaves/2))))) %>%
    ungroup %>% unnest_longer(left) %>% group_by(name1) %>% mutate(powerset = seq_len(n())) %>% ungroup %>%
    unnest_longer(left)
  # table of all possible combinations of target, power set of left pops, and right, where right is right
  # for every left pop of power set
  out = leftpops %>% group_by(name1, powerset) %>% mutate(numleft = n()) %>% filter(numleft > 1) %>%
    left_join(rightpops, by=c('name1', 'left'='name2')) %>% group_by(name1, powerset, right) %>%
    mutate(numoginset = n()) %>% group_by(name1, powerset) %>% filter(numoginset == numleft) %>%
    select(-numoginset) %>% ungroup %>% rename(target = name1)

  # group so that each qpadm model is one line
  out %<>% group_by(target, powerset) %>%
    summarize(left = list(sort(unique(left))), right = list(sort(unique(right)))) %>% ungroup %>% select(-powerset)

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
  if(!nested) out %<>% select(target, left, right) %>%
    mutate(left = map_chr(left, ~paste(., collapse=',')), right = map_chr(right, ~paste(., collapse=',')))

  out

}

#' @export
decompose_graph = function(graph) {
  # splits an admixture graph into trees
  if(!'igraph' %in% class(graph)) {
    graph %<>% graph_from_edgelist
  }
  edges = graph %>% simplify_graph %>% as_edgelist %>%
    as_tibble(.name_repair = ~c('from', 'to'))

  admixmat = edges %>% mutate(i = 1:n()) %>% group_by(to) %>% mutate(cnt = n(), j = 1:n()) %>% ungroup %>%
    filter(cnt == 2) %>% select(to, i, j) %>% spread(to, i) %>% select(-j) %>% as.matrix
  nadmix = ncol(admixmat)
  if(nadmix == 0) return(list(graph))
  #mutate(itree = map(itree, ~delete_vertices(., setdiff(get_leafnames(.), get_leafnames(tree)))))
  indices = expand.grid(replicate(ncol(admixmat), 1:2, simplify = FALSE)) %>% as.matrix
  trees = map(1:nrow(indices), ~slice(edges, -admixmat[cbind(indices[.,], 1:nadmix)]))
  trees %>% map(as.matrix) %>% map(graph_from_edgelist) %>%
    map(~delete_vertices(., setdiff(get_leafnames(.), get_leafnames(graph)))) %>%
    map(simplify_graph) %>% enframe(value = 'graph')
}

#' @export
decomposed_tree_neighbors = function(graph, desimplify = TRUE) {
  graph %>% simplify_graph %>% decompose_graph %$% graph %>% map(tree_neighbors) %>%
    bind_rows %>% rename(graph = itree) %>%
    mutate(isoclass = isomorphism_classes2(graph)) %>%
    filter(!duplicated(isoclass)) %>% select(-isoclass)
}


#' @export
tree_neighbors = function(tree) {
  # returns nested data frame with all trees within edit distance 1
  root = V(tree)[1]
  gen1 = neighbors(tree, root, mode = 'out')
  nodes = difference(V(tree), c(root, gen1))
  map(nodes, ~tree_neighbors_single(tree, .)) %>% bind_rows
}


tree_neighbors_single = function(tree, node) {
  # returns nested data frame with all trees within edit distance 1,
  # that are created by cutting all edges leading to node
  if(class(node) != 'character') node = names(node)
  downstream = subcomponent(tree, node, mode = 'out')
  root = V(tree)[1]
  outgroup = V(tree)[2]
  parent = neighbors(tree, node, mode = 'in')
  sibling = neighbors(tree, parent, mode = 'out')
  from = names(difference(V(tree), c(downstream, root, outgroup, parent, sibling)))
  tibble(from) %>% mutate(to = node, itree = map2(from, to, ~replace_edge(tree, .x, .y)))
}

replace_edge = function(tree, from, to) {
  # removes edge leading to node 'to', and adds edge beginning above node 'from'
  parent = names(neighbors(tree, to, mode = 'in'))
  grandparent = names(neighbors(tree, parent, mode = 'in'))
  sibling = setdiff(names(neighbors(tree, parent, mode = 'out')), to)
  newgrandparent = names(neighbors(tree, from, mode = 'in'))
  tree %>%
    add_edges(c(grandparent, sibling)) %>%
    delete_vertices(parent) %>%
    add_vertices(1, name = parent) %>%
    delete_edges(paste(newgrandparent, from, sep = '|')) %>%
    add_edges(c(parent, to, newgrandparent, parent, parent, from))
}

#' @export
graph_plusone = function(graph, desimplify = TRUE) {
  # returns all graphs with one more edge
  graph %<>% simplify_graph
  newedges = find_newedges(graph)
  fn = ~insert_admix_igraph(graph, .x, .y)
  if(desimplify) fn %<>% compose(desimplify_graph, .dir = 'forward')
  newedges %>% mutate(graph = map2(from, to, fn))
}

#' @export
graph_minusone = function(graph, desimplify = TRUE) {
  # returns all graphs with one admixture edge removed
  graph %<>% simplify_graph
  admixedges = find_admixedges(graph)
  fn = ~delete_edges(graph, paste(.x, .y, sep = '|')) %>%
    delete_vertices(setdiff(get_leafnames(.), get_leafnames(graph))) %>%
    simplify_graph
  if(desimplify) fn %<>% compose(desimplify_graph, .dir = 'forward')
  admixedges %>% mutate(graph = map2(from, to, fn))
}

#' @export
graph_minusplus = function(graph, desimplify = TRUE) {
  graph %>% graph_minusone %>% mutate(graph2 = map(graph, ~graph_plusone(., desimplify=desimplify))) %>%
    rename(from1 = from, to1 = to) %>% select(-graph) %>% unnest(graph2) %>%
    rename(from2 = from, to2 = to) %>% mutate(isoclass = isomorphism_classes2(graph)) %>%
    filter(!duplicated(isoclass)) %>% select(-isoclass)
}

#' @export
graph_flipadmix = function(graph, desimplify = TRUE) {
  graph %<>% simplify_graph
  admixedges = graph %>% find_admixedges
  el = igraph::as_edgelist(graph)
  fn = ~flipadmix(el, .x, .y)
  out = admixedges %>%
    mutate(graph = map2(from, to, fn)) %>%
    filter(!map_lgl(graph, is.null))
  if(nrow(out) == 0) return()
  if(desimplify) out %<>% mutate(graph = map(graph, desimplify_graph))
  out
}

flipadmix = function(edges, from, to) {
  edges %<>% as.matrix
  row = which(edges[,1] == from & edges[,2] == to)
  edges[row,] = c(to, from)
  g = igraph::graph_from_edgelist(edges)
  if(!is.dag(g) || max(degree(g, mode = 'in')) > 2 || max(degree(g, mode = 'out')) > 2) g = NULL
  g
}

#' @export
graph_addleaf = function(graph, pop) {
  graph %>%
    find_normedges(exclude_first = TRUE) %>%
    mutate(graph = map2(from, to, ~insert_leaf(graph, pop, .x, .y)))
}

#' @export
find_newedges = function(graph) {
  # returns a two column data frame with edges that can be inserted into the graph

  edges = matrix(NA, 0, 2)
  firstgen = neighbors(graph, V(graph)[1], mode='out')
  admix = V(graph)[(degree(graph, mode='in') > 1)]
  admixchildren = do.call(c, igraph::ego(graph, 1, admix, mode='out'))
  admixparents = do.call(c, igraph::ego(graph, 1, admix, mode='in'))
  admixsiblings = do.call(c, igraph::ego(graph, 1, admixparents, mode='out'))
  exclude = c(admix, firstgen)
  if(!is.null(admixchildren)) exclude = c(exclude, admixchildren)
  cand_from = igraph::difference(V(graph)[-1], exclude)
  for(fromnode in names(cand_from)) {
    upstream = subcomponent(graph, fromnode, mode='in')
    downstream1 = neighbors(graph, fromnode, mode='out')
    siblings = neighbors(graph, upstream[2], mode='out') # upstream[2] should be parent of fromnode
    uncles = neighbors(graph, neighbors(graph, upstream[2], mode='in'), mode='out')
    exclude = c(upstream, downstream1, siblings, uncles)
    if(!is.null(admixsiblings)) exclude = c(exclude, admixsiblings)
    cand_to = names(igraph::difference(cand_from, exclude))
    if(length(cand_to) > 0) edges %<>% rbind(cbind(fromnode, cand_to))
  }
  edges %>% as_tibble(.name_repair = ~c('from', 'to'))
}

#' @export
find_admixedges = function(graph) {
  # returns a two column data frame with edges that can be removed from the graph

  edges = matrix(NA, 0, 2)
  admixnodes = which(degree(graph, mode = 'in') == 2)
  if(length(admixnodes) == 0) return(tibble(from = character(), to = character()))
  graph %>% adjacent_vertices(admixnodes, mode = 'in') %>%
    map(names) %>% enframe('to', 'from') %>%
    unnest_longer(from) %>% select(2:1)
}

#' @export
find_normedges = function(graph, exclude_first = FALSE) {
  # returns a two column data frame with non-admixture edges
  el = graph %>% as_edgelist
  if(exclude_first) el = el[-1,]
  el %>% as_tibble(.name_repair = ~c('from', 'to')) %>%
    group_by(to) %>% filter(n() == 1) %>% ungroup
}

#' @export
delete_admix = function(graph, from, to, desimplify = TRUE) {
  # returns graph with admixture edge deleted
  # does not conserve internal node names
  # assumes input graph is desimplified
  parents = neighbors(graph, from, mode = 'in')
  if(length(parents) == 1) {
    del = from
    if(length(neighbors(graph, parents, mode = 'in')) == 2) del = c(del, names(parents))
    graph %<>% delete_vertices(del)
  } else {
    graph %<>% delete_edges(paste(from, to, sep = '|'))
  }
  graph %<>% simplify_graph
  if(desimplify) graph %<>% desimplify_graph
  graph
}

#' @export
delete_leaf = function(graph, leaf, desimplify = TRUE) {
  # deletes leaf and all internal nodes leading to no other leaves

  leaves = get_leafnames(graph)
  graph %<>%
    subcomponent(leaf, mode = 'in') %>%
    keep(~length(intersect(leaves, names(subcomponent(graph, .x, mode = 'out')))) == 1) %>%
    delete_vertices(graph, .) %>%
    simplify_graph()
  if(desimplify) graph %<>% desimplify_graph
}

insert_leaf = function(graph, leaf, from, to) {
  # inserts new leaf at edge from -> to
  nn = paste(from, to, sep = '_')
  graph %>%
    add_vertices(2, attr = list('name' = c(nn, leaf))) %>%
    add_edges(c(from, nn, nn, to, nn, leaf)) %>%
    delete_edges(paste(from, to, sep = '|'))
}


generate_all_trees = function(leaves) {

  stopifnot(!'R' %in% leaves)
  if(any(str_detect(leaves, '\\|'))) stop('Leaves cannot have "|" in name!')
  init = graph_from_edgelist(matrix(c('R', leaves[1]), 1))
  add_leaves_rec(init, leaves[-1])
}

generate_all_graphs = function(leaves, nadmix = 0, sure = FALSE, verbose = TRUE) {
  #numtot = admixtools:::numtreesadmix(length(leaves), nadmix)
  #if(numtot > 1000 && !sure) stop(paste0('If you really want to generate ', numtot, ' graphs, set sure to TRUE'))
  nleaves = length(leaves)
  if(nleaves > 5 | nadmix > 2) stop('If you really want to generate that many graphs, set sure to TRUE')
  if(verbose) alert_info(paste0('Generating ', numtrees(nleaves),' trees...\n'))
  trees = generate_all_trees(leaves)
  if(verbose) alert_info(paste0('Adding all possible admixutre edges...\n'))
  graphs = flatten(map(trees, ~add_edges_rec(., nadmix)))
  if(verbose) alert_info(paste0('Identifying isomorpisms in ', length(graphs),' graphs...\n'))
  iso = isomorphism_classes2(graphs)
  graphs[!duplicated(iso)]
}

add_leaves_rec = function(tree, leaves) {
  el = tree %>% E %>% as_ids %>% str_split('\\|')
  newtrees = map(el, ~{
    from = .[1]
    to = .[2]
    nam = paste(from, to, sep='_')
    tree %>%
      add_vertices(2, name = c(nam, leaves[1])) %>%
      add_edges(edges = c(nam, leaves[1])) %>%
      add_edges(edges = c(from, nam)) %>%
      add_edges(edges = c(nam, to)) %>%
      delete_edges(paste(from, to, sep = '|'))
  })
  if(length(leaves) == 1) return(newtrees)
  flatten(map(newtrees, ~add_leaves_rec(., leaves[-1])))
}

add_all_edges = function(graph) {
  eg = E(graph)
  nume = length(eg)
  newe = expand_grid(from = seq_len(nume), to = seq_len(nume)) %>% filter(from != to)
  newe %>% mutate(newg = map2(from, to, ~insert_edge(graph, eg[.x], eg[.y]))) %$%
    newg %>% keep(~is.dag(.))
}

add_edges_rec = function(graph, nadmix) {
  if(nadmix == 0) return(list(graph))
  newgraphs = add_all_edges(graph)
  flatten(map(newgraphs, ~add_edges_rec(., nadmix-1)))
}


#' Return all graphs created from permuting a subclade
#'
#' generates new graphs from basegraph as follows:
#' 1. generates all possible trees using `addpops`` (which are not in basegraph)
#' 2. attaches trees to connection_edge, which is defined by two nodes in basegraph
#' 3. adds edges originating above each edge in `source_node`, to each node above `addpops``
#'
#' @export
#' @param basegraph an admixture graph as igraph object. (convert from edge list using `igraph::graph_from_edgelist`)
#' @param addpops a vector of population labels which are not in `basegraph`. These populations should form a clade. All possible trees will be generated and those trees will be attached to `basegraph`.
#' @param connection_edge edge in `basegraph` where the tree made from `addpops` should be attached
#' @param source_nodes nodes in `basegraph`. edges above these nodes will be added and attached to all terminal edges leading to `addpops`
#' @examples
#' \dontrun{
  #' graphlist = graphmod_pavel(example_igraph, addpops = c('pop1', 'pop2', 'pop3'),
  #'                            connection_edge = c('N2N0', 'N1N'), source_nodes = c('Denisova.DG', 'N2N2'))
  #' results = tibble(graph = graphlist) %>%
  #'   mutate(res = map(graph, ~qpgraph(example_f2_blocks, .))) %>%
  #'   unnest_wider(res) %>%
  #'   mutate(worstz = map_dbl(f3, ~max(abs(.$z))))
#' }
graphmod_pavel = function(basegraph, addpops, connection_edge, source_nodes) {

  c1 = connection_edge[1]
  c2 = connection_edge[2]
  cn = 'connection_node'
  tn = 'tree_R'

  trees = generate_all_trees(addpops) %>% map(~{
    leaves = get_leaves(.)
    internal = difference(V(.), leaves)
    set_vertex_attr(., 'name', index = internal, value = paste0('tree_', names(internal)))
    })
  graphs = trees %<>% map(~{
    igraph::union(basegraph, .) %>%
      add_vertices(1, name = cn) %>%
      add_edges(c(c1, cn, cn, c2, cn, tn)) %>%
      delete_edges(paste(c1, c2, sep='|'))
    })
  e = expand_grid(source_nodes, addpops)
  graphs %>% map(~insert_edges(., e$source_nodes, e$addpops))
}


split_multifurcations = function(graph) {
  # graph is a four column edge list matrix
  # returns graph as tibble with ''

  stopifnot(class(graph)[1] == 'matrix')
  multi = names(which(table(graph[,1]) > 2))
  newedges = matrix(0, 0, 4)

  for(i in seq_len(length(multi))) {
    rows = which(graph[,1] == multi[i])
    torows = which(graph[,2] == multi[i])
    newnam = paste0(graph[rows, 1], '_multi', i, '_', seq_along(rows))
    graph[rows, 1] = newnam
    graph[torows, 2] = newnam[1]
    newedges = rbind(newedges, cbind(newnam[-length(newnam)], newnam[-1], '0', '1e-9'))
  }
  if(ncol(graph) == 2) graph = cbind(graph, NA, NA)
  graph %>% rbind(newedges) %>%
    as_tibble(.name_repair = ~c('from', 'to', 'lower', 'upper')) %>%
    type_convert(col_types = cols())
}



#' Returns a signature of a graph consisting of the left and right descendent leaf nodes of each internal node (sorted and concatenated)
#'
#' Can be used to determine how often internal nodes occur in a list of other well fitting models
#' @export
#' @examples
#' \dontrun{
#' sigs = example_winners %>% mutate(sig = map(igraph, node_signature)) %$% sig %>% unlist %>% table %>% c
#' node_signature(example_winners$igraph[[1]])
#' }
node_signature = function(graph) {
  graph %<>% simplify_graph
  inner = setdiff(names(V(graph)), get_leafnames(graph))
  inner2 = inner %>%
    map(~neighbors(graph, ., mode = 'out')) %>%
    map(names)
  leaves = get_leafnames(graph)
  map(inner2, ~map(., ~names(subcomponent(graph, ., mode = 'out'))) %>%
        map(~intersect(., leaves)) %>%
        map(sort) %>%
        map(~paste(., collapse=':')) %>%
        flatten_chr) %>%
    map(sort) %>%
    map(~paste(., collapse=' ')) %>%
    flatten_chr %>%
    set_names(inner)
}

#' Count how often each node in graph occurs in other graphs
#'
#' @examples
#' \dontrun{
#' sigs = example_winners %>% mutate(sig = map(igraph, node_signature)) %$% sig %>% unlist %>% table %>% c
#' node_signature(example_winners$igraph[[1]])
#' }
node_counts = function(graph, graphlist) {
  counts = map(graphlist, node_signature) %>%
    map(unique) %>%
    unlist %>%
    table %>%
    c %>%
    enframe(name = 'signature', value = 'count')
  graph %>%
    node_signature %>%
    enframe(name = 'node', value = 'signature') %>%
    left_join(counts, by = 'signature')
}


match_graphs = function(graph1, graph2) {

  v1 = names(V(graph1))
  v2 = names(V(graph2))
  l1 = get_leafnames(graph1)
  l2 = get_leafnames(graph2)
  stopifnot(all.equal(sort(l1), sort(l2)))

  vmap = cbind(l1, l1, l1)
  for(i in rep(seq_along(v1), 2)) {
    if(v1[i] %in% vmap[,1]) next
    parents = neighbors(graph1, v1[i], mode = 'in')
    children = neighbors(graph1, v1[i], mode = 'out')
    np = length(parents)
    nc = length(children)
    cl = sort(intersect(l1, names(children)))
    indeg = unname(sort(degree(graph1, children, mode = 'in')))
    if(length(cl) > 0) {
      cands = names(neighbors(graph2, cl[1], mode = 'in'))
      for(cand in cands) {
        print(cand)
        cancchildren = neighbors(graph2, cand, mode = 'out')
        cl2 = sort(intersect(l1, names(cancchildren)))
        if(isTRUE(all.equal(cl, cl2))) {
          if(!v1[i] %in% vmap[,1] && !cand %in% vmap[,2]) {
            vmap = rbind(vmap, c(v1[i], cand, v1[i]))
            break
          }
        }
      }
      if(v1[i] %in% vmap[,1]) next
    }
    matchedparents = intersect(vmap[,1], names(parents))
    if(length(matchedparents) == 0) next
    mp2 = vmap[vmap[,1] == matchedparents[1], 2]
    stopifnot(mp2 %in% names(V(graph2)))
    if(length(mp2) == 1 && !is.na(mp2)) {
      cands = names(neighbors(graph2, mp2, mode = 'out'))
    }
    for(cand in cands) {
      print(cand)
      cancchildren = neighbors(graph2, cand, mode = 'out')
      if(length(cl) > 0) {
        cl2 = sort(intersect(l1, names(cancchildren)))
        if(isTRUE(all.equal(cl, cl2))) {
          if(!v1[i] %in% vmap[,1] && !cand %in% vmap[,2]) {
            vmap = rbind(vmap, c(v1[i], cand, v1[i]))
            break
          }
        }
      } else {
        if(!v1[i] %in% vmap[,1] && !cand %in% vmap[,2]) {
          vmap = rbind(vmap, c(v1[i], cand, v1[i]))
          break
        }
        # indeg2 = unname(sort(degree(graph2, names(cancchildren), mode = 'in')))
        # if(isTRUE(all.equal(indeg, indeg2))) {
        #   vmap = rbind(vmap, c(v1[i], cand, v1[i]))
        #   break
        # }
      }
    }
  }

  # add non-matched nodes
  rem1 = setdiff(v1, vmap[,1])
  rem2 = setdiff(v2, vmap[,2])
  vmap = rbind(vmap, cbind(rem1, NA, paste0('newnode1_', 1:length(rem1))))
  vmap = rbind(vmap, cbind(NA, rem2, paste0('newnode2_', 1:length(rem2))))

  # edge lists with new names
  el1 = as_edgelist(graph1)
  el2 = as_edgelist(graph2)
  el1[,1] = vmap[match(el1[,1], vmap[,1]), 3]
  el1[,2] = vmap[match(el1[,2], vmap[,1]), 3]
  el2[,1] = vmap[match(el2[,1], vmap[,2]), 3]
  el2[,2] = vmap[match(el2[,2], vmap[,2]), 3]
  ee1 = paste(el1[,1], el1[,2])
  ee2 = paste(el2[,1], el2[,2])

  comb = rbind(el1, el2)
  comb = comb[!duplicated(comb),]
  in1 = paste(comb[,1], comb[,2]) %in% ee1
  in2 = paste(comb[,1], comb[,2]) %in% ee2
  comb = cbind(comb, ifelse(in1 & in2, 'both', ifelse(in1, 'g1', 'g2')))
  comb
}














