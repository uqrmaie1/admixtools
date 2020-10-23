
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

#' Count number of admixture nodes
#'
#' @export
#' @param graph An admixture graph
#' @return Number of admixture nodes
numadmix = function(graph) {
  sum(degree(graph, mode='in') == 2)
}

is_valid = function(graph) {
  igraph::is_simple(graph) &&
    igraph::is_connected(graph) &&
    igraph::is_dag(graph) &&
    sum(igraph::degree(graph, mode = 'in') == 0) == 1
}

is_simplified = function(graph) {
  deg = degree(graph)
  sum(deg == 2) == 1 && max(deg) == 3 && igraph::is_simple(graph)
}

#' Convert data frame graph to igraph
#'
#' @export
#' @param edges An admixture graph as an edge list data frame
#' @return An `igraph` object
edges_to_igraph = function(edges) {
  edges[,1:2] %>% as.matrix %>% igraph::graph_from_edgelist()
}

#' Get the population names of a graph
#'
#' @export
#' @param graph An admixture graph
#' @return Population names
get_leafnames = function(graph) {
  graph %>% V %>% {names(which(degree(graph, ., mode='out') == 0))}
}

get_leaves = function(graph) {
  graph %>% V %>% {.[degree(graph, ., mode='out') == 0]}
}

get_leaves2 = function(graph) {
  # uses 'subcomponent' to return leaves in a consistent order
  root = get_root(graph)
  graph %>% subcomponent(root, mode='out') %>% igraph::intersection(get_leaves(graph))
}


get_root = function(graph) {
  root = V(graph)[igraph::degree(graph, mode = 'in') == 0]
  if(length(root) != 1) stop(paste0('Root problem ', root))
  root
}

#' Get the root name
#'
#' @export
#' @param graph An admixture graph
#' @return Root name
get_rootname = function(graph) {
  names(get_root(graph))
}

#' Get the outgroup from a graph (if it exists)
#'
#' @export
#' @param graph An admixture graph
#' @return Outgroup name
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
  } else if(length(children) == 2) {
    leaves = get_leafnames(graph)
    c1 = distances(graph, children[1], leaves) %>% paste(collapse='')
    c2 = distances(graph, children[2], leaves) %>% paste(collapse='')
    if(c1 < c2) children %<>% rev
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

  root = get_root(graph)
  nammap = unify_vertex_names_rec(g, root, ormap, sep1 = sep1, sep2 = sep2)
  nammap[lv] = names(leaves)
  names(nammap)[match(lv, names(nammap))] = names(leaves)
  nammap['R'] = names(root)
  if(!keep_unique) {
    changed = setdiff(names(nammap), c('R', names(leaves)))
    nammap[changed] = paste0('n', as.numeric(as.factor(nammap[changed])))
  }
  nammap %<>% map_chr(digest::digest) %>% paste0('x', .)
  #set_vertex_attr(graph, 'name', names(nammap), nammap)
  set_vertex_attr(graph, 'name', V(graph), nammap)
}


#' Generate a random admixture graph
#'
#' This function randomly generates an admixture graph for a given set of leaf nodes
#' @export
#' @param leaves Names of the leaf nodes, or a number specifying how many leaf nodes there should be
#' @param numadmix Number of admixture events
#' @param simple Should edges leading to admixture nodes consist of separate admix edges and normal edges
#' @param outpop Outgroup population
#' @examples
#' rand_graph = random_admixturegraph(10, numadmix = 5)
#' plot_graph(rand_graph)
random_admixturegraph = function(leaves, numadmix = 0, simple = FALSE, outpop = NULL) {
  # makes a random admixture graph
  # returns an 'igraph' graph object
  # 'leaves' can be a number of leaf nodes, or a character vector of leaf names

  stopifnot(class(leaves)[1] %in% c('numeric', 'character'))
  if(length(leaves) == 1) leaves = paste0('l', 1:leaves)
  #if(is.null(outpop)) outpop = sample(leaves, 1)
  graph = leaves %>%
    setdiff(outpop) %>%
    sample %>%
    random_newick(outpop = outpop) %>%
    newick_to_edges %>%
    graph_from_edgelist %>%
    insert_admix_multi(numadmix)
  #newick = random_newick(sample(setdiff(leaves, outpop)), outpop = outpop)
  #graph = graph_from_edgelist(newick_to_edges(newick)) %>% insert_admix_multi(numadmix)
  if(!simple) graph %<>% desimplify_graph
  graph
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

  if(is_simplified(graph)) return(graph)
  convmat = FALSE
  if(class(graph) == 'matrix') {
    convmat = TRUE
    graph = graph_from_edgelist(graph)
  }
  graph %<>% igraph::simplify()
  i = 0
  repeat({
    i = i+1
    redundant = which(degree(graph, mode='in') == 1 & degree(graph, mode='out') == 1)
    if(length(redundant) == 0) break
    newfrom = names((neighbors(graph, redundant[1], mode='in')))
    newto = names((neighbors(graph, redundant[1], mode='out')))
    graph %<>% igraph::delete_vertices(redundant[1])
    graph %<>% igraph::add_edges(c(newfrom, newto))
    graph %<>% igraph::simplify()
    if(i > 100) browser()
  })
  graph %<>% X_to_H
  if(convmat) graph = igraph::as_edgelist(graph)
  graph
}

X_to_H = function(graph) {

  crosses = names(which(degree(graph, mode = 'in') == 2 & degree(graph, mode = 'out') == 2))
  for(i in seq_along(crosses)) {
    nam = newnodenam(crosses[i], names(V(graph)))
    parents = names(neighbors(graph, crosses[i], mode = 'in'))
    graph %<>% igraph::add_vertices(1, name = nam) %>%
      add_edges(c(parents[1], nam, parents[2], nam, nam, crosses[i])) %>%
      delete_edges(paste0(parents, '|', crosses[i]))
  }
  graph
}

H_to_X = function(graph) {
  # keeps child names

  adm = names(which(degree(graph, mode = 'in') == 2 & degree(graph, mode = 'out') == 1))
  children = names(neighbors(graph, adm, mode = 'out'))
  ind = which(degree(graph, children, mode = 'out') == 2)
  for(i in seq_along(ind)) {
    parents = names(neighbors(graph, adm[i], mode = 'in'))
    graph %<>% delete_vertices(adm[i]) %>%
      add_edges(c(parents[1], children[i], parents[2], children[i]))
  }
  graph
}

#' Add two nodes before each admixture node
#'
#' This is used to revert simplify_graph.
#' @export
#' @param graph An admixture graph
#' @examples
#' simple = simplify_graph(example_igraph)
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
  newnam = newnodenam(paste0(admixnamesrep, c('a', 'b')), names(V(graph)))
  graph = igraph::add_vertices(graph, length(admix)*2, name = newnam)
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
  if(!'igraph' %in% class(graph)) {
    edges = graph
    graph = edges %>% edges_to_igraph %>% simplify_graph
    edges %<>% filter(type == 'admix') %>%
      left_join(edges %>% transmute(fromfrom = from, from = to), by = 'from')
    pickedges = TRUE
  } else pickedges = FALSE
  nodes = frommap = tomap = names(V(graph))
  names(frommap) = names(tomap) = nodes
  fromnodes = tonodes = c()
  deleted = kept = matrix(NA, 0, 2)
  if(!is_valid(graph)) stop('not valid!')
  root = get_rootname(graph)
  leaves = get_leafnames(graph)
  repeat({
    admix = setdiff(names(which(degree(graph, mode='in') == 2)), root)
    if(length(admix) == 0) break
    adm = admix[1]
    parents = sample(setdiff(sample(names(neighbors(graph, adm, mode='in'))), root))
    if(pickedges) {
      p1 = edges %>% filter(to == adm, from %in% parents) %>% slice_min(weight, with_ties = FALSE) %>% pull(from)
      if(length(p1) == 0) p1 = edges %>% filter(to == adm, fromfrom %in% parents) %>% slice_max(weight, with_ties = FALSE) %>% pull(fromfrom)
    } else {
      p1 = parents[1]
    }
    fromnode0 = names(neighbors(graph, p1, mode = 'in'))
    i = 0
    while(length(fromnode0) > 1) {
      i = i+1
      adm = p1
      parents = sample(setdiff(names(neighbors(graph, adm, mode='in')), root))
      if(pickedges) {
        p1 = edges %>% filter(to == adm, from %in% parents) %>% slice_min(weight, with_ties = FALSE) %>% pull(from)
        if(length(p1) == 0) p1 = edges %>% filter(to == adm, fromfrom %in% parents) %>% slice_max(weight, with_ties = FALSE) %>% pull(fromfrom)
      } else {
        p1 = parents[1]
      }
      if(length(parents) == 0 || !p1 %in% names(V(graph)) || root %in% parents) browser()
      fromnode0 = names(neighbors(graph, p1, mode = 'in'))
      if(i > 100) browser()
    }
    # stopifnot(!any(c(parents[1], admix[1]) %in% unlist(admixedges)))
    # too strict; probably won't be able to always reintroduce admixedges at appropriate locations after SPR
    fromnode = setdiff(names(neighbors(graph, p1, mode = 'out')), adm)
    tonode0 = setdiff(names(neighbors(graph, adm, mode = 'in')), p1)
    tonode = names(neighbors(graph, adm, mode = 'out'))
    fromnodes %<>% c(fromnode)
    tonodes %<>% c(tonode)
    # fromnodes %<>% c(fromnode0)
    # tonodes %<>% c(tonode0)
    #if(adm %in% fromnodes) fromnodes[fromnodes == adm] = fromnode
    #if(adm %in% tonodes) tonodes[tonodes == adm] = tonode
    #if(adm == 'pPygmy3' || p1 == 'pPygmy3') browser()
    graphnew = igraph::delete_vertices(graph, adm)
    graphnew = igraph::add_edges(graphnew, c(tonode0, tonode))
    if(tonode %in% admix) frommap[tonode0] = adm
    #if(delete_nodes) {
      graphnew = igraph::delete_vertices(graphnew, p1)
      if(p1 %in% tonodes) tonodes[tonodes == p1] = tonode
      if(p1 %in% fromnodes) fromnodes[fromnodes == p1] = fromnode
      if(length(c(fromnode0, fromnode)) %% 2 != 0) browser()
      graphnew = igraph::add_edges(graphnew, c(fromnode0, fromnode)) %>% simplify_graph()
    #}
    if(is_valid(graphnew) && numadmix(graphnew) < numadmix(graph)) graph = graphnew
    else browser()
    deleted = rbind(deleted, c(frommap[p1], adm))
    if(frommap[p1] != p1) {
      deleted = rbind(deleted, c(p1, frommap[p1]))
    }
    kept = rbind(kept, c(setdiff(parents, p1), adm))
  })
  tree = graph
  if(max(degree(tree, mode = 'in')) > 1) stop('Failed to make tree!')
  namedList(tree, fromnodes, tonodes, deleted, kept)
}

shorten_admixed_leaves = function(graph) {

  # keep terminal branch only for leaves which are non-admixed

  leaves = get_leafnames(graph)
  leaves %<>% intersect(names(which(degree(graph, leaves, mode = 'in') == 1)))
  parents = map_chr(leaves, ~names(neighbors(graph, ., mode = 'in')))
  ind = which(degree(graph, parents, mode = 'in') == 2)
  graph %>%
    delete_vertices(leaves[ind]) %>%
    set_vertex_attr('name', parents[ind], leaves[ind])
}


newnodenam = function(newnam, current) {
  # this function takes a vector of proposed new node names (newnam), checks if they already exist,
  # and if so, returns unique newnams with newnam as prefix
  for(i in seq_along(newnam)) {
    mod = newnam[i]
    while(mod %in% current) {
      mod %<>% paste0(sample(letters, 1))
    }
    newnam[i] = mod
    current %<>% c(mod)
  }
  newnam
}


#' Generate a random binary graph
#' @export
#' @param n The number of terminal nodes, or a vector of population labels.
#' @param start Prefix.
#' @param end Postfix.
#' @param outpop Outgroup (if \code{n} is a vector of labels).
#' @return Tree in newick format.
#' @examples
#' random_newick(5)
#' random_newick(c('a', 'b', 'c', 'd')) # toplogy random, but pop order fixed
#' random_newick(sample(c('a', 'b', 'c', 'd'))) # toplogy and pop order random
random_newick = function(n, start = '', end = '', outpop = NULL) {
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

random_newick_named = function(names, start = '', end = '') {
  # recursive function which returns a labelled random, binary tree in newick format with named leaves in order of input
  n = length(names)
  if(n == 1) return(names)
  n1 = sample(n-1,1)
  return(paste0(start, '(',random_newick_named(names[1:n1]),',',random_newick_named(names[(n1+1):n]),')', end))
}

#' Turn a newick format tree to a matrix of edges
#' @export
#' @param newick Tree in newick format.
#' @param node Root label of the tree.
#' @param edgemat Used for recursive function calls.
#' @return Tree as two column matrix of edges (adjacency list)
#' @examples
#' newick = random_newick(c('a', 'b', 'c', 'd'))
#' newick
#' newick_to_edges(newick)
newick_to_edges = function(newick, node = 'R', edgemat = matrix(NA,0,2)) {
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


#' Insert admixture edges into graph
#'
#' @param graph An admixture graph
#' @param fromnodes List of nodes. New edges will originate above these nodes.
#' @param tonodes List of nodes. New edges will end above these nodes.
#' @param substitute_missing If `TRUE`, an attempt will be made to insert random other edges
#'   if some of the provided edges could not be inserted.
#' @param allow_below_admix Allow insertion of edges which begin or end directly underneath an admixture node
#' @return Adxmiture graph with inserted edges
insert_admix_old = function(graph, fromnodes, tonodes, substitute_missing = FALSE,
                               allow_below_admix = FALSE) {
  # inserts edges fromnodes -> tonodes into graph
  # inverse of 'split_graph'
  # stopifnot(all(c(sapply(admixedges, `[`, c(2, 4))) %in% names(V(graph))))
  # some nodes will get lost when splitting; insert random amix edges instead

  stopifnot(length(fromnodes) == length(tonodes))
  leaves = sort(get_leafnames(graph))
  desimplify = !is_simplified(graph)
  if(desimplify) graph %<>% simplify_graph

  miss = 0
  maxprev = max(degree(graph))
  for(i in rev(seq_len(length(fromnodes)))) {
    fromnode = fromnodes[i]
    tonode = tonodes[i]
    if(!fromnode %in% names(V(graph)) ||
       !tonode %in% names(V(graph)) ||
       is.finite(igraph::distances(graph, tonode, fromnode, mode = 'out'))) {
      miss = miss+1
      next
    }
    toedge2parent = neighbors(graph, tonode, mode='in')
    toedge2sib = setdiff(names(neighbors(graph, toedge2parent, mode='out')), tonode)
    if(!allow_below_admix) {
      if(length(toedge2sib) == 0 ||
         degree(graph, fromnode, mode='in') > 1 ||
         degree(graph, toedge2sib, mode='in') > 1) {
        miss = miss+1
        next
      }
    }
    fromnode_parent = names(neighbors(graph, fromnode, mode='in'))
    tonode_parent = names(neighbors(graph, tonode, mode='in'))
    if(length(fromnode_parent) != 1 || length(tonode_parent) != 1) {
      miss = miss+1
      next
    }
    graph = igraph::delete_edges(graph, c(paste(fromnode_parent, fromnode, sep='|'),
                                          paste(tonode_parent, tonode, sep='|')))
    newnam = newnodenam(c(fromnode_parent, 'admix'), names(V(graph)))
    new1 = newnam[1]
    new2 = newnam[2]
    graph = igraph::add_vertices(graph, 2, name = newnam)
    newedges = c(fromnode_parent, new1, new1, fromnode, tonode_parent, new2, new2, tonode, new1, new2)
    graph = igraph::add_edges(graph, newedges)
    stopifnot(max(degree(graph)) <= maxprev)
  }
  if (miss > 0) {
    if(substitute_missing) graph = insert_admix_random(graph, miss)
    else warning('did not insert or substitute admixedge')
  }
  if(desimplify) graph %<>% simplify_graph %>% desimplify_graph
  if(!isTRUE(all.equal(leaves, sort(get_leafnames(graph))))) browser()
  graph
}


insert_admix_random = function(graph, nadmix) {
  # graph should be simplified
  for(i in seq_len(nadmix)) {
    root = get_root(graph)
    firstgen = neighbors(graph, root, mode='out')
    admix = V(graph)[(degree(graph, mode='in') > 1)]
    admixchildren = do.call(c, igraph::ego(graph, 1, admix, mode='out'))
    admixparents = do.call(c, igraph::ego(graph, 1, admix, mode='in'))
    admixsiblings = do.call(c, igraph::ego(graph, 1, admixparents, mode='out'))
    exclude = c(admix, firstgen, root)
    if(!is.null(admixchildren)) exclude = c(exclude, admixchildren)
    cand_from = sample(igraph::difference(V(graph), exclude))
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
    graphnew = insert_admix_old(graph, fromnode, tonode, substitute_missing = FALSE)
    stopifnot(numadmix(graphnew) > numadmix(graph))
    stopifnot(igraph::is_simple(graphnew))
    stopifnot(igraph::is_dag(graphnew))
    graph = graphnew
  }
  graph
}


#' Insert a single edge into graph
#'
#' @export
#' @param graph An admixture graph
#' @param source_from Parent node of the source edge
#' @param source_to Child node of the source edge
#' @param dest_from Parent node of the destination edge
#' @param dest_to Child node of the destination edge
#' @param substitute Should another edge be inserted, if the one specified doesn't work?
#' @return Adxmiture graph with inserted edge
#' @seealso \code{\link{delete_admix}}, \code{\link{insert_admix_multi}}
insert_admix = function(graph, source_from = NULL, source_to = NULL, dest_from = NULL, dest_to = NULL, substitute = FALSE) {
  # assumes graph is simplified

  leaves = sort(get_leafnames(graph))
  nodes = names(V(graph))
  if(is.null(source_to) || is.null(dest_to) || length(setdiff(c(source_to, dest_to), nodes)) > 0 && substitute) {
    e = graph %>% find_newedges %>% slice_sample
    source_from = e$source_from
    dest_from = e$dest_from
    source_to = e$source_to
    dest_to = e$dest_to
  } else if(is.null(source_from) || is.null(dest_from) || length(setdiff(c(source_from, dest_from), nodes)) > 0 && substitute) {
    getfrom = function(x) sample(names(neighbors(graph, x, mode = 'in')))[1]
    source_from = getfrom(source_to)
    dest_from = getfrom(dest_to)
  }
  nam = newnodenam(c(source_from, 'admix'), names(V(graph)))
  newgraph = graph %>% add_vertices(2, attr = list(name = nam)) %>%
    add_edges(c(source_from, nam[1], nam[1], source_to, dest_from, nam[2], nam[2], dest_to, nam[1], nam[2])) %>%
    delete_edges(paste(c(source_from, dest_from), c(source_to, dest_to), sep = '|'))
  leaves2 = sort(get_leafnames(newgraph))
  if(!is_valid(newgraph) || !isTRUE(all.equal(leaves2, leaves))) {
    if(substitute) graph %<>% insert_admix(substitute = FALSE)
    else stop("Inserting edge failed!")
  } else graph = newgraph
  graph
}


#' Insert admixture edges into graph
#'
#' @export
#' @param graph An admixture graph
#' @param from List of nodes. New edges will originate above these nodes.
#' @param to List of nodes. New edges will end above these nodes.
#' @param substitute Should another edge be inserted, if the one specified doesn't work?
#' @return Admixture graph with inserted edges
#' @seealso \code{\link{insert_admix}} \code{\link{delete_admix}}
insert_admix_multi = function(graph, n = 1, source_from = NULL, source_to = NULL, dest_from = NULL, dest_to = NULL, substitute = FALSE) {

  # if(all(map_lgl(list(source_from, source_to, dest_from, dest_to), is.null))) {
  #   f = function(x, y) insert_admix(x)
  # } else {
  #   if(is.null(source_from) && is.null(dest_from)) {
  #     # implement this so that insert_admix_old can be replaced
  #     getfrom = ~sample(names(neighbors(graph, ., mode = 'in')))[1]
  #     source_from = map_chr(source_to, getfrom)
  #     dest_from = map_chr(dest_to, getfrom)
  #   }
  #   n = length(source_from)
  #   f = function(x, y) insert_admix(x, source_from[y], source_to[y], dest_from[y], dest_to[y])
  # }
  if(!all(map_lgl(list(source_from, source_to, dest_from, dest_to), is.null))) n = length(source_to)
  f = function(x, y) insert_admix(x, source_from[y], source_to[y], dest_from[y], dest_to[y], substitute = substitute)
  seq_len(n) %>% reduce(f, .init = graph)
}

subtree_prune_and_regraft = function(graph, only_leaves = FALSE, fix_outgroup = TRUE) {
  # cuts of a parts of the tree and attaches it at a random location
  # root -> outgroup stays fixed
  root = get_root(graph)
  stopifnot(degree(graph, root, mode='in') == 0)
  stopifnot(max(degree(graph, V(graph), mode='in')) == 1)

  if(fix_outgroup) firstgen = neighbors(graph, root, mode='out')
  else firstgen = c()
  i = 0
  repeat({
    i = i+1
    excl = c(firstgen, root)
    if(only_leaves) excl = c(excl, V(graph)[degree(graph, mode='out') > 0])
    cutnode = sample(igraph::difference(V(graph), excl), 1)
    cutnodes = subcomponent(graph, cutnode, mode='out')
    cutparent = neighbors(graph, cutnode, mode='in')
    cutgrandparent = neighbors(graph, cutparent, mode='in')
    cutsibling = igraph::difference(neighbors(graph, cutparent, mode='out'), cutnode)
    hostnodes = igraph::difference(V(graph), c(root, cutnodes, cutparent, cutsibling, firstgen))
    if(length(hostnodes) > 0) break
    if(i > 100) browser()
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
  stopifnot(igraph::is_dag(graph))
  graph
}


admixturegraph_prune_and_regraft = function(graph, only_leaves = FALSE, fix_outgroup = TRUE) {
  # 1. remove admixture edges (one randomly selected from each admixture node)
  # 2. runs subtree_prune_and_regraft on resulting tree
  # 3. add admixture edges back on
  if(numadmix(graph) == 0) return(subtree_prune_and_regraft(graph, only_leaves = only_leaves, fix_outgroup = fix_outgroup))
  o = graph
  desimplify = !is_simplified(graph)
  if(desimplify) graph %<>% simplify_graph
  if(!is_valid(graph)) browser()
  spl = split_graph(graph)
  graph = subtree_prune_and_regraft(spl$tree, only_leaves = only_leaves, fix_outgroup = fix_outgroup)
  #graph = insert_admix_old(graph, spl$fromnodes, spl$tonodes, substitute_missing = TRUE)
  graph = insert_admix_multi(graph, source_to = spl$fromnodes, dest_to = spl$tonodes, substitute = TRUE)
  stopifnot(is_valid(graph))
  if(desimplify) graph = desimplify_graph(graph)
  graph
}

#' Modify a graph by regrafting a leaf
#'
#' @export
#' @param graph An admixture graph
#' @return A new admixture graph
spr_leaves = function(graph, fix_outgroup = TRUE)
  admixturegraph_prune_and_regraft(graph, only_leaves = TRUE, fix_outgroup = fix_outgroup)

#' Modify a graph by regrafting a subcomponent
#'
#' @export
#' @param graph An admixture graph
#' @return A new admixture graph
spr_all = function(graph, fix_outgroup = TRUE)
  admixturegraph_prune_and_regraft(graph, only_leaves = FALSE, fix_outgroup = fix_outgroup)

#' Modify a graph by moving an admixture edge
#'
#' @export
#' @param graph An admixture graph
#' @param fix_outgroup Keep outgroup in place
#' @return A new admixture graph
move_admixedge_once = function(graph, fix_outgroup = TRUE) {
  # selects random admixture edge, and moves it to next closest possible spot
  # if not possible, select other admix node
  # if none possible, print warning
  desimplify = !is_simplified(graph)
  if(desimplify) graph %<>% simplify_graph
  admix = sample(names(V(graph)[(degree(graph, mode='in') > 1)]))
  nadmix = length(admix)
  root = get_rootname(graph)
  if(fix_outgroup) firstgen = names(neighbors(graph, root, mode='out'))
  else firstgen = ''
  for(i in seq_len(nadmix)) {
    parents = sample(names(neighbors(graph, admix[i], mode='in')))
    for(j in 1:2) {
      sibling = setdiff(names(neighbors(graph, parents[j], mode='out')), admix[i])
      grandparent = names(neighbors(graph, parents[j], mode='in'))
      if(length(sibling) == 0 || length(grandparent) == 0) next
      for(n in names(subcomponent(graph, parents[j]))) {
        # subcomponent will return nodes ordered by distance to n
        if(!n %in% root &&
           !n %in% firstgen &&
           !n %in% parents &&
           !n %in% grandparent &&
           !n %in% admix &&
           !n %in% c(names(neighbors(graph, parents[1], mode='out')),
                     names(neighbors(graph, parents[2], mode='out'))) &
           !n %in% names(subcomponent(graph, admix[i], mode='out'))) {
          newgrandparent = names(neighbors(graph, n, mode='in'))
          if(length(c(grandparent, sibling, newgrandparent, parents[j], parents[j], n)) %% 2 != 0) browser()
          graph = igraph::add_edges(graph, c(grandparent, sibling, newgrandparent, parents[j], parents[j], n))
          graph = igraph::delete_edges(graph, c(paste(newgrandparent, n, sep='|'),
                                              paste(grandparent, parents[j], sep='|'),
                                              paste(parents[j], sibling, sep='|')))
          graph %<>% simplify_graph
          stopifnot(igraph::is_dag(graph))
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

#' Modify a graph by permuting leaf nodes
#'
#' @export
#' @param graph An admixture graph
#' @param fix_outgroup Keep outgroup in place
#' @return A new admixture graph
permute_leaves = function(graph, fix_outgroup = TRUE) {

  leaves = V(graph)[degree(graph, v = V(graph), mode = c('out')) == 0]
  if(fix_outgroup) leaves = leaves[-1]
  nam = names(leaves)
  igraph::set_vertex_attr(graph, 'name', leaves, sample(nam))
}

#' Modify a graph by swapping two leaf nodes
#'
#' @export
#' @param graph An admixture graph
#' @param fix_outgroup Keep outgroup in place
#' @return A new admixture graph
swap_leaves = function(graph, fix_outgroup = TRUE) {

  leaves = V(graph)[degree(graph, v = V(graph), mode = c('out')) == 0]
  if(fix_outgroup) leaves = leaves[-1]
  leaves = sample(leaves)[1:2]
  nam = names(leaves)
  igraph::set_vertex_attr(graph, 'name', leaves, sample(nam))
}

#' Modify a graph flipping the direction of an admixture edge
#'
#' @export
#' @param graph An admixture graph
#' @param fix_outgroup Keep outgroup in place (has no effect here)
#' @return A new admixture graph
flipadmix_random = function(graph, fix_outgroup = TRUE) {

  desimplify = !is_simplified(graph)
  if(desimplify) graph %<>% simplify_graph
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

#' Modify a graph by applying n mutation functions
#'
#' @export
#' @param graph An admixture graph
#' @param n Number of functions to apply
#' @param funs List of function from which to choose
#' @param fix_outgroup Keep outgroup in place
#' @return A new admixture graph
mutate_n = function(graph, n = 2, funs = list(spr_all, spr_leaves, swap_leaves, move_admixedge_once,
                                              flipadmix_random, mutate_n), fix_outgroup = TRUE) {

  funsn = sample(funs, n, replace = TRUE) %>%
    map(~function(x) .x(x, fix_outgroup=fix_outgroup))
  purrr::compose(!!!funsn)(graph)
}


get_mutfuns = function(mutfuns, probs, fix_outgroup = TRUE) {
  # return a list of list of graph mutation functions
  # outer list is for each generation, inner list is for each graph, and is named
  # probs is a matrix numgen x numgraphs with probabilities

  if(is.null(names(mutfuns))) names(mutfuns) = paste0('fun', seq_along(mutfuns))
  if(length(unique(names(mutfuns))) < length(mutfuns)) stop("Names in 'mutfuns' are not unique!")
  map(seq_len(nrow(probs)), ~{
    nam = sample(names(mutfuns), ncol(probs), replace = TRUE, prob = probs[.,])
    map(mutfuns, ~function(graph) .(graph, fix_outgroup = fix_outgroup)) %>%
      set_names(nam)
  })
}


add_generation = function(models, numgraphs, numsel, qpgfun, mutfuns, opt_worst_residual = FALSE, parallel = TRUE, verbose = TRUE) {

  space = paste0(paste(rep(' ', 50), collapse=''), '\r')
  if(verbose) alert_info(paste0('Selecting winners...', space))
  lastgen = max(models$generation)
  oldmodels = models %>% filter(generation == lastgen, !is.na(score))
  numsel = min(nrow(oldmodels), numsel)
  numeach = ceiling(numgraphs/numsel)
  winners = oldmodels %>% mutate(prob = score_to_prob(score)) %>% sample_n(numsel-1, weight=prob) %>% select(-prob)
  optvar = if(opt_worst_residual) 'worst_residual' else 'score'
  winners %<>% bind_rows(oldmodels %>% slice_min(.data[[optvar]], with_ties = FALSE))
  sq = (numsel+1):numgraphs
  newmodels = tibble(generation = lastgen+1, index = sq,
                     igraph = rep(winners$igraph, numeach)[sq],
                     oldscore = rep(winners$score, numeach)[sq],
                     oldindex = rep(winners$index, numeach)[sq])
  mutations = sample(mutfuns, numgraphs-numsel, replace = TRUE)
  if(parallel) {
    map = furrr::future_map
    imap = furrr::future_imap
  }
  if(verbose) alert_info(paste0('Generating new graphs...', space))
  newmodels %<>% mutate(mutation = names(mutations),
                        igraph = imap(igraph, ~tryCatch(exec(mutations[[.y]], .x), error = function(e) .x)))
  if(verbose) alert_info(paste0('Evaluating graphs...', space))
  #if(lastgen == 10) browser()
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


evolve_topology = function(init, numgraphs, numgen, numsel, qpgfun, mutlist, repnum,
                           opt_worst_residual = FALSE, parallel = TRUE,
                           keep = 'all', store_intermediate = NULL, stop_after = NULL, verbose = TRUE) {
  out = init
  stop_at = Sys.time() + stop_after
  for(i in seq_len(numgen)) {
    if(!is.null(stop_after) && Sys.time() > stop_at) break
    newmodels = add_generation(out, numgraphs, numsel, qpgfun, mutlist[[i]], opt_worst_residual = opt_worst_residual,
                               parallel = parallel, verbose = verbose)
    if(!is.null(store_intermediate)) {
      fl = paste0(store_intermediate, '_rep', repnum, '_gen', i, '.rds')
      saveRDS(newmodels, fl)
    }
    optvar = if(opt_worst_residual) 'worst_residual' else 'score'
    if(keep == 'all') out %<>% bind_rows(newmodels)
    else if(keep == 'best') out %<>% group_by(generation) %>% slice_min(.data[[optvar]], with_ties = FALSE) %>%
      ungroup %>% bind_rows(newmodels)
    else out = newmodels
    if(verbose) {
      #best = newmodels %>% filter(index > numsel) %>% top_n(numsel, -jitter(score, amount = 1e-9)) %$% score %>% sort
      best = newmodels %>% slice_min(.data[[optvar]], n = numsel, with_ties = FALSE) %>% pull(.data[[optvar]]) %>% sort
      alert_success(paste0('Generation ', i, '  Best ', optvar,'s: ', paste(round(best, 3), collapse=', '),
                    paste(rep(' ', 30), collapse=''), '\n'))
    }
  }
  if(verbose) cat('\n')
  if(keep == 'best') out %<>% group_by(generation) %>% slice_min(.data[[optvar]], with_ties = FALSE) %>% ungroup
  out
}


optimize_admixturegraph_single = function(pops, precomp, mutlist, repnum, numgraphs = 50, numgen = 5,
                                          numsel = 5, numadmix = 0, numstart = 1, outpop = NA, initgraphs = NULL,
                                          opt_worst_residual = FALSE,
                                          parallel = TRUE, stop_after = NULL, store_intermediate = NULL,
                                          keep = 'all', verbose = TRUE, ...) {

  kp = c('edges', 'score', 'f3')
  if(opt_worst_residual) {
    kp = c(kp, 'worst_residual')
    qpgfun = function(graph) qpgraph(data = precomp, graph = graph,
                                     numstart = numstart, return_f4 = opt_worst_residual, verbose = FALSE, ...)[kp]
  } else {
    qpgfun = function(graph) qpgraph(data = NULL, graph = graph, f3precomp = precomp,
                                     numstart = numstart, return_f4 = opt_worst_residual, verbose = FALSE, ...)[kp]
  }

  # qpgfun = possibly(qpgfun, otherwise = NULL)
  space = paste0(paste(rep(' ', 50), collapse=''), '\r')
  if(verbose) alert_info(paste0('Generate new graphs...', space))
  if(is.null(initgraphs)) initgraphs = replicate(numgraphs, random_admixturegraph(pops, numadmix, outpop = outpop),
                                                simplify = FALSE)
  else initgraphs = initgraphs[round(seq(1, length(initgraphs), numgraphs))]
  if(verbose) alert_info(paste0('Evaluate graphs...', space))
  if(parallel) map = furrr::future_map
  init = tibble(generation=0, index = seq_len(numgraphs),
                igraph = initgraphs, mutation = 'random_admixturegraph') %>%
    mutate(out = map(igraph, qpgfun), isn = map_lgl(out, is.null))
  if(all(init$isn)) stop('All NULL!')

  init %<>% select(-isn) %>% unnest_wider(out) %>% mutate(oldscore = score, oldindex = index)
  if(verbose) {
    optvar = if(opt_worst_residual) 'worst_residual' else 'score'
    best = init %>% filter(!is.na(score)) %>% slice_min(.data[[optvar]], n = min(numsel, nrow(.data)), with_ties = FALSE) %>% pull(.data[[optvar]])
    alert_success(paste0('Generation 0  Best ', optvar, 's: ', paste(round(best), collapse=', '), space, '\n'))
  }

  evolve_topology(init, numgraphs, numgen, numsel, qpgfun, mutlist, repnum, opt_worst_residual = opt_worst_residual,
                  parallel = parallel,
                  keep = keep, store_intermediate = store_intermediate, stop_after = stop_after, verbose = verbose)
}



#' Find well fitting admixture graphs
#'
#' This function generates and evaluates admixture graphs in `numgen` iterations across `numrep` independent repeats
#' to find well fitting admixturegraphs. It uses the function \code{\link[furrr]{future_map}}
#' to parallelize across the independent repeats. The function \code{\link[future]{plan}} can be called
#' to specify the details of the parallelization. This can be used to parallelize across cores or across nodes on
#' a compute cluster. Setting `numadmix` to 0 will search for well fitting trees, which is much faster than searching
#' for admixture graphs with many admixture nodes.
#' @export
#' @param data Input data in one of three forms:
#' \enumerate{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}} (fastest option)
#' \item A directory which contains pre-computed f2-statistics
#' \item The prefix of genotype files (slowest option)
#' }
#' @param pops Populations for which to fit admixture graphs (default all)
#' @param outpop An outgroup population which will split at the root from all other populations in all tested graphs. If one of the populations is know to be an outgroup, designating it as `outpop` will greatly reduce the search space compared to including it and not designating it as `outpop`.
#' @param numrep Number of independent repetitions (each repetition can be run in parallel)
#' @param numgraphs Number of graphs in each generation
#' @param numgen Number of generations
#' @param numsel Number of graphs which are selected in each generation. Should be less than `numgraphs`.
#' @param numadmix Number of admixture events within each graph
#' @param numstart Number of random initializations in each call to `qpgraph`. Defaults to 1, to speed up the graph optimization.
#' @param keep Which models should be returned. One of `all`, `best`, `last`
#' \itemize{
#' \item `all` (default): Return all evaluated graphs
#' \item `best`: Return only the best fitting graph from each repeat and each generation
#' \item `last`: Return all graphs from the last generation
#' }
#' @param initgraphs Optional graph or list of igraphs to start with. If `NULL`, optimization will start with random graphs.
#' @param mutfuns Functions used to modify graphs. Defaults to the following:
#' \itemize{
#' \item \code{\link{spr_leaves}}: Subtree prune and regraft leaves. Cuts a leaf node and attaches it
#' to a random other edge in the graph.
#' \item \code{\link{spr_all}}: Subtree prune and regraft. Cuts any edge and attaches the new orphan node
#' to a random other edge in the graph, keeping the number of admixture nodes constant.
#' \item \code{\link{swap_leaves}}: Swaps two leaf nodes.
#' \item \code{\link{move_admixedge_once}}: Moves an admixture edge to a nearby location.
#' \item \code{\link{flipadmix_random}}: Flips the direction of an admixture edge (if possible).
#' \item \code{\link{mutate_n}}: Apply `n` of the mutation functions in this list to a graph (defaults to 2).
#' }
#' See examples for how to make new mutation functions.
#' @param mutprobs Relative frequencies of each mutation function.
#' \itemize{
#' \item `NULL` (default) means each mutation function is picked with equal probability
#' \item A numeric vector of length equal to `mutfuns` defines the relative frequency of each mutation function
#' \item A matrix of dimensions `numgen` x `length(mutfuns)` defines the relative frequency of each mutation function in each generation
#' }
#' @param opt_worst_residual Optimize for lowest worst residual instead of best score. `FALSE` by default, because the likelihood score is generally a better indicator of the quality of the model fit. Optimizing for the lowest worst residual is also slower (because f4-statistics need to be computed).
#' @param store_intermediate Path and prefix of files for intermediate results to `.rds`. Can be useful if `find_graphs` doesn't finish sucessfully.
#' @param parallel Parallelize over repeats (if `numrep > 1`) or graphs (if `numrep == 1`) by replacing \code{\link[purrr]{map}} with \code{\link[furrr]{future_map}}. Will only be effective if \code{\link[future]{plan}} has been set.
#' @param stop_after Stop optimization after `stop_after` seconds (and after finishing the current generation).
#' @param verbose Print progress updates
#' @param ... Additional arguments passed to \code{\link{qpgraph}}
#' @return A nested data frame with one model per line
#' @seealso \code{\link{qpgraph}}
#' @examples
#' \dontrun{
#' find_graphs(example_f2_blocks, numrep = 200, numgraphs = 100,
#'             numgen = 20, numsel = 5, numadmix = 3)
#' }
#' \dontrun{
#' # Making new mutation functions by modifying or combining existing ones:
#' newfun1 = function(graph, ...) mutate_n(graph, 3, ...)
#' newfun2 = function(graph, ...) flipadmix_random(spr_leaves(graph, ...), ...)
#' find_graphs(f2_blocks, mutfuns = namedList(spr_leaves, newfun1, newfun2), mutprobs = c(0.2, 0.3, 0.5))
#' }
find_graphs = function(data, pops = NULL, outpop = NULL, numrep = 1, numgraphs = 50,
                       numgen = 5, numsel = 5, numadmix = 0, numstart = 1, keep = c('all', 'best', 'last'), initgraphs = NULL,
                       mutfuns = namedList(spr_leaves, spr_all, swap_leaves, move_admixedge_once, flipadmix_random, mutate_n),
                       mutprobs = NULL, opt_worst_residual = FALSE, store_intermediate = NULL, parallel = TRUE, stop_after = NULL, verbose = TRUE, ...) {

  keep = rlang::arg_match(keep)
  if(numsel >= numgraphs || numsel < 1) stop("'numsel' has to be smaller than 'numgraphs' and greater than 0!")
  if(!is.null(pops) && !is.null(initgraphs)) stop("You can't provide 'pops' and 'initgraphs' at the same time!")

  if(!is.null(initgraphs)) {
    if('data.frame' %in% class(initgraphs) || 'matrix' %in% class(initgraphs)) {
      initgraphs = graph_from_edgelist(as.matrix(initgraphs[,1:2])) %>% list
    } else if('character' %in% class(initgraphs)) {
      initgraphs = graph_from_edgelist(as.matrix(parse_qpgraph_graphfile(initgraphs)[,1:2])) %>% list
    } else if('igraph' %in% class(initgraphs)) {
      initgraphs %<>% list
    } else if('list' %in% class(initgraphs)) {
    } else stop("'initgraphs' should be a single graph or a list of 'igraph' objects")
  }

  if(is.null(pops)) {
    if(!is.null(initgraphs)) {
      pops = get_leafnames(initgraphs[[1]])
    } else if(is.array(data)) {
      pops = dimnames(data)[[1]]
    } else stop('Please provide population names!')
  }

  if(!opt_worst_residual) precomp = qpgraph_precompute_f3(data, pops, f3basepop = outpop, ...)
  else precomp = get_f2(data, pops, ...)

  if(class(mutfuns[[1]]) == 'character') mutfuns %<>% rlang::set_names() %>% map(get)
  if(is.null(mutprobs)) {
    if(numadmix == 0 && is.null(initgraphs)) mutfuns[c('move_admixedge_once', 'flipadmix_random')] = NULL
    mutprobs = matrix(1, numgen, length(mutfuns)) %>% set_colnames(names(mutfuns))
  } else if(!is.matrix(mutprobs)) {
    if(length(mutprobs) != length(mutfuns)) stop("'mutfuns' and 'mutprobs' don't match")
    mutprobs = t(replicate(numgen, mutprobs)) %>% set_colnames(names(mutfuns))
  } else {
    stopifnot(nrow(mutprobs) == numgen)
    stopifnot(ncol(mutprobs) == length(mutfuns))
    if(!isTRUE(all.equal(sort(colnames(mutprobs)), sort(names(mutfuns))))) mutprobs %>% set_colnames(names(mutfuns))
  }
  mutlist = get_mutfuns(mutfuns, mutprobs, fix_outgroup = !is.null(outpop))

  oa = function(i) optimize_admixturegraph_single(pops, precomp, mutlist = mutlist, repnum = i,
                                                 numgen = numgen, numsel = numsel,
                                                 numgraphs = numgraphs, numadmix = numadmix, numstart = numstart,
                                                 outpop = outpop, initgraphs = initgraphs, opt_worst_residual = opt_worst_residual,
                                                 parallel = parallel && numrep == 1,
                                                 stop_after = stop_after, store_intermediate = store_intermediate,
                                                 keep = keep, verbose = verbose && numrep == 1, ...)
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

  paths = all_simple_paths(graph, get_root(graph), leaves, mode = 'out') %>%
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
    summarize(x12 = any(topo == 0) & !any(topo == 1),
              x21 = any(topo == 0) & any(topo == 1),
              x13 = any(topo == 1) & !any(topo == 0),
              x31 = any(topo == 1) & any(topo == 0),
              x23 = any(topo == 2) & !any(topo == 0),
              x32 = any(topo == 2) & any(topo == 0)) %>%
    group_by(name1, name2, name3) %>%
    #summarize(x13 = any(x13), x23 = any(x23), x31 = any(x31), x32 = any(x32), x12 = any(x12), x21 = any(x21),
    summarize(across(starts_with('x'), any),
              toposet = paste0(x12+0, x21+0, x13+0, x31+0, x23+0, x32+0, collapse='')) %>%
    ungroup %>%
    mutate(tlr = !x13 & !x31 & (x23 | x32) & (x12 | x21))

  tripletopo
}


#' Summarize triples across graphs
#'
#' This function summarizes topologies of population triples across graphs. If the list of graphs comes from \code{\link{find_graphs}}, only one graph from each run should be included, as graphs within the same run are not independent of each other.
#'
#' @export
#' @param graphs A list of graphs
#' @return A data frame with one line for each population triple and columns summarizing the relationship of each triple across graphs.
#' \itemize{
#' \item `clade12`: Fraction of graphs in which `pop1` and `pop2` form a clade
#' \item `x12`: Fraction of graphs in which `pop1` admixes into `pop2` at the exclusion of `pop3`
#' \item `toptopo`: A binary string representation of the most common topology. Digits represent `x12`, `x21`, `x13`, `x31`, `x23`, `x32`
#' \item `toptopocnt`: The number of graphs in which `toptopo` was observed
#' \item `topos`: The counts of all topologies
#' }
summarize_triples = function(graphs) {
  # results is output from 'find_graphs'
  # takes at most one graph from each independent run

  tibble(graph = graphs) %>%
    rowwise %>%
    mutate(topo = list(summarize_graph(graph))) %>%
    unnest(topo) %>%
    group_by(name1, name2, name3, toposet) %>%
    mutate(cnt = n()) %>%
    group_by(name1, name2, name3) %>%
    mutate(clade12 = !(x13|x31|x23|x32),
           clade13 = !(x12|x21|x23|x32),
           clade23 = !(x12|x21|x13|x31)) %>%
    summarize(across(starts_with(c('x','clade')), mean),
              #clade = mean(substr(toposet, 1, 4) == '0000'),
              toptopo = toposet[cnt = max(cnt)][1],
              toptopocnt = max(cnt),
              topos = list(setNames(cnt, toposet))) %>%
    ungroup %>%
    rename(pop1 = name1, pop2 = name2, pop3 = name3)
}


#' Find identical graphs
#'
#' @param igraphlist A list with admixture graphs
#' @return An integer vector with isomorphism classes.
#' Graphs with the same number have identical topology (but may have different labels).
#' @export
#' @seealso \code{\link{isomorphism_classes2}}
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

#' Find identical graphs
#'
#' @param igraphlist A list with admixture graphs
#' @return An integer vector with isomorphism classes.
#' Graphs with the same number have identical topology and leaf labels (but may have different internal labels).
#' @export
#' @seealso \code{\link{isomorphism_classes}}
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
#' @param graph An admixture graph as igraph object
#' @param add_outgroup Should the graph outgroup be added to the qpAdm right populations?
#' Technically this shouldn't be an informative outgroup for qpAdm.
#' @param nested Should a nested data frame be returned (`TRUE`), or should populations be concatenated
#' into strings (`FALSE`)?
#' @param abbr Maximum number of characters to print for each population. The default (-1) doesn't abbreviate the names.
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

#' Find all trees which are part of the admixture graph
#'
#' @export
#' @param graph Admixture graph in `igraph` format
#' @return A data frame with columns `name` and `graph`
#' @examples
#' \dontrun{
#' trees = graph_splittrees(example_igraph)
#' # now evaluate the trees
#' trees %>%
#'   rowwise %>%
#'   mutate(res = list(qpgraph(example_f2_blocks, graph))) %>%
#'   unnest_wider(res)
#' }
graph_splittrees = function(graph) {
  # splits an admixture graph into trees
  if(!'igraph' %in% class(graph)) {
    graph %<>% graph_from_edgelist
  }
  edges = graph %>% simplify_graph %>% as_edgelist %>%
    as_tibble(.name_repair = ~c('from', 'to'))
  leaves = get_leafnames(graph)

  admixmat = edges %>% mutate(i = 1:n()) %>% group_by(to) %>% mutate(cnt = n(), j = 1:n()) %>% ungroup %>%
    filter(cnt == 2) %>% select(to, i, j) %>% spread(to, i) %>% select(-j) %>% as.matrix
  nadmix = ncol(admixmat)
  if(nadmix == 0) return(list(graph))
  #mutate(itree = map(itree, ~delete_vertices(., setdiff(get_leafnames(.), get_leafnames(tree)))))
  indices = expand.grid(replicate(ncol(admixmat), 1:2, simplify = FALSE)) %>% as.matrix
  am = map(1:nrow(indices), ~slice(edges, -admixmat[cbind(indices[.,], 1:nadmix)])) %>% map(as.matrix)
  trees = am %>% map(graph_from_edgelist)

  trees %>%
    map(~{
      tr = .
      lf = get_leafnames(tr)
      while(length(setdiff(lf, leaves)) > 0) {
        tr = igraph::delete_vertices(tr, setdiff(lf, leaves)) %>% simplify_graph
        lf = get_leafnames(tr)
      }
      tr
      }) %>%
    enframe(value = 'graph') #%>% mutate(edges = am) %>% rowwise %>% mutate(evec = list(sort(paste(edges[,1], ' ' , edges[,2])))) %>% ungroup
}

#' Find all trees within SPR distance of 1 of all graph component trees
#'
#' Returns all trees which can be reached through one iteration of subtree-prune-and-regraft on any graph component tree
#' @export
#' @param graph An admixture graph
#' @return A data frame with all trees
decomposed_tree_neighbors = function(graph) {

  graph %>% simplify_graph %>% graph_splittrees %$% graph %>% map(tree_neighbors) %>%
    bind_rows %>% rename(graph = itree) %>%
    mutate(isoclass = isomorphism_classes2(graph)) %>%
    filter(!duplicated(isoclass)) %>% select(-isoclass)
}


#' Find all trees within SPR distance of 1
#'
#' Returns all trees which can be reached through one iteration of subtree-prune-and-regraft
#' @export
#' @param tree A tree in `igraph` format
#' @return A data frame with all trees
tree_neighbors = function(tree) {
  # returns nested data frame with all trees within edit distance 1
  root = get_root(tree)
  gen1 = neighbors(tree, root, mode = 'out')
  nodes = difference(V(tree), c(root, gen1))
  map(nodes, ~tree_neighbors_single(tree, .)) %>% bind_rows
}


tree_neighbors_single = function(tree, node) {
  # returns nested data frame with all trees within edit distance 1,
  # that are created by cutting all edges leading to node
  if(class(node) != 'character') node = names(node)
  downstream = subcomponent(tree, node, mode = 'out')
  root = get_root(tree)
  outgroup = neighbors(tree, root, mode = 'out') %>% igraph::intersection(get_leaves(tree))
  parent = neighbors(tree, node, mode = 'in')
  sibling = neighbors(tree, parent, mode = 'out')
  from = names(difference(V(tree), c(downstream, root, outgroup, parent, sibling)))
  tibble(from) %>% mutate(to = node, itree = map2(from, to, ~replace_treeedge(tree, .x, .y)))
}

replace_treeedge = function(tree, from, to) {
  # removes edge leading to node 'to', and adds edge beginning above node 'from'
  parent = names(neighbors(tree, to, mode = 'in'))
  grandparent = names(neighbors(tree, parent, mode = 'in'))
  sibling = setdiff(names(neighbors(tree, parent, mode = 'out')), to)
  newgrandparent = names(neighbors(tree, from, mode = 'in'))
  tree %>%
    add_edges(c(grandparent, sibling)) %>%
    igraph::delete_vertices(parent) %>%
    add_vertices(1, name = parent) %>%
    delete_edges(paste(newgrandparent, from, sep = '|')) %>%
    add_edges(c(parent, to, newgrandparent, parent, parent, from))
}

graph_plusone_old = function(graph) {
  # returns all graphs with one more edge
  desimplify = !is_simplified(graph)
  if(desimplify) graph %<>% simplify_graph

  newedges = find_newedges_cautiously(graph)
  fn = ~insert_admix_old(graph, .x, .y)
  if(desimplify) fn %<>% compose(desimplify_graph, .dir = 'forward')
  newedges %>% mutate(graph = map2(from, to, fn))
}

#' Find all graphs which result from adding one admixture edge
#'
#' @export
#' @param graph Admixture graph in `igraph` format
#' @param ntry Specify this to return only a subset of all possible graphs with one more edge
#' @return A data frame with columns `from`, `to`, and `graph`
#' @examples
#' \dontrun{
#' newgraphs = graph_plusone(example_igraph)
#' # now evaluate the new graphs
#' newgraphs %>%
#'   rowwise %>%
#'   mutate(res = list(qpgraph(example_f2_blocks, graph))) %>%
#'   unnest_wider(res)
#' }
graph_plusone = function(graph, ntry = Inf) {
  # returns all graphs with one more edge

  desimplify = !is_simplified(graph)
  if(desimplify) graph %<>% simplify_graph
  out = graph %>% find_newedges %>% slice_sample(n = ntry) %>% rowwise %>%
    mutate(g = list(insert_admix(graph, source_from, source_to, dest_from, dest_to)))
  if(desimplify) out %<>% mutate(g = list(desimplify_graph(g)))
  out %>% ungroup %>% rename(graph = g)
}

#' Find all graphs which result from removing one admixture edge
#'
#' @export
#' @param graph Admixture graph in `igraph` format
#' @return A data frame with columns `from`, `to`, and `graph`
#' @examples
#' \dontrun{
#' newgraphs = graph_minusone(example_igraph)
#' # now evaluate the new graphs
#' newgraphs %>%
#'   rowwise %>%
#'   mutate(res = list(qpgraph(example_f2_blocks, graph))) %>%
#'   unnest_wider(res)
#' }
graph_minusone = function(graph, ntry = Inf) {
  # returns all graphs with one admixture edge removed
  desimplify = !is_simplified(graph)
  if(desimplify) graph %<>% simplify_graph
  # fn = ~delete_edges(graph, paste(.x, .y, sep = '|')) %>%
  #   igraph::delete_vertices(setdiff(get_leafnames(.), get_leafnames(graph))) %>%
  #   simplify_graph
  fn = ~delete_admix(graph, .x, .y)
  if(desimplify) fn %<>% compose(desimplify_graph, .dir = 'forward')
  graph %>% find_admixedges %>% slice_sample(n = ntry) %>% mutate(graph = map2(from, to, fn))
}

#' Find all graphs which result from adding and removing one admixture edge
#'
#' @export
#' @param graph Admixture graph in `igraph` format
#' @return A data frame with columns `source_from`, `source_to`, `dest_from`, `dest_to`, and `graph`
#' @examples
#' \dontrun{
#' newgraphs = graph_minusplus(example_igraph)
#' # now evaluate the new graphs
#' newgraphs %>%
#'   rowwise %>%
#'   mutate(res = list(qpgraph(example_f2_blocks, graph))) %>%
#'   unnest_wider(res)
#' }
graph_minusplus = function(graph) {
  graph %>% graph_minusone %>% mutate(graph2 = map(graph, graph_plusone)) %>%
    rename(source_from = from, source_to = to) %>% select(-graph) %>% unnest(graph2) %>%
    rename(dest_from = from, dest_to = to) %>% mutate(isoclass = isomorphism_classes2(graph)) %>%
    filter(!duplicated(isoclass)) %>% select(-isoclass)
}

graph_plusn = function(graphlist, n = 1, ntry = Inf) {
  if(n == 0) return(graphlist)
  map(graphlist, ~graph_plusone(., ntry = ntry)$graph) %>%
    flatten %>% graph_plusn(n-1, ntry = 1)
}

graph_minusn = function(graphlist, n = 1, ntry = Inf) {
  if(n == 0) return(graphlist)
  map(graphlist, ~graph_minusone(., ntry = ntry)$graph) %>%
    flatten %>% graph_minusn(n-1, ntry = 1)
}

#' Find all valid graphs which result from flipping one admixture edge
#'
#' @export
#' @param graph Admixture graph in `igraph` format
#' @return A data frame with columns `from`, `to`, and `graph`
#' @examples
#' \dontrun{
#' newgraphs = graph_flipadmix(example_igraph)
#' # now evaluate the new graphs
#' newgraphs %>%
#'   rowwise %>%
#'   mutate(res = list(qpgraph(example_f2_blocks, graph))) %>%
#'   unnest_wider(res)
#' }
graph_flipadmix = function(graph) {

  desimplify = !is_simplified(graph)
  if(desimplify) graph %<>% simplify_graph
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

flipadmix = function(graph, from, to) {

  if('igraph' %in% class(graph)) edges = as_edgelist(graph)
  else edges = graph
  edges %<>% as.matrix
  row = which(edges[,1] == from & edges[,2] == to)
  edges[row,] = c(to, from)
  g = igraph::graph_from_edgelist(edges)
  if(!is.dag(g) ||
     max(degree(g, mode = 'in')) > 2 ||
     max(degree(g, mode = 'out')) > 2 ||
     max(degree(g, mode = 'all')) > 3) g = NULL
  g
}

#' Add a population to an admixture graph
#'
#' @export
#' @param graph An admixture graph
#' @param pop Population to add to the graph
#' @return Admixture graph with the added population
#' @seealso \code{\link{insert_leaf}} adds pop at a specific position
graph_addleaf = function(graph, pop) {
  graph %>%
    find_normedges(exclude_first = TRUE) %>%
    mutate(graph = map2(from, to, ~insert_leaf(graph, pop, .x, .y)))
}



#' Find possible new edges
#'
#' @param graph An admixture graph
#' @return A data frame with columns `from` and `to`. New edges which begin above `from`
#'   and end above `to` could be inserted
#' @export
#' @seealso \code{\link{find_normedges}} \code{\link{find_admixedges}}
find_newedges = function(graph, fix_outpop = TRUE, all = TRUE) {
  # edge pairs are defined by 4 vertex names
  # todo: implement remove_redundant

  if(!all) return(find_newedges_cautiously(graph))
  dmat = igraph::distances(graph, mode = 'out')
  edges = graph %>% as_edgelist %>% as_tibble(.name_repair = ~c('source_from', 'source_to'))
  if(fix_outpop) {
    root = get_rootname(graph)
    outpop = get_outpop(graph)
    if(!is.null(outpop)) edges %<>% filter(source_from != root | source_to != outpop)
  }
  edges %>%
    expand_grid(rename(edges, dest_from=source_from, dest_to=source_to)) %>%
    filter(source_from != dest_from, source_to != dest_to, source_to != dest_from, dest_to != source_from) %>%
    rowwise %>% filter(!is.finite(dmat[dest_to, source_from])) %>% ungroup
}

find_newedges_cautiously = function(graph) {
  # returns a two column data frame with edges that can be inserted into the graph

  edges = matrix(NA, 0, 2)
  root = get_root(graph)
  firstgen = neighbors(graph, root, mode='out')
  admix = V(graph)[(degree(graph, mode='in') > 1)]
  admixchildren = do.call(c, igraph::ego(graph, 1, admix, mode='out'))
  admixparents = do.call(c, igraph::ego(graph, 1, admix, mode='in'))
  admixsiblings = do.call(c, igraph::ego(graph, 1, admixparents, mode='out'))
  exclude = c(admix, firstgen, root)
  if(!is.null(admixchildren)) exclude = c(exclude, admixchildren)
  cand_from = igraph::difference(V(graph), exclude)
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


#' Find admixture edges
#'
#' @export
#' @param graph An admixture graph
#' @return A data frame with columns `from` and `to` with admixture edges
#' @seealso \code{\link{find_normedges}} \code{\link{find_newedges}}
find_admixedges = function(graph) {
  # returns a two column data frame with edges that can be removed from the graph

  edges = matrix(NA, 0, 2)
  admixnodes = which(degree(graph, mode = 'in') == 2)
  if(length(admixnodes) == 0) return(tibble(from = character(), to = character()))
  graph %>% adjacent_vertices(admixnodes, mode = 'in') %>%
    map(names) %>% enframe('to', 'from') %>%
    unnest_longer(from) %>% select(2:1)
}

#' Find drift edges
#'
#' @export
#' @param graph An admixture graph
#' @param exclude_first Do not return edge from root to outgroup
#' @return A data frame with columns `from` and `to` with drift edges
#' @seealso \code{\link{find_newedges}} \code{\link{find_admixedges}}
find_normedges = function(graph, exclude_first = FALSE) {
  # returns a two column data frame with non-admixture edges
  el = graph %>% as_edgelist
  if(exclude_first) el = el[-1,]
  el %>% as_tibble(.name_repair = ~c('from', 'to')) %>%
    group_by(to) %>% filter(n() == 1) %>% ungroup
}

#' Delete an admixture edge
#'
#' @export
#' @param graph An admixture graph
#' @param from Edge source node
#' @param to Edge target node
#' @return Admixture graph with one deleted edge
#' @seealso \code{\link{insert_admix}}
delete_admix = function(graph, from = NULL, to = NULL) {
  # returns graph with admixture edge deleted
  # does not conserve internal node names

  leaves = sort(get_leafnames(graph))
  desimplify = !is_simplified(graph)
  if(is.null(from)) {
    if(desimplify) graph %<>% simplify_graph
    if(!is_valid(graph)) browser()
    admix = graph %>% find_admixedges %>% slice_sample
    if(nrow(admix) == 0) stop("No admix edges found!")
    from = admix$from
    to = admix$to
  }
  parents = neighbors(graph, from, mode = 'in')
  if(length(parents) == 1 && degree(graph, from) == 2) {
    del = from
    if(length(neighbors(graph, parents, mode = 'in')) == 2) del = c(del, names(parents))
    browser()
    graph %<>% igraph::delete_vertices(del)
  } else {
    graph %<>% delete_edges(paste(from, to, sep = '|'))
    if(!is_valid(graph)) browser()
  }
  graph %<>% simplify_graph
  newleaves = setdiff(get_leafnames(graph), leaves)
  while(length(newleaves) > 0) {
    g = graph %>% igraph::delete_vertices(newleaves) %>% simplify_graph()
    if(!is_valid(g)) browser()
    newleaves = setdiff(get_leafnames(g), leaves)
    if(length(setdiff(leaves, get_leafnames(g))) > 0) browser()
    graph = g
  }
  if(desimplify) graph %<>% desimplify_graph
  if(length(setdiff(leaves, get_leafnames(graph))) > 0) browser()
  if(!is_valid(graph)) browser()
  if(!isTRUE(all.equal(sort(get_leafnames(graph)), leaves))) stop('xxxx')
  graph
}

#' Remove population from graph
#'
#' @export
#' @param graph An admixture graph
#' @param leaf Population to be removed
#' @return Admixture graph with removed population
#' @seealso \code{\link{insert_leaf}}
delete_leaf = function(graph, leaf) {
  # deletes leaf and all internal nodes leading to no other leaves
  desimplify = !is_simplified(graph)
  if(desimplify) graph %<>% simplify_graph
  leaves = get_leafnames(graph)
  graph %<>%
    subcomponent(leaf, mode = 'in') %>%
    keep(~length(intersect(leaves, names(igraph::subcomponent(graph, .x, mode = 'out')))) == 1) %>%
    igraph::delete_vertices(graph, .) %>%
    simplify_graph()
  if(desimplify) graph %<>% desimplify_graph
  graph
}

#' Add population to graph
#'
#' @export
#' @param graph An admixture graph
#' @param leaf Population to be added
#' @param from Source node of edge onto which `leaf` should be added
#' @param to Target node of edge onto which `leaf` should be added
#' @return Admixture graph with added population
#' @seealso \code{\link{delete_leaf}}, \code{\link{graph_addleaf}} to add `leaf` at any position
insert_leaf = function(graph, leaf, from, to) {
  # inserts new leaf at edge from -> to
  nn = paste(from, to, sep = '_')
  graph %>%
    add_vertices(2, attr = list('name' = c(nn, leaf))) %>%
    add_edges(c(from, nn, nn, to, nn, leaf)) %>%
    delete_edges(paste(from, to, sep = '|'))
}

#' Generate all trees
#'
#' This functions generates all possible trees with for a given set of leaf nodes.
#'
#' @export
#' @param leaves The leaf nodes
#' @return A list of trees in `igraph` format
generate_all_trees = function(leaves) {

  stopifnot(!'R' %in% leaves)
  if(any(str_detect(leaves, '\\|'))) stop('Leaves cannot have "|" in name!')
  init = graph_from_edgelist(matrix(c('R', leaves[1]), 1))
  add_leaves_rec(init, leaves[-1])
}

#' Generate all graphs
#'
#' This functions generates all possible admixture graphs with a set number of admixture events for a given set of leaf nodes.
#'
#' @export
#' @param leaves The leaf nodes
#' @param nadmix The number of admixture nodes
#' @param sure Confirm large number of graphs
#' @param verbose Print progress updates
#' @return A list of graphs in `igraph` format
generate_all_graphs = function(leaves, nadmix = 0, sure = FALSE, verbose = TRUE) {
  #if(numtot > 1000 && !sure) stop(paste0('If you really want to generate ', numtot, ' graphs, set sure to TRUE'))
  nleaves = length(leaves)
  if(nleaves > 5 || nadmix > 2 && !sure) stop("If you really want to generate that many graphs, set 'sure' to TRUE")
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


add_edges_rec = function(graph, nadmix) {
  if(nadmix == 0) return(list(graph))
  newgraphs = graph %>% graph_plusone %>% pull(graph)
  flatten(map(newgraphs, ~add_edges_rec(., nadmix-1)))
}


#' Return all graphs created from permuting a subclade
#'
#' generates new graphs from basegraph as follows:
#' 1. generates all possible trees using `addpops` (which are not in basegraph)
#' 2. attaches trees to connection_edge, which is defined by two nodes in basegraph
#' 3. adds edges originating above each edge in `source_node`, to each node above `addpops`
#'
#' @export
#' @param basegraph an admixture graph as igraph object. (convert from edge list using `igraph::graph_from_edgelist`)
#' @param addpops a vector of population labels which are not in `basegraph`. These populations should form a clade. All possible trees will be generated and those trees will be attached to `basegraph`.
#' @param connection_edge edge in `basegraph` where the tree made from `addpops` should be attached
#' @param source_nodes nodes in `basegraph`. edges above these nodes will be added and attached to all terminal edges leading to `addpops`
#' @examples
#' \dontrun{
  #' graphlist = graphmod_pavel(example_igraph, addpops = c('pop1', 'pop2', 'pop3'),
  #'                            connection_edge = c('N2N0', 'N1N'),
  #'                            source_nodes = c('Denisova.DG', 'N2N2'))
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
  # continue here: update this to work with new insert_admix
  graphs %>% map(~insert_edges(., e$source_nodes, e$addpops))
}


split_multifurcations = function(graph) {
  # input is graph as edge list matrix
  # output is a four column matrix with 'lower' and 'upper'

  # probably doesn't work for nodes with more than 2 or 3 incoming edges

  stopifnot(class(graph)[1] == 'matrix')
  multi = names(which(table(graph[,1]) > 2))
  if(ncol(graph) == 2) graph = cbind(graph, NA, NA)

  for(i in seq_len(length(multi))) {
    rows = which(graph[,1] == multi[i])
    torows = which(graph[,2] == multi[i])
    newnam = paste0(graph[rows, 1], '_multi', i, '_', seq_along(rows))
    graph[rows, 1] = newnam
    graph[torows, 2] = newnam[1]
    graph %<>% rbind(cbind(newnam[-length(newnam)], newnam[-1], '0', '1e-9'))
  }

  leaves = setdiff(graph[,2], graph[,1])
  leafadmix = intersect(leaves, names(which(table(graph[,2]) >= 2)))
  for(i in seq_len(length(leafadmix))) {
    rows = which(graph[,2] == leafadmix[i])
    newnam = paste0(leafadmix[i], '_pre', i)
    graph[rows, 2] = newnam
    graph %<>% rbind(cbind(newnam, leafadmix[i], NA, NA))
  }

  multiin = names(which(table(graph[,2]) > 2))
  for(i in seq_len(length(multiin))) {
    rows = which(graph[,1] == multiin[i])
    torows = which(graph[,2] == multiin[i])
    newnam = paste0(graph[rows, 1], '_multi', i, '_', seq_len(length(torows)-2))
    #graph[rows, 1] = newnam[1]
    graph[torows[-1], 2] = rep(newnam, each=2)
    graph %<>% rbind(cbind(newnam, graph[rows, 1], '0', '1e-9'))
  }

  graph %>%
    as_tibble(.name_repair = ~c('from', 'to', 'lower', 'upper')) %>%
    type_convert(col_types = cols())
}



#' Returns a signature of a graph consisting of the left and right descendent leaf nodes of each internal node (sorted and concatenated)
#'
#' Can be used to determine how often internal nodes occur in a list of other well fitting models
#' @export
#' @param graph An admixture graph
#' @return A graph signature as character vector
#' @examples
#' \dontrun{
#' sigs = example_winners %>% mutate(sig = map(igraph, node_signature)) %$%
#'          sig %>% unlist %>% table %>% c
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
#' @param graph An admixture graph
#' @param graphlist List of graphs
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






leafdistances = function(graph) {

  graph %<>% simplify_graph
  leaves = get_leafnames(graph)
  dist = igraph::distances(graph, mode = 'out') %>% as_tibble(rownames = 'from') %>% pivot_longer(-from, names_to = 'to', values_to = 'dist')
  dist2 = dist %>% filter(is.finite(dist), from != to)
  dist2 %>% left_join(dist2, by = 'from') %>% filter(to.x != to.y, to.x %in% leaves, to.y %in% leaves) %>% mutate(pp = paste(to.x, to.y)) %>% group_by(pp) %>% top_n(1, -jitter(dist.y)) %>% ungroup %>% transmute(from = to.y, to = to.x, dist = dist.y)

}

reconstruct_from_leafdist = function(leafdist) {

  leaves = union(leafdist$from, leafdist$to)
  pref = admixtools:::shortest_unique_prefixes(leaves)

  parents = leafdist %>% group_by(to) %>% filter(dist == 1) %>% mutate(node = paste0(to, '_p1'), dist = 1, d = list(union(from, to))) %>% rowwise %>% mutate(dnum = length(d), desc = paste(sort(d), collapse = ' '), from = to) %>% ungroup %>% select(node, dist, from, d, dnum, desc) %>% distinct %>% group_by(desc) %>% top_n(1, node) %>% ungroup

  # merge with desc from parent generation
  leafdist %>% filter(from %in% parents$from) %>% group_by(to) %>% filter(dist == 2) %>% mutate(node = paste0(to, '_p2'), dist = 2, desc = paste(sort(union(from, to)), collapse = ' '), dnum = length(union(from, to)), from = to) %>% ungroup %>% select(node, dist, from, dnum, desc) %>% distinct %>% group_by(desc) %>% top_n(1, node) %>% ungroup



  parents = leafdist %>% filter(dist == 1) %>% group_by(from) %>% mutate(d = list(union(from, to))) %>% rowwise %>% mutate(dnum = length(d), desc = paste(sort(admixtools:::shortest_unique_prefixes(d)), collapse = '_')) %>% ungroup %>% select(from, dist, d, dnum, desc) %>% distinct

  p2 = leafdist %>% filter(dist == 2) %>% filter(from %in% parents$from) %>% group_by(from) %>% mutate(d = list(union(from, to))) %>% rowwise %>% mutate(dnum = length(d), desc = paste(sort(admixtools:::shortest_unique_prefixes(d)), collapse = '_')) %>% ungroup %>% select(from, dist, d, dnum, desc) %>% distinct

  p2 %>% left_join(parents %>% select(from, d), by = 'from') %>% rowwise %>% mutate(d = list(union(to, d)), dnum = length(d), desc = paste(sort(admixtools:::shortest_unique_prefixes(d)), collapse = '_')) %>% ungroup


  graph = graph.empty()
  graph %<>% add_vertices(length(leaves), name = leaves)
  for(i in seq_along(leaves)) {
    l = leaves[i]
    ld = leafdist %>% filter(from == l)
    reachable = l
    cat(paste0(i,'\r'))
    for(j in seq_len(max(ld$dist))) {
      nn = ld %>% filter(dist == j)
      reachable = union(reachable, nn$to)
      oldnam = if(j == 1) l else newnam
      newnam = admixtools:::shortest_unique_prefixes(reachable) %>% sort %>% paste(collapse = '_') %>% paste0('_', .)
      if(!newnam %in% names(V(graph))) graph %<>% add_vertices(1, name = newnam)
      if(oldnam != newnam) graph %<>% add_edges(c(newnam, oldnam))
      if(nrow(nn) == 0) break
    }
  }


}



triples_to_tree = function() {
  # takes choose(npop, 3) population triples, and assembles a tree

}



f4_to_triples = function(f2_blocks, outgroup) {
  # infers possible 3 population trees (rooted; 4 pops with outpop) from fstats

  pops = dimnames(f2_blocks)[[1]]
  pp = setdiff(pops, outgroup)
  f4stats = f4(example_f2_blocks, pop1 = outgroup, pop2 = pp, pop3 = pp, pop4 = pp) %>%
    filter(pop2 != pop3, pop2 != pop4, pop3 < pop4)
  outprobs = f4stats %>% rowwise %>% mutate(pp = list(sort(c(pop2, pop3, pop4))), p1 = pp[[1]], p2 = pp[[2]], p3 = pp[[3]], out = ifelse(z < 0, pop4, pop3), oz = abs(z)) %>% group_by(p1, p2, p3, out) %>% summarize(oz = mean(oz))

  tripletrees = outprobs %>% top_n(1, oz)

}

tripletrees_to_pairlists = function() {
  # takes one tree for each triple and uses them to characterize each pop pair

  dat = expand_grid(pop1 = pp, pop2 = pp) %>% filter(pop1 < pop2) %>% rowwise %>% mutate(left = list(NULL), right = list(NULL), out = list(NULL))

  t1 = tripletrees %>% rename(l = p1, r = p2, x = p3) %>% group_by(l, r) %>%
    summarize(o = list(x[out == x]), left = list(x[out == l]), right = list(x[out == r]))
  t2 = tripletrees %>% rename(l = p1, r = p3, x = p2) %>% group_by(l, r) %>%
    summarize(o = list(x[out == x]), left = list(x[out == l]), right = list(x[out == r]))
  t3 = tripletrees %>% rename(l = p2, r = p3, x = p1) %>% group_by(l, r) %>%
    summarize(o = list(x[out == x]), left = list(x[out == l]), right = list(x[out == r]))

  t1 = tripletrees %>% rename(l = p1, r = p2, x = p3)
  t2 = tripletrees %>% rename(l = p1, r = p3, x = p2)
  t3 = tripletrees %>% rename(l = p2, r = p3, x = p1)
  tt = bind_rows(t1, t2, t3) %>% group_by(l, r) %>%
    summarize(o = list(x[out == x]), left = list(x[out == l]), right = list(x[out == r])) %>%
    rowwise %>% mutate(lo = length(o), lleft = length(left), lright = length(right))

}


tree_to_triplesig = function(tree) {
  # takes an igraph tree and returns a data frame with columns lr, o, present
  # output sorted by lr, o

  leaves = sort(get_leafnames(tree))
  internal = setdiff(names(igraph::V(tree)), leaves)
  mat = map(internal, ~leaves %in% names(igraph::subcomponent(tree, ., mode = 'out'))) %>%
    do.call(rbind, .) %>% set_colnames(leaves) %>% set_rownames(internal)
  cmb = combn(leaves,2)
  (mat[,cmb[1,]] & mat[,cmb[2,]]) %>%
    set_colnames(paste(cmb[1,], cmb[2,], sep = ' ')) %>%
    as_tibble(rownames = 'internal') %>%
    mutate(num = 1:n()) %>%
    pivot_longer(-c(internal, num), names_to = 'lr', values_to = 'reach') %>%
    filter(reach) %>%
    group_by(lr) %>%
    slice_max(num) %>%
    ungroup %>%
    expand_grid(o = leaves) %>%
    separate(lr, c('l', 'r'), sep = ' ', remove = F) %>%
    filter(o != l, o != r) %>%
    rowwise %>%
    mutate(inside = is.finite(igraph::distances(tree, internal, o, mode = 'out')[,1])) %>%
    ungroup %>%
    transmute(lr, o, inside) %>%
    arrange(lr, o)
}

graph_to_triplesig = function(graph) {

  graph %>%
    simplify_graph %>%
    graph_splittrees %>% pluck(2) %>%
    map(tree_to_triplesig) %>%
    bind_rows(.id = 'tree') %>%
    group_by(lr, o) %>%
    summarize(inside = any(inside)) %>%
    ungroup
}
# triples don't uniquely identify graphs


eval_plusnadmix = function(graph, qpgfun, n = 1, ntry = Inf, verbose = TRUE) {

  newgraphs = tibble(graph = graph_plusn(list(graph), n = n, ntry = ntry))
  num = nrow(newgraphs)
  #newgraphs %<>% slice_sample(n = ntry)
  if(verbose) alert_info(paste0('Evaluating ', nrow(newgraphs), ' graphs...\n'))
  if(nrow(newgraphs) == 0) return(tibble())
  newgraphs %>%
    mutate(res = furrr::future_map(graph, qpgfun, .progress = verbose)) %>%
    unnest_wider(res) %>% arrange(score)
  #%>% select(source_from, source_to, dest_from, dest_to, graph, score)
}

eval_minusnadmix = function(graph, qpgfun, n = 1, ntry = Inf, verbose = TRUE) {

  newgraphs = tibble(graph = graph_minusn(list(graph), n = n, ntry = ntry))
  if(verbose) alert_info(paste0('Evaluating ', nrow(newgraphs), ' graphs...\n'))
  if(nrow(newgraphs) == 0) return(tibble())
  newgraphs %>%
    mutate(res = furrr::future_map(graph, qpgfun, .progress = verbose)) %>%
    unnest_wider(res) %>% arrange(score)
}

eval_plusminusn = function(graph, qpgfun, n = 1, ntry = Inf, verbose = TRUE) {

  if(verbose) alert_info(paste0('Adding ', n, ' edge(s)...\n'))
  plus = eval_plusnadmix(graph, qpgfun, n = n, ntry = ntry, verbose = verbose)
  if(verbose) alert_info(paste0('Best score: ', round(plus$score[[1]], 3), '\n'))

  if(verbose) alert_info(paste0('Removing ', n, ' edge(s)...\n'))
  plus %>% slice_min(score, with_ties = FALSE) %>% pull(graph) %>% pluck(1) %>%
    eval_minusnadmix(qpgfun, n = n, ntry = Inf, verbose = verbose)
}


eval_plusonepop = function(graph, pop, qpgfun, ntry = Inf, verbose = TRUE) {

  root = get_rootname(graph)
  outgroup = intersect(names(neighbors(graph, root)), get_leafnames(graph))
  newgraphs = graph %>% as_edgelist() %>% as_tibble(.name_repair = ~c('from', 'to'))
  if(length(outgroup) > 0) newgraphs %<>% filter(from != root | to != outgroup)
  newgraphs %<>% mutate(graph = map2(from, to, ~insert_leaf(graph, pop, .x, .y)))

  if(verbose) alert_info(paste0('Found ',nrow(newgraphs),' graphs. Evaluating ', min(nrow(newgraphs), ntry), '...\n'))
  newgraphs %>%
    slice_sample(n = ntry) %>%
    mutate(res = furrr::future_map(graph, qpgfun, .progress = verbose)) %>%
    unnest_wider(res) %>% arrange(score) %>% select(from, to, graph, score)
}

evaluate_moreadmix = function(graph, qpgfun, maxadmix, ntry = Inf, verbose = TRUE) {
# adds admixture events to graph one by one, always picking the best one out of ntry randomly chosen possible ones

  nadm = numadmix(graph)
  if(nadm >= maxadmix) stop(paste0('Graph already has ', nadm, ' admixture events!'))

  for(i in seq_len(maxadmix - nadm)) {
    if(verbose) alert_info(paste0('Testing graphs with ',nadm+i,' admixture events...\n'))
    newgraphs = graph %>% eval_plusoneadmix(qpgfun, ntry = ntry)
    # if(nrow(newgraphs) > 0) {
      if(verbose) alert_info(paste0('Best score: ', round(min(newgraphs$score), 3),'\n'))
      graph = newgraphs %>% slice_min(score) %>% pull(graph) %>% pluck(1)
    # } else {
    #   warning("Can't add more admixture events!")
    #   break
    # }
  }
  graph
}


rearrange_negdrift = function(graph, from, to) {
  # graph is simplified, edge is regular edge

  sib = graph %>% neighbors(from, mode = 'out') %>% names %>% setdiff(to)
  children = graph %>% neighbors(to, mode = 'out') %>% names %>% sample
  if(length(children) == 0 || length(sib) == 0) {
    #warning('is leaf!')
    return(graph)
  }
  if(length(c(to, sib, from, children[1])) %% 2 != 0) browser()
  graph %<>%
    add_edges(c(to, sib, from, children[1])) %>%
    delete_edges(paste(c(from, to), c(sib, children[1]), sep = '|'))
  if(is_valid(graph)) return(graph)
}


find_and_rearrange_negdrift = function(graph, qpgfun) {
  # finds most negative drift edge and returns rearranged graph

  leaves = get_leafnames(graph)
  e = qpgfun(graph)$edges %>% filter(type == 'edge', !to %in% leaves) %>% slice_min(weight, with_ties = FALSE)
  if(nrow(e) == 0) stop('No negative edge found!')
  rearrange_negdrift(graph, e$from, e$to)
}


rearrange_negadmix1 = function(graph, from, to) {
  # rearranges negative admix edge by flipping its direction if possible
  newgraph = flipadmix(graph, from, to)
  if(is_valid(newgraph)) return(newgraph)
}

rearrange_negadmix2 = function(graph, from, to) {
  # rearranges negative admix edge by flipping its counterpart's direction if possible

  parent_pos = graph %>% neighbors(to, mode = 'in') %>% names %>% setdiff(from, .)
  newgraph = flipadmix(graph, parent_pos, to)
  if(is_valid(newgraph)) return(newgraph)
}


rearrange_negadmix3 = function(graph, from, to) {
  # rearranges negative admix edge by
  # attaching admixed node to edge above parent of pos weight; neg edge is attached to parent node if possible

  leaves = sort(get_leafnames(graph))
  parent_neg = from
  parent_pos = graph %>% neighbors(to, mode = 'in') %>% names %>% setdiff(parent_neg)
  grandparent_pos = graph %>% neighbors(parent_pos, mode = 'in') %>% names
  if(length(grandparent_pos) != 1) return(NULL)
  newgraph = graph %>%
    add_edges(c(grandparent_pos, to, to, parent_pos)) %>%
    delete_edges(paste(c(grandparent_pos, parent_pos, parent_neg), c(parent_pos, to, to), sep = '|'))
  if(!is_valid(newgraph)) browser()
  if(length(parent_neg) != 1 || length(parent_pos) != 1) stop('aaa')
  newgraph2 = newgraph %>% add_edges(c(parent_neg, parent_pos)) %>% simplify_graph
  if(is_valid(newgraph2)) newgraph = newgraph2
  else newgraph %<>% find_newedges %>% slice_sample %>%
    mutate(g = list(insert_admix(newgraph, source_from, source_to, dest_from, dest_to))) %>%
    pull(g) %>% pluck(1) %>% simplify_graph
  leaves2 = sort(get_leafnames(newgraph))
  if(!isTRUE(all.equal(leaves, leaves2)) || !is_valid(newgraph)) return(NULL)
  nold = numadmix(graph)
  nnew = numadmix(newgraph)
  if(nnew < nold) newgraph %<>% insert_admix_multi(nold - nnew)
  if(is_valid(newgraph)) return(newgraph)
}



find_graphs2 = function(graph, qpgfun, numgen = 1, numgraphs = 10,
                        mutfuns = namedList(spr_leaves, spr_all, swap_leaves, move_admixedge_once, flipadmix_random, mutate_n), verbose = TRUE) {

  #wfuns = namedList(rearrange_negadmix1, rearrange_negadmix2, rearrange_negadmix3)
  wfuns = namedList(rearrange_negadmix3, replace_admix_with_random)
  dfuns = namedList(rearrange_negdrift)
  allfuns = c(wfuns, dfuns, mutfuns)

  models = tibble(g = list(simplify_graph(graph))) %>%
    mutate(res = map(g, qpgfun), hash = map_chr(g, graph_hash), lasthash = graph_hash(graph)) %>%
    unnest_wider(res) %>%
    mutate(edges = map(edges, ~filter(., weight < 0))) %>%
    transmute(gen2 = 0, hash, g, edges, score, expanded = FALSE)
  best = besthist = models$score
  bestmut = 'none'
  gimp = 0
  st = data.tree::Node$new('R')
  st$AddChild(models$hash,
              gen2 = 0,
              score = models$score,
              imp = 0,
              scm_1 = NA,
              scm_5 = NA,
              cntm_1 = NA,
              cntm_5 = NA)
  stdat = st_to_dat(st)
  tm = Sys.time()

  for(i in seq_len(numgen)) {

    if(verbose) {
      alert = if(!is.na(besthist[max(1,i-1)]) && besthist[max(1,i-1)] == best) alert_success else alert_info
      msg = paste0(i, ': sc ', round(besthist[max(1,i-1)], 3), '\tbest ', round(best, 3),'\tnew ', if(i==1) 1 else nrow(newmod), '\ttot ', sum(!is.na(models$score)), ' ', sum(stdat$totalCount == 1), ' ', sum(stdat$score[stdat$totalCount == 1] < best*2),'\t')
    }
    newgraph = pick_graph(stdat, verbose = verbose > 2)
    if(verbose > 1) {
      msg %<>% paste0(newgraph %$% paste0(' ', round(score, 2), ' ', round(score/best, 2),'\tg ',gen2, ' ', gimp, '\tcnt ',cntm_1,' ',cntm_5, '\ttime '), round(as.numeric(Sys.time()-tm), 1), '\t', bestmut)
      tm = Sys.time()
    }
    if(verbose) alert(paste0(msg, '\n'))
    sel = models %>% filter(hash == newgraph$hash[[1]])

    if(any(duplicated(models$hash))) browser()
    models %<>% rows_update(sel %>% mutate(expanded = TRUE), by = 'hash')
    graph = sel$g[[1]]

    if(gimp > 0 && gimp %% 5 == 0) {

      graph = models %>% slice_min(score, with_ties = FALSE) %>% pull(g) %>% pluck(1)
      newmod = eval_plusminusn(graph, qpgfun, n = sample(1:2, 1, prob = c(100,1)), ntry = numgraphs*10) %>%
        mutate(hash = map_chr(graph, graph_hash)) %>%
        slice_min(score, with_ties = FALSE) %>%
        transmute(hash, lasthash = sel$hash[[1]], g = graph, gen2 = i, edges, score, mutfun = 'plusminusn')
      if(newmod$hash == newmod$lasthash) {
        gimp = gimp + 1
        next
      }

    } else {

      e = qpgfun(graph, constrained = FALSE)$edges %>% filter(weight < 0)
      if(nrow(e) == 0) {
        newmod = tibble()
      } else {
        newmod = e %>% slice_min(weight, n = floor(numgraphs/2)) %>%
          mutate(mutfun = ifelse(type == 'admix', sample(names(wfuns), 1), sample(names(dfuns), 1)),
                 fun = allfuns[mutfun],
                 g = pmap(list(fun, from, to), function(x, y, z) x(graph, y, z))) %>%
          rowwise %>% filter(!is.null(g)) %>% ungroup
      }
      remaining = numgraphs - nrow(newmod)

      if(remaining > 0) {
        swape = swap_negdrift(graph, e$from, e$to)
        if(!is.null(swape)) newmod %<>% bind_rows(tibble(g = swape, mutfun = 'swap_negdrift') %>% slice_sample(n = ceiling(remaining)/2))
        remaining = numgraphs - nrow(newmod)
      }

      if(remaining > 0) {
        randmut = tibble(fun = mutfuns, mutfun = names(fun)) %>%
          slice_sample(n = remaining, replace = TRUE) %>%
          rowwise %>% mutate(g = list(fun(graph))) %>% ungroup %>% select(-fun)
        newmod %<>% bind_rows(randmut)
      }

      newmod %<>%
        transmute(g, hash = map_chr(g, graph_hash), lasthash = sel$hash[[1]], mutfun, expanded = FALSE) %>%
        filter(!duplicated(hash), !hash %in% models$hash)

      if(nrow(newmod) == 0) {
        alert_danger('No new models!\n')
        next
      }
      newmod %<>%
        # slice_head(n = min(max(numgraphs-3,0), nrow(.))) %>%
        # bind_rows(slice_sample(newmod)) %>%
        # bind_rows(slice(newmod, pmax(1,floor(nrow(newmod)/c(2,10))))) %>%
        # filter(!duplicated(hash)) %>%
        mutate(res = furrr:::future_map(g, qpgfun)) %>% unnest_wider(res) %>%
        mutate(edges = map(edges, ~filter(., weight < 0))) %>%
        transmute(gen2 = i, hash, lasthash, g, edges, score, mutfun, expanded = FALSE) %>%
        arrange(score)
    }

    if(any(is.na(newmod$score)) || nrow(newmod) == 0) browser()
    if(is.na(newmod$score[1])) browser()
    besthist[i] = newmod$score[1]
    bestmut = newmod$mutfun[1]
    if(besthist[i] < best) {
      best = besthist[i]
      gimp = 0
    } else gimp = gimp + 1
    newmod %<>% filter(!hash %in% models$hash)
    models %<>% bind_rows(newmod, .)

    for(j in seq_len(nrow(newmod))) {
      node = data.tree::FindNode(st, newmod$lasthash[j])
      node$AddChild(newmod$hash[j],
                    gen2 = i,
                    score = newmod$score[j],
                    imp = as.numeric(best != newmod$score[1]),
                    scm_1 = node$parent$score,
                    scm_5 = node$parent$parent$parent$parent$parent$score,
                    cntm_1 = node$parent$totalCount,
                    cntm_5 = node$parent$parent$parent$parent$parent$totalCount)
    }
    stdat = st %>% st_to_dat

  }

  #if(verbose) alert_info(paste0(i, ': ', round(models$score[[1]], 3), '\n'))
  namedList(models, stdat)
}

st_to_dat = function(st) {

  fp_inf = function(node) node$scbest = if(node$isLeaf) node$score else min(node$score, sapply(node$children, fp_inf), na.rm=T)
  fp_1 = function(node, n = 1) node$scp_1 = if(node$isLeaf || n == 0) node$score else min(map_dbl(node$children, function(node) if(is.null(node$score)) NA else node$score), na.rm=T)
  fp_5 = function(node, n = 5) node$scp_5 = if(node$isLeaf || n == 0) node$score else min(sapply(node$children, function(node) fp_5(node, n-1)), na.rm=T)
  fimp = function(node) node$gimp = if(is.null(node$imp) || is.na(node$imp) || node$imp == 0) 0 else node$imp + fimp(node$parent)

  st$Do(fp_inf)
  st$Do(fp_1)
  st$Do(fp_5)
  st$Do(fimp)

  as.data.frame(st, NULL, T, 'name', 'score', 'scbest', 'gimp', 'scp_1', 'scp_5', 'scm_1', 'scm_5', 'height', 'level', 'totalCount', 'cntm_1', 'cntm_5', 'gen2') %>%
    as_tibble %>% rename(hash = name) %>% select(-1)
}

graph_potential = function(stdat) {

  stdat %>% pivot_longer(contains('_'), names_to = c('type', 'gen'), names_pattern = '([A-Za-z]+)_(\\d+)', values_to = 'v') %>% pivot_wider(names_from = type, values_from = v) %>% mutate(relscore = score/min(score, na.rm=T), relm = scm/score, relp = score/scp, relmcnt = relm/cntm, relpcnt = relp/totalCount) %>% pivot_wider(names_from = gen, values_from = c(scp, scm, cntm, relm, relp, relmcnt, relpcnt)) %>% mutate(potential = 1/(relscore * relmcnt_1 * relmcnt_5))

}


pick_graph = function(stdat, minprob = 0.5, cnt1 = 10, cnt5 = 30, pen1 = 1.2, pen5 = 2, gpen = 0.2, verbose = FALSE) {
  # selects new graph from search tree data and returns one row of stdat

  sel = stdat %>% slice_min(score, with_ties = FALSE)
  lastimp = max(stdat$gen2, na.rm = TRUE) - sel$gen2
  #if(sel$totalCount[[1]] == 1 || runif(1) < minprob) {
  if(sel$totalCount[[1]] == 1 || lastimp < 3 || runif(1) > minprob) {
    return(sel)
  }

  out = stdat %>% filter(totalCount == 1) %>%
    mutate(g1 = replace_na(ifelse(cntm_1 > cnt1, cntm_1/cnt1, 1), 1),
           g5 = replace_na(ifelse(cntm_5 > cnt5, cntm_5/cnt5, 1), 1),
           imp = replace_na((gimp+1)^gpen, 1),
           potential = 1000 / (score * g1 * g5 * imp), potential = replace_na(potential, 1/max(score, na.rm=T)))
  maxpot = max(out$potential)
  if(maxpot == Inf) browser()
  out %<>% slice_sample(n = 1, weight_by = potential^10)
  if(verbose) alert_warning(out %$% paste0('\tcnt ', sum(stdat$totalCount == 1), '\tpen ', paste0(round(g1, 2), ' ', round(g5, 2), ' ', round(imp, 2), '\tpot ', round(potential, 2), ' ', round(maxpot, 2), '\n')))
  out
}



swap_negdrift = function(graph, from, to) {
  # from, to are vectors of edges to be swapped

  out = tibble(from, to) %>%
    expand_grid(rename(., from2 = from, to2 = to)) %>%
    filter(from != from2, to != to2, from != to2, to != from2) %>%
    rowwise %>%
    mutate(g = graph %>%
             delete_edges(paste(c(from, from2), c(to, to2), sep='|')) %>%
             add_edges(c(from, to2, from2, to)) %>% list)
  if(nrow(out) > 0) out %>% filter(is_valid(g)) %>% pull(g)
}


replace_admix_with_random = function(graph, from = NULL, to = NULL) {

  if(!is_valid(graph) || !all(c(from, to) %in% names(V(graph)))) browser()
  leaves = sort(get_leafnames(graph))
  g = graph %>% delete_admix(from, to)
  if(!is_valid(g)) browser()
  tryCatch(all.equal(leaves, sort(get_leafnames(g))), error = function() browser())
  g %<>% insert_admix
  if(!is_valid(g)) browser()
  g
}


graph_hash = function(graph) {
  # graph identity is determined by topology and leaf names, but not by internal names

  #graph %>% unify_vertex_names() %>% as_edgelist %>% digest::digest()
  graph %<>% simplify_graph
  leaves = get_leafnames(graph) %>% sort
  internal = setdiff(names(V(graph)), leaves)
  graph %>%
    igraph::distances(internal, leaves, mode = 'out') %>%
    apply(1, function(x) paste(str_replace(as.character(x), 'Inf', '.'), collapse = '')) %>%
    unname %>% sort %>% c(leaves) %>% digest::digest()
}


edge_targets = function(graph) {

  leaves = graph %>% get_leafnames()
  nodes = names(V(graph))
  internal = setdiff(nodes, leaves)
  dst = graph %>% distances(nodes, leaves, mode = 'out')
  reachlist = map(nodes, ~names(which(is.finite(dst[.,])))) %>% set_names(nodes)
  graph %>%
    as_edgelist %>%
    as_tibble(.name_repair = ~c('from', 'to')) %>%
    add_count(to) %>%
    rowwise %>%
    mutate(targets = reachlist[to], tgt = paste0((leaves %in% targets)+0, collapse='')) %>%
    group_by(tgt) %>%
    mutate(identifiable = n > 1 | n() == 1) %>%
    ungroup

}


vs_to_es = function(vs) {
  # turns igraph vertex sequence into a vector of edge names
  if('igraph.vs' %in% class(vs)) nam = names(vs)
  else nam = vs
  paste0(head(nam, -1), '|', nam[-1])
}


path_intersections = function(graph) {

  pp = graph %>% pair_paths %>% rowwise %>% mutate(id = paste(sort(c(from, to)), collapse = ' ')) %>% ungroup
  expand_grid(pp, pp %>% rename_with(~paste0(., 2))) %>%
    filter(id != id2) %>%
    rowwise %>%
    mutate(is = list(intersect(path, path2)), islen = length(is), isstr = paste0(is, collapse = ' '),
           admtot = length(union(admnodes, admnodes2)), id4 = paste(sort(c(id, id2)), collapse = ' ')) %>%
    filter(islen > 0) %>%
    ungroup

}

resolve_paths = function(graph) {

  leaves = graph %>% get_leafnames
  adm = names(which(degree(graph, mode = 'in') > 1))
  unresolved_adm = adm
  pis = path_intersections(graph)

  for(i in 1:10) {
    pis %<>% rowwise %>%
      mutate(unresolved = list(intersect(unresolved_adm, union(admnodes, admnodes2))), unum = length(unresolved)) %>%
      ungroup

    resolved0 = pis %>% filter(unum == 0, islen > 0) %>% pull(isstr) %>% unique
    resolved1 = pis %>% filter(unum == 1, islen > 0) %>% pull(isstr) %>% unique

    resolvedadm = pis %>% filter(isstr %in% intersect(resolved0, resolved1), unum == 1) %$%
      unique(unlist(c(admnodes, admnodes2)))

    unresolved_adm %<>% setdiff(resolvedadm)
  }

}

paths_through = function(graph, leaves, node) {

  nl = graph %>% neighbors(node, mode = 'out') %>% pluck(1) %>% names
  nr = graph %>% neighbors(node, mode = 'out') %>% pluck(2) %>% names

  paths = all_simple_paths(graph, node, leaves, mode = 'out') %>%
    map(names)

  expand_grid(up = paths %>% keep(~.[[2]] == nl), down = paths %>% keep(~.[[2]] == nr))
}

pair_paths = function(graph) {

  leaves = graph %>% get_leafnames
  adm = names(which(degree(graph, mode = 'in') > 1))
  internal = setdiff(names(V(graph)), leaves)
  nonadm = setdiff(internal, adm)
  map(nonadm, ~paths_through(graph, leaves, .)) %>% set_names(nonadm) %>% bind_rows(.id = 'pivot') %>% rowwise %>% mutate(from = tail(up, 1), to = tail(down, 1), path = list(c(vs_to_es(up), vs_to_es(down))), admnodes = list(intersect(adm, c(up, down))), admnum = length(admnodes)) %>% ungroup

}

# continue here: use path intersections to determine edges which are not identifiable

