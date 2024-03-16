
# number of all binary trees with a given number of leaf nodes (== doublefactorial(2*n-3))
numtrees = function(n) factorial(2*n-2)/(2^(n-1)*factorial(n-1))

# Catalan number
catalan = function(n) factorial(2*n)/factorial(n+1)/factorial(n)

# number of possible DAGs
numdags = function(n) {
  if(n <= 1) return(1)
  sum(sapply(1:n, function(k) (-1)^(k+1) * choose(n, k) * 2^(k*(n-k)) * numdags(n-k)))
}

# number of trees with x admixture events
numtreesadmix = function(n, nadmix) numtrees(n) * numadmixplacements(2*n-2, nadmix)

# number of unique ways to add 'nadmix' undirected edges; each added edge increases the number of edges by 3
numadmixplacements = function(numedges, nadmix) {
  if(nadmix == 0) return(1)
  choose(numedges, 2) * numadmixplacements(numedges+3, nadmix-1)
}


doublefactorial = function(n) {
  if(n %% 2 == 0) {
    k = n/2
    out = 2^k*factorial(k)
  } else {
    k = (n+1)/2
    out = factorial(2*k-1)/(2^(k-1)*factorial(k-1))
  }
  out
}

#' Count number of admixture nodes
#'
#' @export
#' @param graph An admixture graph
#' @return Number of admixture nodes
numadmix = function(graph) {
  sum(degree(graph, mode='in') == 2)
}

#' Test if an admixture graph is valid
#'
#' @export
#' @param graph An admixture graph
#' @return `TRUE` if graph is valid, otherwise `FALSE`
is_valid = function(graph) {

  indegree = igraph::degree(graph, mode = 'in')
  outdegree = igraph::degree(graph, mode = 'out')

  igraph::is_simple(graph) &&
    igraph::is_connected(graph) &&
    igraph::is_dag(graph) &&
    sum(indegree == 0) == 1 &&
    all(indegree <= 2) &&
    all(outdegree <= 2) &&
    !any(outdegree == 0 & indegree != 1) &&
    !any(indegree == 0 & outdegree < 2)
}

is_simplified = function(graph) {
  indeg = degree(graph, mode = 'in')
  outdeg = degree(graph, mode = 'out')
  sum(indeg == 1 & outdeg == 1) == 0 && max(indeg + outdeg) == 3 && igraph::is_simple(graph)
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
  #if(length(root) != 1) stop(paste0('Root problem ', paste0(names(root), collapse = ' ')))
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
    c1 = igraph::distances(graph, children[1], leaves) %>% paste(collapse='')
    c2 = igraph::distances(graph, children[2], leaves) %>% paste(collapse='')
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
  # this function changes the vertex names of inner vertices,
  # so that all isomorphic graphs with the same leaf nodes have equally labelled inner nodes
  leaves = get_leaves(graph)
  root = get_root(graph)
  lv = shortest_unique_prefixes(names(leaves))
  g = set_vertex_attr(graph, 'name', leaves, lv)
  ormap = names(V(g))
  names(ormap) = ormap

  nammap = unify_vertex_names_rec(g, root, ormap, sep1 = sep1, sep2 = sep2)
  nammap[lv] = names(leaves)
  names(nammap)[match(lv, names(nammap))] = names(leaves)
  nammap[names(root)] = 'R'
  if(!keep_unique) {
    changed = setdiff(names(nammap), c(names(root), names(leaves)))
    nammap[changed] = paste0('n', as.numeric(as.factor(nammap[changed])))
  }
  nammap %<>% map_chr(digest::digest) %>% paste0('x', .)
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
random_admixturegraph = function(leaves, numadmix = 0, simple = TRUE, outpop = NULL,
                                 nonzero_f4 = NULL, admix_constraints = NULL, event_order = NULL, ntry = 100) {
  # makes a random admixture graph
  # returns an 'igraph' graph object
  # 'leaves' can be a number of leaf nodes, or a character vector of leaf names

  stopifnot(class(leaves)[1] %in% c('numeric', 'integer', 'character'))
  if(length(leaves) == 1) {
    if(leaves > length(LETTERS)) leaves = paste0('l', seq_len(leaves))
    else leaves = LETTERS[seq_len(leaves)]
  }
  for(i in seq_len(ntry)) {
    graph = leaves %>%
      setdiff(outpop) %>%
      sample %>%
      random_newick(outpop = outpop) %>%
      newick_to_edges %>%
      graph_from_edgelist %>%
      insert_admix_n(numadmix, fix_outgroup = !is.null(outpop))
    if(all(map_lgl(list(nonzero_f4, admix_constraints, event_order), is.null)) ||
       satisfies_constraints(graph, nonzero_f4 = nonzero_f4, admix_constraints = admix_constraints,
                             event_order = event_order)) break
    if(i == ntry) warning('Constraints not satisfied!')
  }

  if(!simple) graph %<>% desimplify_graph
  graph
}


simplify_graph_old = function(graph) {
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
    if(i > 100) stop('sg old error')
  })
  graph %<>% X_to_H
  if(convmat) graph = igraph::as_edgelist(graph)
  graph
}

simplify_graph_old = function(graph) {
  # removes redundant nodes

  if(is_simplified(graph)) return(graph)
  convmat = FALSE
  if(class(graph) == 'matrix') {
    convmat = TRUE
    graph = graph_from_edgelist(graph)
  }
  graph %<>% igraph::simplify()
  repeat({
    redundant = which(degree(graph, mode='in') == 1 & degree(graph, mode='out') == 1)
    if(length(redundant) == 0) break
    newfrom = names(neighbors(graph, redundant[1], mode='in'))
    newto = names(neighbors(graph, redundant[1], mode='out'))
    graph %<>% igraph::delete_vertices(redundant[1]) %>% igraph::add_edges(c(newfrom, newto))
    #%>% igraph::simplify()
  })
  graph %<>% igraph::simplify() %>% X_to_H
  if(convmat) graph = igraph::as_edgelist(graph)
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

  indegree = degree(graph, mode='in')
  outdegree = degree(graph, mode='out')
  convmat = FALSE
  if(class(graph) == 'matrix') {
    convmat = TRUE
    graph = graph_from_edgelist(graph)
  }
  graph %<>% igraph::simplify()
  redundant = which(indegree == 1 & outdegree == 1)

  if(length(redundant) != 0) {
    newfrom = redundant %>% igraph::adjacent_vertices(graph, ., mode='in') %>% unlist
    #newto = redundant %>% igraph::adjacent_vertices(graph, ., mode='out') %>% unlist
    newto = redundant %>% map_dbl(~nextnonredundant(graph, .))
    graph %<>% igraph::add_edges(interleave(newfrom, newto)) %>% igraph::delete_vertices(names(redundant))
  }
  graph %<>% igraph::simplify()
  if(any(indegree > 1 & outdegree > 1)) graph %<>% X_to_H
  if(!is_simplified(graph)) graph %<>% simplify_graph()
  if(convmat) graph = igraph::as_edgelist(graph)
  graph
}


nextnonredundant = function(graph, node) {
  n = neighbors(graph, node, mode = 'out')[1]
  if(degree(graph, n, mode = 'in') != 1 || degree(graph, n, mode = 'out') != 1) return(n)
  return(nextnonredundant(graph, n))
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
  ograph = graph
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
  if(max(degree(graph, mode='in')) > 2 || max(degree(graph, mode='out')) > 2) {
    stop('Desimplify failed')
  }
  if(convmat) graph = igraph::as_edgelist(graph)
  graph
}


X_to_H = function(graph) {

  crosses = names(which(degree(graph, mode = 'in') > 1 & degree(graph, mode = 'out') > 1))
  for(i in seq_along(crosses)) {
    nam = newnodenam(crosses[i], names(V(graph)))
    parents = names(neighbors(graph, crosses[i], mode = 'in'))
    graph %<>% igraph::add_vertices(1, name = nam) %>%
      add_edges(c(interleave(parents, rep(nam, length(parents))), nam, crosses[i])) %>%
      delete_edges(paste0(parents, '|', crosses[i]))
  }
  graph
}

H_to_X = function(graph) {
  # keeps child names

  adm = names(which(degree(graph, mode = 'in') > 1 & degree(graph, mode = 'out') == 1))
  children = names(neighbors(graph, adm, mode = 'out'))
  ind = which(degree(graph, children, mode = 'out') == 2)
  for(i in seq_along(ind)) {
    parents = names(neighbors(graph, adm[i], mode = 'in'))
    graph %<>% delete_vertices(adm[i]) %>%
      add_edges(interleave(parents, rep(children[i], length(parents))))
  }
  graph
}

merge_nested_admix = function(graph) {

  graph %<>% X_to_H
  leaves = get_leafnames(graph)

  while(TRUE) {
    adm = names(which(degree(graph, mode = 'in') > 1)) %>% setdiff(leaves)
    admchildren = map_chr(adm, ~names(neighbors(graph, ., mode = 'out')))
    both = which(admchildren %in% names(which(degree(graph, admchildren, mode = 'in') > 1)))
    if(length(both) == 0) break
    parent = adm[both[1]]
    child = admchildren[both[1]]
    grandparents = names(neighbors(graph, parent, mode = 'in'))
    graph %<>%
      delete_vertices(parent) %>%
      add_edges(interleave(grandparents, rep(child, length(grandparents))))
  }
  graph #%<>% H_to_X
}


split_graph = function(graph) {
  # removes admixture nodes from an admixturegraph
  # returns a tree and a list of edges that can be used to reconstruct the original admixturegraph
  if(!'igraph' %in% class(graph)) {
    edges = graph
    graph = edges %>% edges_to_igraph
    edges %<>% filter(type == 'admix') %>%
      left_join(edges %>% transmute(fromfrom = from, from = to), by = 'from')
    pickedges = TRUE
  } else pickedges = FALSE
  graph %<>% simplify_graph()
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
      if(length(parents) == 0 || !p1 %in% names(V(graph)) || root %in% parents) stop('split graph error')
      fromnode0 = names(neighbors(graph, p1, mode = 'in'))
      if(i > 100) stop('split graph error 2')
    }
    # stopifnot(!any(c(parents[1], admix[1]) %in% unlist(admixedges)))
    # too strict; probably won't be able to always reintroduce admixedges at appropriate locations after SPR

    fromnode = setdiff(names(neighbors(graph, p1, mode = 'out')), adm)
    tonode0 = setdiff(names(neighbors(graph, adm, mode = 'in')), p1)
    tonode = names(neighbors(graph, adm, mode = 'out'))
    fromnodes %<>% c(fromnode)
    tonodes %<>% c(tonode)

    graphnew = igraph::delete_vertices(graph, adm)
    graphnew = igraph::add_edges(graphnew, c(tonode0, tonode))
    if(tonode %in% admix) frommap[tonode0] = adm
    #if(delete_nodes) {
      graphnew = igraph::delete_vertices(graphnew, p1)
      if(p1 %in% tonodes) tonodes[tonodes == p1] = tonode
      if(length(fromnode) == 0) browser()
      if(p1 %in% fromnodes) fromnodes[fromnodes == p1] = fromnode
      graphnew2 = igraph::add_edges(graphnew, c(fromnode0, fromnode)) %>% simplify_graph()
      if(!is_simplified(graphnew2)) browser()
      graphnew = graphnew2
    #}
    if(is_valid(graphnew) && numadmix(graphnew) < numadmix(graph)) graph = graphnew
    else stop('Graph invalid')
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

restore_admixed_leaves = function(graph) {

  # keep terminal branch only for leaves which are non-admixed

  leaves = get_leafnames(graph)
  adm = names(which(igraph::degree(graph, mode = 'in') > 1)) %>%
    intersect(leaves)
  if(length(adm) == 0) return(graph)
  nam = paste0('adm_before_', adm)

  graph %>%
    igraph::set_vertex_attr('name', adm, nam) %>%
    igraph::add_vertices(length(adm), attr = list(name = adm)) %>%
    igraph::add_edges(interleave(nam, adm))
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
  if(!isTRUE(all.equal(leaves, sort(get_leafnames(graph))))) stop('insert_admix_old failed!')
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
#' @seealso \code{\link{delete_admix}}, \code{\link{insert_admix_n}}
insert_admix = function(graph, source_from = NULL, source_to = NULL, dest_from = NULL, dest_to = NULL,
                        substitute = FALSE, fix_outgroup = TRUE) {
  # assumes graph is simplified

  if(length(source_to) > 1) {
    n = length(source_to)
    f = function(x, y) insert_admix(x, source_from[y], source_to[y], dest_from[y], dest_to[y], substitute = substitute)
    graph = source_to %>% length %>% seq_len %>% reduce(f, .init = graph)
    return(graph)
  }

  leaves = sort(get_leafnames(graph))
  nodes = names(V(graph))
  if(is.null(source_to) || is.null(dest_to) || length(setdiff(c(source_to, dest_to), nodes)) > 0 && substitute) {
    e = graph %>% find_newedges(fix_outgroup = fix_outgroup) %>% slice_sample(n=1)
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
    else browser()
    #else stop("Inserting edge failed!")
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
insert_admix_n = function(graph, n = 1, fix_outgroup = TRUE) {

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
  f = function(x, y) insert_admix(x, substitute = TRUE, fix_outgroup = fix_outgroup)
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
    if(length(hostnodes) > 0 && length(cutsibling) != 0) break
    if(i > 100) stop('spr error')
  })

  hostnode = sample(names(hostnodes), 1)
  hostparent = names(neighbors(graph, hostnode, mode='in'))

  newnam = paste(hostnode, sample(letters, 1), sep='_')
  while(newnam %in% names(V(graph))) newnam = paste0(newnam, sample(letters, 1))

  tryCatch({
  graph %<>% igraph::add_vertices(1, name=newnam) %>%
    igraph::add_edges(c(names(cutgrandparent), names(cutsibling),
                hostparent, newnam,
                newnam, hostnode,
                newnam, names(cutnode))) %>%
    igraph::delete_vertices(names(cutparent)) %>%
    igraph::delete_edges(paste(hostparent, hostnode, sep='|'))
  }, error = function(e) stop('spr error'))
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
  if(!is_valid(graph)) stop('apr error')
  spl = split_graph(graph)
  graph = subtree_prune_and_regraft(spl$tree, only_leaves = only_leaves, fix_outgroup = fix_outgroup)
  #graph = insert_admix_old(graph, spl$fromnodes, spl$tonodes, substitute_missing = TRUE)
  graph = insert_admix(graph, source_to = spl$fromnodes, dest_to = spl$tonodes, substitute = TRUE)
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
          if(length(c(grandparent, sibling, newgrandparent, parents[j], parents[j], n)) %% 2 != 0) stop('move admix error')
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
  #if(nadmix > 0) warning('no suitable attachment points found for admixture edge')
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


swap_admix = function(graph) {

  admixedges = graph %>% find_admixedges %>% sample_frac(1)
  sw = admixedges %>% filter(!duplicated(to)) %>% slice(1:2) %>% mutate(e = paste(from, to, sep='|'))
  graph %>%
    igraph::delete_edges(sw$e) %>%
    igraph::add_edges(c(sw$from[1], sw$to[2], sw$from[2], sw$to[1]))
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
    newg = admixedges %>% slice(i) %>% mutate(graph = map2(from, to, ~flipadmix2(el, .x, .y))) %$% graph[[1]]
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
    map = function(...) furrr::future_map(..., .options = furrr::furrr_options(seed = TRUE))
    imap = function(...) furrr::future_imap(..., .options = furrr::furrr_options(seed = TRUE))
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
                                     numstart = numstart, return_fstats = opt_worst_residual, verbose = FALSE, ...)[kp]
  } else {
    # qpgfun = function(graph) qpgraph(data = NULL, graph = graph, f3precomp = precomp,
    #                                  numstart = numstart, return_fstats = opt_worst_residual, verbose = FALSE, ...)[kp]
    qpgfun = function(graph) qpgraph(data = precomp, graph = graph,
                                     numstart = numstart, return_fstats = opt_worst_residual, verbose = FALSE, ...)[kp]
  }

  # qpgfun = possibly(qpgfun, otherwise = NULL)
  space = paste0(paste(rep(' ', 50), collapse=''), '\r')
  if(verbose) alert_info(paste0('Generate new graphs...', space))
  if(is.null(initgraphs)) initgraphs = replicate(numgraphs, random_admixturegraph(pops, numadmix, outpop = outpop),
                                                simplify = FALSE)
  else initgraphs = initgraphs[round(seq(1, length(initgraphs), numgraphs))]
  if(verbose) alert_info(paste0('Evaluate graphs...', space))
  if(parallel) map = function(...) furrr::future_map(..., .options = furrr::furrr_options(seed = TRUE))
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
#' @param data Input data in one of three forms:
#' \enumerate{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{f2_from_geno}}
#' \item A directory which contains pre-computed f2-statistics
#' \item The prefix of genotype files
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
#' @param store_intermediate Path and prefix of files for intermediate results to `.rds`. Can be useful if `find_graphs_old` doesn't finish sucessfully.
#' @param parallel Parallelize over repeats (if `numrep > 1`) or graphs (if `numrep == 1`) by replacing \code{\link[purrr]{map}} with \code{\link[furrr]{future_map}}. Will only be effective if \code{\link[future]{plan}} has been set.
#' @param stop_after Stop optimization after `stop_after` seconds (and after finishing the current generation).
#' @param verbose Print progress updates
#' @param ... Additional arguments passed to \code{\link{qpgraph}}
#' @return A nested data frame with one model per line
#' @seealso \code{\link{qpgraph}}
#' @examples
#' \dontrun{
#' find_graphs_old(example_f2_blocks, numrep = 200, numgraphs = 100,
#'             numgen = 20, numsel = 5, numadmix = 3)
#' }
#' \dontrun{
#' # Making new mutation functions by modifying or combining existing ones:
#' newfun1 = function(graph, ...) mutate_n(graph, 3, ...)
#' newfun2 = function(graph, ...) flipadmix_random(spr_leaves(graph, ...), ...)
#' find_graphs_old(f2_blocks, mutfuns = namedList(spr_leaves, newfun1, newfun2), mutprobs = c(0.2, 0.3, 0.5))
#' }
find_graphs_old = function(data, pops = NULL, outpop = NULL, numrep = 1, numgraphs = 50,
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

  # if(!opt_worst_residual) precomp = qpgraph_precompute_f3(data, pops, f3basepop = outpop, ...)
  # else precomp = get_f2(data, pops, ...)
  precomp = get_f2(data, pops, ...)

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
    #mutate(tlr = !x13 & !x31 & (x23 | x32) & (x12 | x21))
    mutate(tlr = !x13 & (x23 | x32) & (x12 | x21))

  tripletopo
}


#' Summarize triples across graphs
#'
#' This function summarizes topologies of population triples across graphs.
#'
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
  # returns integer vector with the same length as 'igraphlist', which assigns each graph to a class
  # only considers topology

  numgraph = length(igraphlist)
  if(numgraph == 0) return(numeric())
  if(numgraph == 1) return(1)

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
  # runs `simplify_graph` on all graphs before comparing them

  numgraph = length(igraphlist)
  if(numgraph == 0) return(numeric())
  if(numgraph == 1) return(1)

  hashes = map_chr(igraphlist, graph_hash)
  factor(hashes, levels = unique(hashes)) %>% as.numeric
}

#' Return all valid qpAdm models for an admixturegraph
#'
#' For large admixturegraph, there may be a large number of valid qpAdm models!
#'
#' @param graph An admixture graph as igraph object
#' @param add_outgroup Should the graph outgroup be added to the qpAdm right populations?
#' Technically this shouldn't be an informative outgroup for qpAdm.
#' @param nested Should a nested data frame be returned (`TRUE`), or should populations be concatenated
#' into strings (`FALSE`)?
#' @param abbr Maximum number of characters to print for each population. The default (-1) doesn't abbreviate the names.
#' @examples
#' \dontrun{
#' qpadm_models_old(igraph2, add_outgroup = TRUE)
#' qpadm_models_old(igraph2, add_outgroup = TRUE) %>% slice(1) %$% list(target, left, right)
#' }
qpadm_models_old = function(graph, add_outgroup=FALSE, nested = TRUE, abbr = -1) {
  # don't do this for large models

  outgroup = get_outpop(graph)
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

identify_edge = function(graph, from, to) {
  # this function returns the original from node for given a desimplified graph, and an admixture edge a simplified graph
  nam = names(V(graph))
  if(!from %in% nam || !to %in% nam) return(from)
  shortest_paths(graph, from, to, mode = 'out')$vpath[[1]] %>% names %>% tail(2) %>% head(1)
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
graph_splittrees = function(graph, return_admix = FALSE, simplify = TRUE) {
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

  if(return_admix) {
  admedges = am %>% map(~as_tibble(.) %>% filter(to %in% colnames(admixmat)) %>%
                          rowwise %>% mutate(from = identify_edge(graph, from, to)) %>% ungroup)
  }
  out = trees %>% map(~prune_extra_leaves(., leaves))
  if(simplify) out %<>% map(simplify_graph)
  out %<>% enframe(value = 'graph')
  if(return_admix) out %<>% mutate(admedges) #%>% rowwise %>% mutate(evec = list(sort(paste(edges[,1], ' ' , edges[,2])))) %>% ungroup
  out
}

prune_extra_leaves = function(graph, keep) {
  leaves = get_leafnames(graph)
  while(length(setdiff(leaves, keep)) > 0) {
    graph = igraph::delete_vertices(graph, setdiff(leaves, keep))
    leaves = get_leafnames(graph)
  }
  graph
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
  out = graph %>% find_newedges
  if(is.finite(ntry)) out %<>% slice_sample(n = ntry, replace = TRUE)
  out %<>% rowwise %>%
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
  graph %>% find_admixedges %>% slice_sample(n = min(100, ntry), replace = TRUE) %>%
    mutate(graph = map2(from, to, fn))
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
  graph %>% graph_minusone %>%
    mutate(graph2 = map(graph, graph_plusone)) %>%
    select(-graph) %>% unnest(graph2) %>%
    mutate(isoclass = isomorphism_classes2(graph)) %>%
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
  # if(!is.dag(g) ||
  #    max(degree(g, mode = 'in')) > 2 ||
  #    max(degree(g, mode = 'out')) > 2 ||
  #    max(degree(g, mode = 'all')) > 3) g = NULL
  if(!is_valid(g)) g = NULL
  g
}

flipadmix2 = function(graph, from, to) {

  if('igraph' %in% class(graph)) edges = as_edgelist(graph)
  else edges = graph
  edges %<>% as.matrix
  row1 = which(edges[,1] == from & edges[,2] == to)
  row2 = which(edges[,2] == from)[1]
  edges[row1,] = c(to, from)
  g = igraph::graph_from_edgelist(edges)
  if(!is.na(row2)) edges[row2,] = rev(edges[row2,])
  g2 = igraph::graph_from_edgelist(edges)
  if(is_valid(g2)) return(g2)
  if(is_valid(g)) return(g)
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
find_newedges = function(graph, fix_outgroup = TRUE, all = TRUE) {
  # edge pairs are defined by 4 vertex names
  # todo: implement remove_redundant

  if(!all) return(find_newedges_cautiously(graph))
  dmat = igraph::distances(graph, mode = 'out')
  edges = graph %>% as_edgelist %>% as_tibble(.name_repair = ~c('source_from', 'source_to'))
  if(fix_outgroup) {
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
  ograph = graph
  if(is.null(from)) {
    if(desimplify) graph %<>% simplify_graph
    if(!is_valid(graph)) stop('delete admix error')
    admix = graph %>% find_admixedges
    if(nrow(admix) == 0) {
      stop("No admix edges found!")
    } else {
      admix %<>% slice_sample(n=1)
    }
    from = admix$from
    to = admix$to
  }
  parents = neighbors(graph, from, mode = 'in')
  if(length(parents) == 1 && degree(graph, from) == 2) {
    del = from
    if(length(neighbors(graph, parents, mode = 'in')) == 2) del = c(del, names(parents))
    graph %<>% igraph::delete_vertices(del)
  } else {
    graph %<>% delete_edges(paste(from, to, sep = '|'))
  }
  graph %<>% simplify_graph
  newleaves = setdiff(get_leafnames(graph), leaves)
  while(length(newleaves) > 0) {
    g = graph %>% igraph::delete_vertices(newleaves) %>% simplify_graph()
    newleaves = setdiff(get_leafnames(g), leaves)
    graph = g
  }
  if(degree(graph, get_root(graph), mode = 'out') == 1) {
    graph %<>% igraph::delete_vertices(get_root(graph))
  }
  if(desimplify) graph %<>% desimplify_graph
  if(!isTRUE(all.equal(sort(get_leafnames(graph)), leaves)) || !is_valid(graph)) stop('xxxx')
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
    igraph::subcomponent(leaf, mode = 'in') %>%
    keep(~length(intersect(leaves, names(igraph::subcomponent(graph, .x, mode = 'out')))) == 1) %>%
    igraph::delete_vertices(graph, .) %>%
    simplify_graph()
  if(desimplify) graph %<>% desimplify_graph
  graph
}

delete_leaves = function(graph, leaves) {
  # assumes simplified graph

  keepleaves = graph %>% get_leafnames %>% setdiff(leaves)
  keep_leaves(graph, keepleaves)
}

keep_leaves = function(graph, leaves) {
  # assumes simplified graph

  keepv = map(leaves, ~igraph::subcomponent(graph, ., mode = 'in')) %>% unlist %>% names %>% unique
  delv = setdiff(names(V(graph)), keepv)
  out = graph %>% igraph::delete_vertices(delv)
  out %>% simplify_graph
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
  add_leaves_rec(init, leaves[-1]) %>% map(~delete_vertices(.,'R'))
}

#' Generate all graphs
#'
#' This functions generates all possible admixture graphs with a set number of admixture events for a given set of leaf nodes. It's pretty slow, and may not terminate in reasonable time for more than 5 leaves and 2 admixture events. The function is similar to the \code{\link[admixturegraph]{all_graphs}} function in the \code{admixturegraph} package, but there are a few differences:
#' \itemize{
#' \item The function does not return graphs with fewer than `nadmix` admixture events
#' \item The function does not return most graphs which are unidentifiable and would have equal fits as simpler identifiable graphs (for example it does not return graphs where a node is expanded to a loop)
#' \item The function does not return duplicated graphs, as identified by the \code{\link{graph_hash}} function
#' \item The function generates unique graphs which are missing in the output of \code{\link[admixturegraph]{all_graphs}}
#' }
#' @export
#' @param leaves The leaf nodes
#' @param nadmix The number of admixture nodes
#' @param verbose Print progress updates
#' @return A list of graphs in `igraph` format
#' @seealso \code{\link[admixturegraph]{all_graphs}}, \code{\link{generate_all_trees}}, \code{\link{graph_hash}}
#' @examples
#' \dontrun{
#' graphs = generate_all_graphs(letters[1:4], 1)
#' }
generate_all_graphs = function(leaves, nadmix = 0, verbose = TRUE) {
  nleaves = length(leaves)
  if(verbose) alert_info(paste0('Generating ', numtrees(nleaves),' trees...\n'))
  trees = generate_all_trees(leaves)
  if(verbose) alert_info(paste0('Adding all possible admixutre edges...\n'))
  graphs = flatten(map(trees, ~add_edges_rec(., nadmix)))
  if(verbose) alert_info(paste0('Identifying isomorpisms in ', length(graphs),' graphs...\n'))
  hashes = map_chr(graphs, graph_hash)
  keep = which(!duplicated(hashes))
  if(verbose) alert_info(paste0('Returning ', length(keep),' unique graphs...\n'))
  graphs[keep]
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


permutations = function(vec) {

  if(length(vec) == 1) return(vec)
  out = list()
  for(i in seq_len(length(vec))) {
    new = map(permutations(vec[-i]), ~c(vec[i], .))
    out = c(out, new)
  }
  out
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

#' Split nodes with more than two edges
#'
#' @export
#' @param graph An admixture graph
#' @return A new admixture graph in which no node has more than two incoming or outgoing edges
split_multifurcations = function(graph) {
  # input is graph as edge list matrix
  # output is a four column matrix with 'lower' and 'upper'

  # probably doesn't work for nodes with more than 2 or 3 incoming edges

  ig = class(graph)[1] == 'igraph'
  if(ig) graph %<>% igraph::as_edgelist()

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

  graph %<>%
    as_tibble(.name_repair = ~c('from', 'to', 'lower', 'upper')) %>%
    type_convert(col_types = cols())
  if(ig) graph %<>% edges_to_igraph()
  graph
}



#' Returns a signature of a graph consisting of the left and right descendent leaf nodes of each internal node (sorted and concatenated)
#'
#' Can be used to determine how often internal nodes occur in a list of other well fitting models
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
  dist = igraph::distances(graph, mode = 'out') %>% as_tibble(rownames = 'from') %>%
    pivot_longer(-from, names_to = 'to', values_to = 'dist')
  dist2 = dist %>% filter(is.finite(dist), from != to)
  dist2 %>% left_join(dist2, by = 'from') %>% filter(to.x != to.y, to.x %in% leaves, to.y %in% leaves) %>% mutate(pp = paste(to.x, to.y)) %>% group_by(pp) %>% top_n(1, -jitter(dist.y)) %>% ungroup %>% transmute(from = to.y, to = to.x, dist = dist.y)

}

reconstruct_from_leafdist = function(leafdist) {

  leaves = union(leafdist$from, leafdist$to)
  pref = shortest_unique_prefixes(leaves)

  parents = leafdist %>% group_by(to) %>% filter(dist == 1) %>% mutate(node = paste0(to, '_p1'), dist = 1, d = list(union(from, to))) %>% rowwise %>% mutate(dnum = length(d), desc = paste(sort(d), collapse = ' '), from = to) %>% ungroup %>% select(node, dist, from, d, dnum, desc) %>% distinct %>% group_by(desc) %>% top_n(1, node) %>% ungroup

  # merge with desc from parent generation
  leafdist %>% filter(from %in% parents$from) %>% group_by(to) %>% filter(dist == 2) %>% mutate(node = paste0(to, '_p2'), dist = 2, desc = paste(sort(union(from, to)), collapse = ' '), dnum = length(union(from, to)), from = to) %>% ungroup %>% select(node, dist, from, dnum, desc) %>% distinct %>% group_by(desc) %>% top_n(1, node) %>% ungroup



  parents = leafdist %>% filter(dist == 1) %>% group_by(from) %>% mutate(d = list(union(from, to))) %>% rowwise %>% mutate(dnum = length(d), desc = paste(sort(shortest_unique_prefixes(d)), collapse = '_')) %>% ungroup %>% select(from, dist, d, dnum, desc) %>% distinct

  p2 = leafdist %>% filter(dist == 2) %>% filter(from %in% parents$from) %>% group_by(from) %>% mutate(d = list(union(from, to))) %>% rowwise %>% mutate(dnum = length(d), desc = paste(sort(shortest_unique_prefixes(d)), collapse = '_')) %>% ungroup %>% select(from, dist, d, dnum, desc) %>% distinct

  p2 %>% left_join(parents %>% select(from, d), by = 'from') %>% rowwise %>% mutate(d = list(union(to, d)), dnum = length(d), desc = paste(sort(shortest_unique_prefixes(d)), collapse = '_')) %>% ungroup


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
      newnam = shortest_unique_prefixes(reachable) %>% sort %>% paste(collapse = '_') %>% paste0('_', .)
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


eval_plusnadmix = function(graph, qpgfun, n = 1, ntry = Inf, nonzero_f4 = NULL, admix_constraints = NULL,
                           event_order = NULL, verbose = TRUE) {

  newgraphs = tibble(graph = graph_plusn(list(graph), n = n, ntry = ntry))
  num = nrow(newgraphs)
  #newgraphs %<>% slice_sample(n = ntry)
  if(verbose) alert_info(paste0('Evaluating ', nrow(newgraphs), ' graphs...\n'))
  if(nrow(newgraphs) == 0) return(tibble())
  newgraphs %<>%
    rowwise %>%
    filter(satisfies_constraints(graph, nonzero_f4 = nonzero_f4, admix_constraints = admix_constraints,
                                 event_order = event_order)) %>% ungroup
  if(nrow(newgraphs) == 0) return(tibble())
  newgraphs %>%
    mutate(res = furrr::future_map(graph, qpgfun, .progress = verbose, .options = furrr::furrr_options(seed = TRUE))) %>%
    unnest_wider(res) %>% arrange(score)
  #%>% select(source_from, source_to, dest_from, dest_to, graph, score)
}

eval_minusnadmix = function(graph, qpgfun, n = 1, ntry = Inf, nonzero_f4 = NULL,
                            admix_constraints = NULL, event_order = NULL, verbose = TRUE) {

  newgraphs = tibble(graph = graph_minusn(list(graph), n = n, ntry = ntry))
  if(verbose) alert_info(paste0('Evaluating ', nrow(newgraphs), ' graphs...\n'))
  if(nrow(newgraphs) == 0) return(tibble())
  newgraphs %<>%
    rowwise %>%
    filter(satisfies_constraints(graph, nonzero_f4 = nonzero_f4, admix_constraints = admix_constraints,
                                 event_order = event_order)) %>% ungroup
  if(nrow(newgraphs) == 0) return(tibble())
  newgraphs %>%
    mutate(res = furrr::future_map(graph, qpgfun, .progress = verbose, .options = furrr::furrr_options(seed = TRUE))) %>%
    unnest_wider(res) %>% arrange(score)
}

eval_plusminusn = function(graph, qpgfun, n = 1, ntry = Inf, nonzero_f4 = NULL,
                           admix_constraints = NULL, event_order = NULL, verbose = TRUE) {

  if(verbose) alert_info(paste0('Adding ', n, ' edge(s)...\n'))
  plus = eval_plusnadmix(graph, qpgfun, n = n, ntry = ntry, nonzero_f4 = nonzero_f4,
                         admix_constraints = admix_constraints, event_order = event_order, verbose = verbose)
  if(verbose) alert_info(paste0('Best score: ', round(plus$score[[1]], 3), '\n'))

  if(verbose) alert_info(paste0('Removing ', n, ' edge(s)...\n'))
  plus %>% slice_min(score, with_ties = FALSE) %>% pull(graph) %>% pluck(1) %>%
    eval_minusnadmix(qpgfun, n = n, ntry = Inf, nonzero_f4 = nonzero_f4,
                     admix_constraints = admix_constraints, event_order = event_order, verbose = verbose)
}


eval_plusonepop = function(graph, pop, qpgfun, ntry = Inf, verbose = TRUE) {

  root = get_rootname(graph)
  outgroup = intersect(names(neighbors(graph, root)), get_leafnames(graph))
  newgraphs = graph %>% as_edgelist() %>% as_tibble(.name_repair = ~c('from', 'to'))
  if(length(outgroup) > 0) newgraphs %<>% filter(from != root | to != outgroup)
  newgraphs %<>% mutate(graph = map2(from, to, ~insert_leaf(graph, pop, .x, .y)))

  if(verbose) alert_info(paste0('Found ',nrow(newgraphs),' graphs. Evaluating ', min(nrow(newgraphs), ntry), '...\n'))
  newgraphs %>%
    slice_sample(n = min(100, ntry), replace = TRUE) %>%
    mutate(res = furrr::future_map(graph, qpgfun, .progress = verbose, .options = furrr::furrr_options(seed = TRUE))) %>%
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
  # assumes graph is simplified, edge is regular edge

  sib = graph %>% neighbors(from, mode = 'out') %>% names %>% setdiff(to)
  children = graph %>% neighbors(to, mode = 'out') %>% names %>% sample
  if(length(children) == 0 || length(sib) == 0) {
    return(graph)
  }
  og = get_outpop(graph)
  gnew = graph %>%
    add_edges(c(to, sib, from, children[1])) %>%
    delete_edges(paste(c(from, to), c(sib, children[1]), sep = '|'))
  if(is_valid(gnew) && (is.null(og) || isTRUE(og == get_outpop(gnew)))) return(gnew)
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
  if(length(grandparent_pos) != 1) return()
  newgraph = graph %>%
    add_edges(c(grandparent_pos, to, to, parent_pos)) %>%
    delete_edges(paste(c(grandparent_pos, parent_pos, parent_neg), c(parent_pos, to, to), sep = '|'))
  if(!is_valid(newgraph)) {
    return()
  }
  if(length(parent_neg) != 1 || length(parent_pos) != 1) stop('aaa')
  newgraph2 = newgraph %>% add_edges(c(parent_neg, parent_pos)) %>% simplify_graph
  if(is_valid(newgraph2)) newgraph = newgraph2
  else newgraph %<>% find_newedges %>% slice_sample(n=1) %>%
    mutate(g = list(insert_admix(newgraph, source_from, source_to, dest_from, dest_to))) %>%
    pull(g) %>% pluck(1) %>% simplify_graph
  leaves2 = sort(get_leafnames(newgraph))
  if(!isTRUE(all.equal(leaves, leaves2)) || !is_valid(newgraph)) return()
  nold = numadmix(graph)
  nnew = numadmix(newgraph)
  if(nnew < nold) newgraph %<>% insert_admix_n(nold - nnew)
  if(is_valid(newgraph)) return(newgraph)
}


#' Find well fitting admixture graphs
#'
#' This function generates and evaluates admixture graphs in `numgen` iterations
#' to find well fitting admixturegraphs.
#' @export
#' @param data Input data in one of three forms:
#' \enumerate{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{f2_from_geno}}
#' \item A directory which contains pre-computed f2-statistics
#' \item The prefix of genotype files
#' }
#' @param numadmix Number of admixture events within each graph. (Only relevant if `initgraph = NULL`)
#' @param outpop Name of the outgroup population
#' @param numgraphs Number of graphs in each generation
#' @param stop_gen Total number of generations after which to stop
#' @param stop_gen2 Number of generations without improvement after which to stop
#' @param stop_sec Number of seconds after which to stop
#' @param stop_score Stop once this score has been reached
#' @param initgraph Graph to start with. If it is specified, `numadmix` and `outpop` will be inferred from this graph.
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
#' @param opt_worst_residual Optimize for lowest worst residual instead of best score. `FALSE` by default, because the likelihood score is generally a better indicator of the quality of the model fit, and because optimizing for the lowest worst residual is slower (because f4-statistics need to be computed).
#' @param return_searchtree Return the search tree in addition to the models. Output will be a list with three items: models, search tree, search tree as data frame
#' @param plusminus_generations If the best score does not improve after `plusminus_generations` generations, another approach to improving the score will be attempted: A number of graphs with on additional admixture edge will be generated and evaluated. The resulting graph with the best score will be picked, and new graphs will be created by removing any one admixture edge (bringing the number back to what it was originally). The graph with the lowest score will then be selected. This often makes it possible to break out of local optima, but is slower than regular graph modifications.
#' If the current number of admixture events is lower than `max_numadmix`, the last step (removing an admixture edge) will be skipped.
#' @param admix_constraints A data frame with constraints on the number of admixture events for each population.
#' See \code{\link{satisfies_numadmix}}
#' As soon as one graph happens to satisfy these constraints, all subsequently generated graphs will be required to also satisfy them.
#' @param event_constraints A data frame with constraints on the order of events in an admixture graph.
#' See \code{\link{satisfies_eventorder}}
#' As soon as one graph happens to satisfy these constraints, all subsequently generated graphs will be required to also satisfy them.
#' @param reject_f4z If this is a number greater than zero, all f4-statistics with `abs(z) > reject_f4z` will be used to constrain the search space of admixture graphs: Any graphs in which f4-statistics greater than `reject_f4z` are expected to be zero will not be evaluated.
#' @param max_admix Maximum number of admixture edges. By default, this number is equal to `numadmix`, or to the number of admixture edges in `initgraph`, so the number of admixture edges stays constant. Setting this to a higher number will lead to more admixture edges being added occasionally (see `plusminus_generations`). Graphs with additional admixture edges will only be accepted if they improve the score by 5% or more.
#' @param verbose Print progress updates
#' @param ... Additional arguments passed to \code{\link{qpgraph}}
#' @return A nested data frame with one model per line
#' @seealso \code{\link{qpgraph}}, \code{\link{find_graphs_old}}
#' @examples
#' \dontrun{
#' res = find_graphs(example_f2_blocks, numadmix = 2)
#' res %>% slice_min(score)
#' }
#' \dontrun{
#' # Start with a graph with 0 admixture events, increase up to 3, and stop after 10 generations of no improvement
#' pops = dimnames(example_f2_blocks)[[1]]
#' initgraph = random_admixturegraph(pops, 0, outpop = 'Chimp.REF')
#' res = find_graphs(example_f2_blocks, initgraph = initgraph, stop_gen2 = 10, max_admix = 3)
#' res %>% slice_min(score)
#' }
find_graphs = function(data, numadmix = 0, outpop = NULL, stop_gen = 100, stop_gen2 = 15, stop_score = 0, stop_sec = NULL,
                       initgraph = NULL, numgraphs = 10,
                       mutfuns = namedList(spr_leaves, spr_all, swap_leaves, move_admixedge_once, flipadmix_random,
                                            place_root_random, mutate_n),
                       opt_worst_residual = FALSE, plusminus_generations = 5, return_searchtree = FALSE,
                       admix_constraints = NULL, event_constraints = NULL, reject_f4z = 0, max_admix = numadmix,
                       verbose = TRUE, ...) {

  nodups = TRUE
  traceback_gen = Inf

  f2_blocks = get_f2(data)
  if(is.null(initgraph)) {
    pops = dimnames(f2_blocks)[[1]]
    graph = random_admixturegraph(pops, numadmix, outpop = outpop)
  } else {
    graph = initgraph
    numadmix = numadmix(graph)
  }
  ell = list(...)
  if(!all(names(ell) %in% names(formals(qpgraph)))) {
    notused = setdiff(names(ell), names(formals(qpgraph)))
    stop(paste0("The following arguments are not recognized: '", paste0(notused, collapse = "', '"), "'"))
  }
  if('f3basepop' %in% names(ell)) {
    f3basepop = ell$f3basepop
  } else {
    f3basepop = get_leafnames(graph)[1]
  }
  if('numstart' %in% names(ell)) {
    numstart = ell$numstart
  } else {
    numstart = 1
  }
  if(!isFALSE(opt_worst_residual)) {
    qpgfun = function(graph, ...) {
      res = qpgraph(f2_blocks, graph, numstart = numstart, return_fstats = opt_worst_residual, f3basepop = f3basepop, ...)
      res$score = res$worst_residual
      res$worst_residual = res$f4 = NULL
      res
      }
  } else {
    if('f3precomp' %in% names(ell)) {
      f3precomp = ell$f3precomp
    } else {
      f3precomp = qpgraph_precompute_f3(f2_blocks, get_leafnames(graph), f3basepop = f3basepop)
    }
    qpgfun = function(graph, ...) qpgraph(NULL, graph, numstart = numstart, f3precomp = f3precomp, ...)
  }
  stop_at = Sys.time() + stop_sec
  wfuns = namedList(rearrange_negadmix3, replace_admix_with_random)
  dfuns = namedList(rearrange_negdrift)
  allfuns = c(wfuns, dfuns, mutfuns)

  models = tibble(g = list(simplify_graph(graph))) %>%
    mutate(res = map(g, qpgfun), hash = map_chr(g, graph_hash), lasthash = graph_hash(graph)) %>%
    unnest_wider(res) %>%
    #mutate(edges = map(edges, ~filter(., weight < 0))) %>%
    transmute(gen2 = 0, hash, g, edges, score)
  best = besthist = models$score
  bestmut = 'none'
  gimp = numtraceback = 0
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
  nonzero_f4 = if(reject_f4z > 0) f4(f2_blocks) %>% filter(abs(z) > reject_f4z) else NULL
  nzf4 = admixc = eventc = lasttracebacklevel = NULL
  lasttracebackscore = Inf

  for(i in seq_len(stop_gen)) {

    if(best <= stop_score || gimp > stop_gen2 || !is.null(stop_sec) && Sys.time() > stop_at) break
    if(verbose) {
      sc = besthist[max(1,i-1)]
      alert = if(isTRUE(sc == best)) alert_success else alert_info
      #msg = paste0(i, ': sc ', round(besthist[max(1,i-1)], 3), '\tbest ', round(best, 3),'\tnew ', if(i==1) 1 else nrow(newmod), '\ttot ', sum(!is.na(models$score)), ' ', sum(stdat$totalCount == 1), ' ', sum(stdat$score[stdat$totalCount == 1] < best*2),'\t')
      msg = paste0(i, ': score ', round(besthist[max(1,i-1)], 3), '\tbest score ', round(best, 3),'\t')
      na = if(i == 1 || is.null(newmod) || nrow(newmod) == 0) numadmix(graph) else numadmix(newmod$g[[1]])
      if(numadmix != max_admix && i > 1) msg %<>% paste0('numadmix ', na, '\t')
      msg %<>% paste0(bestmut, '\ttot ', sum(!is.na(models$score)))
    }

    if(gimp >= traceback_gen) {
      if(is.null(lasttracebacklevel)) lasttracebacklevel = newgraph$level#-traceback_gen
      #if(is.null(lasttracebacklevel)) lasttracebacklevel = 2
      else lasttracebacklevel = lasttracebacklevel - 1
      numtraceback = numtraceback + 1
      #data.tree::Prune(st, function(x) x$level < lasttracebacklevel)
      if(lasttracebacklevel > 0) {
        st$Do(function(x) x$closed = x$level >= lasttracebacklevel)
        stdat = st %>% st_to_dat
        gimp = 0
        lasttracebackscore = stdat %>% filter(level <= lasttracebacklevel) %>% slice_min(score) %>% pluck('score', 1)
        if(is.null(lasttracebackscore)) lasttracebackscore = Inf
        if(verbose) alert_info(paste0('Traceback to level ', lasttracebacklevel, ' with score ', round(lasttracebackscore, 3), '\n'))
      }
    }
    #newgraph = pick_graph(stdat, verbose = verbose > 2)
    newgraph = stdat %>% filter(!isTRUE(closed)) %>% slice_min(score, with_ties = FALSE)
    if(nrow(newgraph) == 0) newgraph = stdat %>% slice_min(score, with_ties = FALSE)

    if(verbose > 1) {
      msg %<>% paste0(newgraph %$% paste0(' ', round(score, 2), ' ', round(score/best, 2),'\tg ',gen2, ' ', gimp, '\tcnt ',cntm_1,' ',cntm_5, '\ttime '), round(as.numeric(Sys.time()-tm), 1), '\t', bestmut)
      tm = Sys.time()
    }
    if(verbose) alert(paste0(msg, '\n'))
    sel = models %>% filter(hash == newgraph$hash[[1]])

    graph = sel$g[[1]]
    if(!is.null(outpop) && is.null(get_outpop(graph))) stop('fg error')

    if(is.null(nzf4) && reject_f4z > 0 && satisfies_nonzerof4(graph, nonzero_f4)) {
      nzf4 = nonzero_f4
      if(verbose) alert_info('All zero f4 ok!\n')
    }
    if(is.null(admixc) && !is.null(admix_constraints) && satisfies_numadmix(graph, admix_constraints)) {
      admixc = admix_constraints
      if(verbose) alert_info('Admix constraints ok!\n')
    }
    if(is.null(eventc) && !is.null(event_constraints) && satisfies_eventorder(graph, event_constraints)) {
      eventc = event_constraints
      if(verbose) alert_info('Event constraints ok!\n')
    }

    if(gimp > 0 && gimp %% plusminus_generations == 0) {

      sel = models %>% slice_min(score, with_ties = FALSE)
      graph = sel$g[[1]]
      tryCatch({
      newmod = eval_plusnadmix(graph, qpgfun, n = 1, ntry = numgraphs*10,
                               nonzero_f4 = nzf4, admix_constraints = admixc,
                               event_order = eventc, verbose = verbose)
      if(nrow(newmod) == 0) next
      newmod %<>% slice_min(score, with_ties = FALSE)
      }, error = function(e) stop('fg error 2'))
      mf = 'plusnadmix'
      excess_admix = numadmix(newmod$graph[[1]]) - max_admix
      if(excess_admix > 0) {
        mf = 'plusminusn'
        newmod = eval_minusnadmix(newmod$graph[[1]], qpgfun, n = excess_admix, ntry = numgraphs*10,
                                  nonzero_f4 = nzf4, admix_constraints = admixc,
                                  event_order = eventc, verbose = verbose)
        if(nrow(newmod) == 0) next
        newmod %<>% slice_min(score, with_ties = FALSE)
      }
      newmod %<>% mutate(hash = map_chr(graph, graph_hash)) %>%
        transmute(hash, lasthash = sel$hash[[1]], g = graph, gen2 = i, edges, score, mutfun = mf)
      if(newmod$score > min(models$score, na.rm=F)*0.95) {
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
                 g = pmap(list(fun, from, to), function(x, y, z) possibly(x, NULL)(graph, y, z))) %>%
          rowwise %>% filter(!is.null(g)) %>% ungroup
      }
      remaining = numgraphs - nrow(newmod)

      if(remaining > 0) {
        swape = swap_negdrift(graph, e$from, e$to)
        if(!is.null(swape) && length(swape) > 0) {
          newswap = tibble(g = swape, mutfun = 'swap_negdrift') %>%
            slice_sample(n = min(length(swape), ceiling(remaining/2)))
          newmod %<>% bind_rows(newswap)
        }
        remaining = numgraphs - nrow(newmod)
      }
      if(remaining > 0) {

      tryCatch({
        randmut = tibble(mutfun = names(mutfuns), fun = map(mutfuns, ~possibly(., NULL))) %>%
          slice_sample(n = remaining, replace = TRUE) %>%
          rowwise %>% mutate(g = list(fun(graph))) %>% filter(!is.null(g)) %>% ungroup %>% select(-fun)
        newmod %<>% bind_rows(randmut)
        }, error = function(e) if(verbose) alert_warning("Could not apply mutation function"))
      }

      newmod %<>% rowwise %>% filter(satisfies_constraints(g, nzf4, admixc, eventc)) %>% ungroup
      newmod %<>%
        transmute(g, hash = map_chr(g, graph_hash), lasthash = sel$hash[[1]], mutfun)
      if(nodups) newmod %<>% filter(!duplicated(hash), !hash %in% models$hash)
      tryCatch({
      if(!is.null(outpop) && nrow(newmod) > 0) newmod %<>% rowwise %>% mutate(og = list(get_outpop(g))) %>%
        filter(!is.null(og) && og == outpop) %>% ungroup
      }, error = function(e) browser())
      if(nrow(newmod) == 0) {
        if(verbose) alert_danger('No new models!\n')
        next
      }
      newmod %<>%
        mutate(res = furrr::future_map(g, qpgfun, .options = furrr::furrr_options(seed = TRUE))) %>%
        unnest_wider(res) %>%
        transmute(gen2 = i, hash, lasthash, g, edges, score, mutfun) %>%
        arrange(score)
    }

    besthist[i] = newmod$score[1]
    bestmut = newmod$mutfun[1]
    if(besthist[i] < best) {
      best = besthist[i]
      gimp = 0
    } else gimp = gimp + 1
    if(nodups) newmod %<>% filter(!hash %in% models$hash)
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

  mout = models %>% transmute(generation = gen2, graph = unname(g), edges, score, mutation = mutfun, hash, lasthash)
  if(return_searchtree) return(list(models = mout, tree = st, dat = stdat))
  mout
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

  as.data.frame(st, NULL, T, 'name', 'score', 'scbest', 'gimp', 'scp_1', 'scp_5', 'scm_1', 'scm_5', 'height', 'level', 'totalCount', 'cntm_1', 'cntm_5', 'gen2', 'closed') %>%
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
  if(maxpot == Inf) stop('maxpot inf')
  out %<>% slice_sample(n = 1, weight_by = potential^10)
  if(verbose) alert_warning(out %$% paste0('\tcnt ', sum(stdat$totalCount == 1), '\tpen ', paste0(round(g1, 2), ' ', round(g5, 2), ' ', round(imp, 2), '\tpot ', round(potential, 2), ' ', round(maxpot, 2), '\n')))
  out
}

pick_graph2 = function(stdat, lasttracebacklevel, tracebackafter = 4, verbose = FALSE) {

  sel = stdat %>% slice_min(score, with_ties = FALSE)
  lastimp = max(stdat$gen2, na.rm = TRUE) - sel$gen2
  if(sel$totalCount[[1]] == 1 || lastimp < tracebackafter) {
    return(sel)
  }
  if(is.null(lasttracebacklevel)) lasttracebacklevel = max(stdat$level)
  sel = stdat %>% filter(level < lasttracebacklevel, height == 1) %>% slice_max(level) %>% slice_min(score, with_ties = FALSE)
  if(nrow(sel) == 0) stop('pick graph error')
  if(verbose) print(sel)
  sel
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

  leaves = sort(get_leafnames(graph))
  g = graph %>% delete_admix(from, to)
  g %<>% insert_admix
  g
}

#' Get unique hash of an admixture graph
#'
#' This can be used to check if two graphs are identical. For two graphs to be identical, the topology and leaf nodes have to match, but internal node names do not matter.
#' @export
#' @param graph Admixture graph in igraph format
#' @return A hash string of the admixture graph
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
    group_by(from, to, from2, to2) %>%
    mutate(uniq = list(names(which(table(unlist(is))==1)))) %>%
    rowwise %>%
    mutate(allunique = all(is %in% uniq)) %>%
    ungroup

}

# unidentifiable_admixture = function(graph) {
#
#   stopifnot(is_valid(graph) && is_simplified(graph))
#
#   leaves = graph %>% get_leafnames
#   adm = names(which(degree(graph, mode = 'in') > 1))
#   unresolved_adm = adm
#   pis = path_intersections(graph)
#
#   for(i in 1:10) {
#     if(length(unresolved_adm) == 0) break
#
#     pis %<>% rowwise %>%
#       mutate(unresolved = list(intersect(unresolved_adm, union(admnodes, admnodes2))), unum = length(unresolved)) %>%
#       ungroup
#
#     resolved0 = pis %>% filter(unum == 0, islen > 0) %>% pull(isstr) %>% unique
#     resolved1 = pis %>% filter(unum == 1, islen > 0, allunique) %>% pull(isstr) %>% unique
#
#     resolved_adm = pis %>% filter(isstr %in% intersect(resolved0, resolved1), unum == 1, allunique) %$%
#       unique(unlist(c(admnodes, admnodes2)))
#
#     if(length(resolved_adm) == 0) break
#
#     unresolved_adm %<>% setdiff(resolved_adm)
#   }
#   unresolved_adm
# }

paths_through = function(graph, leaves, node) {

  nl = graph %>% neighbors(node, mode = 'out') %>% pluck(1) %>% names
  nr = graph %>% neighbors(node, mode = 'out') %>% pluck(2) %>% names

  paths = all_simple_paths(graph, node, leaves, mode = 'out') %>%
    map(names)

  expand_grid(up = paths %>% keep(~.[[2]] == nl), down = paths %>% keep(~.[[2]] == nr)) %>%
    bind_rows(rename(., up=down, down=up))
}

pair_paths = function(graph) {

  leaves = graph %>% get_leafnames
  adm = names(which(degree(graph, mode = 'in') > 1))
  internal = setdiff(names(V(graph)), leaves)
  nonadm = setdiff(internal, adm)
  map(nonadm, ~paths_through(graph, leaves, .)) %>% set_names(nonadm) %>% bind_rows(.id = 'pivot') %>% rowwise %>% mutate(from = tail(up, 1), to = tail(down, 1), path = list(c(vs_to_es(up), vs_to_es(down))), admnodes = list(intersect(adm, c(up, down))), admnum = length(admnodes)) %>% ungroup

}

# continue here: use path intersections to determine edges which are not identifiable


tabulate_subgraphs = function(graphlist, pops) {
  # for a list of graphs and set of populations in that graph
  # returns data frame counting all subgraphs consisting of those populations

  #delete_leaves = function(graph, leaves) reduce(leaves, delete_leaf, .init = graph)
  graphlist %>%
    tibble(.name_repair = ~'graph') %>%
    rowwise %>%
    mutate(subgraph = list(keep_leaves(graph, pops)), hash = graph_hash(subgraph)) %>%
    ungroup %>%
    select(-graph) %>%
    add_count(hash) %>%
    filter(!duplicated(hash)) %>%
    arrange(-n)
}

tabulate_subgraphs_n = function(graphlist, pops, n, outpop = NULL, verbose = TRUE) {

  cmb = t(combn(setdiff(pops, outpop), n))
  if(verbose) alert_info(paste0('Tabulating ', nrow(cmb), ' combinations of ', n, ' populations...'))
  #pb = progress::progress_bar$new(total = nrow(cmb))
  as_tibble(cmb, .name_repair = ~paste0('X', seq_len(n))) %>%
    mutate(i = 1:n()) %>% rowwise %>%
    mutate(subgraphs = list(tabulate_subgraphs(graphlist, c(outpop, cmb[i,])))) %>%
    unnest(subgraphs) %>% arrange(-n)
}


count_admix = function(graph) {

  leaves = get_leafnames(graph)
  admnodes = names(which(degree(graph, mode = 'in') > 1))
  map_dbl(leaves, ~length(intersect(names(subcomponent(graph, .x, mode = 'in')), admnodes))) %>%
    set_names(leaves) %>% enframe('pop', 'nadmix') %>% arrange(pop)
}

unadmixed = function(graph) {

  graph %>% count_admix %>% filter(nadmix == 0) %>% pull(pop)
}


#' Assign proxy populations to admixed populations
#'
#' @export
#' @param graph An admixture graph in igraph format, or as edge list with column `weight`
#' @return A data frame with columns `pop`, `nadmix`, `proxy`, `nproxies` (and `weight`)
#' @examples
#' \dontrun{
#' summarize_proxies(example_igraph)
#' }
summarize_proxies = function(graph) {

  if('data.frame' %in% class(graph)) {
    edges = graph
    graph %<>% edges_to_igraph()
  } else {
    edges = NULL
  }
  leaves = get_leafnames(graph)
  admnodes = names(which(degree(graph, mode = 'in') > 1))
  admparents = map(admnodes, ~names(neighbors(graph, ., mode = 'in'))) %>% do.call(rbind, .)

  unadm_desc = function(graph, v, leaves) {
    all_simple_paths(graph, v, leaves, mode = 'out') %>%
      discard(~any(names(.)[-1] %in% admnodes)) %>%
      set_names(map_chr(., ~tail(names(.), 1))) %>% names %>% table %>%
      enframe %>% filter(value == 1) %>% pull(name)
  }
  vnam = names(V(graph))
  ud = map(vnam, ~unadm_desc(graph, ., leaves)) %>% set_names(vnam) #%>% keep(~length(.) > 0)

  # very inefficient; change this if function is slow
  find_keys = function(graph, v, stopv) {
    root = names(get_root(graph))
    for(n in names(all_simple_paths(graph, v, root, mode = 'in')[[1]])) {
      #if(stopv %in% names(subcomponent(graph, n, mode = 'out'))) return(NULL)
      if(length(ud[[n]]) > 0) return(ud[[n]])
      if(degree(graph, n, mode = 'in') > 1) {
        par = names(neighbors(graph, n, mode = 'in'))
        return(c(find_keys(graph, par[1], par[2]), find_keys(graph, par[2], par[1])))
      }
    }
  }
  admkeys = tibble(admnodes,
                   left = map(1:nrow(admparents), ~find_keys(graph, admparents[.,1], admparents[.,2])),
                   right = map(1:nrow(admparents), ~find_keys(graph, admparents[.,2], admparents[.,1]))) %>%
    rowwise %>%
    mutate(is = list(intersect(left, right)), left = list(setdiff(left, is)), right = list(setdiff(right, is)),
           left = list(ifelse(length(left) == 0, NA, left)), right = list(ifelse(length(right) == 0, NA, right))) %>%
    select(-is) %>% ungroup
  if(!is.null(edges)) {
    admkeys %<>% left_join(edges %>% filter(type == 'admix', !duplicated(to)) %>%
                             transmute(admnodes = to, wl = weight), by = 'admnodes') %>%
      mutate(wr = 1-wl)
  } else {
    admkeys %<>% mutate(wl = 0, wr = 0)
  }
  lastadm = map(leaves, ~intersect(names(subcomponent(graph, .x, mode = 'in')), admnodes)) %>%
    set_names(leaves) %>% discard(~length(.) == 0) %>% map(~head(., 1)) %>% unlist %>% enframe('pop', 'admnodes')

  proxies = lastadm %>%
    left_join(admkeys %>% transmute(admnodes, proxy = left, weight = wl), by = 'admnodes') %>%
    unnest(proxy) %>%
    bind_rows(lastadm %>%
                left_join(admkeys %>% transmute(admnodes, proxy = right, weight = wr), by = 'admnodes') %>%
                unnest(proxy)) %>%
    select(-admnodes) %>% add_count(pop, name = 'nproxies')
  if(is.null(edges)) proxies %<>% select(-weight)
  graph %>% count_admix %>% filter(nadmix > 0) %>% left_join(proxies, by = 'pop')
}

#' List proxy populations in graphlist
#'
#' @export
#' @param graphlist A list of admixture graphs
#' @return A data frame with columns `pop`, `proxy`, `nadmix`, `nproxies`, `n`, `frac`
#' @examples
#' \dontrun{
#' summarize_proxies_list(graphlist)
#' }
summarize_proxies_list = function(graphlist) {
  # need to update for list of fits

  graphlist %>% map(summarize_proxies) %>% bind_rows(.id = 'g') %>%
    group_by(pop) %>% add_count(proxy) %>%
    group_by(pop, proxy) %>%
    summarize(nadmix = mean(nadmix), nproxies = mean(nproxies), n = mean(n), .groups = 'drop') %>%
    mutate(frac = n/length(graphlist)) %>% arrange(-n) %>% ungroup
}



proxypred = function(sim, graphlist, auc = FALSE) {

  pops = get_leafnames(sim)

  out = expand_grid(pop = pops, proxy = pops) %>% filter(pop != proxy) %>%
    left_join(summarize_proxies(sim), by = c('pop', 'proxy')) %>%
    transmute(pop, proxy, obs = !is.na(nadmix)) %>%
    left_join(summarize_graphlist(graphlist), by = c('pop', 'proxy')) %>%
    transmute(pop, proxy, obs, pred = replace_na(frac, 0))
  if(auc) out = out %$% pROC::auc(obs, pred) %>% c
  out
}


node_events = function(graph, leaves, node) {
  # for a node, list associated population splits and merges

  children = names(neighbors(graph, node, mode = 'out'))
  if(length(children) == 0) {
    #out = tibble(type = 'leaf', pop1 = node, pop2 = NA)
    out = tibble()
  } else if(length(children) == 1) {
    # parents = neighbors(graph, node, mode = 'in')
    # desc1 = subcomponent(graph, parents[1], mode = 'out') %>% names %>% intersect(leaves)
    # desc2 = subcomponent(graph, parents[2], mode = 'out') %>% names %>% intersect(leaves)
    # out = expand_grid(type = 'merge', pop1 = setdiff(desc1, desc2), pop2 = setdiff(desc2, desc1))
    out = tibble()
  } else if(length(children) == 2) {
    desc1 = subcomponent(graph, children[1], mode = 'out') %>% names %>% intersect(leaves)
    desc2 = subcomponent(graph, children[2], mode = 'out') %>% names %>% intersect(leaves)
    out = expand_grid(type = 'split', pop1 = setdiff(desc1, desc2), pop2 = setdiff(desc2, desc1))
    if(children[1] %in% leaves) out %<>% bind_rows(tibble(type = 'leaf', pop1 = children[1], pop2 = NA))
    if(children[2] %in% leaves) out %<>% bind_rows(tibble(type = 'leaf', pop1 = children[2], pop2 = NA))
  } else {
    stop('More than two children!')
  }
  out
}


all_node_events = function(graph) {
  # for each node in a graph, list associated population splits and merges

  leaves = get_leafnames(graph)
  graph %>% V %>% names %>% set_names %>% map(~node_events(graph, leaves, .)) %>% bind_rows(.id = 'node')

}

#' List population split events in a graph
#'
#' @export
#' @param graph An admixture graph
#' @return A data frame with columns `earlier1`, `earlier2`, `later1`, `later2`
#' @examples
#' \dontrun{
#' summarize_eventorder(example_igraph)
#' }
summarize_eventorder = function(graph, unique_only = TRUE) {

  node_pairs = graph %>% igraph::distances(mode = 'out') %>% as_tibble(rownames = 'from') %>%
    pivot_longer(-from, names_to = 'to', values_to = 'order') %>%
    mutate(order = is.finite(order)+0) %>% filter(order > 0, from != to) %>%
    bind_rows(rename(., from = to, to = from) %>% mutate(order = -order))

  events = all_node_events(graph) %>%
    bind_rows(rename(., pop2 = pop1, pop1 = pop2) %>% filter(type != 'leaf'))

  out = expand_grid(rename_with(events, ~paste0(., '_1')), rename_with(events, ~paste0(., '_2'))) %>%
    filter(node_1 != node_2) %>% left_join(node_pairs, by = c('node_1'='from', 'node_2'='to')) %>%
    mutate(order = replace_na(order, 0))
  if(unique_only) {
    out %<>% group_by(type_1, pop1_1, pop2_1, type_2, pop1_2, pop2_2) %>%
      filter(all(order != -1) & any(order == 1)) %>% ungroup %>%
      transmute(earlier1 = pmin(pop1_1, pop2_1, na.rm=T),
                earlier2 = pmax(pop1_1, pop2_1),
                later1 = pmin(pop1_2, pop2_2, na.rm=T),
                later2 = pmax(pop1_2, pop2_2), type1 = type_1, type2 = type_2) %>% distinct
  }
  out
}

#' List population split events in a list of graphs
#'
#' @export
#' @param graphlist A list of admixture graphs
#' @return A data frame with columns `earlier1`, `earlier2`, `later1`, `later2`, `n`, `frac`
#' @examples
#' \dontrun{
#' summarize_eventorder(graphlist)
#' }
summarize_eventorder_list = function(graphlist) {

  tot = length(graphlist)
  graphlist %>% map(summarize_eventorder) %>% bind_rows(.id = 'g') %>%
    count(earlier1, earlier2, later1, later2) %>% mutate(frac = n/tot)
}


# unresolvable_params = function(groebner, vrs) {
#
#   resolved = c()
#   remaining = groebner %>% map(~intersect(vrs, unlist(map(., ~head(names(.), -1))))) %>% keep(~length(.) > 0)
#
#   repeat {
#   newresolved = remaining %>% keep(~length(.) == 1) %>% unlist %>% unique
#   resolved %<>% c(newresolved)
#   remaining %<>% keep(~length(.) > 1) %>% map(~setdiff(., resolved))
#   if(length(newresolved) == 0 || length(remaining) == 0) break
#   }
#   remaining
# }


root_paths = function(graph) {

  leaves = get_leafnames(graph)
  root = get_rootname(graph)
  adm = names(which(degree(graph, mode = 'in') > 1))
  parents = adm %>% rlang::set_names() %>% map(~names(neighbors(graph, ., mode = 'in')))
  admedges = parents %>% imap(~paste(.x, .y, sep = '|'))

  paths = all_simple_paths(graph, root, leaves, mode = 'out') %>% map(names) %>%
    tibble(.name_repair = ~'path') %>% rowwise %>%
    mutate(pop = tail(path, 1),
           edges = list(setdiff(vs_to_es(path), admedges)),
           admix = list(intersect(path, adm)),
           lr = list(map(admedges[admix], 1) %in% edges + 0),
           edges = list(setdiff(edges, unlist(admedges))),
           admedges = list(map2_chr(admedges[admix], 2-lr, ~.x[.y]))) %>%
    ungroup %>% add_count(pop)
  paths
}

path_pairs = function(graph) {

  # returns all pairs of paths from root to all leaves

  paths = root_paths(graph) %>% select(-admedges)
  expand_grid(rename_with(paths, ~paste0(.,'1')),
              rename_with(paths, ~paste0(.,'2'))) %>%
    filter(pop1 != pop2) %>% rowwise %>%
    mutate(edges = list(setdiff(union(edges1, edges2), intersect(edges1, edges2))),
           admix = list(c(admix1, admix2)), lr = list(c(lr1, lr2))) %>%
    ungroup
}

path_triples = function(graph) {

  paths = root_paths(graph) #%>% select(-edges)
  expand_grid(rename_with(paths, ~paste0(.,'1')),
              rename_with(paths, ~paste0(.,'2')),
              rename_with(paths, ~paste0(.,'3'))) %>%
    filter(pop1 != pop2, pop1 != pop3, pop2 != pop3) %>% rowwise %>%
    mutate(admix12 = list(c(admix1, admix2)),
           admix13 = list(c(admix1, admix3)),
           admix23 = list(c(admix2, admix3)),
           lr12 = list(c(lr1, lr2)),
           lr13 = list(c(lr1, lr3)),
           lr23 = list(c(lr2, lr3)),
           last12 = tail(intersect(path1, path2), 1),
           last13 = tail(intersect(path1, path3), 1),
           last23 = tail(intersect(path2, path3), 1),
           m121 = match(last12, path1),
           m131 = match(last13, path1),
           og = case_when(m121 < m131 ~ pop2,
                          m121 < match(last23, path2) ~ pop1,
                          m131 < match(last12, path1) ~ pop3)) %>%
    ungroup %>% select(-m121, m131)
}

#' Find well fitting admixture graphs
#'
#' This function generates and evaluates admixture graphs in `numgen` iterations
#' to find well fitting admixturegraphs.
#' @export
#' @param graph An admixture graph
#' @param substitute Should edge names be represented by shorter symbols?
#' @param nam Symbols used to shorten edge names
#' @return A list with two data frames: `equations` holds the equtions for all f2-statistics; `coding` has the mapping from edge names to edge symbols, which is used when `substitute = TRUE`
graph_equations = function(graph, substitute = TRUE, nam = c('a', 'e', 'f'), return_everything = FALSE) {

  pp = path_pairs(graph)
  leaves = get_leafnames(graph)
  e = pp$edges %>% unlist %>% unique
  a = names(which(degree(graph, mode = 'in') > 1))
  f = apply(combn(sort(leaves), 2), 2, paste, collapse = ' ')
  if(substitute) {
    avars = paste0(nam[1], seq_along(a))[seq_along(a)] %>% rlang::set_names(a)
    evars = paste0(nam[2], seq_along(e)) %>% rlang::set_names(e)
    fvars = paste0(nam[3], seq_len(choose(length(leaves), 2))) %>% rlang::set_names(f)
  } else {
    evars = e %>% paste0('`',. ,'`') %>% rlang::set_names(e)
    avars = (a %>% paste0('`',. ,'`'))[seq_along(a)] %>% rlang::set_names(a)
    fvars = f %>% paste0('`',. ,'`') %>% rlang::set_names(f)
  }
  coding = map2(list(avars, evars, fvars), nam, ~enframe(.x, 'edge', 'symbol') %>% mutate(type = .y)) %>% bind_rows

  eq = pp %>% filter(pop1 < pop2) %>%
    rowwise %>%
    mutate(mul = list(ifelse(lr == 0, avars[admix], paste0('(1 - ',avars[admix],')'))), add = list(paste0(evars[edges])),
           eq = list(paste0(paste0(mul, collapse = ' * '), ' * (', paste0(add, collapse = ' + '), ')') %>% str_replace_all('^ \\* ', '')),
           res = prod(ifelse(lr == 0, 1/4, 3/4)) * length(add),  resold = (1/4)^length(mul) * length(add),
           fvar = fvars[paste(pop1, pop2)]) %>%
    group_by(pop1, pop2) %>%
    summarize(equation = paste0(eq, collapse = ' + '), fvar = fvar[1], res = sum(res)) %>% ungroup %>%
    mutate(eq2 = paste0('(', equation, ')*2^10 - ', res*2^10)) %>% suppressMessages
  if(!return_everything) eq %<>% select(pop1, pop2, equation)

  list(equations = eq, coding = coding)
}

graph_to_function1 = function(graph, ge = NULL) {

  if(is.null(ge)) ge = graph_equations(graph)
  body1 = c('a = x[seq_len(na)]; e = x[(na+1):length(x)]; ')
  body2 = map_chr(ge$equations$equation, ~str_replace_all(., '([ae])([0-9]+)', '\\1\\[\\2\\]')) %>%
    paste(collapse = ', ') %>% paste('c(', ., ')')
  body = rlang::parse_expr(paste0('{', body1, body2, '}'))

  args = list(x = NULL, na = 0)
  eval(call("function", as.pairlist(args), body), env = parent.frame())
}

graph_to_function2 = function(graph, ge = NULL) {

  ge = graph_equations(graph)
  xx = map_chr(ge$equations$equation, ~str_replace_all(., '([ae])([0-9]+)', '\\1\\[\\2\\]'))
  Rcpp::cppFunction(paste0('NumericVector cppt(NumericVector x, int na) {x.push_front(0); NumericVector out(',length(ge$equations$equation),'); NumericVector a = x[Range(0,na)]; NumericVector e = x[Range(na-1,x.length())] = x[Range(na,x.length())]; ',paste0('out[', seq_along(xx)-1, '] = ', xx, collapse = '; '),'; return out;}'), plugins = 'cpp11')

}

graph_to_function3 = function(graph, ge = NULL) {

  if(is.null(ge)) ge = graph_equations(graph)
  na = numadmix(graph)
  ne = length(E(graph)) - 2*na

  pp = path_pairs(graph) %>% select(pop1, pop2, edges, admix, lr) %>% mutate(i = 1:n()) %>% filter(pop1 < pop2)
  pp2 = pp %>% select(-pop1, -pop2, -edges) %>% unnest(admix) %>% mutate(lr = pp %>% unnest(lr) %>% pull(lr))

  function(x, na) {

    a = x[seq_len(na)]
    e = x[(na+1):length(x)]

    adat = ge$coding %>% filter(type == 'a') %>% transmute(admix = edge, aval = x[seq_len(na)])
    edat = ge$coding %>% filter(type == 'e') %>% transmute(edges = edge, eval = x[(na+1):length(x)])
    admixvals = pp2 %>% left_join(adat, by = 'admix') %>%
      mutate(aval = ifelse(lr == 0, aval, 1-aval)) %>%
      group_by(i) %>% summarize(aval = prod(aval), .groups = 'drop')

    pp %>% left_join(admixvals, by = 'i') %>% mutate(aval = replace_na(aval, 1)) %>% unnest(edges) %>% left_join(edat, by = 'edges') %>% group_by(pop1, pop2) %>% summarize(f2pred = sum(aval * eval), .groups = 'drop') %>% ungroup %>% pull(f2pred)

  }
}

#' Make a function representing a graph
#'
#' This function takes an igraph object and turns it into a function that takes edge weights as input,
#' and outputs the expected f2-statistics.
#' @export
#' @param graph An admixture graph
#' @param admix_default The default weights for admixture edges
#' @param drift_default The default weights for drift edges
#' @param random_defaults Set default weights randomly for each edge between 0 and 1
#' @return A function mapping edge weights to f2-statistics
#' @examples
#' \dontrun{
#' mygraph = graph_f2_function(example_igraph)
#' mygraph(N3N8 = 0.1, `N2N1|Vindija.DG` = 0.4)
#' }
graph_f2_function = function(graph, admix_default = 0.5, drift_default = 0.01, random_defaults = FALSE) {

  ge = graph_equations(graph, substitute = FALSE)
  body1 = ge$equations$equation %>% paste(collapse = ', ') %>% paste('dat <- c(', ., '); ') %>% rlang::parse_expr()
  body2 = rlang::quo(d2 <- ge$coding %>% filter(type == 'f') %>% select(edge) %>% separate(edge, c('pop1', 'pop2'), ' '))
  body3 = rlang::expr(mutate(d2, f2 = dat))
  body = rlang::expr({!!body1; !!body2; !!body3})

  if(random_defaults) {
    na = ge$coding %>% filter(type == 'a') %>% nrow
    ne = ge$coding %>% filter(type == 'e') %>% nrow
    admix_default = runif(na)
    drift_default = runif(ne)
  }
  args = ge$coding %>% filter(type != 'f') %>%
    transmute(edge, val = ifelse(type == 'a', admix_default, drift_default)) %>%
    deframe %>% as.list
  rlang::eval_tidy(call("function", as.pairlist(args), body), env = parent.frame())
}



graph_jacobian = function(graph, ge = NULL, args = NULL) {

  if(is.null(ge)) ge = graph_equations(graph)
  if(is.null(args)) {
    na = numadmix(graph)
    ne = length(E(graph)) - 2*na
    args = c(runif(na), rep(1,ne))
  }
  symb = ge$coding %>% filter(type != 'f') %>% pull(symbol)

  jac = ge$equations$equation %>%
    map(~paste0('~', .) %>% as.formula %>% deriv(symb, symb) %>%
          exec(!!!args) %>% attr('gradient')) %>%
    do.call(rbind, .)
  rownames(jac) = paste(ge$equations$pop1, ge$equations$pop2)
  jac
}

graph_rank = function(graph) {

  jac = graph_jacobian(graph)
  qr(jac)$rank
}


#' Find all unidentifiable edges
#'
#' This function generates and evaluates admixture graphs in `numgen` iterations
#' to find well fitting admixturegraphs.
#' @export
#' @param graph An admixture graph
#' @return A data frame with all unidentifiable graph parameters
unidentifiable_edges = function(graph) {

  ge = graph_equations(graph)
  jac = graph_jacobian(graph, ge)

  dep = which(qr(jac)$rank == sapply(1:ncol(jac), function (i) qr(jac[,-i])$rank))
  adm = names(which(degree(graph, mode = 'in') > 1))
  parents = adm %>% rlang::set_names() %>% map(~names(neighbors(graph, ., mode = 'in')))
  out = ge$coding %>% slice(dep)
  oa = out %>% filter(type == 'a') %>% rename(to = edge) %>% rowwise %>%  mutate(from = parents[to]) %>%
    unnest(from) %>% select(from, to, type)
  oe = out %>% filter(type == 'e') %>% separate(edge, c('from', 'to'), '\\|') %>% select(-symbol)
  bind_rows(oa, oe) %>% mutate(type = ifelse(type == 'a', 'admix', 'edge'))
}

identifiable_sets = function(graph, jac, edge = NULL, n = 2) {

  rnk = qr(jac)$rank
  edges = colnames(jac)
  dep = which(rnk == sapply(1:ncol(jac), function (i) qr(jac[,-i])$rank))
  if(n == 1) return(t(t(edges[setdiff(1:ncol(jac), dep)])))
  if(n > length(dep)) return(NULL)
  if(is.null(edge)) cmb = combn(dep, n)
  else cmb = rbind(edge, combn(dep, n-1))
  #depn = which(sapply(1:ncol(cmb), function (i) qr(jac[,-cmb[,i]])$rank < rnk && all(map_lgl(seq_len(n), ~qr(jac[,-cmb[-.,i]])$rank == rnk))))
  depn = which(sapply(1:ncol(cmb), function (i) qr(jac[,-cmb[,i]])$rank < rnk))
  map(depn, ~edges[cmb[,.]]) %>% do.call('rbind', .)
}

identifiable_comb = function(graph, edge, jac = NULL, verbose = TRUE) {
  # returns the smallest sets in which an edge becomes identifiable

  stopifnot(edge %in% c(attr(E(graph), 'vnames'), names(V(graph))))
  if(is.null(jac)) jac = graph_jacobian(graph)
  edges = colnames(jac)

  for(i in seq_len(length(E(graph)))) {
    if(verbose) alert_info(paste0(i, '...\r'))
    us = identifiable_sets(graph, jac, edge = which(edges == edge), n = i)
    if(edge %in% unlist(us)) break
  }
  us[map_lgl(seq_len(nrow(us)), ~edge %in% us[.,]),,drop=FALSE]
}

predicted_f2 = function(graph, a = NULL, e = NULL) {

  ge = graph_equations(graph)
  fun = graph_to_function3(graph, ge)
  na = numadmix(graph)
  if(is.null(a)) a = runif(na)
  if(is.null(e)) e = rep(1e-2, length(E(graph))-2*length(a))
  ge$eq %>% select(pop1, pop2) %>% mutate(f2 = fun(c(a,e), na))
}

predicted_f4 = function(graph, a = NULL, e = NULL) {

  f2 = predicted_f2(graph, a = a, e = e) %>% bind_rows(rename(., pop1=pop2, pop2=pop1))
  pops = get_leafnames(graph)
  expand_grid(pop1 = pops, pop2 = pops, pop3 = pops, pop4 = pops) %>%
    filter(pop1 < pop2, pop1 < pop3, pop1 < pop4, pop3 < pop4, pop2 != pop3, pop2 != pop4) %>%
    left_join(f2 %>% rename(pop1 = pop1, pop3 = pop2, f13 = f2)) %>%
    left_join(f2 %>% rename(pop2 = pop1, pop4 = pop2, f24 = f2)) %>%
    left_join(f2 %>% rename(pop1 = pop1, pop4 = pop2, f14 = f2)) %>%
    left_join(f2 %>% rename(pop2 = pop1, pop3 = pop2, f23 = f2)) %>%
    mutate(f4 = (f14 + f23 - f13 - f24)/2) %>% ungroup %>% suppressMessages

}

graph_to_groebner = function(graph) {

  require(m2r)
  ge = graph_equations(graph, return_everything = TRUE)
  #vrs = eq %>% str_replace_all('[\\*\\+\\(\\)]', ' ') %>% str_replace_all('1-', '') %>% str_replace_all('-', '') %>%
  #  str_squish() %>% str_split(' ') %>% unlist %>% unique
  m2r::stop_m2()
  rng = m2r::ring_.(ge$coding$symbol, coefring = "QQ")

  out = m2r::gb_(ge$equations$eq2 %>% str_replace_all(' ', ''))
  m2r::stop_m2()
  out
}

find_invariants = function(graph, eps = 1e-7) {
  # find f4-statistics which are exactly 0

  pops = get_leafnames(graph)
  ge = graph_equations(graph)
  eq = ge$eq %>% transmute(p1 = pop1, p2 = pop2, eq = paste0('(', eq, ')')) %>%
    bind_rows(rename(., p1 = p2, p2 = p1))

  myenv = new.env()
  ge$coding %>% mutate(val = ifelse(type == 'e', 1, ifelse(type == 'a', runif(n()), NA))) %$%
    map2(symbol, val, ~assign(.x, .y, myenv)) %>% invisible

  expand_grid(pop1 = pops, pop2 = pops, pop3 = pops, pop4 = pops) %>%
    filter(pop1 < pop2, pop1 < pop3, pop1 < pop4, pop3 < pop4, pop2 != pop3, pop2 != pop4) %>%
    left_join(eq %>% rename(pop1 = p1, pop3 = p2, eq13 = eq)) %>%
    left_join(eq %>% rename(pop2 = p1, pop4 = p2, eq24 = eq)) %>%
    left_join(eq %>% rename(pop1 = p1, pop4 = p2, eq14 = eq)) %>%
    left_join(eq %>% rename(pop2 = p1, pop3 = p2, eq23 = eq)) %>%
    rowwise %>%
    mutate(f4 = paste(eq14, ' + ', eq23, ' - ', eq13, ' - ', eq24), res = eval(parse(text = f4), envir=myenv)) %>%
    ungroup %>%
    filter(between(res, -eps, eps)) %>% select(pop1:pop4) %>% suppressMessages
}


#' List clades in a graph
#'
#' @export
#' @param graph An admixture graph
#' @param eps Used to determine what counts as zero
#' @return A data frame with all (non-redundant) zero f4-statistics
#' @examples
#' \dontrun{
#' summarize_zerof4(example_igraph)
#' }
summarize_zerof4 = function(graph, eps = 1e-7) {

  pp = path_pairs(graph) %>% select(pop1, pop2, edges, admix, lr) %>% mutate(i = 1:n()) %>% filter(pop1 < pop2)
  pp2 = pp %>% select(-pop1, -pop2, -edges) %>% unnest(admix) %>% mutate(lr = pp %>% unnest(lr) %>% pull(lr))
  adm = names(which(degree(graph, mode = 'in') > 1))
  admval = runif(length(adm)) %>% set_names(adm)
  #ne = length(E(graph)) - 2*length(adm)
  driftval = runif(length(E(graph))) %>% set_names(attr(E(graph), 'vnames'))

  admixvals = pp2 %>%
    mutate(aval = ifelse(lr == 0, admval[admix], 1-admval[admix])) %>%
    group_by(i) %>% summarize(aval = prod(aval), .groups = 'drop')

  f2 = pp %>% left_join(admixvals, by = 'i') %>% mutate(aval = replace_na(aval, 1)) %>% unnest(edges) %>% mutate(eval = driftval[edges]) %>% group_by(pop1, pop2) %>% summarize(f2 = sum(aval * eval), .groups = 'drop') %>% ungroup %>% bind_rows(rename(., pop1 = pop2, pop2 = pop1))

  pops = get_leafnames(graph)
  expand_grid(pop1 = pops, pop2 = pops, pop3 = pops, pop4 = pops) %>%
    filter(pop1 < pop2, pop1 < pop3, pop1 < pop4, pop3 < pop4, pop2 != pop3, pop2 != pop4) %>%
    left_join(f2 %>% rename(pop1 = pop1, pop3 = pop2, f13 = f2)) %>%
    left_join(f2 %>% rename(pop2 = pop1, pop4 = pop2, f24 = f2)) %>%
    left_join(f2 %>% rename(pop1 = pop1, pop4 = pop2, f14 = f2)) %>%
    left_join(f2 %>% rename(pop2 = pop1, pop3 = pop2, f23 = f2)) %>%
    mutate(f4 = (f14 + f23 - f13 - f24)/2) %>% ungroup %>% suppressMessages %>%
    filter(between(f4, -eps, eps)) %>% select(pop1:pop4)
}

#' List clades in a list of graphs
#'
#' @export
#' @param graphlist A list of admixture graphs
#' @return A data frame with columns `pop1`, `pop2`, `pop3`, `pop4`, `n`, `frac`
#' @examples
#' \dontrun{
#' summarize_zerof4_list(graphlist)
#' }
summarize_zerof4_list = function(graphlist) {

  tot = length(graphlist)
  graphlist %>% map(summarize_zerof4) %>% bind_rows(.id = 'g') %>%
    count(pop1, pop2, pop3, pop4) %>% mutate(frac = n/tot)
}

normalize_zerof4 = function(zerof4) {
  zerof4[,1:4] %>% as_tibble(.name_repair = ~c('pop1', 'pop2', 'pop3', 'pop4')) %>%
    transmute(p1 = pmin(pop1, pop2),
              p2 = pmax(pop1, pop2),
              p3 = pmin(pop3, pop4),
              p4 = pmax(pop3, pop4)) %>%
    transmute(pop1 = if_else(p1 < p3, p1, p3),
              pop2 = if_else(p1 < p3, p2, p4),
              pop3 = if_else(p1 < p3, p3, p1),
              pop4 = if_else(p1 < p3, p4, p2)) %>% distinct
}

#' Test f4 constraints on a graph
#'
#' This function returns `TRUE` if and only if the admixture graph is compatible with
#' the f4-statistics listed in `nonzero_f4` being non-zero
#' @export
#' @param graph An admixture graph
#' @param nonzero_f4 A data frame or matrix with four columns. One row for each f4-statistic which is
#' expected to be non-zero
#' @return `TRUE` if all constraints are satisfied, else `FALSE`
#' @examples
#' \dontrun{
#' # Test whether f4(A,B; C,D) is expected to be non-zero in this graph:
#' nonzero_f4 = matrix(c('A', 'B', 'C', 'D'), 1)
#' satisfies_nonzerof4(random_admixturegraph(5, 2), nonzero_f4)
#' }
satisfies_nonzerof4 = function(graph, nonzero_f4) {

  if(is.null(nonzero_f4) || nrow(nonzero_f4) == 0) return(TRUE)
  unexpected_f4 = nonzero_f4 %>% normalize_zerof4 %>%
    inner_join(summarize_zerof4(graph)) %>% suppressMessages()
  nrow(unexpected_f4) == 0
}

#' Test f4 constraints on a graph
#'
#' This function returns `TRUE` if and only if the admixture graph is compatible with
#' the f4-statistics listed in `nonzero_f4` being non-zero
#' @export
#' @param graph An admixture graph
#' @param zero_f4 A data frame or matrix with four columns. One row for each f4-statistic which is
#' expected to be zero
#' @return `TRUE` if all constraints are satisfied, else `FALSE`
#' @examples
#' \dontrun{
#' # Test whether Chimp and Altai are a clade with respect to all populations X and Y:
#' # (whether f4("Chimp", "Altai"; X, Y) is 0 for all pairs of X and Y)
#' zero_f4 = expand_grid(pop1 = "Chimp", pop2 = "Altai", pop3 = X, pop4 = Y)
#' satisfies_zerof4(random_admixturegraph(5, 2), zero_f4)
#' }
satisfies_zerof4 = function(graph, zero_f4) {

  if(is.null(zero_f4) || nrow(zero_f4) == 0) return(TRUE)
  unexpected_f4 = zero_f4 %>% normalize_zerof4 %>%
    anti_join(summarize_zerof4(graph)) %>% suppressMessages()
  nrow(unexpected_f4) == 0
}



satisfies_oneevent = function(graph, earlier1, earlier2, later1, later2, type = 1) {

  root = get_rootname(graph)

  if(is.na(earlier2)) {
    earlier2 = names(neighbors(graph, earlier1, mode='in'))
    if(earlier2 == root) return(TRUE)
  }
  if(is.na(later2)) {
    later2 = names(neighbors(graph, later1, mode='in'))
    if(later2 == root) return(FALSE)
  }
  pe1 = all_simple_paths(graph, earlier1, root, mode = 'in')
  pe2 = all_simple_paths(graph, earlier2, root, mode = 'in')
  pl1 = all_simple_paths(graph, later1, root, mode = 'in')
  pl2 = all_simple_paths(graph, later2, root, mode = 'in')
  if(length(pe1) == 0 || length(pe2) == 0 || length(pl1) == 0 || length(pl2) == 0) stop('satisfies_oneevent error')

  ppis = expand_grid(expand_grid(pe1, pe2) %>% mutate(i = 1:n()),
                     expand_grid(pl1, pl2) %>% mutate(j = 1:n())) %>% rowwise %>%
    mutate(is1 = list(intersect(pe1, pe2)), is2 = list(intersect(pl1, pl2)),
           diff1 = list(setdiff(is1, is2)), diff2 = list(setdiff(is2, is1)),
           d1len = length(diff1), d2len = length(diff2), cont1 = is1[1] %in% is2) %>% ungroup %>%
    group_by(j) %>% mutate(maxdiff1len = max(d1len)) %>%
    group_by(i) %>% mutate(maxdiff2len = max(d2len)) %>% ungroup

  if(type == 1) {
    m10 = min(ppis$maxdiff1len) == 0
    m20 = min(ppis$maxdiff2len) == 0
    if(m10 && !m20) res = TRUE
    else if(m20 && !m10) res = FALSE
    else res = NA
  } else {
    res = all(ppis$cont1)
  }
  res
}

#' Test f4 constraints on a graph
#'
#' This function returns `TRUE` if and only if the admixture graph is compatible with
#' all event orders listed in eventorder
#' @export
#' @param graph An admixture graph
#' @param eventorder A data frame with columns `earlier1`, `earlier2`, `later1`, `later2`, optionally `type`
#' @param strict What to do in case some events are not determined by the graph.
#' If `strict = TRUE` (the default), the function will only return `TRUE` if there are no ambiguous constraints.
#' Otherwise, `TRUE` will be returned as long as no constraint is directly contradicted.
#' @return `TRUE` if all constraints are satisfied, else `FALSE`
#' @details Each row in `eventorder` represents a constraint that `earlier1` and `earlier2` split earlier than `later1` and `later2`. If `later2` is `NA`, `later2` will be set to the parent node of `later1`. By default (`type = 1`), a constraint will be satisfied as long as there is any lineage in which a split between `earlier1` and `earlier2` is ancestral to a split between `later1` and `later2`. `type = 2` is stricter and requires that the `earlier1`, `earlier2` split is ancestral to the `later1`, `later2` split in all lineages. In graphs with multiple admixture events there can be multiple splits between `earlier1`, `earlier2` and `later1`, `later2`, and many ways in which these splits can relate to each other. The current implementation only covers some of the many possible topological relationships.
#' @examples
#' \dontrun{
#' # Test whether the split between A and B is earlier than the split between C and D,
#' #   and whether the split between C and D is earlier than the terminal branch leading to E
#' constrain_events = tribble(
#'   ~earlier1, ~earlier2, ~later1, ~later2,
#'   'A', 'B', 'C', 'D',
#'   'C', 'D', 'E', NA)
#' satisfies_eventorder(random_admixturegraph(5, 0), eventorder = constrain_events)
#' }
satisfies_eventorder = function(graph, eventorder, strict = TRUE) {

  if(is.null(eventorder)) return(TRUE)
  if(!'type' %in% names(eventorder)) {
    eventorder$type = 1
  }
  status = eventorder %>% rowwise %>%
    mutate(ok = satisfies_oneevent(graph, earlier1, earlier2, later1, later2, type)) %>%
    ungroup %>% pull(ok)
  if(strict) return(isTRUE(all(status)))
  all(na.omit(status))
}

frac_eventorder = function(graph, eventorder, strict = TRUE) {

  if(is.null(eventorder)) return(TRUE)
  status = eventorder %>% rowwise %>%
    mutate(ok = satisfies_oneevent(graph, earlier1, earlier2, later1, later2)) %>%
    ungroup %>% pull(ok)
  if(strict) return(mean(!is.na(status)))
  mean(status, na.rm=T)
}


which_eventorder = function(graph, eventorder, strict = TRUE) {

  eventorder %>% rowwise %>%
    mutate(status = satisfies_oneevent(graph, earlier1, earlier2, later1, later2)) %>%
    ungroup
}



#' List number of admixture events for each population
#'
#' @export
#' @param graph An admixture graph
#' @return A data frame with columns `pop` and `nadmix`
#' @examples
#' \dontrun{
#' summarize_numadmix(example_igraph)
#' }
summarize_numadmix = function(graph) {
  # returns a data frame which lists the maximum number of admixture events for each population

  adm = names(which(degree(graph, mode = 'in') > 1))
  root = get_rootname(graph)
  graph %>% get_leafnames %>% set_names %>%
    map(~all_simple_paths(graph, root, ., mode = 'out') %>% map_dbl(~length(intersect(adm, names(.)))) %>% max) %>%
    unlist %>% enframe('pop', 'nadmix')
}

#' List number of admixture events for each population in a list of graphs
#'
#' @export
#' @param graphlist A list of admixture graphs
#' @return A data frame with columns `pop`, `mean`, `mean_nonzero`, `min`, `max`
#' @examples
#' \dontrun{
#' summarize_numadmix_list(graphlist)
#' }
summarize_numadmix_list = function(graphlist) {

  tot = length(graphlist)
  graphlist %>% map(summarize_numadmix) %>% bind_rows(.id = 'g') %>% group_by(pop) %>%
    summarize(mean = mean(nadmix), mean_nonzero = mean(nadmix != 0),
              min = min(nadmix), max = max(nadmix), .groups = 'drop')
}


#' Test admixture constraints on a graph
#'
#' This function returns `TRUE` if and only if the admixture graph satisfies all constraints on
#' the number of admixture events for the populations in `admix_constraints`
#' @export
#' @param graph An admixture graph
#' @param admix_constraints A data frame with columns `pop`, `min`, `max`
#' @return `TRUE` if all admixture constraints are satisfied, else `FALSE`
#' @examples
#' \dontrun{
#' # At least one admixture event for C, and none for D:
#' constrain_cd = tribble(
#'   ~pop, ~min, ~max,
#'   'C', 1, NA,
#'   'D', NA, 0)
#' satisfies_numadmix(random_admixturegraph(5, 2), constrain_cd)
#' }
satisfies_numadmix = function(graph, admix_constraints) {
  # admix_constraints is data frame with minimum and maximum admixture for each population

  if(is.null(admix_constraints)) return(TRUE)
  if(length(setdiff(c('min', 'max', 'pop'), names(admix_constraints))) > 0)
    stop("'admix_constraints' should have columns 'pop', 'min', 'max'!")
  unexpected_admix = admix_constraints %>%
    left_join(summarize_numadmix(graph), by = 'pop') %>%
    filter(nadmix < min | nadmix > max)
  nrow(unexpected_admix) == 0
}


#' Test constraints on a graph
#'
#' This function tests whether a graph satisfies a number of different types of constraints
#' @export
#' @param graph An admixture graph
#' @param nonzero_f4 A data frame or matrix with four columns. One row for each f4-statistic which is
#' observed to be significantly non-zero
#' @param admix_constraints A data frame with columns `pop`, `min`, `max`
#' @param event_order A data frame with columns `earlier1`, `earlier2`, `later1`, `later2`
#' @return `TRUE` if all admixture constraints are satisfied, else `FALSE`
#' @seealso \code{\link{satisfies_numadmix}}, \code{\link{satisfies_zerof4}}, \code{\link{satisfies_eventorder}}
#' \dontrun{
#' # At least one admixture event for C, and none for D:
#' constrain_cd = tibble(pop = c('C', 'D'), min = c(1, NA), max = c(NA, 0))
#' satisfies_numadmix(random_admixturegraph(5, 2), constrain_cd)
#' }
satisfies_constraints = function(graph, nonzero_f4 = NULL, admix_constraints = NULL, event_order = NULL) {

  satisfies_nonzerof4(graph, nonzero_f4) &&
    satisfies_numadmix(graph, admix_constraints) &&
    satisfies_eventorder(graph, event_order)
}


#' Test if a tree is part of a graph
#'
#' This function tests whether a tree is part of a graph. This is useful for testing whether a Y-chromosome tree is consistent with an autosomal admixture graph. Leaf node names matter, but internal node names are ignored.
#' @export
#' @param tree An admixture graph without admixture event
#' @param graph An admixture graph
#' @return `TRUE` if all admixture constraints are satisfied, else `FALSE`
#' \dontrun{
#' tree = graph_splittrees(example_igraph) %>% pull(graph) %>% pluck(1)
#' tree_in_graph(tree, example_igraph)
#' }
tree_in_graph = function(tree, graph) {
  # returns TRUE if tree is in graph

  treehash = graph_hash(tree)
  splithashes = graph_splittrees(graph) %>% rowwise %>%
    mutate(hash = graph_hash(graph)) %>% pull(hash)
  treehash %in% splithashes
}


is_clade = function(geno, pops1, pops2, i1=1, i2=1, perm=100) {

  g1 = geno[,pops1[-i1]]-geno[,pops1[i1]]
  g2 = geno[,pops2[-i2]]-geno[,pops2[i2]]
  #cancor(g1, g2)$Summary
  obs = cancor(g1, g2)$cor[1]
  expected = c()
  for(i in seq_len(perm)) {
    expected[i] = cancor(g1, g2[sample(seq_len(nrow(g2))),])$cor[1]
  }
  pmax(1/perm, mean(expected > obs))
}


condense_graph = function(graph) {
  # groups pops by shared admixture, return graph of those groups
  # should be redone in a nicer way

  leaves = get_leafnames(graph)
  adm = names(which(degree(graph, mode = 'in') > 1))
  root = get_rootname(graph)
  #dist = igraph::distances(graph, mode = 'out') %>% as_tibble(rownames = 'from') %>%
  #  pivot_longer(-from, names_to = 'to', values_to = 'dist') %>% filter(is.finite(dist), dist > 0)


  grps = delete_vertices(graph, adm) %>% components %$% membership %>%
    split(., .) %>% map(names) %>% map(~intersect(., c(leaves, adm))) %>% compact
  nam = map_chr(grps, ~shortest_paths(graph, .[1], root, mode = 'in') %$%
                  vpath %>% pluck(1) %>% names %>% intersect(c(adm, root)) %>% `[`(1))
  names(grps) = nam
  nam2 = map_chr(grps, ~paste(sort(.), collapse=' '))
  #sg = adm %>% intersect(nam) %>% set_names %>% map(~neighbors(graph, ., mode = 'in') %>% names %>% map_chr(~shortest_paths(graph, ., root, mode = 'in') %$% vpath %>% pluck(1) %>% names %>% intersect(nam) %>% `[`(1)) %>% unique) %>% enframe('to', 'from') %>% select(2:1) %>% unnest(from)
  # sg = adm %>% set_names %>% map(~neighbors(graph, ., mode = 'in') %>% names %>% map_chr(~shortest_paths(graph, ., root, mode = 'in') %$% vpath %>% pluck(1) %>% names %>% intersect(nam) %>% `[`(1)) %>% unique) %>% enframe('to0', 'from0') %>% select(2:1) %>% unnest(from0)

  # nam2 %<>% c(setdiff(adm, names(grps)) %>% set_names %>% map_chr(~shortest_paths(graph, ., leaves, mode = 'out')$vpath %>% suppressWarnings %>% compact %>% pluck(1) %>% names %>% intersect(nam2) %>% `[`(1)))
  # sg %>% mutate(from = nam2[from0], to = nam2[to0]) %>% distinct(from, to)

  sgi = distances(graph, nam, nam, mode = 'out') %>% igraph::graph_from_adjacency_matrix() %>% igraph::simplify()
  while(!all(names(V(sgi)) %in% nam)) {
    node = setdiff(names(V(sgi)), nam)[1]
    parents = sgi %>% neighbors(node, mode = 'in') %>% names
    children = sgi %>% neighbors(node, mode = 'out') %>% names
    newedges = expand_grid(a=parents, b=children) %>% as.matrix %>% t %>% c
    sgi %<>% delete_vertices(node) %>% add_edges(newedges)
  }
  sgi %>% as_edgelist() %>% as_tibble(.name_repair = ~c('from0', 'to0')) %>% mutate(from = nam2[from0], to = nam2[to0]) %>% distinct(from, to)
}

pair_admix = function(graph) {

  graph %>% condense_graph %>% edges_to_igraph() %>% igraph::distances(mode = 'out') %>%
    igraph::graph_from_adjacency_matrix() %>% igraph::simplify() %>%
    as_edgelist %>% as_tibble(.name_repair = ~c('from', 'to')) %>% rowwise %>%
    mutate(from = str_split(from, ' '), to = str_split(to, ' ')) %>% unnest(from) %>% unnest(to) %>%
    transmute(from, to, dir = 1) %>%
    bind_rows(rename(., from = to, to = from) %>% mutate(dir = -dir)) %>%
    complete(from, to, fill = list(dir = 0)) %>% filter(from < to) %>% arrange(from, to)
}


triplet_proportions = function(fit) {

  graph = fit %>% edges_to_igraph()
  triples = path_triples(graph)
  adm = fit %>% filter(type == 'admix') %>%
    transmute(e = paste0(from, '|', to), weight) %>% deframe
  prop = triples %>% rowwise %>%
    mutate(w1 = prod(adm[admedges1]),
           w2 = prod(adm[admedges2]),
           w3 = prod(adm[admedges3]), w = w1*w2*w3) %>%
    group_by(pop1, pop2, pop3, og) %>% summarize(sm = sum(w)) %>% ungroup %>% suppressMessages()
  x = prop %>% select(1:3) %>% distinct
  miss = bind_rows(x %>% mutate(og = pop1), x %>% mutate(og = pop2), x %>% mutate(og = pop3)) %>%
    mutate(sm = 0) %>% anti_join(prop, by = c('pop1', 'pop2', 'pop3', 'og'))
  prop %>% bind_rows(miss) %>% arrange(pop1, pop2, pop3) %>% rename(outgroup = og, proportion = sm)
}


tree_triplets = function(tree) {
  # returns the population triplets that make up a tree
  # currently too slow to run on many trees

  leaves = get_leafnames(tree)
  root = get_rootname(tree)
  igraph::distances(tree, V(tree), leaves, mode = 'out') %>% as_tibble(rownames = 'from') %>% pivot_longer(-from, names_to = 'to', values_to = 'dist') %>% filter(!from %in% c(leaves, root)) %>% group_by(from) %>% summarize(l = list(to[is.finite(dist)]), out = list(to[!is.finite(dist)])) %>% rowwise %>% mutate(l = list(t(combn(l, 2)) %>% as_tibble)) %>% unnest(l) %>% unnest(out) %>% transmute(out, X1 = pmin(V1, V2), X2 = pmax(V1, V2)) %>% distinct

}

tree_isoutgroup = function(tree, outgroup, pop1, pop2) {

  root = get_rootname(tree)
  p1 = shortest_paths(tree, root, pop1)$vpath[[1]]
  p2 = shortest_paths(tree, root, pop2)$vpath[[1]]
  og = shortest_paths(tree, root, outgroup)$vpath[[1]]
  length(intersect(p1, og)) < length(intersect(p1, p2))
}



place_root = function(graph, from, to, outpop = NULL) {
  # places root at edge from -> to
  # may reduce number of admixture edges

  root = get_rootname(graph)
  children = neighbors(graph, root, mode = 'out') %>% names
  oldleaves = get_leafnames(graph)
  newroot = newnodenam('root', names(V(graph)))
  newg = graph %>% delete_vertices(root) %>% add_edges(c(children[1], children[2])) %>%
    add_vertices(1, name = newroot) %>% delete_edges(paste0(from, '|', to)) %>%
    add_edges(c(newroot, from, newroot, to)) %>%
    igraph::as.undirected(mode = 'collapse')
  dist = igraph::distances(newg, newroot)
  out = newg %>% as_edgelist %>% as_tibble(.name_repair = ~c('v1', 'v2')) %>%
    rowwise %>%
    mutate(from = ifelse(dist[,v1] <= dist[,v2], v1, v2), to = ifelse(from == v1, v2, v1)) %>%
    ungroup %>% select(from, to) %>% edges_to_igraph()
  count = 0
  newleaves =  setdiff(get_leafnames(out), oldleaves)
  while(length(newleaves) > 0) {
    g = out %>% igraph::delete_vertices(newleaves) %>% simplify_graph()
    newleaves = setdiff(get_leafnames(g), oldleaves)
    out = g
  }
  out
}

#' Modify a graph by changing the position of the root node
#'
#' @export
#' @param graph An admixture graph
#' @param fix_outgroup Keep outgroup in place
#' @return A new admixture graph
place_root_random = function(graph, fix_outgroup = TRUE) {

  stopifnot(is_simplified(graph) && is_valid(graph))
  outpop = get_outpop(graph)
  if(fix_outgroup && !is.null(outpop)) {
    graph %<>% delete_vertices(c(outpop, get_rootname(graph)))
  }
  randedge = graph %>% delete_vertices(c(get_rootname(.), get_leafnames(.))) %>%
    E %>% attr('vnames') %>% sample(1) %>% str_split('\\|') %>% pluck(1)
  graph %<>% place_root(randedge[1], randedge[2])
  if(fix_outgroup && !is.null(outpop)) {
    newroot = newnodenam('root', names(V(graph)))
    oldroot = get_rootname(graph)
    graph %<>% add_vertices(2, name = c(newroot, outpop)) %>% add_edges(c(newroot, outpop, newroot, oldroot))
  }
  graph
}



adjlist_find_paths = function(a, n, m, path = c()) {
  # Find paths from node index n to m using adjacency list a.
  path = c(path, n)
  if(n == m) return(list(path))
  paths = list()
  for(child in a[[n]]) {
    child = as.numeric(child)
    if (!child %in% path) {
      child_paths = adjlist_find_paths(a, child, m, path)
      for(child_path in child_paths) {
        paths = c(paths, list(child_path))
      }
    }
  }
  paths
}

paths_from_to = function(graph, from, to) {
  # Find paths in graph from vertex source to vertex dest
  adj = get.adjlist(graph, mode = 'out')
  nam = names(V(graph))
  adjlist_find_paths(adj, which(nam == from), which(nam == to)) %>% map(~nam[.])
}

all_paths = function(graph) {
  # returns all paths from root to leaves; list with one list for each leaf node
  root = get_rootname(graph)
  leaves = get_leafnames(graph)
  leaves %>% set_names %>% map(~paths_from_to(graph, root, .))
}


pair_percentages = function(graph) {
  # returns how often each population ordered pair of populations is observed among the internal node descendants

  leaves = get_leafnames(graph)
  internal = setdiff(names(V(graph)), leaves)
  dest = graph %>% distances(internal, leaves, mode = 'out') %>% as_tibble(rownames = 'from') %>%
    pivot_longer(-from, names_to = 'to', values_to = 'order') %>% filter(is.finite(order)) %>% select(-order)
  single = dest %>% count(to)
  dest %>% full_join(dest, by = 'from') %>% count(to.x, to.y) %>%
    left_join(single %>% transmute(to.x = to, n1 = n)) %>%
    left_join(single %>% transmute(to.y = to, n2 = n)) %>%
    mutate(frac1 = n/n1, frac2 = n/n2) %>% arrange(to.x, to.y) %>% suppressMessages()
}

graph_distance = function(graph1, graph2) {

  x1 = pair_percentages(graph1)
  x2 = pair_percentages(graph2)
  mean((c(x1$frac1, x1$frac2) - c(x2$frac1, x2$frac2))^2)
}

#' Pairwise distance estimates for graphs
#'
#' Computes a distance estimate for each graph pair. Each graph is first summarized as a vector which counts for every leaf pair how many internal nodes reach that pair. The distance between two graphs is the Euclidean distance between the vectors of two graphs, and is scaled to fall between 0 and 1.
#' @param graphlist List of graphs
#' @return A data frame with graph distances
graph_distances = function(graphlist) {

  numgraphs = length(graphlist)
  pp = map(graphlist, pair_percentages) %>% map(~c(.$frac1, .$frac2))
  expand_grid(graph1 = seq_len(numgraphs), graph2 = seq_len(numgraphs)) %>% filter(graph1 < graph2) %>%
    rowwise %>% mutate(dist = mean((pp[[graph1]] - pp[[graph2]])^2)) %>% ungroup
}


consistent_with_qpadm = function(graph, left, right, target) {

  pops = get_leafnames(graph)
  internal = setdiff(names(V(graph)), pops)
  graph %>% distances(internal, pops, mode = 'out') %>% as_tibble(rownames = 'from') %>%
    pivot_longer(-from, names_to = 'to', values_to = 'order') %>% filter(is.finite(order)) %>% select(-order) %>% mutate(type = case_when(to %in% left ~ 'left', to %in% right ~ 'right', to %in% target ~ 'target')) %>% group_by(from) %>% filter('right' %in% type & 'target' %in% type & !'left' %in% type) %>% nrow %>% equals(0)

}

#' Partition a list of populations into left and right populations
#'
#' @export
#' @param pops A vector of populations
#' @param allpops Return only models which use all provided populations
#' @param more_right Return only models with more right than left populations
#' @return A data frame with one qpadm model per row
#' @seealso \code{\link{qpadm_models}}, \code{\link{graph_f2_function}}
#' @examples
#' \dontrun{
#' graph_to_qpadm(get_leafnames(example_igraph))
#' }
qpadm_models = function(pops, allpops = TRUE, more_right = TRUE) {

  models = tibble(l = power_set(pops))
  if(!is.logical(allpops)) stop('"allpops" should be TRUE or FALSE!')
  if(allpops) {
    models %<>% rowwise %>% mutate(r = list(setdiff(pops, l)))
  } else {
    models %<>% rowwise %>% mutate(r = list(power_set(setdiff(pops, l)))) %>% unnest(r)
  }
  models %<>% rowwise %>% filter(length(r) > 0) %>% ungroup
  if(more_right) models %<>% rowwise %>% filter(length(r) > length(l)) %>% ungroup
  models
}


#' Get all qpadm models for a graph
#'
#' This function tests which qpadm models should be valid for an admixture graph and a target population. By default, all models returned by \code{\link{qpadm_models}} are tested. For large graphs this will be too slow, and you may want test only some models by providing the `models` argument, or only a single model by providing the `left` and `right` arguments.
#' @details Two validity criteria are tested for each qpadm model: Rank validity and weight validity. Rank validity means that the rank of the f4 statistic matrix for only left and right populations is the same as the rank of the f4 statistic matrix that includes the target population among the left populations. Weight validity means that the estimated admixture weights for all left populations are between 0 and 1.
#' @export
#' @param graph An admixture graph
#' @param target Name of the target population.
#' @param left Left populations (provide this optionally if you want to test only a single qpadm model)
#' @param right Right populations (provide this optionally if you want to test only a single qpadm model)
#' @param models A two column nested data frame with models to be evaluated, one model per row. The first column, `l`, should contain the left populations, the second column, `r`, should contain the right populations. The target population is provided separately in the `target` argument.
#' @param weights Set this to `FALSE` to return only information on the ranks, not the weights, of each qpadm model. The ranks should depend only on the graph topology, while the weights and weight-validity (all weights for left populations between 0 and 1) can depend on the branch lengths of the graph. By default f4-statistics are based on equal branch lengths and admixture weights of 0.5. This can be overridden by providing `f4dat`.
#' @param f4dat A data frame of f4-statistics which can be provided to override the default branch lengths.
#' @param allpops Evaluate only models which use all populations in the admixture graph. See \code{\link{qpadm_models}}
#' @param return_f4 Include f4 statistic matrices in the results (default `FALSE`)
#' @param eps Epsilon value close to zero which is used for determining which f4 matrix elements should be considered non-zero, and which weights are strictly between 0 and 1.
#' @return A data frame with one qpadm model per row and columns `valid_rank` and `valid_weights` indicating whether a model should be valid under the graph.
#' @note An earlier version of this function tried to use the graph topology for identifying valid qpadm models, but this didn't work reliably. Christian Huber had the idea of using the ranks of expected f4 statistic matrices instead.
#' @seealso \code{\link{qpadm_models}}, \code{\link{graph_f2_function}}
#' @examples
#' \dontrun{
#' graph2 = example_igraph %>% simplify_graph() %>%
#'   delete_admix('N3N0', 'N3N4') %>% delete_admix('N3N1', 'N3N8')
#' graph_to_qpadm(graph2, 'Mbuti.DG') %>% filter(valid_rank, valid_weights)
#' }
graph_to_qpadm = function(graph, target, left = NULL, right = NULL, models = NULL,
                          weights = TRUE, f4dat = NULL, allpops = TRUE, more_right = TRUE, return_f4 = FALSE, eps = 1e-10) {

  if(is.null(models)) {
    if(!is.null(left) && !is.null(right)) {
      models = tibble(l = list(left), r = list(right))
    } else {
      models = get_leafnames(graph) %>% setdiff(target) %>%
        qpadm_models(allpops = allpops, more_right = more_right)
    }
  }
  if(is.null(f4dat)) f4dat = graph_f2_function(graph, random_defaults=FALSE)() %>% f2dat_f4dat()
  out = models %>% rowwise %>%
    mutate(target,
           m1 = list(f4dat_f4mat(f4dat, l, r, eps = eps)),
           m2 = list(f4dat_f4mat(f4dat, c(target, l), r, eps = eps)),
           rank1 = qr(m1)$rank,
           rank2 = qr(m2)$rank,
           fullrank = length(l)-1 == rank1,
           stablerank = rank1 == rank2)
  if(weights) {
    out %<>%
      mutate(w = list(cpp_qpadm_weights(m2, diag(prod(dim(m2))), rnk = min(dim(m2)), qpsolve=qpsolve)$weights),
             wmin = min(w), wmax = max(w), targetadmixed = isTRUE(wmin > 0+eps & wmax < 1-eps),
             valid = fullrank & stablerank & targetadmixed) #%>% select(-w)
  }
  if(!return_f4) out %<>% select(-m1, -m2)
  out %>% ungroup
}

qpadm_to_graphs = function(target, left, right, nadmix=0, ntry = 100) {

  pops = c(target, left, right)
  graphs = rerun(ntry, random_admixturegraph(as.numeric(length(pops)), nadmix))
  map(graphs, ~graph_to_qpadm(.x, target, left, right)) %>%
    bind_rows() %>% mutate(graph=graphs) %>% select(-c(l,r,target))
}


f4dat_f4mat = function(f4dat, left, right, eps = NA) {

  mat = f4dat %>%
    filter(pop1 == left[1], pop2 %in% left[-1], pop3 == right[1], pop4 %in% right[-1]) %>%
    select(-pop1,-pop3) %>% pivot_wider(pop2, names_from='pop4', values_from='f4') %>%
    select(-pop2) %>% as.matrix
  mat[abs(mat) < eps] = 0
  mat
}

consistent_with_qpwave = function(graph, left, right) {

  edges = graph %>% delete_vertices(pops) %>% {attr(E(.), 'vnames')}

  out = tibble()
  for(i in 1:3) {
    ecomb = t(combn(edges, i))
    new = map(seq_len(nrow(ecomb)), ~{
      comp = components(delete_edges(graph, ecomb[.,]))
      if(comp$no == 2 && min(comp$csize) > 1) comp$membership %>% enframe %>%
        filter(name %in% pops) %$% split(name, value) %>% enframe %>% rowwise %>% mutate(len = length(value))
    }) %>% bind_rows(.id = 'comp') %>% group_by(comp) %>% filter(min(len) > 1) %>% select(-len) %>% pivot_wider(names_from = name, values_from = value) %>% ungroup %>% transmute(num = i, l = `1`, r = `2`)
    out %<>% bind_rows(new)
  }

}

#' Simulate allele frequncies under an admixture graph
#'
#' This function performs a very crude simulation of allele frequencies
#' under an admixture graph model
#' @export
#' @param graph An admixture graph as igraph object, or as edge list data frame
#' with a column `weight`, as returned by `qpgraph()$edges`
#' @param nsnps Number of SNPs to simulate
#' @param drift_default Default branch lengths. Ignored if `graph` is a data frame with weights
#' @param admix_default Default admixture weights. Ignored if `graph` is a data frame with weights
#' @param leaves_only Return allele frequencies for leaf nodes only
#' @return A data frame with allele frequencies
#' @seealso \code{\link{graph_to_pcs}}
#' @examples
#' \dontrun{
#' afs = graph_to_afs(example_igraph)
#' }
graph_to_afs = function(graph, nsnps = 1e4, drift_default=0.02, admix_default=0.5, leaves_only = FALSE) {

  if(is.data.frame(graph)) {
    weights = graph
    graph = graph %>% select(from, to) %>% edges_to_igraph
  } else {
    weights = graph %>% igraph::as_edgelist() %>%
      as_tibble(.name_repair = ~c('from', 'to')) %>% add_count(to) %>%
      mutate(weight = ifelse(n > 1, admix_default, drift_default)) %>% select(-n)
  }

  rootaf = runif(nsnps)
  nodes = names(igraph::V(graph))
  afs = list()
  aw = function(graph, node) {
    parents = names(igraph::neighbors(graph, node, mode = 'in'))
    children = names(igraph::neighbors(graph, node, mode = 'out'))
    isadmix = length(parents) > 1
    if(!isadmix) {
      if(length(parents) == 0) {
        afs[[node]] <<- rootaf
      } else {
        w = weights %>% filter(from == parents, to == node) %>% pluck('weight', 1)
        afs[[node]] <<- pmin(1,pmax(0,afs[[parents]] + w*runif(nsnps)))
      }
      for(child in children) {
        aw(graph, child)
      }
    } else {
      if(!is.null(afs[[parents[1]]]) && !is.null(afs[[parents[2]]])) {
        w1 = weights %>% filter(from == parents[1], to == node) %>% pluck('weight', 1)
        w2 = weights %>% filter(from == parents[2], to == node) %>% pluck('weight', 1)
        afs[[node]] <<- afs[[parents[1]]]*w1 + afs[[parents[2]]]*w2
        for(child in children) {
          aw(graph, child)
        }
      }
    }
  }
  aw(graph, get_rootname(graph))
  out = afs %>% bind_cols
  if(leaves_only) out = out[get_leafnames(graph)]
  out
}

#' Simulate PCs under an admixture graph
#'
#' This function simulates PCA of allele frequencies under an admixture graph model
#' @export
#' @param graph An admixture graph as igraph object, or as edge list data frame
#' with a column `weight`, as returned by `qpgraph()$edges`
#' @param nsnps Number of SNPs to simulate
#' @param drift_default Default branch lengths. Ignored if `graph` is a data frame with weights
#' @param admix_default Default admixture weights. Ignored if `graph` is a data frame with weights
#' @param leaves_only Return PCs for leaf nodes only
#' @return A data frame with PCs for each population
#' @seealso \code{\link{graph_to_afs}}
#' @examples
#' \dontrun{
#' pcs = graph_to_pcs(example_igraph)
#' pcs %>% ggplot(aes(PC1, PC2, label = pop)) + geom_text() + geom_point()
#' }
graph_to_pcs = function(graph, nsnps = 1e4, drift_default=0.02, admix_default=0.5, leaves_only = TRUE) {

  afs = graph_to_afs(graph, nsnps = nsnps, drift_default = drift_default,
                     admix_default = admix_default, leaves_only = leaves_only)
  prcomp(t(as.matrix(afs)))$x %>% as_tibble(rownames = 'pop')
}


#' List leaf nodes for all internal nodes
#'
#' @export
#' @param graph An admixture graphs
#' @return A data frame with columns `from`, `to`, `id`
#' @examples
#' \dontrun{
#' summarize_descendants(graph)
#' }
summarize_descendants = function(graph) {

  leaves = get_leafnames(graph)
  internal = setdiff(names(V(graph)), leaves)
  graph %>% distances(internal, leaves, mode = 'out') %>% as_tibble(rownames = 'from') %>%
    pivot_longer(-from, names_to = 'to', values_to = 'order') %>% filter(is.finite(order)) %>%
    select(-order) %>% group_by(from) %>% mutate(id = paste(sort(to), collapse = ' ')) %>% ungroup
}

#' List leaf nodes for all internal nodes in a list of graphs
#'
#' @export
#' @param graphlist A list of admixture graphs
#' @param rename If `FALSE` (default), the output will be a data frame indicating how often each node in each graph is observed in all other graphs. If `TRUE`, the output will be a list, where the inner nodes will be renamed to have these percentages as part of their name. \code{\link{plot_graph}} will print the percentages of graphs renamed in this way.
#' @return A data frame with columns `graph`, `from`, `n`, `frac`
#' @examples
#' \dontrun{
#' summarize_descendants_list(graphlist)
#' }
summarize_descendants_list = function(graphlist, rename = FALSE) {

  ngraphs = length(graphlist)
  dat = graphlist %>% map(summarize_descendants) %>% bind_rows(.id = 'graph') %>% select(-to) %>% distinct %>%
    mutate(graph = as.numeric(graph))
  counts = dat %>% select(-from) %>% distinct %>% count(id)
  out = dat %>% left_join(counts, by = 'id') %>% mutate(frac = n/ngraphs)
  if(rename) {
    pops = get_leafnames(graphlist[[1]])
    out = map2(graphlist, split(out, out$graph), ~{
    nam = .y %>% transmute(from, to = paste0(from, ' (', round(frac*100), ')')) %>%
      bind_rows(tibble(from = pops, to = pops)) %>% deframe
    .x %>% as_edgelist() %>% as_tibble(.name_repair = ~c('from', 'to')) %>%
      mutate(from = nam[from], to = nam[to]) %>% edges_to_igraph
  })
  }
  out
}


partition_fit = function(fit) {

  leaves = get_leafnames(graph)
  internal = setdiff(names(V(graph)), leaves)
  dest = graph %>% igraph::distances(internal, leaves, mode = 'out') %>% as_tibble(rownames = 'from') %>%
    pivot_longer(-from, names_to = 'to', values_to = 'order') %>% filter(is.finite(order)) %>% select(-order)
}


summarize_driftpartition = function(edges, all = FALSE) {

  graph = edges %>% edges_to_igraph()
  leaves = get_leafnames(graph)
  internal = setdiff(names(V(graph)), leaves)
  desc = graph %>% igraph::distances(internal, leaves, mode = 'out') %>% as_tibble(rownames = 'from') %>%
    pivot_longer(-from, names_to = 'to', values_to = 'order') %>% filter(is.finite(order)) %>%
    select(-order) %>% group_by(from) %>% mutate(id = paste(sort(to), collapse = ' ')) %>% ungroup

  out = edges %>% filter(type == 'edge') %>% select(from, to, weight) %>% left_join(desc %>% select(-to) %>% distinct, by = c('to' = 'from')) %>% mutate(id = ifelse(is.na(id), to, id)) %>% rowwise %>% mutate(id2 = paste0(sort(setdiff(leaves, str_split(id, ' ')[[1]])), collapse = ' ')) %>% ungroup %>% mutate(s1 = pmin(id, id2), s2 = pmax(id, id2)) %>% group_by(s1, s2) %>% summarize(weight = sum(weight), .groups = 'drop') %>% bind_rows(rename(., s1=s2, s2=s1)) %>% rowwise %>% transmute(s1, s2, c1 = str_count(s1, ' ')+1, c2 = str_count(s2, ' ')+1, weight) %>% ungroup
  if(all) {
    x = tibble(s1 = head(power_set(sort(leaves)), -1)) %>% rowwise %>%
      mutate(s2 = list(sort(setdiff(leaves, s1))),
             c1 = length(s1), c2 = length(s2),
             s1 = paste0(s1, collapse = ' '),
             s2 = paste0(s2, collapse = ' '), weight = 0) %>% ungroup
    out = x %>% anti_join(out, by = c('s1', 's2')) %>% bind_rows(out)
  }
  out %>% arrange(s1, s2)
}




graph_boot_score = function(bootfit, graph = NULL, f2_blocks = NULL) {

  if(is.null(bootfit)) bootfit = qpgraph_resample_snps(f2_blocks, boot = 100, graph = graph, numstart = 5, return_fstats = TRUE)
  if(!'f4' %in% names(bootfit)) stop("Need to run with `return_fstats = TRUE`!")
  x = bootfit$f4 %>% bind_rows(.id = 'id') %>% mutate(id = as.numeric(id)) %>%
    filter(pop1 == pop3, pop2 == pop4) %>% select(-pop3, -pop4) %>% arrange(id, pop1, pop2) %>%
    transmute(id, p = paste(pop1, pop2), diff)
  means = x %>% group_by(p) %>% summarize(mn = mean(diff), var = var(diff), .groups = 'drop')

  covmat = x %>% left_join(x, by = 'id') %>%
    group_by(p.x, p.y) %>% summarize(cv = cov(diff.x, diff.y), .groups = 'drop') %>%
    pivot_wider(p.x, names_from = p.y, values_from = cv) %>% select(-p.x) %>% as.matrix

  stat = means$mn %*% (pracma::pinv(covmat)) %*% means$mn
  stat[1,1]

}


graph_boot_pval = function(bootfit) {

  stat = graph_boot_score(bootfit)
  rnk = bootfit$edges[[1]] %>% edges_to_igraph() %>% graph_rank()
  pchisq(stat, rnk, lower.tail = F)
}

#' Convert agraph to igraph
#'
#' `agraph` is the format used by the `admixturegraph` packge. `igraph` is used by the `admixtools` package
#' @export
#' @param agraph An admixture graph in \code{\link{agraph}} format
#' @return An admixture graph in \code{\link{igraph}} format
#' @examples
#' \dontrun{
#' agraph_to_igraph(agraph)
#' }
agraph_to_igraph = function(agraph) {
  igraph::graph_from_adjacency_matrix(agraph$children)
}

#' Convert igraph to agraph
#'
#' `agraph` is the format used by the `admixturegraph` packge. `igraph` is used by the `admixtools` package
#' @export
#' @param agraph An admixture graph in \code{\link{igraph}} format
#' @return An admixture graph in \code{\link{agraph}} format
#' @examples
#' \dontrun{
#' igraph_to_agraph(example_igraph)
#' }
igraph_to_agraph = function(igraph) {

  leaves = get_leafnames(igraph)
  inner_nodes = setdiff(names(V(igraph)), leaves)
  e = igraph %>% as_edgelist() %>% `[`(,c(2,1)) %>% cbind(NA)
  admixturegraph::agraph(leaves, inner_nodes, e)
}



rename_internal = function(graph) {

  isedge = FALSE
  if(!'igraph' %in% class(graph)) {
    isedge = TRUE
    edges = graph
    graph %<>% edges_to_igraph()
  }

  leaves = get_leafnames(graph)
  nammap = graph %>% distances(names(V(graph)), leaves, mode = 'out') %>% as_tibble(rownames = 'from') %>%
    pivot_longer(-from, names_to = 'to', values_to = 'order') %>% filter(is.finite(order)) %>%
    select(-order) %>% group_by(from) %>% mutate(id = paste(sort(to), collapse = '\n')) %>% ungroup %>%
    select(-to) %>% distinct %>% group_by(id) %>%
    mutate(n = 1:n(), id = ifelse(n == n(), id, paste(id, n, sep = '_'))) %>% select(-n) %>% deframe
  if(isedge) {
    edges %<>% mutate(from = nammap[from], to = nammap[to])
    return(edges)
  }
  graph %>% set_vertex_attr('name', value = nammap[names(V(graph))])
}


label_internal = function(edges) {

  # assigns leaf nodes to internal nodes, if the internal node contributes more to one leaf than to all others

  graph = edges %>% edges_to_igraph()
  leaves = get_leafnames(graph)
  internal = setdiff(names(V(graph)), leaves)
  wvec = edges %>% mutate(weight = ifelse(type == 'admix', weight, 1)) %>% transmute(e = paste(from, to, sep = '|'), weight) %>% deframe

  iweights = imap_dfr(internal, ~{
    node = .x
    igraph::all_simple_paths(graph, node, leaves) %>%
      set_names(map_chr(., ~attr(tail(., 1), 'name'))) %>%
      map(admixtools:::vs_to_es) %>% map(~wvec[.]) %>% map(prod) %>%
      enframe %>% unnest(value) %>% group_by(name) %>%
      summarize(weight = sum(value)) %>% transmute(from = node, to = name, weight)
  })

  int_labels = iweights %>% complete(to, from, fill = list(weight = 0)) %>% group_by(from) %>% arrange(-weight) %>% slice_head(n = 2) %>% mutate(diff = weight[1]-weight[2]) %>% filter(diff > 0) %>% slice_head(n=1) %>% ungroup

}

leaf_internal_ancestry = function(edges) {

  # returns data frame with all leaf nodes with ancestry that has been mapped to
  # another leaf node by 'label_internal'

  graph = edges %>% edges_to_igraph()
  leaves = get_leafnames(graph)
  internal = setdiff(names(V(graph)), leaves)
  root = get_rootname(graph)

  int_labels = label_internal(edges)
  map_dfr(leaves, ~{
    node = .
    x = igraph::all_simple_paths(graph, root, node) %>% unlist %>% names %>% unique
    from = int_labels %>% filter(to != node, from %in% x) %>% pull(to) %>% unique
    tibble(from = from, to = node)
  })
}

#' Count zero-length edges
#'
#' `agraph` is the format used by the `admixturegraph` packge. `igraph` is used by the `admixtools` package
#' @param edges Edges data frame from fitted admixture graph
#' @param epsilon Every edge with length
#' @return The number of edges with length < epsilon
#' @examples
#' \dontrun{
#' fit = qpgraph(example_f2_blocks, example_igraph)
#' count_zero_edges(fit$edges)
#' }
count_zero_edges = function(edges, epsilon = 1e-6) {

  if(!'type' %in% names(edges)) stop("'edges' should be a data frame with a column named 'type'!")
  if(!'weight' %in% names(edges)) stop("'edges' should be a data frame with a column named 'weight'!")

  edges %>%
    filter(weight < epsilon, type == 'edge') %>%
    nrow()
}



