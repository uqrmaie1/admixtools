# Contract tests for three byte identical optimizations in the qpgraph scoring
# path (pull requests 163, 164, 165). Each optimization is meant to leave
# observable output unchanged, so each test below pins the invariant that makes
# it safe. A future change that breaks the invariant fails here rather than
# silently perturbing find_graphs and qpgraph results.
#
#  163  graph_to_pwts and graph_to_weightind resolve each root to leaf path to
#       its edge ids by passing the integer vertex sequence to get_edge_ids
#       instead of vertex names. Locked two ways. An independent oracle
#       recomputes the path edge ids from the unnamed edge list, and a graph
#       level invariance check rebuilds each graph under a different vertex id
#       assignment and requires the path weight matrix to be unchanged.
#  164  inout_degree derives both degrees from the unnamed edge list. It must
#       equal igraph::degree for every vertex in both directions.
#  165  qpgraph computes the root and the admixture edge vector once and shares
#       them with both helpers through optional arguments. Supplying them must
#       equal recomputing them.

suppressMessages(library(igraph))

# A battery of random admixture graphs spanning trees and a few admixture
# counts. random_admixturegraph returns an igraph object directly. Seeds are
# fixed so the battery is reproducible.
make_graph_battery <- function() {
  grid <- expand.grid(nleaf = c(4, 6, 8), nadmix = c(0, 1, 2), seed = 1:4,
                      KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(grid)), function(i) {
    set.seed(grid$seed[i] * 1000 + grid$nleaf[i] * 10 + grid$nadmix[i])
    g <- suppressWarnings(random_admixturegraph(grid$nleaf[i], numadmix = grid$nadmix[i]))
    list(g = g, label = sprintf("nleaf=%d nadmix=%d seed=%d",
                                grid$nleaf[i], grid$nadmix[i], grid$seed[i]))
  })
}

graph_battery <- make_graph_battery()

# The admixture edge vector qpgraph builds and shares with the two helpers.
# Mirror its derivation exactly so the 165 fallback test feeds the helpers what
# qpgraph would feed them.
shared_admixedges <- function(g) {
  admixnodes <- which(igraph::degree(g, mode = 'in') == 2)
  unlist(igraph::incident_edges(g, admixnodes, mode = 'in'))
}

# ---- 164  inout_degree equals igraph::degree -------------------------------

test_that("inout_degree equals igraph::degree in both directions", {
  for (case in graph_battery) {
    io <- admixtools:::inout_degree(case$g)
    # value identical to igraph for every vertex. inout_degree intentionally
    # returns integer where degree returns double, so compare by value.
    expect_equal(io$indeg,  unname(igraph::degree(case$g, mode = 'in')),
                 info = case$label)
    expect_equal(io$outdeg, unname(igraph::degree(case$g, mode = 'out')),
                 info = case$label)
  }
})

test_that("inout_degree returns integer vectors and handles boundary graphs", {
  # the integer return type is part of the contract. Callers consume it only
  # through == and sum, where 1L and 1.0 agree.
  io <- admixtools:::inout_degree(graph_battery[[1]]$g)
  expect_type(io$indeg, "integer")
  expect_type(io$outdeg, "integer")

  # explicit small graph, plus the canonical example fixture
  g <- igraph::graph_from_edgelist(
    matrix(c('R','a', 'R','b', 'a','x', 'a','y'), ncol = 2, byrow = TRUE),
    directed = TRUE)
  io <- admixtools:::inout_degree(g)
  expect_equal(io$indeg,  unname(igraph::degree(g, mode = 'in')))
  expect_equal(io$outdeg, unname(igraph::degree(g, mode = 'out')))

  io2 <- admixtools:::inout_degree(example_igraph)
  expect_equal(io2$indeg,  unname(igraph::degree(example_igraph, mode = 'in')))
  expect_equal(io2$outdeg, unname(igraph::degree(example_igraph, mode = 'out')))
})

# ---- 165  shared root and admixedges equal recomputing them ----------------

test_that("graph_to_pwts is unchanged when root and admixedges are supplied", {
  for (case in graph_battery) {
    g <- case$g
    pops <- get_leafnames(g)                 # exactly what qpgraph passes
    root <- admixtools:::get_root(g)
    ae <- shared_admixedges(g)
    expect_identical(
      admixtools:::graph_to_pwts(g, pops),
      admixtools:::graph_to_pwts(g, pops, root = root, admixedges = ae),
      info = case$label)
  }
})

test_that("graph_to_weightind is unchanged when root and admixedges are supplied", {
  checked <- 0
  for (case in graph_battery) {
    g <- case$g
    if (numadmix(g) == 0) next               # qpgraph calls weightind only when nadmix > 0
    root <- admixtools:::get_root(g)
    ae <- shared_admixedges(g)
    expect_identical(
      admixtools:::graph_to_weightind(g),
      admixtools:::graph_to_weightind(g, root = root, admixedges = ae),
      info = case$label)
    checked <- checked + 1
  }
  # guard against the battery silently degrading to all nadmix == 0 graphs,
  # which would make every iteration skip and pass this test vacuously.
  expect_gt(checked, 0)
})

# ---- 163  path edge resolution and id assignment invariance ----------------

test_that("path edge id resolution matches an independent oracle", {
  # Recompute the edge ids of a path straight from the unnamed edge list, with
  # no call to get_edge_ids, and require the helper's resolution to match.
  oracle_edges <- function(g, vids) {
    elf <- igraph::as_edgelist(g, names = FALSE)
    from <- vids[-length(vids)]; to <- vids[-1]
    vapply(seq_along(from),
           function(k) which(elf[, 1] == from[k] & elf[, 2] == to[k]),
           integer(1))
  }
  checked <- 0
  for (case in graph_battery) {
    g <- case$g
    root <- admixtools:::get_root(g)
    leaves <- get_leafnames(g)
    for (p in igraph::all_simple_paths(g, root, leaves, mode = 'out')) {
      via_helper <- as.vector(igraph::get_edge_ids(g, as.numeric(admixtools:::expand_path(p))))
      via_oracle <- oracle_edges(g, as.numeric(p))
      # get_edge_ids returns doubles, the oracle's which() returns integers;
      # the edge id values are what the helpers consume, so compare by value.
      expect_equal(via_helper, via_oracle, info = case$label)
      checked <- checked + 1
    }
  }
  expect_gt(checked, 0)
})

test_that("graph_to_pwts is invariant to vertex id assignment", {
  # graph_from_edgelist numbers vertices by first appearance, so reversing the
  # edge rows yields the same topology and the same names under a different
  # vertex id assignment. A correct id based resolution gives the same path
  # weight matrix once rows are aligned by edge label.
  pwts_sorted <- function(g, pops) {
    m <- admixtools:::graph_to_pwts(g, pops)
    m[order(rownames(m)), , drop = FALSE]
  }
  saw_relabel <- FALSE
  for (case in graph_battery) {
    g1 <- case$g
    el <- igraph::as_edgelist(g1)
    g2 <- igraph::graph_from_edgelist(el[nrow(el):1, , drop = FALSE], directed = TRUE)
    saw_relabel <- saw_relabel ||
      !identical(igraph::vertex_attr(g1, 'name'), igraph::vertex_attr(g2, 'name'))
    pops <- sort(get_leafnames(g1))          # same name set and order for both graphs
    expect_equal(pwts_sorted(g1, pops), pwts_sorted(g2, pops), info = case$label)
  }
  # confirm the battery actually exercised a changed id assignment
  expect_true(saw_relabel)
})
