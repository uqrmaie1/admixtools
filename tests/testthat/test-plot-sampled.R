# Data-layer tests for plot_graph's two-state node encoding (sampled vs
# unsampled), driven through the shared layout helper graph_to_plotdat().
# No vdiffr; assert the partition the plot is built from, not pixels.

test_that("graph_to_plotdat styles an internal-sampled node as sampled", {
  g  <- make_test_nodes_graph()        # anc internal-sampled; modern internal-unsampled
  pd <- graph_to_plotdat(g)
  expect_true("anc" %in% pd$nodes$name)                                   # sampled glyph set
  expect_false(!is.null(pd$internal) && "anc" %in% pd$internal$name)      # not in internal set
})

test_that("graph_to_plotdat unchanged for a leaf-only graph", {
  g  <- make_minimal_graph()           # no nodes attr -> sampled_set == leaves
  pd <- graph_to_plotdat(g)
  expect_setequal(pd$nodes$name, get_leafnames(edges_to_igraph(g)))
})
