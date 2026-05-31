# Tests for add_sampled_tips(), the tip-materialisation transform that makes
# a graph with sampled internal nodes fittable, and the qpgraph guard that
# refuses an un-transformed sampled-internal-node graph.
#
# These are pure-topology assertions: no f2 data needed, so they run in CI.
# The real-data FIT (that the transformed graph fits non-degenerately) is
# validated separately by the dogfood script against the 14-pop f2 cache.

test_that("add_sampled_tips makes a sampled internal node a leaf tip", {
  g <- make_test_nodes_graph()                 # anc internal, samples=1
  gt <- add_sampled_tips(g)
  expect_true("anc" %in% get_leafnames(edges_to_igraph(gt)))
  expect_true("anc_anc" %in% c(gt$from, gt$to))   # new ancestor created
  # anc's former children now descend from anc_anc, not anc
  expect_true(all(gt$from[gt$to %in% c("v", "modern")] != "anc"))
})

test_that("add_sampled_tips is a no-op when no sampled internal nodes", {
  g <- make_minimal_graph()                    # leaves only, no nodes attr
  gt <- add_sampled_tips(g)
  expect_setequal(c(gt$from, gt$to), c(g$from, g$to))
})

test_that("add_sampled_tips errors on a name collision", {
  g <- make_test_nodes_graph()
  g <- dplyr::bind_rows(g, tibble::tibble(from = "anc_anc", to = "x",
                                          type = "normal", weight = NA_real_))
  g <- set_node_attrs(g, "anc", samples = 1L)
  expect_error(add_sampled_tips(g), class = "admixtools_invalid_graph")
})

# --- Structural breadth: the shapes the dogfood's single fit does not cover.

test_that("add_sampled_tips handles MULTIPLE sampled internal nodes", {
  g <- tibble::tribble(
    ~from, ~to,
    "R", "out", "R", "P", "P", "Q", "Q", "leaf1", "Q", "leaf2", "P", "leaf3"
  )
  # P and Q are both internal; mark both sampled.
  g <- set_node_attrs(g, c("P", "Q"), samples = c(1L, 1L))
  gt <- add_sampled_tips(g)
  lv <- get_leafnames(edges_to_igraph(gt))
  expect_true(all(c("P", "Q") %in% lv))            # both now leaf tips
  expect_true(all(c("P_anc", "Q_anc") %in% c(gt$from, gt$to)))
})

test_that("add_sampled_tips preserves admixture when the sampled internal node is an admix node", {
  # Sampled internal node M has in-degree 2 (admixture) and a descendant.
  g <- tibble::tribble(
    ~from, ~to,
    "R", "out", "R", "A", "R", "B",
    "A", "M", "B", "M",          # M is admixed (in-degree 2)
    "M", "X"                     # M ancestral to X (internal)
  )
  g <- set_node_attrs(g, "M", samples = 1L)
  gt <- add_sampled_tips(g)
  # Both admixture parents must remain on the new ancestor (in-degree 2 kept).
  expect_equal(sum(gt$to == "M_anc"), 2L)
  expect_true("M" %in% get_leafnames(edges_to_igraph(gt)))   # M now a tip
})

test_that("read_lgo(rha20) -> add_sampled_tips yields a fittable topology", {
  # The headline case. No f2 needed; assert structure only.
  g  <- suppressWarnings(read_lgo(path = test_path("fixtures", "rha20.lgo")))
  nt <- graph_nodes(g)
  internal_sampled <- setdiff(nt$name[!is.na(nt$samples)],
                              setdiff(g$to, g$from))
  skip_if(length(internal_sampled) == 0, "rha20 fixture has no internal samples")
  gt <- add_sampled_tips(g)
  lv <- get_leafnames(edges_to_igraph(gt))
  # every previously-internal sampled node is now a fitted leaf
  expect_true(all(internal_sampled %in% lv))
  # result is a valid edge tibble
  expect_s3_class(gt, "tbl_df")
})

# --- The qpgraph guard: refuse a graph that marks internal nodes as sampled.

test_that("qpgraph refuses a graph with a sampled internal node", {
  g <- make_test_nodes_graph()                 # anc internal, samples=1
  expect_error(
    qpgraph(example_f2_blocks, g),             # any f2 arg; guard fires first
    class = "admixtools_internal_samples_need_tips")
})

test_that("qpgraph guard is silent once the flag is cleared (intent B)", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "anc", samples = NA)  # leaf-only intent
  # No guard error. (May still error later for unrelated reasons; we only
  # assert the guard class does NOT fire.)
  err <- tryCatch(qpgraph(example_f2_blocks, g), error = function(e) e)
  expect_false(inherits(err, "admixtools_internal_samples_need_tips"))
})

test_that("qpgraph guard is silent after add_sampled_tips (intent A)", {
  g <- add_sampled_tips(make_test_nodes_graph())
  err <- tryCatch(qpgraph(example_f2_blocks, g), error = function(e) e)
  expect_false(inherits(err, "admixtools_internal_samples_need_tips"))
})

# --- Robustness: orphan rows, samples = 0, and edge-column preservation
# (all three call sites share get_internal_sampled, so they agree).

test_that("add_sampled_tips ignores orphan nodes-tibble rows (no fabricated edge)", {
  g  <- make_test_nodes_graph()
  nt <- graph_nodes(g)
  nt <- dplyr::bind_rows(nt, tibble::tibble(name = "ghost", samples = 1L))  # not in any edge
  attr(g, "nodes") <- nt
  gt <- add_sampled_tips(g)
  expect_false(any(c("ghost", "ghost_anc") %in% c(gt$from, gt$to)))   # orphan not materialised
  expect_true("anc" %in% get_leafnames(edges_to_igraph(gt)))          # real internal still tipped
})

test_that("qpgraph guard ignores orphan sampled rows", {
  g <- tibble::tibble(from = c("R", "R"), to = c("A", "B"),
                      type = "normal", weight = NA_real_)
  attr(g, "nodes") <- tibble::tibble(name = "ghost", samples = 1L)     # orphan only
  err <- tryCatch(qpgraph(example_f2_blocks, g), error = function(e) e)
  expect_false(inherits(err, "admixtools_internal_samples_need_tips"))
})

test_that("add_sampled_tips and the guard treat samples = 0 as unsampled", {
  g <- set_node_attrs(make_test_nodes_graph(), "anc", samples = 0L)   # 0 = not sampled
  gt <- add_sampled_tips(g)
  expect_false("anc_anc" %in% c(gt$from, gt$to))                      # not materialised
  expect_setequal(c(gt$from, gt$to), c(g$from, g$to))                 # no-op
  err <- tryCatch(qpgraph(example_f2_blocks, g), error = function(e) e)
  expect_false(inherits(err, "admixtools_internal_samples_need_tips"))
})

test_that("add_sampled_tips carries every edge column through (e.g. time)", {
  g  <- make_test_nodes_graph()
  nt <- attr(g, "nodes")
  g$time <- seq_len(nrow(g)) * 1.0      # a non-(type/weight) edge column
  attr(g, "nodes") <- nt                # re-attach (column assignment can drop attrs)
  gt <- add_sampled_tips(g)
  expect_true("time" %in% names(gt))                  # column survives
  expect_true(is.numeric(gt$time))                    # type preserved
  expect_equal(sum(!is.na(gt$time)), nrow(g))         # all original times retained
  expect_equal(sum(is.na(gt$time)), 1L)               # the one new tip edge is NA
})
