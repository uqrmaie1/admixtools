# Tests for the nodes-tibble attribute API

test_that("graph_nodes returns empty 7-column tibble when no attr", {
  edges <- tibble::tibble(from = "A", to = "B",
                          type = "normal", weight = NA_real_)
  nt <- graph_nodes(edges)
  expect_s3_class(nt, "tbl_df")
  expect_setequal(names(nt),
    c("name", "samples", "twoN_param", "twoN",
      "time_param", "time", "admix_event_time"))
  expect_equal(nrow(nt), 0)
})

test_that("graph_nodes returns the attached tibble verbatim", {
  edges <- tibble::tibble(from = "A", to = "B",
                          type = "normal", weight = NA_real_)
  attr(edges, "nodes") <- tibble::tibble(name = "B", samples = 1L)
  nt <- graph_nodes(edges)
  expect_equal(nt$name, "B")
  expect_equal(nt$samples, 1L)
})

test_that("graph_nodes errors on igraph input", {
  ig <- igraph::graph_from_edgelist(matrix(c("A", "B"), ncol = 2))
  expect_error(graph_nodes(ig), class = "admixtools_invalid_graph")
})

test_that("set_node_attrs rejects unknown name", {
  g <- make_test_nodes_graph()
  expect_error(
    set_node_attrs(g, "not_a_real_node", samples = 1L),
    class = "admixtools_invalid_graph")
})

test_that("set_node_attrs rejects duplicate name", {
  g <- make_test_nodes_graph()
  expect_error(
    set_node_attrs(g, c("v", "v"), samples = c(1L, 2L)),
    class = "admixtools_invalid_graph")
})

test_that("set_node_attrs rejects length mismatch", {
  g <- make_test_nodes_graph()
  expect_error(
    set_node_attrs(g, c("v", "afr"), samples = c(1L, 2L, 3L)),
    class = "admixtools_invalid_graph")
})

test_that("set_node_attrs recycles scalar to length(name)", {
  g <- make_test_nodes_graph()
  g2 <- set_node_attrs(g, c("v", "afr"), samples = 2L)
  nt <- graph_nodes(g2)
  expect_equal(nt$samples[nt$name == "v"],   2L)
  expect_equal(nt$samples[nt$name == "afr"], 2L)
})

test_that("set_node_attrs takes per-name vector when length matches", {
  g <- make_test_nodes_graph()
  g2 <- set_node_attrs(g, c("v", "afr"), samples = c(3L, 5L))
  nt <- graph_nodes(g2)
  expect_equal(nt$samples[nt$name == "v"],   3L)
  expect_equal(nt$samples[nt$name == "afr"], 5L)
})

test_that("set_node_attrs adds new rows for canonical-but-not-in-nodes names", {
  edges <- tibble::tibble(from = c("R","R"), to = c("A","B"),
                          type = "normal", weight = NA_real_)
  # No nodes attr to start
  g <- set_node_attrs(edges, "A", samples = 1L)
  nt <- graph_nodes(g)
  expect_equal(nt$name, "A")
  expect_equal(nt$samples, 1L)
})

test_that("set_node_attrs coerces numeric samples to integer", {
  g <- make_test_nodes_graph()
  g2 <- set_node_attrs(g, "v", samples = 2)   # numeric, not integer
  nt <- graph_nodes(g2)
  expect_type(nt$samples, "integer")
  expect_equal(nt$samples[nt$name == "v"], 2L)
})

test_that("prune_nodes_attr drops orphan rows", {
  g <- make_test_nodes_graph()
  nt <- graph_nodes(g)
  nt <- dplyr::bind_rows(nt, tibble::tibble(name = "ghost", samples = 1L))
  attr(g, "nodes") <- nt
  g2 <- prune_nodes_attr(g)
  expect_false("ghost" %in% graph_nodes(g2)$name)
})

test_that("prune_nodes_attr preserves in-set rows", {
  g <- make_test_nodes_graph()
  g2 <- prune_nodes_attr(g)
  expect_equal(graph_nodes(g)$name, graph_nodes(g2)$name)
})

test_that("refresh_edge_times rebuilds edges$time from nodes$time", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, c("v", "afr", "eur"), time = c(1.5, 0.0, 0.0))
  # set_node_attrs writes the denormalized edge view itself; verify it
  expect_equal(g$time[g$to == "v"],   1.5)
  expect_equal(g$time[g$to == "afr"], 0.0)
  expect_equal(g$time[g$to == "eur"], 0.0)
})

test_that("refresh_edge_times rebuilds admix_event_time", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "v", admix_event_time = 0.05)
  expect_equal(g$admix_event_time[g$to == "v"], 0.05)
})

test_that("refresh_edge_times is a no-op when nodes tibble is empty", {
  edges <- tibble::tibble(from = "A", to = "B",
                          type = "normal", weight = NA_real_)
  out <- refresh_edge_times(edges)
  expect_identical(out, edges)
})

test_that("validate_edge_tibble(strict = TRUE) errors on orphan rows", {
  g <- make_test_nodes_graph()
  nt <- graph_nodes(g)
  nt <- dplyr::bind_rows(nt, tibble::tibble(name = "ghost", samples = 1L))
  attr(g, "nodes") <- nt
  expect_error(
    validate_edge_tibble(g, strict = TRUE),
    class = "admixtools_invalid_graph")
})

test_that("validate_edge_tibble(strict = FALSE) warns and prunes orphan rows", {
  g <- make_test_nodes_graph()
  nt <- graph_nodes(g)
  nt <- dplyr::bind_rows(nt, tibble::tibble(name = "ghost", samples = 1L))
  attr(g, "nodes") <- nt
  expect_warning(
    g2 <- validate_edge_tibble(g, strict = FALSE),
    class = "admixtools_orphan_nodes_pruned")
  expect_false("ghost" %in% graph_nodes(g2)$name)
})

test_that("validate_edge_tibble(strict = TRUE) tolerates 1e-12 float drift", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "v", time = 1.0)
  # Manually drift one edge by 1e-12 (still inside tolerance)
  g$time[g$to == "v"] <- 1.0 + 1e-12
  expect_no_error(validate_edge_tibble(g, strict = TRUE))
})

test_that("validate_edge_tibble(strict = TRUE) errors on float drift past tolerance", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "v", time = 1.0)
  g$time[g$to == "v"] <- 1.0 + 1e-6  # well outside tolerance
  expect_error(
    validate_edge_tibble(g, strict = TRUE),
    class = "admixtools_invalid_graph")
})

test_that("validate_edge_tibble(strict = TRUE) errors on duplicate node names", {
  g <- make_test_nodes_graph()
  nt <- graph_nodes(g)
  # Force a duplicate via direct attr manipulation
  nt <- dplyr::bind_rows(nt, nt[nt$name == "v", ])
  attr(g, "nodes") <- nt
  expect_error(
    validate_edge_tibble(g, strict = TRUE),
    class = "admixtools_invalid_graph")
})

test_that("validate_edge_tibble(strict = FALSE) warns on duplicate node names", {
  g <- make_test_nodes_graph()
  nt <- graph_nodes(g)
  nt <- dplyr::bind_rows(nt, nt[nt$name == "v", ])
  attr(g, "nodes") <- nt
  expect_warning(
    validate_edge_tibble(g, strict = FALSE),
    class = "admixtools_duplicate_node_name")
})

test_that("as_edge_tibble drops attr cleanly when empty", {
  edges <- tibble::tibble(from = "A", to = "B",
                          type = "normal", weight = NA_real_)
  attr(edges, "nodes") <- tibble::tibble(
    name = character(0), samples = integer(0),
    twoN_param = character(0), twoN = numeric(0),
    time_param = character(0), time = numeric(0),
    admix_event_time = numeric(0)
  )
  out <- expect_no_warning(as_edge_tibble(edges))
  expect_null(attr(out, "nodes"))
})

test_that("as_edge_tibble warns on non-empty drop", {
  g <- make_test_nodes_graph()
  expect_warning(
    as_edge_tibble(g),
    class = "admixtools_dropping_node_attrs")
})

test_that("node_times returns nodes$time as named numeric", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, c("v", "afr", "eur"),
                      time = c(1.5, 0.0, 0.0))
  nt_vec <- node_times(g)
  expect_type(nt_vec, "double")
  expect_named(nt_vec)
  expect_equal(nt_vec[["v"]], 1.5)
})

test_that("node_times returns empty named numeric when nodes tibble empty", {
  edges <- tibble::tibble(from = "A", to = "B",
                          type = "normal", weight = NA_real_)
  nt_vec <- node_times(edges)
  expect_length(nt_vec, 0)
})

test_that("coerce_to_edge_tibble extracts V(g)$samples into nodes tibble", {
  ig <- igraph::graph_from_edgelist(matrix(
    c("R","A", "R","B"), ncol = 2, byrow = TRUE))
  igraph::V(ig)$samples <- c(NA, 1, 1)  # length = vcount(ig) = 3
  edges <- coerce_to_edge_tibble(ig)
  nt <- graph_nodes(edges)
  expect_type(nt$samples, "integer")
  expect_equal(nt$samples[nt$name == "A"], 1L)
  expect_equal(nt$samples[nt$name == "B"], 1L)
})

test_that("coerce_to_edge_tibble warns admixtools_ambiguous_nodes_metadata", {
  ig <- igraph::graph_from_edgelist(matrix(c("R","A"), ncol = 2))
  igraph::V(ig)$samples <- c(NA, 1)
  attr(ig, "nodes") <- tibble::tibble(name = "A", samples = 99L)
  expect_warning(
    edges <- coerce_to_edge_tibble(ig),
    class = "admixtools_ambiguous_nodes_metadata")
  nt <- graph_nodes(edges)
  # vertex_attr wins: samples = 1 (not 99), and type must be integer
  expect_type(nt$samples, "integer")
  expect_equal(nt$samples[nt$name == "A"], 1L)
})

test_that("coerce_to_edge_tibble errors on igraph with no vertex names", {
  ig <- igraph::make_empty_graph(n = 3, directed = TRUE)
  # No names assigned
  expect_error(
    coerce_to_edge_tibble(ig),
    class = "admixtools_invalid_graph")
})

test_that("coerce_to_edge_tibble preserves edge_attr extraction (regression)", {
  # the read direction behavior: edge attrs weight, type, time still extracted
  ig <- igraph::graph_from_edgelist(matrix(c("R","A"), ncol = 2))
  igraph::E(ig)$weight <- 0.5
  igraph::V(ig)$name <- c("R", "A")  # vertex names are required
  edges <- coerce_to_edge_tibble(ig)
  expect_equal(edges$weight, 0.5)
})

test_that("coerce_to_edge_tibble always emits full 7-column nodes tibble", {
  ig <- igraph::graph_from_edgelist(matrix(c("R","A","R","B"), ncol=2, byrow=TRUE))
  # No vertex attrs beyond auto-assigned names
  edges <- coerce_to_edge_tibble(ig)
  nt <- graph_nodes(edges)
  expect_setequal(names(nt),
    c("name","samples","twoN_param","twoN","time_param","time","admix_event_time"))
  expect_equal(nrow(nt), 3L)
  # set_node_attrs must work on this result without crashing
  expect_no_error(set_node_attrs(edges, "A", samples = 1L))
})

test_that("refresh_edge_times does not clobber existing times for nodes outside the nodes tibble", {
  # Use a custom graph so we control exactly which nodes are tracked
  edges <- tibble::tibble(
    from = c("R", "R", "A"),
    to   = c("A", "B", "C"),
    type = "normal", weight = NA_real_
  )
  edges$time <- c(5, 4, 3)
  # Only track node A in the nodes tibble
  g <- set_node_attrs(edges, "A", time = 1.5)
  # A is tracked: its edge time must be updated
  expect_equal(g$time[g$to == "A"], 1.5)
  # B and C are not tracked: their edge times must be unchanged
  expect_equal(g$time[g$to == "B"], 4)
  expect_equal(g$time[g$to == "C"], 3)
})

test_that("refresh_edge_times does not create admix_event_time column on tree graphs", {
  g <- make_test_nodes_graph()   # pure tree, no admix events
  g <- set_node_attrs(g, "v", time = 1.0)
  # time column should exist, admix_event_time should NOT be created
  expect_true("time" %in% names(g))
  expect_false("admix_event_time" %in% names(g))
})

test_that("node_times returns all-NA named numeric when nodes tibble lacks time column", {
  edges <- tibble::tibble(from = "R", to = "A", type = "normal", weight = NA_real_)
  attr(edges, "nodes") <- tibble::tibble(name = "A", samples = 1L)
  nt_vec <- node_times(edges)
  expect_type(nt_vec, "double")
  expect_named(nt_vec)
  expect_true(is.na(nt_vec[["A"]]))
})

test_that("as_edge_tibble warns when nodes tibble has only non-standard columns with data", {
  edges <- tibble::tibble(from = "R", to = "A", type = "normal", weight = NA_real_)
  attr(edges, "nodes") <- tibble::tibble(name = "A", my_custom_col = 42)
  expect_warning(as_edge_tibble(edges), class = "admixtools_dropping_node_attrs")
})

test_that("refresh_edge_times does not clobber edge time when tracked node has NA time", {
  # Bug: if node A is tracked (non-time attr) and node B is tracked with time=t,
  # triggering refresh must not overwrite A's pre-existing edge time with NA.
  edges <- tibble::tibble(
    from = c("R", "R"), to = c("A", "B"),
    type = "normal", weight = NA_real_
  )
  edges$time <- c(5.0, 4.0)
  # Track A for a non-time attribute; its edge time (5) must be preserved.
  edges <- set_node_attrs(edges, "A", twoN = 1000.0)
  # Now track B with an explicit time; this triggers refresh_edge_times.
  edges <- set_node_attrs(edges, "B", time = 3.0)
  expect_equal(edges$time[edges$to == "A"], 5.0)  # must not be clobbered
  expect_equal(edges$time[edges$to == "B"], 3.0)
})

test_that("refresh_edge_times does not clobber admix_event_time for NA-time tracked node", {
  edges <- tibble::tibble(
    from = c("R", "R"), to = c("A", "B"),
    type = "admix", weight = 0.5
  )
  edges$admix_event_time <- c(10.0, 8.0)
  edges <- set_node_attrs(edges, "A", twoN = 500.0)
  edges <- set_node_attrs(edges, "B", admix_event_time = 7.0)
  expect_equal(edges$admix_event_time[edges$to == "A"], 10.0)
  expect_equal(edges$admix_event_time[edges$to == "B"], 7.0)
})

test_that("coerce_to_edge_tibble preserves integer type for samples from igraph", {
  ig <- igraph::graph_from_edgelist(matrix(c("R","A"), ncol = 2))
  igraph::V(ig)$samples <- c(NA_real_, 2.0)
  edges <- coerce_to_edge_tibble(ig)
  expect_type(graph_nodes(edges)$samples, "integer")
})

test_that("coerce_to_edge_tibble skips V(ig)$time when E(ig)$time already extracted", {
  ig <- igraph::graph_from_edgelist(matrix(c("R","A"), ncol = 2))
  igraph::E(ig)$time <- 0.5          # branch length in edge attrs
  igraph::V(ig)$time <- c(1.0, 0.0) # absolute times in vertex attrs
  edges <- coerce_to_edge_tibble(ig)
  # edges$time must hold the edge-level branch length, not the vertex-level abs time
  expect_equal(edges$time, 0.5)
  # nodes$time must be all-NA: vertex-level time was skipped because edge attr is non-NA
  expect_true(all(is.na(graph_nodes(edges)$time)))
  # validate_edge_tibble(strict=TRUE) must not false-positive abort
  # (vacuous: strict block skips when nodes$time is all-NA)
  expect_no_error(validate_edge_tibble(edges, strict = TRUE))
})

test_that("validate_edge_tibble(strict=TRUE) errors when nodes$time set but edges$time absent", {
  edges <- tibble::tibble(from = "R", to = "A", type = "normal", weight = NA_real_)
  # Bypass set_node_attrs to get nodes$time without triggering refresh_edge_times
  attr(edges, "nodes") <- tibble::tibble(
    name = "A", samples = NA_integer_, twoN_param = NA_character_,
    twoN = NA_real_, time_param = NA_character_,
    time = 1.5, admix_event_time = NA_real_
  )
  # edges has no $time column - strict mode must catch this
  expect_error(
    validate_edge_tibble(edges, strict = TRUE),
    class = "admixtools_invalid_graph")
})

# --- regression tests for time-consistency edge cases ---------------------

test_that("validate_edge_tibble(strict=TRUE) errors when nodes$admix_event_time set but edges$admix_event_time absent", {
  # admix_event_time block must mirror the time block
  # and abort when the edge column is missing.
  # Use type="normal" so the admix-pair invariant check does not fire first.
  edges <- tibble::tibble(from = "R", to = "A", type = "normal", weight = NA_real_)
  attr(edges, "nodes") <- tibble::tibble(
    name = "A", samples = NA_integer_, twoN_param = NA_character_,
    twoN = NA_real_, time_param = NA_character_,
    time = NA_real_, admix_event_time = 0.5
  )
  expect_error(
    validate_edge_tibble(edges, strict = TRUE),
    class = "admixtools_invalid_graph")
})

test_that("set_node_attrs time=NA_real_ clears the edge time entry", {
  # targeted write must propagate NA explicitly,
  # not skip it the way refresh_edge_times (full-rebuild) would.
  g <- tibble::tibble(from = "R", to = "A", type = "normal", weight = NA_real_)
  g <- set_node_attrs(g, "A", time = 1.5)
  expect_equal(g$time, 1.5)
  g <- set_node_attrs(g, "A", time = NA_real_)
  expect_true(is.na(g$time))
})

test_that("coerce_to_edge_tibble preserves V(ig)$time when E(ig)$time is all-NA", {
  # skip guard must check any(!is.na(edges[[col]]));
  # an all-NA edge column must NOT suppress extraction of vertex-level times.
  ig <- igraph::graph_from_edgelist(matrix(c("R", "A"), ncol = 2))
  igraph::E(ig)$time <- NA_real_     # edge attr exists but is all-NA
  igraph::V(ig)$time <- c(0.05, 0.0) # vertex-level absolute times
  edges <- coerce_to_edge_tibble(ig)
  nt <- graph_nodes(edges)
  expect_equal(nt$time[nt$name == "R"], 0.05, tolerance = 1e-12)
  expect_equal(nt$time[nt$name == "A"], 0.0,  tolerance = 1e-12)
})

# --- strict mode must flag a present-but-NA edge value --------------------
# A node carrying a time requires its incoming edge to carry that value in
# the denormalized edges$time view. An all-NA edge column is NOT "consistent";
# the old predicate (which required !is.na(x$time)) silently passed it.

test_that("validate_edge_tibble(strict=TRUE) errors when edges$time is NA but nodes$time is set", {
  edges <- tibble::tibble(from = "R", to = "A", type = "normal", weight = NA_real_)
  edges$time <- NA_real_                       # column present, value missing
  attr(edges, "nodes") <- tibble::tibble(
    name             = c("R", "A"),
    samples          = NA_integer_,   twoN_param = NA_character_, twoN = NA_real_,
    time_param       = NA_character_, time       = c(NA_real_, 0.05),  # A has a time
    admix_event_time = NA_real_)
  expect_error(
    validate_edge_tibble(edges, strict = TRUE),
    class = "admixtools_invalid_graph")
})

test_that("validate_edge_tibble(strict=TRUE) errors when edges$admix_event_time is NA but nodes value is set", {
  # type='normal' so the admix-pair invariant check does not fire first.
  edges <- tibble::tibble(from = "R", to = "A", type = "normal", weight = NA_real_)
  edges$admix_event_time <- NA_real_
  attr(edges, "nodes") <- tibble::tibble(
    name             = c("R", "A"),
    samples          = NA_integer_,   twoN_param = NA_character_, twoN = NA_real_,
    time_param       = NA_character_, time       = NA_real_,
    admix_event_time = c(NA_real_, 0.5))
  expect_error(
    validate_edge_tibble(edges, strict = TRUE),
    class = "admixtools_invalid_graph")
})

test_that("strict validation flags a coerced igraph with V time but all-NA E time, and refresh_edge_times reconciles it", {
  # The end-to-end path the bug actually arises on: coerce extracts vertex
  # times into nodes$time but leaves edges$time all-NA. That graph is
  # inconsistent under strict mode; refresh_edge_times is the documented fix.
  ig <- igraph::graph_from_edgelist(matrix(c("R", "A"), ncol = 2))
  igraph::E(ig)$time <- NA_real_
  igraph::V(ig)$time <- c(0.05, 0.0)
  edges <- coerce_to_edge_tibble(ig)
  expect_error(
    validate_edge_tibble(edges, strict = TRUE),
    class = "admixtools_invalid_graph")
  edges2 <- refresh_edge_times(edges)
  expect_no_error(validate_edge_tibble(edges2, strict = TRUE))
})

test_that("validate_edge_tibble(strict=TRUE) still passes a fully consistent denormalized graph", {
  # Guard against over-firing: when every tracked node's time matches its
  # incoming edge, strict mode must NOT error.
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "v", time = 1.0)   # sets nodes$time AND edges$time
  expect_no_error(validate_edge_tibble(g, strict = TRUE))
})

test_that("validate_edge_tibble(strict=TRUE) does not crash on non-finite (Inf) times", {
  # A naive predicate computes Inf - Inf = NaN, NaN > tol = NA, and
  # `if (any(NA))` throws a base-R 'missing value' error. The check must stay
  # in its typed-condition contract: Inf==Inf is consistent (no error); an
  # Inf node time against a finite edge is a genuine mismatch.
  g <- tibble::tibble(from = "R", to = "A", type = "normal", weight = NA_real_)
  g <- set_node_attrs(g, "A", time = Inf)            # nodes$time AND edges$time = Inf
  expect_no_error(validate_edge_tibble(g, strict = TRUE))

  g$time[g$to == "A"] <- 0.5                          # node=Inf, edge=finite -> mismatch
  expect_error(
    validate_edge_tibble(g, strict = TRUE),
    class = "admixtools_invalid_graph")
})
