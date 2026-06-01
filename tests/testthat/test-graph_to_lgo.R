test_that("graph_to_lgo errors when samples named-vector has unknown leaf names", {
  # Typo: NoSuchLeaf doesn't exist in the minimal graph (only X is a leaf)
  expect_error(
    graph_to_lgo(make_minimal_graph(), samples = c(NoSuchLeaf = 5L)),
    class = "legofit_invalid_input"
  )
})

test_that("graph_to_lgo igraph input produces byte-identical output to golden", {
  ig      <- make_minimal_igraph()
  expect_warning(txt <- graph_to_lgo(ig), class = "legofit_unfittable_lgo")
  golden  <- readLines(testthat::test_path("fixtures", "minimal.lgo"))
  expect_equal(strsplit(txt, "\n", fixed = TRUE)[[1]], golden)
})

test_that("validate_via_roundtrip catches corrupted .lgo", {
  edges <- make_minimal_graph()
  txt <- paste(
    readLines(testthat::test_path("fixtures", "minimal.lgo")),
    collapse = "\n"
  )
  # Corrupt: introduce a reference to a non-existent parent "Z"
  corrupted <- gsub("derive A from R", "derive A from Z", txt, fixed = TRUE)
  expect_error(
    validate_via_roundtrip(corrupted, edges),
    class = "legofit_validation_failed"
  )
})

test_that("validate_via_roundtrip accepts an uncorrupted .lgo", {
  edges <- make_minimal_graph()
  txt <- paste(
    readLines(testthat::test_path("fixtures", "minimal.lgo")),
    collapse = "\n"
  )
  expect_true(validate_via_roundtrip(txt, edges))
})

test_that("assemble_lgo produces the canonical .lgo for minimal_graph", {
  g     <- make_minimal_graph()
  edges <- coerce_to_edge_tibble(g)
  validate_edge_tibble(edges)
  params     <- generate_param_names(edges)
  times      <- compute_times(edges, "fix_admix", default_drift_to_time,
                              dates_terminal = 0)
  twoN_decls <- resolve_twoN_decls(edges, NULL)
  leaves      <- setdiff(edges$to, edges$from)
  all_nodes   <- unique(c(edges$from, edges$to))
  samples_vec <- resolve_scalar_or_named(1L, leaves, default = NA_integer_)
  full_samples <- setNames(rep(NA_integer_, length(all_nodes)), all_nodes)
  full_samples[leaves] <- as.integer(samples_vec[leaves])
  txt <- assemble_lgo(edges, params, times, twoN_decls, full_samples)

  expected <- readLines(testthat::test_path("fixtures", "minimal.lgo"))
  expect_equal(strsplit(txt, "\n", fixed = TRUE)[[1]], expected)
})

test_that("compute_times fix_admix gives expected exact values", {
  g <- make_minimal_graph()
  result <- compute_times(g, "fix_admix", default_drift_to_time, 0)
  expect_equal(result$value[["R"]], 0.07, tolerance = 1e-9)
  expect_equal(result$value[["A"]], 0.02, tolerance = 1e-9)
  expect_equal(result$value[["B"]], 0.02, tolerance = 1e-9)
  expect_equal(result$value[["M"]], 0.02, tolerance = 1e-9)
  expect_equal(result$value[["X"]], 0.00, tolerance = 1e-9)
  expect_true(result$free[["R"]])
  # Under Path B (ghost-segment encoding), admix dests get their own
  # `t=` declaration (free, like other internals); the LEGOFIT
  # param-sharing constraint is handled by ghost segments anchoring at
  # T_admix_<dest> instead of by omitting the dest's t=.
  expect_true(result$free[["M"]])
  expect_false(result$free[["X"]])  # leaf still not free
  expect_length(result$omit, 0)     # omit semantics retired in Path B
})

test_that("compute_times init: errors when default drift_to_time gives inconsistent admix times", {
  # With default drift_to_time (identity) on minimal_graph,
  # the two admix branch lengths (0.60, 0.40) propagate to inconsistent
  # root times via the two paths. This is expected behavior.
  g <- make_minimal_graph()
  expect_error(
    compute_times(g, "init", default_drift_to_time, 0),
    class = "legofit_invalid_input"
  )
})

test_that("F-2: additive-inconsistency abort points to free/init alternatives", {
  g <- make_minimal_graph()
  expect_error(
    compute_times(g, "init", default_drift_to_time, 0),
    regexp = "free|init", class = "legofit_invalid_input"
  )
})

test_that("F-4: fix_times emits `time fixed` for internal nodes (recovers absolute twoN)", {
  expect_warning(
    txt <- graph_to_lgo(make_minimal_graph(), time_handling = "free",
                        fix_times = TRUE, validate = FALSE),
    class = "legofit_unfittable_lgo")
  expect_match(txt, "time fixed T_R")      # internal tree-node time is fixed
  expect_false(grepl("time free", txt))    # no free time remains
  expect_type(read_lgo(text = txt), "list")  # still parses
})

test_that("F-4: fix_times works with free mode on an admix graph", {
  # free is the only mode that exports additively-inconsistent admix graphs;
  # fix_times must still pin every time (including the admix-event time) so the
  # absolute scale is identifiable rather than degenerate against free twoN.
  txt <- graph_to_lgo(make_rogers2020_graph(), time_handling = "free",
                      fix_times = TRUE, validate = FALSE)
  expect_match(txt, "time fixed")
  expect_false(grepl("time free", txt))
})

test_that("F-4: fix_times defaults FALSE (internal times stay free)", {
  expect_warning(
    txt <- graph_to_lgo(make_minimal_graph(), time_handling = "free",
                        validate = FALSE),
    class = "legofit_unfittable_lgo")
  expect_match(txt, "time free T_R")
})

test_that("compute_times errors when graph has no leaves", {
  # Cycle: every node appears as both `from` and `to`
  cyclic <- tibble::tribble(
    ~from, ~to, ~type,  ~weight,
    "A",   "B", "normal", 0.1,
    "B",   "A", "normal", 0.1
  )
  expect_error(
    compute_times(cyclic, "fix_admix", default_drift_to_time, 0),
    class = "legofit_invalid_input"
  )
})

test_that("compute_times init: respects time column on admix edges", {
  g <- make_minimal_graph()
  # Set explicit time column: 0 on admix edges avoids inconsistency
  g$time <- c(0.05, 0.05, 0, 0, 0.02)
  result <- compute_times(g, "init", default_drift_to_time, 0)
  # M is NOT in omit under "init" mode
  expect_false("M" %in% result$omit)
  # M is free under "init" mode
  expect_true(result$free[["M"]])
  # Same absolute times as fix_admix (because time col has same values)
  expect_equal(result$value[["R"]], 0.07, tolerance = 1e-9)
  expect_equal(result$value[["M"]], 0.02, tolerance = 1e-9)
})

test_that("compute_times free: non-leaf times are valid starting values (not NA)", {
  # LEGOFIT requires `free` declarations to carry a starting value.
  # We emit topological depth (deeper = larger time) so parent > child
  # holds and the optimizer has a sensible starting point.
  g <- make_minimal_graph()
  result <- compute_times(g, "free", default_drift_to_time, 0)
  leaves <- setdiff(g$to, g$from)
  non_leaves <- setdiff(unique(c(g$from, g$to)), leaves)
  expect_true(all(!is.na(result$value[non_leaves])))
  expect_true(all(result$free[non_leaves]))
  expect_false(any(result$free[leaves]))
  expect_equal(length(result$omit), 0)
  # Specifically: depth(R) > depth(A) > depth(M) > depth(X=leaf=0)
  # because R is grandparent of M and M is grandparent of X.
  expect_gt(result$value[["R"]], result$value[["M"]])
  expect_gt(result$value[["M"]], result$value[["X"]])
})

test_that("strip_outgroup strips outgroup and re-roots", {
  g <- make_minimal_graph_with_outgroup()
  result <- strip_outgroup(g, "Outgroup")
  # 7 input rows - 2 root-outgoing edges = 5 rows
  expect_equal(nrow(result), 5)
  # IngroupRoot is the new root (appears in from, not in to)
  new_root <- setdiff(result$from, result$to)
  expect_equal(new_root, "IngroupRoot")
})

test_that("strip_outgroup errors when outpop not found", {
  g <- make_minimal_graph_with_outgroup()
  expect_error(
    strip_outgroup(g, "NoSuchNode"),
    class = "legofit_invalid_input"
  )
})

test_that("strip_outgroup errors on multifurcating root", {
  bad <- tibble::tribble(
    ~from, ~to,        ~type,  ~weight,
    "R",   "Outgroup", "normal", 0.05,
    "R",   "A",        "normal", 0.05,
    "R",   "B",        "normal", 0.05
  )
  expect_error(
    strip_outgroup(bad, "Outgroup"),
    class = "legofit_invalid_input"
  )
})

test_that("strip_outgroup errors when outpop not adjacent to root", {
  g <- make_minimal_graph()
  # X is a leaf but not adjacent to root R
  expect_error(
    strip_outgroup(g, "X"),
    class = "legofit_invalid_input"
  )
})

test_that("resolve_twoN_decls NULL -> coalescent-unit form", {
  g <- make_minimal_graph()
  result <- resolve_twoN_decls(g, NULL)
  expect_equal(result$declarations, "twoN fixed one=1")
  expect_true(all(result$per_segment == "one"))
})

test_that("resolve_twoN_decls scalar -> fixed shared", {
  g <- make_minimal_graph()
  result <- resolve_twoN_decls(g, 1000)
  expect_match(result$declarations, "twoN fixed shared=1000")
  expect_true(all(result$per_segment == "shared"))
})

test_that("resolve_twoN_decls named vector -> free per-segment", {
  g <- make_minimal_graph()
  all_nodes <- unique(c(g$from, g$to))
  twoN_vals <- setNames(rep(10000, length(all_nodes)), all_nodes)
  result <- resolve_twoN_decls(g, twoN_vals)
  expect_match(result$declarations, "twoN free")
  expect_true(all(startsWith(result$per_segment, "twoN_")))
})

test_that("resolve_twoN_decls errors on unknown node names", {
  g <- make_minimal_graph()
  all_nodes <- unique(c(g$from, g$to))
  twoN_vals <- setNames(rep(10000, length(all_nodes)), all_nodes)
  twoN_vals[["ExtraNode"]] <- 5000
  expect_error(
    resolve_twoN_decls(g, twoN_vals),
    class = "legofit_invalid_input"
  )
})

test_that("generate_param_names matches the parameter naming convention", {
  g <- make_minimal_graph()
  p <- generate_param_names(g)
  expect_equal(p$time["R"],    c(R = "T_R"))
  expect_equal(p$twoN["A"],    c(A = "twoN_A"))
  expect_equal(p$mixFrac["M"], c(M = "m_M"))
  expect_setequal(names(p$mixFrac), "M")  # only admix destinations
})

# ---------------------------------------------------------------------------
# Master plan Tier 0 gap closure
# ---------------------------------------------------------------------------

# U-W12: graph_to_lgo(file=) writes exactly the text it returns when file=NULL.
# The minimal graph emits the unfittable warning (one leaf); that is orthogonal
# to the write contract, so it is suppressed here.
test_that("U-W12: graph_to_lgo(file=) writes content equal to the file=NULL text return", {
  g   <- make_minimal_graph()
  txt <- suppressWarnings(graph_to_lgo(g))
  tmp <- tempfile(fileext = ".lgo")
  on.exit(unlink(tmp))
  suppressWarnings(graph_to_lgo(g, file = tmp))
  expect_identical(readLines(tmp), strsplit(txt, "\n", fixed = TRUE)[[1]])
})
