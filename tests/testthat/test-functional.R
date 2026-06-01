# Tier 1 functional tests for the graph_to_lgo writer.
#
# Unlike test-graph_to_lgo.R (which exercises individual helpers on the
# minimal 5-node fixture), this file runs the integrated pipeline on
# realistic admixtools inputs and exercises corners the minimal graph
# doesn't reach.

# T1.1 ------------------------------------------------------------------
test_that("T1.1: graph_to_lgo() handles admixtools' example_igraph", {
  # example_igraph is the canonical admixtools fixture: 25 nodes, 7 leaves,
  # 2 admixture events. No edge attributes (no `type`, no `weight`); the
  # coerce/derive path must handle that.
  txt <- graph_to_lgo(example_igraph, validate = FALSE)
  expect_type(txt, "character")
  expect_true(nchar(txt) > 0)

  parsed <- read_lgo(text = txt)
  # 25 internal+leaf nodes -> 24 derive edges + 2*2 admix edges - the
  # admixture destinations get 2 incoming admix edges each.
  # Total edges should equal example_igraph's edge count.
  expect_equal(nrow(parsed), igraph::ecount(example_igraph))
})

test_that("T1.1: graph_to_lgo() handles random_sim() topology output", {
  # random_sim() emits weight=1 placeholders on every edge, which means
  # the topology walk will see inconsistent times under "fix_admix" or
  # "init". The right mode for topology-only graphs is "free", which
  # declares all times free with no starting values.
  set.seed(42)
  sim <- random_sim(nadmix = 2, nleaf = 5)
  txt <- graph_to_lgo(sim$edges, time_handling = "free", validate = FALSE)
  expect_type(txt, "character")

  # Round-trip topology preservation
  parsed <- read_lgo(text = txt)
  expect_equal(nrow(parsed), nrow(sim$edges))
  expect_setequal(unique(c(parsed$from, parsed$to)),
                  unique(c(sim$edges$from, sim$edges$to)))
})

# T1.4 ------------------------------------------------------------------
test_that("T1.4: multi-admixture graph round-trips (rha20 / Rogers et al. 2020)", {
  # Real-world published topology: 18 nodes, 4 admixture events. Sourced
  # from rha20.lgo (Rogers et al. 2020 Sci Adv, bundled with LEGOFIT
  # under ISC license). Exercises admix_parent_pairs() with 4 events,
  # not just 1.
  g <- make_rogers2020_graph()
  expect_equal(nrow(g), 21)               # 13 normal + 8 admix
  expect_equal(sum(g$type == "admix"), 8) # 4 events × 2 parents
  expect_equal(sum(g$type == "normal"), 13)

  # Use init mode + the published times as the `time` column on each
  # edge. dates_terminal = c(...) anchors the three leaves (d is sampled
  # at Td=3484.25 in the original LEGOFIT model). The topology walk
  # reconstructs all internal times from these branch lengths.
  txt <- graph_to_lgo(
    g,
    time_handling  = "init",
    dates_terminal = c(x = 0, y = 0, d = 3484.25),
    validate       = TRUE
  )

  # All 4 mixFrac declarations should appear
  expect_match(txt, "mixFrac free m_d=")
  expect_match(txt, "mixFrac free m_y=")
  expect_match(txt, "mixFrac free m_a=")
  expect_match(txt, "mixFrac free m_nd=")

  # All 4 mix statements should appear
  expect_match(txt, "mix d from")
  expect_match(txt, "mix y from")
  expect_match(txt, "mix a from")
  expect_match(txt, "mix nd from")

  # Round-trip topology preservation
  parsed <- read_lgo(text = txt)
  # read_lgo attaches a nodes tibble; the topology check compares
  # from/to/type only, so drop it before the tibble comparison.
  attr(parsed, "nodes") <- NULL
  expect_equal(nrow(parsed), nrow(g))
  g_topo <- g %>% dplyr::select(from, to, type) %>% dplyr::arrange(from, to)
  p_topo <- parsed %>% dplyr::select(from, to, type) %>% dplyr::arrange(from, to)
  expect_equal(p_topo, g_topo)
})

test_that("T1.4: admix_parent_pairs picks correct high/low for each event", {
  # The published mixfracs are all < 0.5, so the "high" parent is the
  # first-listed parent in rha20.lgo's mix statements (d2 for d, y2 for
  # y, a2 for a, nd2 for nd) and the "low" parent is the second.
  g <- make_rogers2020_graph()
  txt <- graph_to_lgo(
    g,
    time_handling  = "init",
    dates_terminal = c(x = 0, y = 0, d = 3484.25),
    validate       = FALSE
  )
  # Under Path B (ghost-segment encoding), the mix statement references
  # ghost names <dest>_<parent>, not the real parent populations.
  # Format: "mix <to> from <to>_<high> + m_<to> * <to>_<low>"
  expect_match(txt, "mix d from d_d2 \\+ m_d \\* d_s")
  expect_match(txt, "mix y from y_y2 \\+ m_y \\* y_n")
  expect_match(txt, "mix a from a_a2 \\+ m_a \\* a_xy2")
  expect_match(txt, "mix nd from nd_nd2 \\+ m_nd \\* nd_s2")
})

test_that("T1.4: rha20 topology walk reconstructs published LEGOFIT times", {
  # The branch lengths in make_rogers2020_graph were derived from the
  # rha20.lgo fitted times. With dates_terminal anchoring the 3 leaves,
  # the topology walk should reconstruct each internal node's published
  # absolute time.
  g <- make_rogers2020_graph()
  txt <- graph_to_lgo(
    g,
    time_handling  = "init",
    dates_terminal = c(x = 0, y = 0, d = 3484.25),
    validate       = FALSE
  )
  # Spot-check a few of the reconstructed times by regex against the
  # output. The published values (with %g formatting) appear verbatim:
  expect_match(txt, "T_xynds=82008\\.2")     # root
  expect_match(txt, "T_xynd=25920")          # fixed Txynd
  expect_match(txt, "T_nd=25416\\.9")        # Tnd
  expect_match(txt, "T_av=16307")            # Tav
  expect_match(txt, "T_xy=774\\.856")        # Txy
  expect_match(txt, "T_v=2511\\.5")          # Tv
})

# T1.1.b ----------------------------------------------------------------
test_that("T1.1: qpgraph(example_igraph) round-trips via time_handling='free'", {
  # qpgraph emits type='edge' (vs random_sim's 'normal'); we accept both.
  # Fitted qpgraph drifts are optimized for f-stat fit, not time-walk
  # additive consistency, so "fix_admix" / "init" would error on the
  # consistency check. "free" mode bypasses the topology walk and emits
  # depth-based starting values, so qpgraph outputs round-trip cleanly.
  skip_if_not(exists("example_f2_blocks"))
  result <- tryCatch(
    qpgraph(example_f2_blocks, example_igraph, return_fstats = FALSE),
    error = function(e) NULL
  )
  skip_if(is.null(result), "qpgraph fit failed")
  skip_if(is.null(result$edges), "qpgraph did not return edges")
  fitted_edges <- result$edges

  # Validate doesn't error — confirms 'edge' type accepted.
  expect_no_error(validate_edge_tibble(fitted_edges))

  # Full export via "free" mode + topology round-trip.
  txt <- graph_to_lgo(fitted_edges, time_handling = "free", validate = TRUE)
  expect_type(txt, "character")
  expect_true(nchar(txt) > 0)

  parsed <- read_lgo(text = txt)
  expect_equal(nrow(parsed), nrow(fitted_edges))
  expect_setequal(unique(c(parsed$from, parsed$to)),
                  unique(c(fitted_edges$from, fitted_edges$to)))
})

# T1.2 ------------------------------------------------------------------
test_that("T1.2: time_handling='init' produces minimal-init.lgo verbatim", {
  g <- make_minimal_graph()
  g$time <- c(0.05, 0.05, 0, 0, 0.02)  # admix edges zeroed for consistency
  expect_warning(txt <- graph_to_lgo(g, time_handling = "init", validate = FALSE),
                 class = "legofit_unfittable_lgo")
  expected <- readLines(testthat::test_path("fixtures", "minimal-init.lgo"))
  expect_equal(strsplit(txt, "\n", fixed = TRUE)[[1]], expected)
})

test_that("T1.2: time_handling='free' produces minimal-free.lgo verbatim", {
  expect_warning(
    txt <- graph_to_lgo(make_minimal_graph(), time_handling = "free", validate = FALSE),
    class = "legofit_unfittable_lgo"
  )
  expected <- readLines(testthat::test_path("fixtures", "minimal-free.lgo"))
  expect_equal(strsplit(txt, "\n", fixed = TRUE)[[1]], expected)
})

# T1.3 ------------------------------------------------------------------
test_that("T1.3: twoN=scalar produces minimal-twoN-scalar.lgo verbatim", {
  expect_warning(txt <- graph_to_lgo(make_minimal_graph(), twoN = 1000, validate = FALSE),
                 class = "legofit_unfittable_lgo")
  expected <- readLines(testthat::test_path("fixtures", "minimal-twoN-scalar.lgo"))
  expect_equal(strsplit(txt, "\n", fixed = TRUE)[[1]], expected)
})

test_that("T1.3: twoN=named-vector produces minimal-twoN-named.lgo verbatim", {
  twoN_named <- c(R = 10000, A = 5000, B = 5000, M = 5000, X = 1000)
  expect_warning(txt <- graph_to_lgo(make_minimal_graph(), twoN = twoN_named, validate = FALSE),
                 class = "legofit_unfittable_lgo")
  expected <- readLines(testthat::test_path("fixtures", "minimal-twoN-named.lgo"))
  expect_equal(strsplit(txt, "\n", fixed = TRUE)[[1]], expected)
})

# T1.5 ------------------------------------------------------------------
test_that("T1.5: outpop is fully stripped from the output", {
  g <- make_minimal_graph_with_outgroup()
  expect_warning(txt <- graph_to_lgo(g, outpop = "Outgroup", validate = TRUE),
                 class = "legofit_unfittable_lgo")

  # No segment, derive, or mix should reference Outgroup
  expect_false(grepl("Outgroup", txt))

  # The new root (IngroupRoot) should appear as a segment
  expect_match(txt, "segment IngroupRoot")

  # Parsed result should not contain Outgroup either
  parsed <- read_lgo(text = txt)
  expect_false("Outgroup" %in% c(parsed$from, parsed$to))
})

# T1.6 ------------------------------------------------------------------
test_that("T1.6: custom drift_to_time function scales output times", {
  # Use init mode with a time column on admix edges (=0 to avoid the
  # natural inconsistency), so non-admix branch lengths run through
  # drift_to_time.
  g <- make_minimal_graph()
  g$time <- c(0.05, 0.05, 0, 0, 0.02)  # admix edges zeroed

  scale <- 1000
  expect_warning(
    graph_to_lgo(
      g,
      time_handling = "init",
      drift_to_time = function(d, s) d * scale,
      validate = FALSE
    ),
    class = "legofit_unfittable_lgo"
  )

  # Note: when `time` column is present and non-NA, compute_times uses
  # it directly and bypasses drift_to_time. So this test is really
  # testing that the call shape is accepted, not that the scaling
  # applies. Below: remove the time column to exercise drift_to_time.
  g2 <- make_minimal_graph()
  # Use only normal edges so the topology walk doesn't error on admix
  g2 <- g2 %>% dplyr::filter(type == "normal")
  # But that breaks the graph topology. Use a pure-tree fixture instead:
  tree <- tibble::tribble(
    ~from, ~to, ~type,    ~weight,
    "R",   "A", "normal", 0.05,
    "A",   "B", "normal", 0.03,
    "B",   "C", "normal", 0.02
  )
  expect_warning(
    txt_tree <- graph_to_lgo(
      tree,
      time_handling = "init",
      drift_to_time = function(d, s) d * scale,
      validate = FALSE
    ),
    class = "legofit_unfittable_lgo"
  )
  # Tree: C=0, B=0.02*1000=20, A=20+0.03*1000=50, R=50+0.05*1000=100
  expect_match(txt_tree, "T_R=100")
  expect_match(txt_tree, "T_A=50")
  expect_match(txt_tree, "T_B=20")
  expect_match(txt_tree, "T_C=0")
})

test_that("T1.6: drift_to_time receives correct segment_vec argument", {
  # Custom function that errors if segment_vec is NULL — proves the
  # second argument is being passed correctly.
  tree <- tibble::tribble(
    ~from, ~to, ~type,    ~weight,
    "R",   "A", "normal", 0.05,
    "A",   "B", "normal", 0.03
  )
  picky_drift_to_time <- function(d, s) {
    if (is.null(s)) stop("segment_vec was NULL")
    if (length(s) != length(d)) stop("length mismatch")
    d
  }
  expect_warning(
    expect_no_error(
      graph_to_lgo(tree, time_handling = "init",
                   drift_to_time = picky_drift_to_time,
                   validate = FALSE)
    ),
    class = "legofit_unfittable_lgo"
  )
})

# T1.7 ------------------------------------------------------------------
test_that("T1.7: read_lgo handles a 100-node graph in under 500ms", {
  skip_on_cran()
  # random_sim with many leaves; use "free" mode since weights are
  # placeholder.
  set.seed(123)
  big <- random_sim(nadmix = 3, nleaf = 50)
  expect_gte(length(unique(c(big$edges$from, big$edges$to))), 50)
  txt <- graph_to_lgo(big$edges, time_handling = "free", validate = FALSE)
  t <- system.time(read_lgo(text = txt))[["elapsed"]]
  expect_lt(t, 0.5)
})

# PF-4 ------------------------------------------------------------------
test_that("PF-4: graph_to_lgo on a 100-node graph stays bounded (format_* at scale)", {
  skip_on_cran()
  # The documented performance bound is only on read_lgo; the format_* writer
  # helpers were never profiled at scale. Time the writer on the same large graph.
  set.seed(123)
  big <- random_sim(nadmix = 3, nleaf = 50)
  expect_gte(length(unique(c(big$edges$from, big$edges$to))), 50)
  t <- system.time(
    graph_to_lgo(big$edges, time_handling = "free", validate = FALSE)
  )[["elapsed"]]
  expect_lt(t, 1.0)
})

# T1.8 ------------------------------------------------------------------
test_that("T1.8: graph_to_lgo output is deterministic across runs", {
  g <- make_minimal_graph()
  # minimal graph is unfittable (one leaf), so each call warns; we test
  # determinism here, not the warning, so suppress the expected noise
  outputs <- suppressWarnings(replicate(10, graph_to_lgo(g)))
  expect_length(unique(outputs), 1)
})

test_that("T1.8: graph_to_lgo output is deterministic for igraph input", {
  ig <- make_minimal_igraph()
  # same expected warning as the determinism test for edge tibble input
  outputs <- suppressWarnings(replicate(10, graph_to_lgo(ig)))
  expect_length(unique(outputs), 1)
})

# T2.3 ------------------------------------------------------------------
# Parse rha20.lgo (Rogers et al. 2020, bundled with LEGOFIT under ISC
# license; vendored in tests/testthat/fixtures/rha20.lgo). Confirms the
# parser handles constructs beyond what graph_to_lgo emits: whitespace
# around `=`, multi-line param declarations, and `time constrained`
# arithmetic expressions evaluated safely against earlier params.
#
# The narrow ghost-detection rule (3-condition AND) is used instead of the
# old greedy rule. rha20's intermediate segments (d2, s2, a2, y2, nd2, etc.)
# do NOT match the <dest>_<parent> naming convention used by our writer, so
# they are preserved as real edges. Result: 13 derive edges + 8 admix edges
# = 21 total (not 13).
test_that("T2.3: read_lgo on rha20.lgo preserves intermediate segments (narrow ghost detection)", {
  result <- read_lgo(path = test_path("fixtures", "rha20.lgo"))
  expect_s3_class(result, "tbl_df")
  expect_setequal(names(result), c("from", "to", "type", "weight"))
  # 13 derive statements (all preserved under narrow rule)
  # + 4 mix statements (each yields 2 admix edges = 8 admix edges)
  # = 21 total edges
  expect_equal(nrow(result), 21)
  expect_setequal(unique(result$type), c("normal", "admix"))
  # 8 admix edges from the 4 mix statements
  admix <- result[result$type == "admix", ]
  expect_equal(nrow(admix), 8)
  # Intermediate segment d2 is preserved (not collapsed to its real parent nd)
  expect_true("d2" %in% c(result$from, result$to))
  # Mixfracs from rha20's published values appear on admix edges
  highs <- admix$weight[admix$weight > 0.5]
  expect_length(highs, 4)
  expect_true(all(highs > 0.95))
})

test_that("T2.3: safe_eval_arith rejects function calls in constrained RHS", {
  # Construct a synthetic .lgo with a malicious constrained RHS
  bad <- "time fixed a=1\ntime constrained b = system(\"echo pwned\")"
  expect_error(
    read_lgo(text = bad),
    class = "legofit_lgo_unsupported"
  )
})

# Review-feedback regression tests --------------------------------------
test_that("compute_node_depths errors on unreachable nodes (cyclic input)", {
  # A→B→A loop with no leaf: every node appears in `from`, so the leaves
  # check inside compute_times will fire first. Construct a graph WITH a
  # leaf but containing a cycle to exercise the depth-walk's unreachable
  # branch via free mode directly:
  cyclic_with_leaf <- tibble::tribble(
    ~from, ~to, ~type,    ~weight,
    "A",   "B", "normal", NA_real_,
    "B",   "A", "normal", NA_real_,  # cycle between A and B
    "A",   "X", "normal", NA_real_   # X is a leaf
  )
  # A is reachable from X (via A→X edge), but B is only reachable via
  # A→B and B is also a parent of A — circular, so depth never settles.
  expect_error(
    compute_node_depths(cyclic_with_leaf),
    class = "legofit_invalid_input"
  )
})

test_that("graph_to_lgo errors when node names collide with ghost-segment names", {
  # Construct a graph where a real node "M_A" exists AND there's an admix
  # event with destination M and parent A — our writer would emit a ghost
  # named "M_A", colliding with the existing node.
  g <- tibble::tribble(
    ~from, ~to,   ~type,    ~weight,
    "R",   "A",   "normal", 0.05,
    "R",   "B",   "normal", 0.05,
    "R",   "M_A", "normal", 0.05,
    "A",   "M",   "admix",  0.60,
    "B",   "M",   "admix",  0.40,
    "M_A", "X1",  "normal", 0.01,
    "M",   "X2",  "normal", 0.01
  )
  expect_error(
    graph_to_lgo(g, validate = FALSE),
    class = "legofit_invalid_input"
  )
})
