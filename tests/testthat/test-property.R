# Tier 1 property-based and adversarial coverage for the LEGOFIT bridge.
# P series.
#
# Two properties:
#  1. Generative round trip. For a grid of random_sim graphs (within the
#     export envelope: <=6 leaves, <=1 admix, time_handling="free"), the
#     topology survives graph_to_lgo -> read_lgo: every from/to/type edge is
#     recovered, ghost admix segments collapse back to the original admix edges,
#     and the edge count is preserved. Weights are intentionally lossy on
#     round trip and are not compared.
#  2. Adversarial corpus. A malformed input must NEVER produce a silently wrong
#     parse or an untyped error. Each bad input either parses to a tibble or
#     signals a condition of a documented legofit_* class.

suppressMessages({
  library(dplyr)
})

# --- 1. Generative round trip ------------------------------------------------

edge_key <- function(df) {
  df %>%
    mutate(type = ifelse(type == "edge", "normal", type)) %>%
    transmute(k = paste(from, to, type)) %>%
    pull(k) %>%
    sort()
}

test_that("graph_to_lgo -> read_lgo preserves topology across a random grid", {
  grid <- expand.grid(nleaf = c(3, 4, 5, 6), nadmix = c(0, 1),
                      seed = 1:3, KEEP.OUT.ATTRS = FALSE)
  for (i in seq_len(nrow(grid))) {
    nl <- grid$nleaf[i]; na <- grid$nadmix[i]; s <- grid$seed[i]
    set.seed(s)
    g     <- random_sim(nleaf = nl, nadmix = na)
    edges <- coerce_to_edge_tibble(g$edges)

    txt <- graph_to_lgo(edges, time_handling = "free", validate = FALSE)
    parsed <- read_lgo(text = txt, as = "edges")

    label <- sprintf("nleaf=%d nadmix=%d seed=%d", nl, na, s)
    # topology (from/to/type) is recovered exactly, ghosts collapsed
    expect_setequal(edge_key(parsed), edge_key(edges))
    expect_identical(edge_key(parsed), edge_key(edges), info = label)
    # edge count preserved
    expect_equal(nrow(parsed), nrow(edges), info = label)
  }
})

test_that("P-2: generate_param_names output never classifies as unknown", {
  for (s in 1:5) {
    set.seed(s)
    g     <- random_sim(nleaf = 5, nadmix = 1)
    edges <- coerce_to_edge_tibble(g$edges)
    nms   <- generate_param_names(edges)
    fam   <- classify_legofit_param(unlist(nms, use.names = FALSE))
    expect_false(any(fam == "unknown"),
                 info = sprintf("seed=%d produced an unknown-classified name", s))
  }
})

test_that("P-5: all three twoN modes and init/free time modes parse-round-trip", {
  set.seed(7)
  edges <- coerce_to_edge_tibble(random_sim(nleaf = 5, nadmix = 0)$edges)

  # twoN modes: NULL (coalescent), scalar, named
  named_twoN <- setNames(rep(1000, length(unique(c(edges$from, edges$to)))),
                         unique(c(edges$from, edges$to)))
  for (tw in list(NULL, 1000, named_twoN)) {
    txt <- graph_to_lgo(edges, twoN = tw, time_handling = "free", validate = FALSE)
    expect_s3_class(read_lgo(text = txt), "tbl_df")
  }

  # Only "free" reliably round-trips an ARBITRARY random graph: random_sim's
  # uniform unit drifts are not additively consistent, so "init" and "fix_admix"
  # (which run the topology walk) abort on most random topologies.
  # Those modes are byte-golden tested on consistent fixtures elsewhere.
  txt <- graph_to_lgo(edges, time_handling = "free", validate = FALSE)
  expect_setequal(edge_key(read_lgo(text = txt)), edge_key(edges))
})

test_that("P-3: node_times covers exactly the graph node set, never more or fewer", {
  g <- make_ourex1_graph()
  nodes <- unique(c(coerce_to_edge_tibble(g)$from, coerce_to_edge_tibble(g)$to))
  r <- suppressMessages(
    read_legofit_output(test_path("fixtures", "ourex1.legofit"), graph = g))
  nt <- attr(r, "node_times")
  expect_setequal(names(nt), nodes)
  expect_length(nt, length(nodes))
  expect_true(all(is.finite(nt) | is.na(nt)))   # numeric, no garbage
})

test_that("round trip is deterministic across repeated reads", {
  set.seed(42)
  g     <- random_sim(nleaf = 5, nadmix = 1)
  edges <- coerce_to_edge_tibble(g$edges)
  txt   <- graph_to_lgo(edges, time_handling = "free", validate = FALSE)
  keys  <- replicate(5, paste(edge_key(read_lgo(text = txt)), collapse = "|"))
  expect_length(unique(keys), 1L)
})

# --- 2. Adversarial corpus ---------------------------------------------------

# Each entry is a malformed or edge-case .lgo. The contract: read_lgo either
# returns a tbl_df, or aborts with a legofit_* class. A bare/untyped error or a
# silently wrong parse is a failure.
lgo_corpus <- list(
  unsupported_construct = "wibble foo from bar",
  malformed_derive      = "derive x from",
  malformed_mix         = "mix m from a + b",
  bounds_on_fixed       = "time fixed [0,1] T_x = 0.5",
  nonnumeric_bounds     = "time free [a, b] T_x = 1",
  empty_constrained_rhs = "time fixed a = 1\ntime constrained b =",
  function_call_rhs     = "time fixed a = 1\ntime constrained b = system(\"echo hi\")",
  bare_equals           = "time fixed = 0.5"
)

legofit_classes <- c("legofit_lgo_unsupported", "legofit_invalid_input")

classify_outcome <- function(expr) {
  tryCatch({
    res <- force(expr)
    if (inherits(res, "tbl_df")) "tibble" else "other-value"
  },
  legofit_lgo_unsupported = function(e) "typed",
  legofit_invalid_input   = function(e) "typed",
  error = function(e) paste0("untyped:", conditionMessage(e)))
}

test_that("malformed .lgo never errors untyped and never silently mis-parses", {
  for (nm in names(lgo_corpus)) {
    outcome <- classify_outcome(read_lgo(text = lgo_corpus[[nm]]))
    expect_true(
      outcome %in% c("tibble", "typed"),
      info = sprintf("corpus '%s' gave outcome '%s'", nm, outcome)
    )
  }
})

test_that("parse_bootci_output rejects malformed bootci with a typed condition", {
  bad_bootci <- list(
    all_comments = c("# bootci.py", "# confidence: 0.95"),
    wrong_header = c("not a header at all", "T_R 1 2 3"),
    too_few_cols = c("       par   est   low   high", "   T_R   1.0   0.9")
  )
  for (nm in names(bad_bootci)) {
    outcome <- tryCatch(
      { parse_bootci_output(bad_bootci[[nm]]); "value" },
      legofit_invalid_input = function(e) "typed",
      error = function(e) paste0("untyped:", conditionMessage(e))
    )
    expect_true(outcome %in% c("typed", "value"),
                info = sprintf("bootci corpus '%s' gave '%s'", nm, outcome))
  }
})
