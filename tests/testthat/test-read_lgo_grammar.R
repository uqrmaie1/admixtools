# Tests for Phase 2 read_lgo grammar additions:
#   Step 1: line continuation
#   Step 2: bounded `free` declarations
#   Step 3: as = "igraph" return type

# ---------------------------------------------------------------------------
# Step 1 — line continuation (tests 14, 15, 18)
# ---------------------------------------------------------------------------

test_that("T14: multi-line constrained expression evaluates correctly", {
  lgo <- "
time fixed Ta=2
time fixed Tb=3
time fixed Tc=1
time constrained Tfoo = 0.5 * (Ta + Tb +
                Tc)
segment foo t=Tfoo twoN=one
derive foo from bar
segment bar t=Ta twoN=one
"
  result <- read_lgo(text = lgo)
  expect_s3_class(result, "tbl_df")
  # Tfoo = 0.5 * (2 + 3 + 1) = 3 — not in the edge tibble directly,
  # but the parse must succeed (no error) to get here.
  expect_equal(nrow(result), 1L)
})

test_that("T15a: continuation triggers on +", {
  lgo <- "
time fixed Ta=1
time fixed Tb=2
time constrained Tc = Ta +
Tb
segment x t=Tc twoN=one
derive x from y
segment y t=Ta twoN=one
"
  expect_no_error(read_lgo(text = lgo))
})

test_that("T15b: continuation triggers on -", {
  lgo <- "
time fixed Ta=3
time fixed Tb=1
time constrained Tc = Ta -
Tb
segment x t=Tc twoN=one
derive x from y
segment y t=Ta twoN=one
"
  expect_no_error(read_lgo(text = lgo))
})

test_that("T15c: continuation triggers on *", {
  lgo <- "
time fixed Ta=2
time fixed Tb=4
time constrained Tc = Ta *
Tb
segment x t=Tc twoN=one
derive x from y
segment y t=Ta twoN=one
"
  expect_no_error(read_lgo(text = lgo))
})

test_that("T15d: continuation triggers on /", {
  lgo <- "
time fixed Ta=8
time fixed Tb=2
time constrained Tc = Ta /
Tb
segment x t=Tc twoN=one
derive x from y
segment y t=Ta twoN=one
"
  expect_no_error(read_lgo(text = lgo))
})

test_that("T15e: ^ does NOT trigger continuation", {
  # A line ending with ^ would not be joined; the ^ would be treated as
  # the start of an unsupported construct or a value, causing an error
  # downstream. We only verify that ^ alone does NOT cause silent joining.
  # We test this by verifying join_continuations leaves ^ lines alone.
  lgo_lines <- c("time fixed Ta=2^", "3")  # Not continuation; two separate lines
  # join_continuations should NOT join these
  joined <- admixtools:::join_continuations(lgo_lines)
  expect_equal(length(joined), 2L)
})

test_that("T15f: = does NOT trigger continuation", {
  lgo_lines <- c("time fixed Ta=", "2")  # should NOT join
  joined <- admixtools:::join_continuations(lgo_lines)
  expect_equal(length(joined), 2L)
})

test_that("T18: EOF mid-continuation aborts with legofit_lgo_unsupported", {
  expect_error(
    admixtools:::join_continuations(c("time fixed Ta = 2 +")),
    class = "legofit_lgo_unsupported"
  )
})

# ---------------------------------------------------------------------------
# Step 2 — bounded free declarations (tests 16, 17, 19)
# ---------------------------------------------------------------------------

test_that("T16: bounded free captures bounds before name", {
  lgo <- "
time free [0.1, 1.0] T_R = 0.5
twoN fixed one=1
segment R t=T_R twoN=one
derive A from R
segment A t=T_R twoN=one
"
  result <- read_lgo(text = lgo)
  bounds <- attr(result, "param_bounds")
  expect_type(bounds, "list")
  expect_named(bounds, "T_R")
  expect_equal(bounds$T_R, c(lo = 0.1, hi = 1.0))
})

test_that("T16b: bounds tolerate arbitrary internal whitespace", {
  lgo <- "
time free [ 0.1 , 1.0 ] T_R = 0.5
twoN fixed one=1
segment R t=T_R twoN=one
derive A from R
segment A t=T_R twoN=one
"
  result <- read_lgo(text = lgo)
  bounds <- attr(result, "param_bounds")
  expect_equal(bounds$T_R, c(lo = 0.1, hi = 1.0))
})

test_that("T17: bounds on fixed declaration errors with legofit_lgo_unsupported", {
  # per parse.c:237-242 — only `free` accepts bounds
  lgo <- "time fixed [0.1, 1.0] T_R=0"
  expect_error(
    read_lgo(text = lgo),
    class = "legofit_lgo_unsupported"
  )
})

test_that("T17b: bounds on constrained declaration errors", {
  lgo <- "
time fixed Ta=1
time constrained [0.1, 1.0] T_R = Ta
"
  expect_error(
    read_lgo(text = lgo),
    class = "legofit_lgo_unsupported"
  )
})

test_that("T19: malformed bounded free (unclosed bracket) errors", {
  # The PARAM_RE requires a closing ] before the name.
  # Unclosed bracket means the regex won't match, triggering abort.
  lgo <- "time free [0.1, T_R=0.5"
  expect_error(
    read_lgo(text = lgo),
    class = "legofit_lgo_unsupported"
  )
})

# ---------------------------------------------------------------------------
# Step 3 — as = "igraph" return type (test 20)
# ---------------------------------------------------------------------------

test_that("T20: as='igraph' returns igraph with same topology as as='edges'", {
  lgo <- "
time fixed T_x=0
time fixed T_y=0
time free T_R=2
twoN fixed one=1
segment R t=T_R twoN=one
segment x t=T_x twoN=one samples=1
segment y t=T_y twoN=one samples=1
derive x from R
derive y from R
"
  edges_result <- read_lgo(text = lgo, as = "edges")
  igraph_result <- read_lgo(text = lgo, as = "igraph")

  expect_true(igraph::is_igraph(igraph_result))
  expect_equal(igraph::ecount(igraph_result), nrow(edges_result))

  # Node names match (order may differ)
  expect_setequal(
    igraph::V(igraph_result)$name,
    unique(c(edges_result$from, edges_result$to))
  )

  # Edge attributes present
  expect_true("type"   %in% igraph::edge_attr_names(igraph_result))
  expect_true("weight" %in% igraph::edge_attr_names(igraph_result))
})

test_that("T20b: as='igraph' preserves param_bounds as graph attribute", {
  lgo <- "
time free [0.1, 5.0] T_R = 2
twoN fixed one=1
segment R t=T_R twoN=one
derive A from R
segment A t=T_R twoN=one
"
  ig <- read_lgo(text = lgo, as = "igraph")
  bounds <- igraph::graph_attr(ig, "param_bounds")
  expect_type(bounds, "list")
  expect_equal(bounds$T_R, c(lo = 0.1, hi = 5.0))
})

# ---------------------------------------------------------------------------
# Regression tests for code-review fixes
# ---------------------------------------------------------------------------

# Fix 5: parse() syntax error in constrained RHS → legofit_lgo_unsupported
test_that("Fix5: invalid constrained expression aborts legofit_lgo_unsupported", {
  lgo <- "
time free T_A = 1.0
time constrained T_B = 1 +* 2
segment A twoN=one samples=1
segment B twoN=one samples=1
derive B from A
"
  expect_error(read_lgo(text = lgo), class = "legofit_lgo_unsupported")
})

# Fix 3: regex-metachar in dest name does not cause false ghost collapse
test_that("Fix3: admix dest with period does not falsely collapse real segments", {
  # Dest is "M.east" (period in name), real segment is "M_east_root" — the
  # unescaped dot would previously match M_east via the ghost pattern.
  lgo <- "
time  free T_M.east   = 1.0
time  free T_root     = 2.0
time  free T_Mxeast   = 0.8
mixFrac free m_M.east = 0.3
twoN  fixed one = 1
segment root      t=T_root   twoN=one samples=1
segment M.east    t=T_M.east twoN=one samples=1
segment Mxeast    t=T_Mxeast twoN=one samples=1
derive M.east  from root
derive Mxeast  from root
mix    M.east  from root + m_M.east * Mxeast
"
  result <- read_lgo(text = lgo)
  # Mxeast should appear as a real edge, not be collapsed as a ghost
  expect_true("Mxeast" %in% c(result$from, result$to))
})

# Fix 9: as="igraph" on an edge-less .lgo returns an igraph not a tibble
test_that("Fix9: empty-edge lgo with as='igraph' returns igraph", {
  lgo <- "
time free T_A = 1.0
twoN fixed one = 1
segment A t=T_A twoN=one samples=1
"
  result <- read_lgo(text = lgo, as = "igraph")
  expect_true(igraph::is_igraph(result))
})

# ---------------------------------------------------------------------------
# Master plan Tier 0 gap closure
# ---------------------------------------------------------------------------

# U-G4: bounds present but non-numeric. PARAM_RE matches the bracket, so the
# numeric coercion of lo/hi is what must reject it (parse.c:211-235 requires
# numeric bounds). Distinct from T17 (bounds on a non-free declaration).
test_that("U-G4: non-numeric bounds on free declaration abort legofit_lgo_unsupported", {
  lgo <- "
time free [a, b] T_R = 0.5
twoN fixed one = 1
segment R t=T_R twoN=one
"
  expect_error(
    read_lgo(text = lgo),
    class = "legofit_lgo_unsupported"
  )
})

# U-G5: the generic `param` keyword (parse.c:693) is a fourth declaration type
# beyond time/twoN/mixFrac. The pre-PR parser could not handle it and would
# abort as an unsupported construct; PARAM_RE now accepts it. The declared
# param is inert (not referenced by a segment) so the edge tibble is unaffected.
test_that("U-G5: generic `param` declaration parses and edges still build", {
  lgo <- "
param free pmix = 0.5
time fixed T_x = 0
time fixed T_y = 0
time free  T_R = 2
twoN fixed one = 1
segment R t=T_R twoN=one
segment x t=T_x twoN=one samples=1
segment y t=T_y twoN=one samples=1
derive x from R
derive y from R
"
  result <- read_lgo(text = lgo)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2L)
  expect_setequal(c(result$from, result$to), c("R", "x", "y"))
})
