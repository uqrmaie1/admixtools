# Tests for classify_legofit_param() (Step 4)
# and the write-time admix_<X> collision guard in generate_param_names.

# ---------------------------------------------------------------------------
# classify_legofit_param — vectorized pattern matching
# ---------------------------------------------------------------------------

test_that("classify_legofit_param: admix_time (T_admix_*) before time (T_*)", {
  # Order of pattern application matters:
  # T_admix_M must classify as admix_time, NOT as time for node "admix_M"
  expect_equal(classify_legofit_param("T_admix_M"), "admix_time")
  expect_equal(classify_legofit_param("T_admix_XY"), "admix_time")
})

test_that("classify_legofit_param: time (T_*)", {
  expect_equal(classify_legofit_param("T_R"),   "time")
  expect_equal(classify_legofit_param("T_A"),   "time")
  expect_equal(classify_legofit_param("T_xyz"), "time")
  # Names with dots and colons (legal per parse.c:248-254)
  expect_equal(classify_legofit_param("T_A.B"), "time")
  expect_equal(classify_legofit_param("T_A:B"), "time")
})

test_that("classify_legofit_param: twoN (twoN_*)", {
  expect_equal(classify_legofit_param("twoN_R"),   "twoN")
  expect_equal(classify_legofit_param("twoN_xyz"), "twoN")
})

test_that("classify_legofit_param: twoN sentinels (one, shared)", {
  expect_equal(classify_legofit_param("one"),    "twoN_one")
  expect_equal(classify_legofit_param("shared"), "twoN_shared")
})

test_that("classify_legofit_param: mixFrac (m_*)", {
  expect_equal(classify_legofit_param("m_M"),   "mixFrac")
  expect_equal(classify_legofit_param("m_XY"),  "mixFrac")
})

test_that("classify_legofit_param: unknown patterns", {
  expect_equal(classify_legofit_param("bogus"),  "unknown")
  expect_equal(classify_legofit_param("zero"),   "unknown")
  expect_equal(classify_legofit_param("T"),      "unknown")  # no underscore suffix
  expect_equal(classify_legofit_param("twoN"),   "unknown")  # no underscore suffix
})

test_that("classify_legofit_param: vectorized", {
  names <- c("T_admix_M", "T_R", "twoN_R", "one", "shared", "m_M", "BOGUS")
  expected <- c("admix_time", "time", "twoN", "twoN_one", "twoN_shared", "mixFrac", "unknown")
  expect_equal(classify_legofit_param(names), expected)
})

test_that("classify_legofit_param: zero-length input", {
  expect_equal(classify_legofit_param(character(0)), character(0))
})

# ---------------------------------------------------------------------------
# write-time admix_<X> collision guard in generate_param_names
# ---------------------------------------------------------------------------

test_that("collision guard: node admix_M where M is an admix dest aborts", {
  bad <- tibble::tribble(
    ~from,      ~to,       ~type,    ~weight,
    "R",        "admix_M", "normal",   0.05,
    "R",        "A",       "normal",   0.05,
    "A",        "M",       "admix",    0.60,
    "admix_M",  "M",       "admix",    0.40,
    "M",        "X",       "normal",   0.02
  )
  expect_error(
    graph_to_lgo(bad, validate = FALSE),
    class = "legofit_invalid_input"
  )
})

test_that("collision guard: node admix_X where X is NOT an admix dest is OK", {
  # admix_X exists as a node, but X is not an admix destination — no ambiguity
  ok <- tibble::tribble(
    ~from,     ~to,      ~type,    ~weight,
    "R",       "admix_X", "normal",  0.05,
    "R",       "Y",       "normal",  0.05,
    "admix_X", "Z",       "normal",  0.02,
    "Y",       "Z",       "normal",  0.02
  )
  # This should pass the collision guard (admix_X node exists but Z is not an
  # admix dest because there are no admix-typed edges).
  # graph_to_lgo will fail the validate_edge_tibble (Z has 2 incoming normal
  # edges, which is fine actually? no, each `to` with type normal can have 1).
  # Let's just test the collision guard passes for this case.
  # Re-build a simpler example to isolate the guard.
  ok2 <- tibble::tribble(
    ~from,     ~to,      ~type,    ~weight,
    "R",       "admix_X", "normal",  0.05,
    "admix_X", "Y",       "normal",  0.05
  )
  # No admix-typed edges → admix_dests is empty → no collision → guard passes
  expect_no_error(admixtools:::generate_param_names(ok2))
})
