# F-1 structural identifiability flag.
# attr(result, "identifiability") classifies each parameter without ever
# consulting bootstrap CI width.

write_legofit <- function(free_lines, fixed_lines = "       one = 1") {
  f <- tempfile(fileext = ".legofit")
  writeLines(c(
    "Initial parameter values", "Fixed:", fixed_lines,
    "Free:", free_lines,
    "DiffEv reached_goal. cost=1e-16 spread=3e-07",
    "Fitted parameter values", "Free:", free_lines
  ), f)
  f
}
idf <- function(r) {
  t <- attr(r, "identifiability")
  setNames(t$identifiability, t$parameter)
}

test_that("F-1: admix-event times are flagged structural_none, mixFrac is not", {
  f <- write_legofit(c("     T_R = 2", "     T_admix_m = 1", "      m_m = 0.3"))
  r <- suppressMessages(read_legofit_output(f))
  cls <- idf(r)
  expect_equal(cls[["T_admix_m"]], "structural_none")   # single-lineage segment
  expect_equal(cls[["m_m"]],       "identifiable")       # mixFrac always identifiable
  expect_equal(cls[["T_R"]],       "identifiable")       # free time, twoN fixed -> fine
  expect_equal(cls[["one"]],       "fixed")
})

test_that("F-1: free times + free twoN trigger scale_degenerate (F-4)", {
  f <- write_legofit(
    free_lines  = c("     T_R = 2", "  twoN_R = 1000", "  twoN_A = 500"),
    fixed_lines = "       T_x = 0")
  cls <- idf(suppressMessages(read_legofit_output(f)))
  expect_equal(cls[["T_R"]],    "scale_degenerate")
  expect_equal(cls[["twoN_R"]], "scale_degenerate")
  expect_equal(cls[["T_x"]],    "fixed")
})

test_that("F-1: default one=1 (twoN fixed) keeps tree times identifiable", {
  f <- write_legofit(c("     T_R = 2", "      T_S = 1"))
  cls <- idf(suppressMessages(read_legofit_output(f)))
  expect_equal(unname(cls[c("T_R", "T_S")]), c("identifiable", "identifiable"))
  expect_false(any(cls == "scale_degenerate"))
})

test_that("F-1b: deep no-leaf-child split is weak_unconstrained; shallow stays identifiable", {
  # Root R's children are both internal (no attached leaf) -> stricter thresholds
  # (well<=4, unconstrained>6). T_R at coalescent depth 10 (one=1) -> unconstrained.
  # A and B each have leaf children and are shallow (D=2) -> identifiable.
  g <- tibble::tribble(
    ~from, ~to, ~type,     ~weight,
    "R", "A", "normal", NA_real_,
    "R", "B", "normal", NA_real_,
    "A", "x", "normal", NA_real_,
    "A", "y", "normal", NA_real_,
    "B", "z", "normal", NA_real_,
    "B", "w", "normal", NA_real_)
  f <- write_legofit(c("     T_R = 10", "      T_A = 2", "      T_B = 2"))
  cls <- idf(suppressMessages(read_legofit_output(f, graph = g)))
  expect_equal(cls[["T_R"]], "weak_unconstrained")
  expect_equal(cls[["T_A"]], "identifiable")
  expect_equal(cls[["T_B"]], "identifiable")
})

test_that("F-1b: a leaf-child split holds deeper (weak_identified at D=7, not unconstrained)", {
  # R has a directly attached leaf child (x) -> looser thresholds (well<=5,
  # weak 5<D<=9). T_R at depth 7 -> weak_identified, not unconstrained.
  g <- tibble::tribble(
    ~from, ~to, ~type,     ~weight,
    "R", "x", "normal", NA_real_,
    "R", "A", "normal", NA_real_,
    "A", "y", "normal", NA_real_,
    "A", "z", "normal", NA_real_)
  f <- write_legofit(c("     T_R = 7", "      T_A = 2"))
  cls <- idf(suppressMessages(read_legofit_output(f, graph = g)))
  expect_equal(cls[["T_R"]], "weak_identified")
  expect_equal(cls[["T_A"]], "identifiable")
})

test_that("F-1: identifiability attaches on both graph and graph=NULL paths", {
  g <- make_ourex1_graph()
  r_graph <- suppressMessages(
    read_legofit_output(test_path("fixtures", "ourex1.legofit"), graph = g))
  r_raw   <- suppressMessages(
    read_legofit_output(test_path("fixtures", "ourex1.legofit")))
  expect_s3_class(attr(r_graph, "identifiability"), "tbl_df")
  expect_s3_class(attr(r_raw,   "identifiability"), "tbl_df")
  expect_setequal(names(attr(r_graph, "identifiability")),
                  c("parameter", "family", "identifiability", "reason"))
})
