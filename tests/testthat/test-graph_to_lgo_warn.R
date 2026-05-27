test_that("graph_to_lgo warns when fewer than 2 leaves receive samples", {
  # make_minimal_graph() has only one leaf (X), so scalar samples = 1 places
  # `samples=1` on a single segment. LEGOFIT requires >= 2 sampled segments
  # to produce observable site patterns, so emit a structured warning.
  expect_warning(
    graph_to_lgo(make_minimal_graph(), validate = FALSE),
    class = "legofit_unfittable_lgo"
  )
})

test_that("legofit_unfittable_lgo warning names the sampled segments", {
  w <- tryCatch(
    graph_to_lgo(make_minimal_graph(), validate = FALSE),
    legofit_unfittable_lgo = function(c) c
  )
  expect_s3_class(w, "legofit_unfittable_lgo")
  msg <- paste(w$message, collapse = "\n")
  expect_match(msg, "\\bX\\b")
  expect_match(msg, "Sampled segments \\(1\\)")
})

test_that("graph_to_lgo counts only positive samples (samples=0 is unfittable)", {
  # Two leaves would normally be fittable, but samples = 0 on every leaf
  # yields zero observable site patterns. The count must exclude zeros, so
  # the warning still fires and reports no sampled segments.
  e <- tibble::tibble(from = c("R", "R"), to = c("A", "B"))
  w <- tryCatch(
    graph_to_lgo(e, samples = c(A = 0, B = 0), validate = FALSE),
    legofit_unfittable_lgo = function(c) c
  )
  expect_s3_class(w, "legofit_unfittable_lgo")
  expect_match(paste(w$message, collapse = "\n"),
               "Sampled segments (0): none", fixed = TRUE)
})

test_that("graph_to_lgo is silent when at least 2 leaves are sampled", {
  e <- tibble::tibble(from = c("R", "R"), to = c("A", "B"))
  expect_warning(
    graph_to_lgo(e, validate = FALSE),
    class = "legofit_unfittable_lgo",
    regexp = NA
  )
})
