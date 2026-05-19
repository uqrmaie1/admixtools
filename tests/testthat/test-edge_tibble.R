test_that("validate_edge_tibble errors on missing required columns", {
  bad <- tibble::tibble(from = "A", to = "B")
  expect_error(
    validate_edge_tibble(bad),
    class = "legofit_invalid_input"
  )
})

test_that("validate_edge_tibble errors on invalid type values", {
  bad <- make_minimal_graph()
  bad$type <- "bogus"
  expect_error(
    validate_edge_tibble(bad),
    class = "legofit_invalid_input"
  )
})

test_that("validate_edge_tibble errors when admix invariant violated", {
  bad <- tibble::tibble(
    from = "A", to = "M", type = "admix", weight = 0.6
  )
  expect_error(
    validate_edge_tibble(bad),
    class = "legofit_invalid_input"
  )
})

test_that("coerce_to_edge_tibble accepts a data.frame", {
  df <- as.data.frame(make_minimal_graph())
  result <- coerce_to_edge_tibble(df)
  expect_s3_class(result, "tbl_df")
  expect_setequal(names(result), c("from", "to", "type", "weight"))
})

test_that("coerce_to_edge_tibble accepts an igraph with attributes", {
  ig <- make_minimal_igraph()
  result <- coerce_to_edge_tibble(ig)
  expect_s3_class(result, "tbl_df")
  expect_true("type" %in% names(result))
  expect_true("weight" %in% names(result))
  expect_equal(nrow(result), 5)
})
