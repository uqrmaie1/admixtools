test_that("make_minimal_graph returns expected shape", {
  g <- make_minimal_graph()
  expect_s3_class(g, "tbl_df")
  expect_setequal(names(g), c("from", "to", "type", "weight"))
  expect_equal(nrow(g), 5)
  expect_setequal(unique(g$type), c("normal", "admix"))
})

test_that("make_minimal_igraph returns an igraph", {
  ig <- make_minimal_igraph()
  expect_s3_class(ig, "igraph")
  expect_equal(igraph::ecount(ig), 5)
})
