test_that("assemble_lgo -> read_lgo preserves topology", {
  g     <- make_minimal_graph()
  edges <- coerce_to_edge_tibble(g)
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

  parsed <- read_lgo(text = txt, as = "edges")
  # read_lgo attaches a nodes tibble; the topology check compares
  # from/to/type only, so drop it before the tibble comparison.
  attr(parsed, "nodes") <- NULL

  g_topology <- g %>% dplyr::select(from, to, type) %>%
    dplyr::arrange(from, to)
  p_topology <- parsed %>% dplyr::select(from, to, type) %>%
    dplyr::arrange(from, to)
  expect_equal(p_topology, g_topology)
})

test_that("read_lgo parses the golden file without error", {
  result <- read_lgo(path = testthat::test_path("fixtures", "minimal.lgo"))
  expect_s3_class(result, "tbl_df")
  expect_setequal(names(result), c("from", "to", "type", "weight"))
  expect_equal(nrow(result), 5)
})

test_that("read_lgo errors when both path and text supplied", {
  expect_error(
    read_lgo(path = "x.lgo", text = "some text"),
    class = "legofit_invalid_input"
  )
})

test_that("read_lgo errors when neither path nor text supplied", {
  expect_error(
    read_lgo(),
    class = "legofit_invalid_input"
  )
})
