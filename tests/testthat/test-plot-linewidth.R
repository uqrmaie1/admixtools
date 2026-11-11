# ggplot2 >= 3.4.0 deprecates `size` for line widths in favour of `linewidth`
# (#175). testthat 3e reports every lifecycle deprecation, so the old
# `size = 1` on the highlight layer fails here; the width must stay at 1.

test_that("plot_graph draws the unidentifiable-edge highlight with linewidth, not size", {
  expect_no_warning(
    p <- plot_graph(example_igraph, highlight_unidentifiable = TRUE),
    message = "`size` aesthetic", class = "lifecycle_warning_deprecated"
  )
  hl <- ggplot2::layer_data(p, length(p$layers))
  expect_true(all(hl$colour == "red"))
  expect_true(all(hl$linewidth == 1))
})
