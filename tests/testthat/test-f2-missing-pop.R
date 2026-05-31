# The f2-cache missing-population error is a typed condition
# (admixtools_pop_missing_in_f2_cache), thrown by get_f2() when a requested
# population is absent from the f2 data. These exercise the class contract so
# a future message edit or refactor cannot silently drop the typed class.

test_that("get_f2 errors with admixtools_pop_missing_in_f2_cache for an absent population", {
  f2b <- array(0.1, dim = c(2, 2, 1),
               dimnames = list(c("A", "B"), c("A", "B"), "l1"))
  err <- tryCatch(get_f2(f2b, pops = c("A", "ZZZ")), error = function(e) e)
  expect_s3_class(err, "admixtools_pop_missing_in_f2_cache")
  expect_match(conditionMessage(err), "ZZZ")   # message names the missing pop
})

test_that("get_f2 does not error when all requested populations are present", {
  f2b <- array(0.1, dim = c(2, 2, 1),
               dimnames = list(c("A", "B"), c("A", "B"), "l1"))
  expect_no_error(get_f2(f2b, pops = c("A", "B")))
})
