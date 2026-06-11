test_that("run_legofit errors clearly when the binary is absent", {
  expect_error(
    run_legofit(make_minimal_graph(), patterns = tempfile(),
                bin = "/nonexistent/legofit"),
    class = "legofit_binary_not_found")
})

test_that("run_legofit errors when the pattern file is missing", {
  # A fake bin that exists so the first guard passes.
  fake <- tempfile(); file.create(fake); Sys.chmod(fake, "0755")
  on.exit(unlink(fake))
  expect_error(
    run_legofit(make_minimal_graph(), patterns = "/nonexistent.opf", bin = fake),
    class = "legofit_invalid_input")
})
