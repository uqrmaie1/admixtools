test_that("run_legofit errors clearly when the binary is absent", {
  expect_error(
    run_legofit(make_minimal_graph(), patterns = tempfile(),
                bin = "/nonexistent/legofit"),
    class = "legofit_binary_not_found")
})

test_that("a bare bin name is resolved on PATH, not from the working dir", {
  # A bare command name must go through Sys.which(), even if a like-named
  # file sits in the working directory; otherwise the guard would pass on a
  # non-executable local file that system2() would never actually run.
  wd <- tempfile(); dir.create(wd)
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  file.create("definitely_not_legofit_xyz")
  expect_error(
    run_legofit(make_minimal_graph(), patterns = tempfile(),
                bin = "definitely_not_legofit_xyz"),
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

test_that("run_legofit errors clearly on a non-existent .lgo path", {
  # A length-1 character graph_or_lgo is treated as a path; a missing one
  # should report file-not-found rather than failing inside graph_to_lgo().
  fake <- tempfile(); file.create(fake); Sys.chmod(fake, "0755")
  pat  <- tempfile(); file.create(pat)
  on.exit(unlink(c(fake, pat)))
  expect_error(
    run_legofit("/nonexistent.lgo", patterns = pat, bin = fake),
    class = "legofit_invalid_input")
})

test_that("run_legofit aborts (not hangs) on an unfittable graph in the default time mode", {
  # graph_to_lgo()'s default fix_admix would loop forever on a graph whose
  # branch lengths are NA - an admixture graph with no `time` column, exactly
  # what read_lgo() produces. run_legofit() predicts the NA branch lengths and
  # aborts with guidance toward time_handling = "free" before shelling out.
  fake <- tempfile(); file.create(fake); Sys.chmod(fake, "0755")
  pat  <- tempfile(); file.create(pat)
  on.exit(unlink(c(fake, pat)))
  g <- read_lgo(system.file("extdata/legofit", "rha20.lgo", package = "admixtools"))

  expect_error(
    run_legofit(g, patterns = pat, bin = fake),
    class = "legofit_invalid_input")
})

test_that("time_handling = 'free' skips the unfittable-graph guard", {
  # Under `free` the topology walk is bypassed, so the guard must let the call
  # through; execution then reaches the binary, which here exits non-zero ->
  # legofit_run_failed. (Reaching the run at all is the point: the guard would
  # have aborted earlier with legofit_invalid_input.) Gated off Windows, whose
  # shebang-script exec semantics differ.
  skip_on_os("windows")
  fake <- tempfile(); writeLines(c("#!/bin/sh", "exit 3"), fake); Sys.chmod(fake, "0755")
  pat  <- tempfile(); file.create(pat)
  on.exit(unlink(c(fake, pat)))
  g <- read_lgo(system.file("extdata/legofit", "rha20.lgo", package = "admixtools"))

  expect_error(
    suppressWarnings(
      run_legofit(g, patterns = pat, bin = fake, time_handling = "free")),
    class = "legofit_run_failed")
})

test_that("run_legofit fits the demo graph end-to-end (real binary)", {
  # Happy path: only runs when a real legofit binary is available. Point it
  # at a build via the LEGOFIT_BIN env var, or install legofit on PATH;
  # otherwise this is skipped so the suite stays green without the tool.
  bin <- Sys.getenv("LEGOFIT_BIN", unset = Sys.which("legofit"))
  skip_if(!nzchar(bin) || !file.exists(bin), "legofit binary not available")

  opf <- system.file("extdata/legofit", "demo.opf", package = "admixtools")
  rds <- system.file("extdata/legofit", "demo_fit.rds", package = "admixtools")
  skip_if(!nzchar(opf) || !nzchar(rds), "demo fixtures not installed")

  # The 3-leaf no-admix graph behind the cached demo fixtures.
  demo_graph <- tibble::tribble(
    ~from,  ~to,   ~type,    ~weight,  ~time,
    "xyz",  "xy",  "normal", NA_real_, 1.5,
    "xyz",  "z",   "normal", NA_real_, 2,
    "xy",   "x",   "normal", NA_real_, 0.5,
    "xy",   "y",   "normal", NA_real_, 0.5)

  # Deterministic flags so the fit reproduces the cached result.
  fit <- run_legofit(demo_graph, patterns = opf, bin = bin,
                     args = c("--threads", "1", "-1", "-d", "0"))

  expect_s3_class(fit, "tbl_df")
  expect_equal(nrow(fit), nrow(demo_graph))
  # The read_legofit_output() metadata attributes ride along.
  expect_true(all(c("node_times", "fit_convergence", "identifiability") %in%
                    names(attributes(fit))))
  expect_identical(attr(fit, "fit_convergence")$status, "reached_goal")

  # The live fit matches the cached fixture the vignette ships. The DE
  # optimizer is not bit-reproducible across runs (the split-time optimum
  # sits on a near-flat cost surface; observed run-to-run spread ~6e-3), so
  # this is a "same optimum" check, not an exact-reproduction one.
  cached <- readRDS(rds)
  expect_identical(names(fit), names(cached))
  expect_equal(fit$time, cached$time, tolerance = 2e-2)

  # A .lgo *path* input returns the same edge-tibble shape: read_lgo() recovers
  # the topology under the hood, so graph= need not be supplied. (Previously a
  # path with no graph= yielded the raw param table instead.)
  lgo <- system.file("extdata/legofit", "demo.lgo", package = "admixtools")
  fit_path <- run_legofit(lgo, patterns = opf, bin = bin,
                          args = c("--threads", "1", "-1", "-d", "0"))
  expect_s3_class(fit_path, "tbl_df")
  expect_true(all(c("from", "to", "time") %in% names(fit_path)))
  expect_false(is.null(attr(fit_path, "node_times")))
})
