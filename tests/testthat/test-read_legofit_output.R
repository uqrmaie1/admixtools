# Tests for read_legofit_output() Phase B (LLD §3, §7.2 cases 1-9)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

lgo_path     <- function(...) test_path("fixtures", ...)
ourex1_path  <- function() lgo_path("ourex1.legofit")
ourex1_graph <- function() make_ourex1_graph()

# ---------------------------------------------------------------------------
# extract_param_section (Step 6 helper)
# ---------------------------------------------------------------------------

test_that("extract_param_section: Fitted block returns free lines only", {
  lines <- readLines(ourex1_path(), encoding = "UTF-8")
  sec   <- extract_param_section(lines, "Fitted parameter values")
  expect_equal(sec$fixed, character(0))
  expect_length(sec$free, 2)   # T_xyz and T_xy
})

test_that("extract_param_section: Initial block returns both fixed and free", {
  lines <- readLines(ourex1_path(), encoding = "UTF-8")
  sec   <- extract_param_section(lines, "Initial parameter values")
  expect_length(sec$fixed, 4)  # T_z, T_x, T_y, one
  expect_length(sec$free,  2)  # T_xyz, T_xy
})

test_that("extract_param_section: missing header returns empty lists", {
  sec <- extract_param_section(c("some text", "other text"), "Nonexistent header")
  expect_equal(sec$fixed, character(0))
  expect_equal(sec$free,  character(0))
})

# ---------------------------------------------------------------------------
# parse_param_lines (Step 6 helper)
# ---------------------------------------------------------------------------

test_that("parse_param_lines: returns tibble with name/value", {
  lines <- c("     T_xyz = 2", "      T_xy = 0.5")
  result <- parse_param_lines(lines)
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("name", "value"))
  expect_equal(result$name,  c("T_xyz", "T_xy"))
  expect_equal(result$value, c(2, 0.5))
})

test_that("parse_param_lines: handles no-space-around-= form", {
  result <- parse_param_lines(c("T_R=0.5"))
  expect_equal(result$name,  "T_R")
  expect_equal(result$value, 0.5)
})

test_that("parse_param_lines: empty input returns zero-row tibble", {
  result <- parse_param_lines(character(0))
  expect_equal(nrow(result), 0L)
  expect_named(result, c("name", "value"))
})

# ---------------------------------------------------------------------------
# extract_convergence_status (Step 6 helper)
# ---------------------------------------------------------------------------

test_that("extract_convergence_status: reached_goal parsed correctly", {
  lines <- readLines(ourex1_path(), encoding = "UTF-8")
  conv  <- extract_convergence_status(lines)
  expect_equal(conv$status, "reached_goal")
  expect_true(is.numeric(conv$cost))
  expect_true(is.numeric(conv$spread))
  expect_true(conv$cost   < 1e-10)
  expect_true(conv$spread < 1e-5)
})

test_that("extract_convergence_status: NA when no DiffEv line", {
  conv <- extract_convergence_status(c("some", "lines", "without", "difdev"))
  expect_true(is.na(conv$status))
  expect_true(is.na(conv$cost))
  expect_true(is.na(conv$spread))
})

# ---------------------------------------------------------------------------
# inform_convergence_status (Step 6 helper)
# ---------------------------------------------------------------------------

test_that("inform_convergence_status: no message on reached_goal", {
  expect_no_message(
    inform_convergence_status(list(status = "reached_goal", cost = 1e-16, spread = 1e-7))
  )
})

test_that("inform_convergence_status: legofit_fit_incomplete on finished_iterations", {
  expect_message(
    inform_convergence_status(list(status = "finished_iterations", cost = 4.32e-3, spread = 1.05e-1)),
    class = "legofit_fit_incomplete"
  )
})

# ---------------------------------------------------------------------------
# read_legofit_output — end-to-end (LLD §7.2 cases 1-9)
# ---------------------------------------------------------------------------

# Case 1: smoke — parses without error
test_that("T1: smoke — read_legofit_output parses ourex1.legofit", {
  expect_no_error(read_legofit_output(ourex1_path()))
})

# Case 2: returns expected columns when graph supplied
test_that("T2: returns expected columns with graph", {
  result <- read_legofit_output(ourex1_path(), graph = ourex1_graph())
  expected_cols <- c("from", "to", "type", "weight", "time",
                     "admix_event_time", "twoN", "admix_prop", "admix_event")
  expect_true(all(expected_cols %in% names(result)))
})

# Case 3: fitted times correct — T_xy = 0.5, T_xyz = 2 (from ourex1)
test_that("T3: fitted times attach to correct edges", {
  result <- read_legofit_output(ourex1_path(), graph = ourex1_graph())
  # T_xy = 0.5 is the time of segment xy; edges with to=xy should have time~0.5
  xy_edges <- result[result$to == "xy", ]
  expect_true(nrow(xy_edges) > 0)
  expect_true(all(abs(xy_edges$time - 0.5) < 1e-4))
})

# Case 4: node_times attribute contains root time
test_that("T4: node_times attribute contains all nodes including root", {
  result <- read_legofit_output(ourex1_path(), graph = ourex1_graph())
  nt <- attr(result, "node_times")
  expect_type(nt, "double")
  expect_true("xyz" %in% names(nt))
  expect_equal(nt[["xyz"]], 2, tolerance = 1e-4)
  expect_equal(nt[["xy"]],  0.5, tolerance = 1e-4)
})

# Case 5: no admix graph — admix_event_time column is all NA
test_that("T5: admix_event_time is NA for no-admix graph", {
  result <- read_legofit_output(ourex1_path(), graph = ourex1_graph())
  expect_true(all(is.na(result$admix_event_time)))
})

# Case 6: graph = NULL returns raw param table with exactly name/value/family
test_that("T6: graph=NULL returns raw parameter table", {
  result <- read_legofit_output(ourex1_path(), graph = NULL)
  expect_s3_class(result, "tbl_df")
  # Public API: exactly these three columns (no internal source/free cols)
  expect_setequal(names(result), c("name", "value", "family"))
  # Should have T_xyz, T_xy among the free-fitted params
  expect_true("T_xyz" %in% result$name)
  expect_true("T_xy"  %in% result$name)
})

# Case 7: include_fixed = FALSE excludes fixed params
test_that("T7: include_fixed=FALSE excludes fixed params", {
  full    <- read_legofit_output(ourex1_path(), graph = NULL)
  no_fixed <- read_legofit_output(ourex1_path(), graph = NULL, include_fixed = FALSE)
  # Full result should include fixed params (T_x, T_y, T_z, one)
  expect_true("T_z" %in% full$name)
  # No-fixed result should only have free-fitted params
  expect_false("T_z" %in% no_fixed$name)
})

# Case 8: param mismatch inform fires for extra params
test_that("T8: legofit_param_mismatch inform fires for extra param in file", {
  # Build a synthetic .legofit content with an extra param T_BOGUS
  # that doesn't correspond to any node in ourex1_graph()
  content <- readLines(ourex1_path(), encoding = "UTF-8")
  # Inject T_BOGUS = 99 into the Fitted Free: block
  fitted_idx <- which(grepl("^Fitted parameter values", content))
  free_idx   <- which(grepl("^Free:", content))
  free_idx   <- free_idx[free_idx > fitted_idx[[1]]][[1]]
  content_aug <- append(content, "    T_BOGUS = 99", after = free_idx)
  tmp <- tempfile(fileext = ".legofit")
  writeLines(content_aug, tmp)
  on.exit(unlink(tmp))
  expect_message(
    read_legofit_output(tmp, graph = ourex1_graph()),
    class = "legofit_param_mismatch"
  )
})

# Case 9: malformed input (no Fitted block) errors
test_that("T9: legofit_invalid_input when no Fitted parameter values block", {
  tmp <- tempfile(fileext = ".legofit")
  writeLines(c("some header", "Initial parameter values", "Free:", "  T_R = 0.5"),
             tmp)
  on.exit(unlink(tmp))
  expect_error(
    read_legofit_output(tmp),
    class = "legofit_invalid_input"
  )
})

# fit_convergence attribute is attached
test_that("fit_convergence attribute attached with correct status", {
  result <- read_legofit_output(ourex1_path(), graph = ourex1_graph())
  conv <- attr(result, "fit_convergence")
  expect_type(conv, "list")
  expect_equal(conv$status, "reached_goal")
})

# T_admix read-time collision check
test_that("read-time T_admix collision aborts legofit_invalid_input", {
  # Synthesize a .legofit with T_admix_M AND supply a graph that has
  # node admix_M AND M as admix dest.
  tmp <- tempfile(fileext = ".legofit")
  writeLines(c(
    "Initial parameter values",
    "Fixed:",
    "       zero = 0",
    "Free:",
    "     T_admix_M = 0.1",
    "DiffEv reached_goal. cost=1e-10 spread=1e-7",
    "Fitted parameter values",
    "Free:",
    "     T_admix_M = 0.1"
  ), tmp)
  on.exit(unlink(tmp))
  # Graph with admix dest M AND a node named admix_M
  bad_graph <- tibble::tribble(
    ~from,      ~to,      ~type,    ~weight,
    "R",        "admix_M", "normal",  0.05,
    "R",        "A",       "normal",  0.05,
    "A",        "M",       "admix",   0.60,
    "admix_M",  "M",       "admix",   0.40,
    "M",        "X",       "normal",  0.02
  )
  expect_error(
    read_legofit_output(tmp, graph = bad_graph),
    class = "legofit_invalid_input"
  )
})

# ---------------------------------------------------------------------------
# Regression tests for code-review fixes
# ---------------------------------------------------------------------------

# Fix 1: seq() descending bug — header-is-last-line should give a clean error
# not a "subscript out of bounds" crash.  The minimal trigger is a file with NO
# Initial block (so fixed_initial is also 0 rows) and "Fitted parameter values"
# as the very last line.
test_that("Fix1: truncated legofit file (header is last line) aborts cleanly, not subscript OOB", {
  tmp <- tempfile(fileext = ".legofit")
  writeLines(c(
    "DiffEv reached_goal. cost=1e-10 spread=1e-7",
    "Fitted parameter values"   # last line; no Free: subsection follows
  ), tmp)
  on.exit(unlink(tmp))
  # Should abort with legofit_invalid_input (no Free block found), not crash
  expect_error(read_legofit_output(tmp), class = "legofit_invalid_input")
})

# Fix 4: include_fixed=FALSE must not emit spurious twoN mismatch for coalescent-unit models
test_that("Fix4: include_fixed=FALSE does not spuriously warn about missing twoN on coalescent-unit model", {
  mismatch_fired <- FALSE
  withCallingHandlers(
    read_legofit_output(ourex1_path(), graph = ourex1_graph(), include_fixed = FALSE),
    legofit_param_mismatch = function(m) {
      mismatch_fired <<- TRUE
      tryInvokeRestart("muffleMessage")
    }
  )
  expect_false(mismatch_fired)
})

# Fix 10: separate files get independent mismatch warnings per session
test_that("Fix10: mismatch warn fires for each file independently (path-keyed frequency_id)", {
  content <- readLines(ourex1_path(), encoding = "UTF-8")
  fitted_idx <- which(grepl("^Fitted parameter values", content))
  free_idx   <- which(grepl("^Free:", content))
  free_idx   <- free_idx[free_idx > fitted_idx[[1]]][[1]]

  # File 1: inject T_BOGUS1
  tmp1 <- tempfile(fileext = ".legofit")
  writeLines(append(content, "    T_BOGUS1 = 99", after = free_idx), tmp1)
  on.exit(unlink(tmp1), add = TRUE)

  # File 2: inject T_BOGUS2
  tmp2 <- tempfile(fileext = ".legofit")
  writeLines(append(content, "    T_BOGUS2 = 88", after = free_idx), tmp2)
  on.exit(unlink(tmp2), add = TRUE)

  # Both files should independently trigger legofit_param_mismatch
  expect_message(
    read_legofit_output(tmp1, graph = ourex1_graph()),
    class = "legofit_param_mismatch"
  )
  expect_message(
    read_legofit_output(tmp2, graph = ourex1_graph()),
    class = "legofit_param_mismatch"
  )
})

# ---------------------------------------------------------------------------
# Regression tests for second code-review pass (5 findings)
# ---------------------------------------------------------------------------

# Finding 1: convergence inform must also fire in graph=NULL path
test_that("Rv1: finished_iterations inform fires in graph=NULL path", {
  tmp <- tempfile(fileext = ".legofit")
  on.exit(unlink(tmp))
  content <- readLines(ourex1_path(), encoding = "UTF-8")
  content <- gsub("reached_goal", "finished_iterations", content)
  writeLines(content, tmp)
  expect_message(
    read_legofit_output(tmp, graph = NULL),
    class = "legofit_fit_incomplete"
  )
})

# Finding 3: graph=NULL returns exactly name/value/family, not source/free
test_that("Rv3: graph=NULL result has no internal source/free columns", {
  result <- read_legofit_output(ourex1_path(), graph = NULL)
  expect_false("source" %in% names(result))
  expect_false("free"   %in% names(result))
})

# Finding 4: empty constrained RHS gives targeted error, not subscript OOB
test_that("Rv4: empty constrained RHS aborts with readable message", {
  lgo <- "
time free T_A = 1.0
time constrained T_B =
twoN fixed one = 1
segment A t=T_A twoN=one samples=1
segment B t=T_B twoN=one samples=1
derive B from A
"
  # Expect legofit_lgo_unsupported with a message that names the parameter
  err <- tryCatch(read_lgo(text = lgo), error = function(e) e)
  expect_s3_class(err, "legofit_lgo_unsupported")
  expect_match(conditionMessage(err), "T_B", fixed = TRUE)
})
