# Minimal valid .lgo (one admixture) whose `M` and `M_A` segment lines can
# be overridden to exercise the sampled-internal-node detection. `M` is a
# derive parent (X derives from M); `M_A` is a mix parent only.
minimal_lgo_lines <- function(M   = "segment M t=T_M twoN=one",
                              M_A = "segment M_A t=T_admix_M twoN=one") {
  c(
    "time fixed T_X=0", "time free T_R=0.07", "time free T_A=0.02",
    "time free T_B=0.02", "time free T_M=0.02", "time free T_admix_M=0.02",
    "twoN fixed one=1", "mixFrac free m_M=0.4",
    "segment R t=T_R twoN=one", "segment A t=T_A twoN=one",
    "segment B t=T_B twoN=one", M, "segment X t=T_X twoN=one samples=1",
    M_A, "segment M_B t=T_admix_M twoN=one",
    "derive A from R", "derive B from R", "derive X from M",
    "derive M_A from A", "derive M_B from B",
    "mix M from M_A + m_M * M_B"
  )
}

test_that("read_lgo warns when sampled segments appear as derive parents", {
  # rha20.lgo contains segments v and a that both carry `samples=1` and
  # also appear as derive parents (n derives from v; v derives from a).
  # The edge tibble cannot represent per-segment sampling, so warn.
  expect_warning(
    read_lgo(testthat::test_path("fixtures", "rha20.lgo")),
    class = "legofit_lossy_round_trip"
  )
})

test_that("legofit_lossy_round_trip warning names the affected internal nodes", {
  w <- tryCatch(
    read_lgo(testthat::test_path("fixtures", "rha20.lgo")),
    legofit_lossy_round_trip = function(c) c
  )
  expect_s3_class(w, "legofit_lossy_round_trip")
  msg <- paste(w$message, collapse = "\n")
  expect_match(msg, "\\bv\\b")
  expect_match(msg, "\\ba\\b")
})

test_that("read_lgo is silent on .lgo with only leaf-sampled segments", {
  # minimal.lgo has `samples=1` only on X, which is a leaf (never a
  # derive parent). No lossy round trip, so no warning.
  expect_warning(
    read_lgo(testthat::test_path("fixtures", "minimal.lgo")),
    class = "legofit_lossy_round_trip",
    regexp = NA
  )
})

test_that("read_lgo warns when a sampled segment is only a mix parent", {
  # M_A feeds the admixture `mix M from M_A + ...` but is never a derive
  # parent. Its samples tag is still lost in the edge tibble, so warn.
  w <- tryCatch(
    read_lgo(text = minimal_lgo_lines(
      M_A = "segment M_A t=T_admix_M twoN=one samples=1")),
    legofit_lossy_round_trip = function(c) c
  )
  expect_s3_class(w, "legofit_lossy_round_trip")
  expect_match(paste(w$message, collapse = "\n"), "M_A", fixed = TRUE)
})

test_that("read_lgo treats samples=0 on an internal node as unsampled", {
  # `samples=0` means the segment is not sampled, so nothing is lost even
  # though M is a derive parent. No lossy round trip warning.
  expect_warning(
    read_lgo(text = minimal_lgo_lines(M = "segment M t=T_M twoN=one samples=0")),
    class = "legofit_lossy_round_trip",
    regexp = NA
  )
})

test_that("read_lgo detects samples= regardless of whitespace or casing", {
  # Some writers emit `Samples = 1` with spaces and mixed case. M is a
  # derive parent, so the lossy round trip warning must still fire.
  w <- tryCatch(
    read_lgo(text = minimal_lgo_lines(M = "segment M t=T_M twoN=one Samples = 1")),
    legofit_lossy_round_trip = function(c) c
  )
  expect_s3_class(w, "legofit_lossy_round_trip")
  expect_match(paste(w$message, collapse = "\n"), "Affected nodes: M", fixed = TRUE)
})
