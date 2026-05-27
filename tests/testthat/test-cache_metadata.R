# Tests for cache_metadata.json (the schema-versioned JSON sidecar written by
# extract_f2() alongside .f2_cache_id) and the read_f2_cache_metadata() accessor.
# Re-uses the fixture builders in helper-pfile.R.

test_that("extract_f2 writes cache_metadata.json with expected fields", {
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    outdir = file.path(getwd(), "f2_out")
    suppressMessages(suppressWarnings(
      extract_f2(pref, outdir, pops = fix$fid_pops, verbose = FALSE)))

    meta_path = file.path(outdir, "cache_metadata.json")
    expect_true(file.exists(meta_path))

    meta = read_f2_cache_metadata(outdir)
    expect_equal(meta$schema_version, 1L)
    expect_equal(meta$admixtools_version,
                 as.character(utils::packageVersion("admixtools")))
    expect_setequal(meta$pops, fix$fid_pops)
    expect_gt(meta$n_snps, 0)
    expect_gt(meta$n_blocks, 0)
    expect_equal(meta$blgsize, 0.05)
    expect_false(meta$qpfstats)
    expect_match(meta$cache_id, "^sha256:[0-9a-f]{64}$")
    # built_at is ISO 8601 with millisecond precision in UTC.
    expect_match(meta$built_at,
                 "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}\\.\\d{3}Z$")
  })
})


test_that("pops always serializes as a JSON array (I() preserves array shape)", {
  # Non-R orchestrators read cache_metadata.json with a strict schema:
  # `pops` is documented as an array. extract_f2 itself requires >=2 pops,
  # but the auto_unbox = TRUE call in the production code path would still
  # collapse a length-1 character vector to a JSON string scalar without the
  # I() wrap — silently breaking the array contract for any future code path
  # that emits a single pop. Two checks:
  #   1. The extract_f2 output (2+ pops) parses as a JSON array
  #      (simplifyVector = FALSE makes the array-vs-scalar distinction visible).
  #   2. A direct jsonlite call documents the I() contract: even length-1
  #      stays as an array.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    outdir = file.path(getwd(), "f2_out")
    suppressMessages(suppressWarnings(
      extract_f2(pref, outdir, pops = fix$fid_pops, verbose = FALSE)))

    raw = jsonlite::fromJSON(file.path(outdir, "cache_metadata.json"),
                             simplifyVector = FALSE)
    expect_true(is.list(raw$pops))

    # Single-element pops still emits a JSON array thanks to I().
    json1 = jsonlite::toJSON(list(pops = I("Pop_A")), auto_unbox = TRUE)
    expect_match(as.character(json1), '"pops":\\["Pop_A"\\]')
  })
})


test_that("cache_metadata.json cache_id matches the .f2_cache_id sidecar", {
  # The two sidecars must be consistent: a non-R orchestrator reading
  # cache_metadata.json must see the same hash the cache-key probe
  # (.f2_cache_id) emits. They come from a single compute_f2_cache_id()
  # call but are written by separate writeLines — this guards against a
  # future refactor that recomputes either one independently.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    outdir = file.path(getwd(), "f2_out")
    suppressMessages(suppressWarnings(
      extract_f2(pref, outdir, pops = fix$fid_pops, verbose = FALSE)))

    sidecar_value = readLines(file.path(outdir, ".f2_cache_id"), warn = FALSE)[1]
    meta = read_f2_cache_metadata(outdir)
    expect_equal(meta$cache_id, sidecar_value)
  })
})


test_that("compute_f2_cache_id falls back to cache_metadata.json when sidecar missing", {
  # extract_f2 writes cache_metadata.json FIRST (atomically: tempfile +
  # file.rename) and .f2_cache_id second (also atomic). A SIGKILL between the
  # two leaves cache_metadata.json present and .f2_cache_id absent —
  # cache_metadata.json's POSIX-atomic rename guarantees no truncation at
  # that final path. compute_f2_cache_id's mode-1 must recover the hash from
  # cache_metadata.json in that case, otherwise an orchestrator probe would
  # mis-report a populated cache as missing and trigger a needless rebuild.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    outdir = file.path(getwd(), "f2_out")
    suppressMessages(suppressWarnings(
      extract_f2(pref, outdir, pops = fix$fid_pops, verbose = FALSE)))

    sidecar = file.path(outdir, ".f2_cache_id")
    expected = readLines(sidecar, warn = FALSE)[1]

    # Simulate the post-crash state: cache_metadata.json present, sidecar gone.
    file.remove(sidecar)
    expect_false(file.exists(sidecar))
    expect_true(file.exists(file.path(outdir, "cache_metadata.json")))

    recovered = compute_f2_cache_id(outdir)
    expect_equal(recovered, expected)
  })
})


test_that("compute_f2_cache_id errors when both sidecars are missing", {
  # The error message must mention cache_metadata.json so orchestrators
  # debugging a missing-cache report know to check both sidecars.
  withr::with_tempdir({
    expect_error(compute_f2_cache_id(getwd()),
                 "cache_metadata.json")
  })
})


test_that("compute_f2_cache_id warns + errors clearly when cache_metadata.json is malformed", {
  # A truncated cache_metadata.json (from a mid-write crash) is a different
  # signal from a missing one. Orchestrators should see a warning naming the
  # parse failure, not a generic "no sidecar" error that sends them looking
  # in the wrong place.
  withr::with_tempdir({
    dir.create("outdir")
    # Write malformed JSON (missing close brace) and no .f2_cache_id sidecar.
    writeLines('{"schema_version": 1, "cache_id": "sha256:',
               file.path("outdir", "cache_metadata.json"))
    expect_warning(
      expect_error(compute_f2_cache_id("outdir"), "cache_metadata.json"),
      "cache_metadata.json present but unparseable")
  })
})


test_that("read_f2_cache_metadata errors clearly when the file is missing", {
  # Orchestrators pointed at a directory built by an older admixtools (or
  # by something else entirely) must get an actionable error, not a silent
  # NULL or a parser stack trace.
  withr::with_tempdir({
    expect_error(read_f2_cache_metadata(getwd()),
                 "cache_metadata.json not found")
  })
})


test_that("compute_f2_cache_id rejects a malformed cache_id field in JSON", {
  # The JSON parses cleanly but the cache_id field doesn't match the
  # sha256:<64hex> contract (truncated hash from a bad tool, manual edit,
  # or older incompatible schema). compute_f2_cache_id must fall through
  # to the same stop() that fires on a missing field — not silently return
  # a bogus hash that an orchestrator would then use as a cache key.
  withr::with_tempdir({
    dir.create("outdir")
    # Well-formed JSON, malformed cache_id (too short / wrong charset).
    writeLines('{"schema_version":1,"cache_id":"sha256:not-a-valid-hash"}',
               file.path("outdir", "cache_metadata.json"))
    expect_error(compute_f2_cache_id("outdir"), "cache_metadata.json")
  })
})


test_that("cache_metadata.json is written atomically (no truncated file on disk)", {
  # `.write_atomic` uses tempfile + file.rename so the final path either
  # holds the previous version or fully-committed new bytes — never a
  # partial write. We can't simulate a SIGKILL in-process, but we can
  # verify the file ends up well-formed (the rename committed in full)
  # and that no .tmp sidecar leaks into the outdir after a successful
  # extract_f2 run.
  withr::with_tempdir({
    fix = build_pfile_fixture(getwd(), with_fid = TRUE)
    pref = fix$bed_pref
    testthat::skip_if(is.na(pref), "BED fixture unavailable")

    outdir = file.path(getwd(), "f2_out")
    suppressMessages(suppressWarnings(
      extract_f2(pref, outdir, pops = fix$fid_pops, verbose = FALSE)))

    # File parses (committed in full).
    expect_silent(jsonlite::fromJSON(file.path(outdir, "cache_metadata.json"),
                                     simplifyVector = TRUE))
    # No leftover tempfiles (file.rename succeeded; cleanup happened).
    leftovers = list.files(outdir, pattern = "\\.tmp$", all.files = TRUE)
    expect_length(leftovers, 0)
  })
})
