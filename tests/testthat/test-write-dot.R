# Tests for write_dot (PR #122).
#
# Headline contract: the function produces a Graphviz-parseable dot file
# whose contents reflect the args. New args added by this PR are
# `fontsize`, `color`, `hide_weights`, `highlight_unidentifiable`,
# `nodesep`, `ranksep`, `fix_names`; pre-existing args `size1`, `size2`,
# `title`, `dot2pdf` retain their original positional slots so callers
# like `write_dot(g, file, 8, 12, "title")` are unaffected.

# Convenience: render to a tempfile and read back as a single string
# for substring assertions.
.write_dot_str = function(graph, ...) {
  f = tempfile(fileext = ".dot")
  on.exit(unlink(f), add = TRUE)
  suppressWarnings(write_dot(graph, outfile = f, ...))
  paste(readLines(f), collapse = "\n")
}

.example_igraph = function() {
  data("example_igraph", package = "admixtools", envir = environment())
  get("example_igraph", envir = environment())
}

test_that("write_dot emits a digraph with size/nodesep/ranksep/title from args", {
  out = .write_dot_str(.example_igraph(),
                       size1 = 8, size2 = 12, title = "MyGraph",
                       nodesep = 0.4, ranksep = 0.7)
  expect_match(out, "digraph G \\{", fixed = FALSE)
  expect_match(out, 'label = "MyGraph"', fixed = TRUE)
  expect_match(out, 'size = "8,12"',     fixed = TRUE)
  expect_match(out, 'nodesep = "0.4"',   fixed = TRUE)
  expect_match(out, 'ranksep = "0.7"',   fixed = TRUE)
})

test_that("write_dot preserves the pre-PR positional argument contract", {
  # Pre-PR signature: write_dot(graph, outfile, size1, size2, title, dot2pdf).
  # The positional call below must still assign 8 -> size1, 12 -> size2,
  # "positional" -> title (i.e. new args must NOT have been inserted
  # between outfile and size1, which would shift these).
  f = tempfile(fileext = ".dot")
  on.exit(unlink(f), add = TRUE)
  suppressWarnings(write_dot(.example_igraph(), f, 8, 12, "positional"))
  txt = paste(readLines(f), collapse = "\n")
  expect_match(txt, 'size = "8,12"',          fixed = TRUE)
  expect_match(txt, 'label = "positional"',   fixed = TRUE)
})

test_that("write_dot: fontsize is injected into every label", {
  out = .write_dot_str(.example_igraph(), fontsize = 22)
  # Edge labels and node declarations both carry the fontsize.
  expect_match(out, 'fontsize = "22"', fixed = TRUE)
  # Default fontsize is 14; switching to 22 should remove all 14s.
  expect_false(grepl('fontsize = "14"', out, fixed = TRUE))
})

test_that("write_dot: color = FALSE produces only Black colors (no ggplot bridging)", {
  out = .write_dot_str(.example_igraph(), color = FALSE)
  # Every `color = "..."` in the output must be "Black".
  matches = regmatches(out, gregexpr('color = "[^"]+"', out))[[1]]
  expect_true(length(matches) > 0, info = "expected at least one color attribute")
  expect_true(all(matches == 'color = "Black"'),
              info = "color=FALSE must not produce any non-Black color attrs")
})

test_that("write_dot: color = TRUE produces at least one non-Black color (ggplot-bridged)", {
  # Allow ggplot-bridging warnings (e.g. y-coordinate collisions on small
  # example graphs); those don't affect the contract that some non-Black
  # color is emitted.
  out = .write_dot_str(.example_igraph(), color = TRUE)
  matches = regmatches(out, gregexpr('color = "[^"]+"', out))[[1]]
  expect_true(any(matches != 'color = "Black"'),
              info = "color=TRUE should produce some non-Black colors via plot_graph bridging")
})

test_that("write_dot: hide_weights = TRUE blanks non-admix edge labels", {
  out_default = .write_dot_str(.example_igraph(), hide_weights = FALSE)
  out_hidden  = .write_dot_str(.example_igraph(), hide_weights = TRUE)
  # Default: non-admix edges carry numeric labels (round(weight * 1000));
  # at example_igraph's typical weights those numbers are >= 1.
  expect_match(out_default, 'label = "[0-9]', fixed = FALSE)
  # Hidden: non-admix edges carry a single-space label.
  expect_match(out_hidden, 'label = " "', fixed = TRUE)
})

test_that("write_dot: fix_names = TRUE replaces dots and dashes in node names with underscores", {
  # Construct a tiny igraph with `.` and `-` in node names so we can
  # observe the fix_names rewrite. (Edge-list input requires a `type`
  # column the function computes only on the igraph branch, so we go
  # igraph -> dot here.)
  edges = tibble::tibble(from = c("R", "R", "Pop.A", "Pop-B"),
                         to   = c("Pop.A", "Pop-B", "L1", "L2"))
  g = admixtools:::edges_to_igraph(edges)

  out_default = suppressWarnings(.write_dot_str(g, color = FALSE,
                                                fix_names = FALSE))
  out_fixed   = suppressWarnings(.write_dot_str(g, color = FALSE,
                                                fix_names = TRUE))

  # fix_names=TRUE normalizes dots/dashes to underscores before the
  # dot-stripping step; result is Pop_A / Pop_B in node IDs.
  expect_match(out_fixed, "Pop_A", fixed = TRUE)
  expect_match(out_fixed, "Pop_B", fixed = TRUE)

  # Default path strips [\\.-] entirely from node IDs in the emitted
  # dot text (existing behavior), producing PopA / PopB. Neither path
  # should leave dotted-name strings in the output.
  expect_false(grepl("Pop\\.A", out_default))
  expect_false(grepl("Pop-B",   out_default, fixed = TRUE))
})

test_that("write_dot: graphviz `dot` parses the emitted file (if dot is on PATH)", {
  # Strongest cross-check: the output isn't just substring-matchable but
  # is actually valid Graphviz dot syntax. Skip if `dot` isn't installed.
  dot_bin = Sys.which("dot")
  testthat::skip_if(dot_bin == "", "graphviz `dot` not on PATH")

  f = tempfile(fileext = ".dot")
  on.exit(unlink(f), add = TRUE)
  suppressWarnings(write_dot(.example_igraph(), outfile = f, color = FALSE))

  # `dot -Tsvg -o /dev/null` parses + lays out the graph without writing
  # output; exit 0 = valid syntax, non-zero = parse error.
  status = suppressWarnings(system2(dot_bin,
                                    args = c("-Tsvg", shQuote(f),
                                             "-o", "/dev/null"),
                                    stdout = FALSE, stderr = FALSE))
  expect_equal(status, 0L,
               info = "graphviz `dot` failed to parse write_dot's output")
})
