# Round trip of per-segment nodes metadata through .lgo (read side)

test_that("read_lgo captures samples on rha20 sampled nodes", {
  res <- read_lgo(path = test_path("fixtures", "rha20.lgo"))
  nt <- graph_nodes(res)
  # v, a, d are sampled in rha20 and survive the narrow-rule round trip
  expect_equal(nt$samples[nt$name == "v"], 1L)
  expect_equal(nt$samples[nt$name == "a"], 1L)
  expect_equal(nt$samples[nt$name == "d"], 1L)
})

test_that("read_lgo captures twoN_param and resolves twoN value", {
  res <- read_lgo(path = test_path("fixtures", "rha20.lgo"))
  nt <- graph_nodes(res)
  # rha20 uses named twoN params (twoNn, one). Whatever v uses, the token
  # is captured and the numeric resolves from the param block.
  vp <- nt$twoN_param[nt$name == "v"]
  expect_false(is.na(vp))
  expect_false(is.na(nt$twoN[nt$name == "v"]))
})

test_that("read_lgo captures time_param for non-conforming t= token", {
  res <- read_lgo(path = test_path("fixtures", "rha20.lgo"))
  nt <- graph_nodes(res)
  # rha20 names times Tv, Ta, etc. (not T_<node>); the token round trips.
  expect_false(is.na(nt$time_param[nt$name == "v"]))
})

test_that("read_lgo does not create nodes rows for synthetic ghosts", {
  # graph_to_lgo's own output names ghosts <dest>_<parent>. The round trip
  # should not surface them as nodes rows.
  g   <- make_minimal_graph()
  txt <- suppressWarnings(graph_to_lgo(g, validate = FALSE))
  res <- suppressWarnings(read_lgo(text = txt))
  nt  <- graph_nodes(res)
  expect_false(any(grepl("_", nt$name)))   # no <dest>_<parent> names
})

test_that("read_lgo reads rha20 without warnings (samples preserved)", {
  # With narrow ghost detection plus per-segment capture, rha20's internal
  # samples (v, a, d) survive, so nothing is lost and read_lgo stays silent.
  expect_no_warning(read_lgo(path = test_path("fixtures", "rha20.lgo")))
})

test_that("read_lgo on inline-literal twoN captures the numeric", {
  lgo <- paste(
    "twoN fixed one = 1",
    "time fixed T_A = 0",
    "segment R t=T_A twoN=one",
    "segment A t=T_A twoN=10000 samples=1",
    "derive A from R",
    "segment B t=T_A twoN=one samples=1",
    "derive B from R",
    sep = "\n")
  res <- suppressWarnings(read_lgo(text = lgo))
  nt  <- graph_nodes(res)
  expect_equal(nt$twoN_param[nt$name == "A"], "10000")
  expect_equal(nt$twoN[nt$name == "A"], 10000)
})

# graph_to_lgo write side. make_test_nodes_graph carries no edge drift
# (weights are NA, like a read_lgo result), so these use time_handling =
# "free", which skips the drift-consistency time walk.

test_that("graph_to_lgo emits samples on an internal sampled node", {
  g <- make_test_nodes_graph()   # anc is internal AND has samples=1
  txt <- graph_to_lgo(g, time_handling = "free", validate = FALSE)
  seg_anc <- grep("^segment anc ", strsplit(txt, "\n", fixed = TRUE)[[1]], value = TRUE)
  expect_match(seg_anc, "samples=1")
})

test_that("graph_to_lgo samples= arg conflict warns, nodes tibble wins", {
  g <- make_test_nodes_graph()   # anc samples=1 in the tibble
  expect_warning(
    txt <- graph_to_lgo(g, samples = 99, time_handling = "free", validate = FALSE),
    class = "admixtools_samples_arg_overridden")
  seg_anc <- grep("^segment anc ", strsplit(txt, "\n", fixed = TRUE)[[1]], value = TRUE)
  expect_match(seg_anc, "samples=1")   # not 99
})

test_that("graph_to_lgo warns on samples override only when samples= is supplied", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "eur", samples = 3L)   # a leaf captured above the samples=1 default
  # Default call: samples= not supplied, so the nodes tibble wins silently
  # (the samples=1 default must not be read as an explicit override).
  expect_no_warning(
    graph_to_lgo(g, time_handling = "free", validate = FALSE),
    class = "admixtools_samples_arg_overridden")
  # An explicit samples= that disagrees still warns, even at the default value.
  expect_warning(
    graph_to_lgo(g, samples = 1L, time_handling = "free", validate = FALSE),
    class = "admixtools_samples_arg_overridden")
})

test_that("graph_to_lgo twoN= arg overriding captured twoN warns; arg wins", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "anc", twoN_param = "N_anc", twoN = 12345)
  expect_warning(
    txt <- graph_to_lgo(g, twoN = 99999, time_handling = "free", validate = FALSE),
    class = "admixtools_nodes_twoN_overridden")
  seg_anc <- grep("^segment anc ", strsplit(txt, "\n", fixed = TRUE)[[1]], value = TRUE)
  expect_match(seg_anc, "twoN=shared")   # arg's generated param wins, not N_anc
})

test_that("graph_to_lgo twoN= arg equal to captured value does not warn", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "anc", twoN_param = "N_anc", twoN = 12345)
  # Plain (unscoped) no-warning: the equal-value path must be fully silent,
  # not merely free of the override class.
  expect_no_warning(
    graph_to_lgo(g, twoN = 12345, time_handling = "free", validate = FALSE))
})

test_that("graph_to_lgo twoN= length-1 named vector is treated as scalar, not keyed by name", {
  # resolve_twoN_decls() treats any length-1 numeric as a scalar regardless of
  # its name, so the override check must too. Keying c(ignored = 500) by the
  # captured node names would yield NA and crash the any() gate.
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "anc", twoN_param = "N_anc", twoN = 12345)
  expect_warning(
    graph_to_lgo(g, twoN = c(ignored = 500), time_handling = "free", validate = FALSE),
    class = "admixtools_nodes_twoN_overridden")
})

test_that("graph_to_lgo emits captured twoN_param token AND declares it", {
  # Option B emits a captured param only with a value to declare, so the
  # segment references N_anc and the Parameters block declares it (valid
  # LEGOFIT, not a dangling reference).
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "anc", twoN_param = "N_anc", twoN = 12345)
  lines <- strsplit(graph_to_lgo(g, time_handling = "free", validate = FALSE),
                    "\n", fixed = TRUE)[[1]]
  expect_match(grep("^segment anc ", lines, value = TRUE), "twoN=N_anc")
  expect_true(any(grepl("^twoN (free|fixed) N_anc=12345$", lines)))
})

test_that("graph_to_lgo emits captured time_param token AND declares it", {
  g <- make_test_nodes_graph()
  g <- set_node_attrs(g, "anc", time_param = "Tanc", time = 3)
  lines <- strsplit(graph_to_lgo(g, time_handling = "free", validate = FALSE),
                    "\n", fixed = TRUE)[[1]]
  expect_match(grep("^segment anc ", lines, value = TRUE), "t=Tanc")
  expect_true(any(grepl("^time (free|fixed) Tanc=3$", lines)))
})

test_that("minimal graph (no nodes tibble) writes unchanged", {
  g <- make_minimal_graph()
  txt <- strsplit(suppressWarnings(graph_to_lgo(g, validate = FALSE)), "\n", fixed = TRUE)[[1]]
  gold <- readLines(test_path("fixtures", "minimal.lgo"))
  expect_equal(txt[nzchar(txt)], gold[nzchar(gold)])
})

# Idempotent fixed point: read -> write -> read is stable on topology
# and sample tags. Not byte-for-byte (graph_to_lgo emits its own ghost
# convention). The write uses time_handling = "free" because a read_lgo
# result carries no edge drift for the consistency time walk.
test_that("rha20 read->write->read is a topology + samples fixed point", {
  r1 <- read_lgo(path = test_path("fixtures", "rha20.lgo"))
  r2 <- suppressWarnings(read_lgo(
    text = graph_to_lgo(r1, time_handling = "free", validate = FALSE)))
  topo <- function(e) { x <- e[order(e$from, e$to), c("from", "to", "type")]
    rownames(x) <- NULL; attr(x, "nodes") <- NULL; x }
  expect_equal(topo(r1), topo(r2))
  s1 <- graph_nodes(r1); s2 <- graph_nodes(r2)
  expect_equal(s1$samples[order(s1$name)], s2$samples[order(s2$name)])
})

test_that("rha20 captured graph passes graph_to_lgo's validate=TRUE self-check", {
  # The default validate=TRUE round-trips the emitted .lgo back through read_lgo
  # and checks topology, so it guards that the narrow ghost rule and
  # graph_to_lgo's <dest>_<parent> ghost naming stay consistent on a captured,
  # admix-bearing nodes-tibble graph (rha20 has 4 mix events). The other write
  # tests use validate = FALSE, so this is the lock on that combination.
  r1 <- read_lgo(path = test_path("fixtures", "rha20.lgo"))
  expect_no_error(
    suppressWarnings(graph_to_lgo(r1, time_handling = "free", validate = TRUE)))
})

# --- Validity of the emitted .lgo (Option B) -------------------------------
# The writer must declare every parameter it references; otherwise LEGOFIT
# rejects the file ("Parameter X is undefined") even though admixtools can
# re-read it. This pure-R check mirrors `legosim --network` and runs in CI.
lgo_undeclared_params <- function(txt) {
  lines <- trimws(strsplit(txt, "\n", fixed = TRUE)[[1]])
  decl <- unlist(lapply(c("time", "twoN", "mixFrac"), function(kw) {
    d <- grep(paste0("^", kw, " "), lines, value = TRUE)
    trimws(sub(paste0("^", kw, " +\\S+ +([^=]+)=.*$"), "\\1", d))
  }))
  segs    <- grep("^segment ", lines, value = TRUE)
  seg_ref <- sub("^.*=", "",
    unlist(regmatches(segs, gregexpr("(?:^|\\s)(?:t|twoN)=\\S+", segs, perl = TRUE))))
  mixl    <- grep("^mix ", lines, value = TRUE)
  mix_ref <- gsub("[+ *]", "",
    unlist(regmatches(mixl, gregexpr("\\+ *\\S+ *\\*", mixl))))
  ref <- c(seg_ref, mix_ref)
  ref <- ref[nzchar(ref) & is.na(suppressWarnings(as.numeric(ref)))]
  setdiff(unique(ref), unique(decl))
}

test_that("graph_to_lgo emits no undefined parameters (rha20 round trip)", {
  r1  <- read_lgo(path = test_path("fixtures", "rha20.lgo"))
  txt <- graph_to_lgo(r1, time_handling = "free", validate = FALSE)
  expect_equal(lgo_undeclared_params(txt), character(0))
})

test_that("graph_to_lgo emits no undefined parameters (fresh nodes graph)", {
  txt <- graph_to_lgo(make_test_nodes_graph(), time_handling = "free",
                      validate = FALSE)
  expect_equal(lgo_undeclared_params(txt), character(0))
})

test_that("round trip re-emits the real twoN value, not a regenerated 1", {
  r1 <- read_lgo(path = test_path("fixtures", "rha20.lgo"))
  nt <- graph_nodes(r1)
  vp <- nt$twoN_param[nt$name == "v"]; vv <- nt$twoN[nt$name == "v"]
  lines <- strsplit(graph_to_lgo(r1, time_handling = "free", validate = FALSE),
                    "\n", fixed = TRUE)[[1]]
  decl <- grep(paste0("^twoN (free|fixed) ", vp, "="), lines, value = TRUE)
  expect_length(decl, 1)
  expect_equal(as.numeric(sub(".*=", "", decl)), vv, tolerance = 1e-4)
})

# --- Read-side samples= parsing (coverage moved off the deleted warn test) --

test_that("read_lgo treats samples=0 as unsampled", {
  lgo <- paste("twoN fixed one=1", "time fixed T0=0",
    "segment R t=T0 twoN=one",
    "segment A t=T0 twoN=one samples=0", "derive A from R",
    "segment B t=T0 twoN=one samples=1", "derive B from R", sep = "\n")
  nt <- graph_nodes(suppressWarnings(read_lgo(text = lgo)))
  expect_true(is.na(nt$samples[nt$name == "A"]))
  expect_equal(nt$samples[nt$name == "B"], 1L)
})

test_that("read_lgo captures samples= despite whitespace and casing", {
  lgo <- paste("twoN fixed one=1", "time fixed T0=0",
    "segment R t=T0 twoN=one",
    "segment A t=T0 twoN=one Samples = 2", "derive A from R",
    "segment B t=T0 twoN=one samples=1", "derive B from R", sep = "\n")
  nt <- graph_nodes(suppressWarnings(read_lgo(text = lgo)))
  expect_equal(nt$samples[nt$name == "A"], 2L)
})

test_that("read_lgo keeps the samples tag on a sampled mix-parent-only segment", {
  # P feeds an admixture and is a real (non <dest>_<parent>) segment, so the
  # narrow rule preserves it and its samples tag survives.
  lgo <- paste("twoN fixed one=1",
    "time fixed T0=0", "time free T1=1", "mixFrac free m=0.3",
    "segment R t=T1 twoN=one",
    "segment P t=T1 twoN=one samples=1",
    "segment Q t=T1 twoN=one",
    "segment C t=T0 twoN=one samples=1",
    "segment X t=T0 twoN=one samples=1",
    "derive P from R", "derive Q from R", "derive X from R",
    "mix C from P + m * Q", sep = "\n")
  nt <- graph_nodes(suppressWarnings(read_lgo(text = lgo)))
  expect_equal(nt$samples[nt$name == "P"], 1L)
})

test_that("graph_to_lgo tolerates a partial / hand-built nodes tibble", {
  # A nodes tibble carrying only name + samples (no twoN/time columns) must not
  # crash the writer: the override degrades gracefully to generated params for
  # the absent dimensions. (Guards the former missing-column crash.)
  g <- make_test_nodes_graph()
  attr(g, "nodes") <- tibble::tibble(name = "anc", samples = 1L)
  expect_no_error(
    suppressWarnings(graph_to_lgo(g, time_handling = "free", validate = FALSE)))
})
