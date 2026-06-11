# Regression tests for calls to functions that had been renamed away or never
# imported, so the call site raised "could not find function" when reached.
# Each test fails on the pre-fix code and passes once the call is repaired.

test_that("joint_sfs(pref=) builds an SFS via anygeno_to_aftable (dropped geno_to_afs)", {
  fix <- build_pfile_fixture(tempdir())
  sfs <- suppressWarnings(joint_sfs(NULL, pref = fix$pfile_pref))
  expect_s3_class(sfs, "tbl_df")
})

test_that("f4blockdat_from_geno(allsnps='qpfs') errors clearly instead of calling a dropped function", {
  expect_error(
    f4blockdat_from_geno("nonexistent", left = c("A", "B"), right = c("C", "D"), allsnps = "qpfs"),
    "qpfstats"
  )
})

test_that("proxypred uses summarize_proxies_list (dropped summarize_graphlist)", {
  g <- example_igraph
  expect_s3_class(admixtools:::proxypred(g, list(g)), "data.frame")
})

test_that("consistent_with_qpadm resolves equals() and returns a logical", {
  g <- example_igraph
  L <- admixtools:::get_leafnames(g)
  res <- admixtools:::consistent_with_qpadm(g, left = L[1], right = L[2:3], target = L[4])
  expect_type(res, "logical")
})

test_that("paths_from_to resolves the modern igraph as_adj_list (dropped get.adjlist)", {
  g <- example_igraph
  L <- admixtools:::get_leafnames(g)
  expect_type(admixtools:::paths_from_to(g, L[1], L[2]), "list")
})

# The graph.empty -> igraph::make_empty_graph rename in reconstruct_from_leafdist, and
# the eval_plusoneadmix -> eval_plusnadmix rename in evaluate_moreadmix, sit in internal
# code paths that need a fitted-graph scoring function and a tree-shaped leafdist to
# exercise, neither cheap to build here. Both renames are guarded instead by R CMD check,
# which flags any reappearance of the dropped name as an undefined global. This test pins
# only that the dropped eval_plusoneadmix name does not return.
test_that("evaluate_moreadmix no longer references the dropped eval_plusoneadmix", {
  src <- paste(deparse(body(admixtools:::evaluate_moreadmix)), collapse = "\n")
  expect_false(grepl("eval_plusoneadmix", src))
})
