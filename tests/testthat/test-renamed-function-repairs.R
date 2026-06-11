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
