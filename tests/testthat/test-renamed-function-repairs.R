# Regression tests for calls to functions that had been renamed away or never
# imported, so the call site raised "could not find function" when reached.
# Each test fails on the pre-fix code and passes once the call is repaired.

test_that("joint_sfs(pref=) builds an SFS via anygeno_to_aftable (dropped geno_to_afs)", {
  fix <- build_pfile_fixture(tempdir())
  sfs <- suppressWarnings(joint_sfs(NULL, pref = fix$pfile_pref))
  expect_s3_class(sfs, "tbl_df")
})
