# Shared synthetic aftable + popcomb fixture for the cpp_aftable_to_dstat*
# kernel tests (test-cpp-aftable-omp.R and test-cpp-aftable-rowmeans.R).
#
# The RNG call sequence is fixed, so a given (npop, nsnp, ncomb, nmodels,
# na_frac, seed) reproduces byte-identical inputs across runs and across the
# two test files. The two files pass different defaults but share this one
# builder, so the dstat-kernel input contract (aftable shape, usesnps row
# count == nmodels, 1-based p-vectors) is defined in a single place.
build_dstat_fixture = function(npop = 12L, nsnp = 1000L, ncomb = 256L,
                               nmodels = 4L, na_frac = 0.05, seed = 20260525L) {
  withr::with_seed(seed, {
    aft = matrix(runif(npop * nsnp), nrow = npop)
    aft[sample.int(length(aft), size = round(na_frac * length(aft)))] = NA_real_
    p1 = sample.int(npop, ncomb, replace = TRUE)
    p2 = sample.int(npop, ncomb, replace = TRUE)
    p3 = sample.int(npop, ncomb, replace = TRUE)
    p4 = sample.int(npop, ncomb, replace = TRUE)
    # nmodels models so the !allsnps path exercises modelvec selection; usesnps
    # is an (nmodels x nsnp) 0/1 mask keeping ~90% of SNPs per model.
    modelvec = sample.int(nmodels, ncomb, replace = TRUE)
    usesnps = matrix(sample(c(0, 1), nmodels * nsnp, replace = TRUE,
                            prob = c(0.1, 0.9)),
                     nrow = nmodels, ncol = nsnp)
  })
  list(aft = aft, p1 = p1, p2 = p2, p3 = p3, p4 = p4,
       modelvec = modelvec, usesnps = usesnps)
}
