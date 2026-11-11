# Tests for the attach-time BLAS note (R/zzz.R). The classification is a pure
# string function so it can be checked against real library paths from each
# platform without depending on the BLAS the test machine happens to use.

test_that("blas_looks_optimized recognizes optimized backends across vendors", {
  optimized <- c(
    "/opt/homebrew/Cellar/openblas/0.3.33/lib/libopenblasp-r0.3.33.dylib",   # OpenBLAS (Homebrew)
    "/usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3",               # OpenBLAS (Debian alt)
    "/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_rt.so.2",               # Intel MKL
    "/usr/lib/x86_64-linux-gnu/blis-openmp/libblis.so.4",                    # BLIS / AMD AOCL
    "/usr/lib/x86_64-linux-gnu/atlas/libblas.so.3",                          # ATLAS (Debian)
    "/usr/lib/libtatlas.so.3",                                               # ATLAS (threaded)
    "/usr/lib/libsatlas.so.3",                                               # ATLAS (serial)
    "/opt/openblas/lib/libgoto2_haswellp-r1.13.so",                          # GotoBLAS2
    "/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/libBLAS.dylib", # Apple Accelerate
    "/opt/cray/pe/libsci/23.09/GNU/12.3/x86_64/lib/libsci_gnu_82_mp.so",     # Cray LibSci
    "/opt/ibm/essl/lib64/libesslsmp.so.1",                                   # IBM ESSL
    "/opt/arm/armpl/lib/libarmpl_lp64_mp.so",                                # Arm Performance Libraries
    "/usr/lib/nvpl/lib/libnvpl_blas_lp64_gomp.so.0",                         # NVIDIA NVPL
    "/usr/lib/libnvblas.so.12",                                              # NVBLAS (GPU interposer)
    "/opt/FJSVxos/lib64/libfjlapackexsve_lp64.so",                           # Fujitsu A64FX SSL2
    "/opt/nec/ve/nlc/lib/libblas_openmp.so",                                 # NEC NLC
    "/opt/oracle/lib/libsunperf.so.5",                                       # Oracle/Sun Performance
    "/usr/lib64/libflexiblas.so.3")                                          # FlexiBLAS dispatcher
  for (p in optimized) expect_true(admixtools:::blas_looks_optimized(p), info = p)
})

test_that("blas_looks_optimized does not recognize the reference BLAS", {
  reference <- c(
    "/Library/Frameworks/R.framework/Resources/lib/libRblas.dylib",   # macOS CRAN default
    "/usr/lib/R/lib/libRblas.so",
    "/usr/lib/x86_64-linux-gnu/blas/libblas.so.3",                    # Debian netlib reference
    "/usr/lib64/blas/reference/libblas.so.3")                         # generic netlib reference
  for (p in reference) expect_false(admixtools:::blas_looks_optimized(p), info = p)
})

test_that("a reference BLAS under a vendor-named parent directory is not misclassified", {
  # the optimized tokens are anchored so a coincidental directory or user name
  # does not suppress the note for a real reference BLAS
  tricky <- c(
    "/opt/atlassian/R/lib/libRblas.so",      # 'atlas' inside 'atlassian'
    "/home/dessler/R/lib/libRblas.so",       # 'essl' inside 'dessler'
    "/home/bliss/R/lib/libRblas.so",         # 'blis' inside 'bliss'
    "/srv/accelerate-team/lib/libRblas.so",  # 'accelerate' without '.framework'
    "/opt/libscience/lib/libRblas.so")       # 'libsci' inside 'libscience'
  for (p in tricky) expect_false(admixtools:::blas_looks_optimized(p), info = p)
})

test_that("classify_blas honors the FLEXIBLAS backend override", {
  fb_path <- "/usr/lib64/libflexiblas.so.3"
  old <- Sys.getenv("FLEXIBLAS", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("FLEXIBLAS") else Sys.setenv(FLEXIBLAS = old), add = TRUE)
  Sys.unsetenv("FLEXIBLAS")
  expect_false(admixtools:::classify_blas(fb_path))            # distro default, treat as optimized
  Sys.setenv(FLEXIBLAS = "netlib")
  expect_true(admixtools:::classify_blas(fb_path))             # explicit reference backend
  Sys.setenv(FLEXIBLAS = "openblas")
  expect_false(admixtools:::classify_blas(fb_path))            # explicit fast backend
})

test_that("classify_blas follows a symlink to the real backend", {
  # update-alternatives (Debian) and the macOS Accelerate recipe both expose the
  # optimized library only through a symlink at a generic name, so classify_blas
  # must resolve the link before deciding.
  d <- tempfile("blas"); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  real <- file.path(d, "libopenblas.so.3"); file.create(real)
  link <- file.path(d, "libRblas.so")
  ok <- tryCatch(file.symlink(real, link), warning = function(w) FALSE, error = function(e) FALSE)
  if (!isTRUE(ok)) skip("symlinks not available on this platform")
  expect_false(admixtools:::classify_blas(link))   # resolves to libopenblas, not reference
  ref <- file.path(d, "libRblas2.so"); file.create(ref)
  expect_true(admixtools:::classify_blas(ref))      # a real file with no optimized name
})

test_that("classify_blas returns NA for input that is not a library path", {
  for (x in list("", "3.11.0", NA_character_, character(0)))
    expect_true(is.na(admixtools:::classify_blas(x)), info = paste(x, collapse = ","))
})

test_that("blas_is_reference returns a single logical or NA", {
  v <- admixtools:::blas_is_reference()
  expect_length(v, 1)
  expect_true(is.logical(v))
  expect_true(is.na(v) || v %in% c(TRUE, FALSE))
})

test_that("the attach note fires on a reference BLAS and is silenced by the option", {
  # Drive both branches without depending on the test machine's BLAS.
  testthat::local_mocked_bindings(blas_is_reference = function() TRUE, .package = "admixtools")
  expect_message(admixtools:::.onAttach("lib", "admixtools"), "reference BLAS")

  old <- options(admixtools.blas_check = FALSE)
  on.exit(options(old), add = TRUE)
  expect_message(admixtools:::.onAttach("lib", "admixtools"), NA)   # opted out, silent even on reference
})

test_that("the attach note stays silent on an optimized BLAS", {
  testthat::local_mocked_bindings(blas_is_reference = function() FALSE, .package = "admixtools")
  expect_message(admixtools:::.onAttach("lib", "admixtools"), NA)
})
