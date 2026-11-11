# BLAS backend check. The qpgraph optimizer and the f-statistic kernels are
# dominated by dense linear algebra, so the BLAS that R is linked against
# matters a lot: R's bundled reference (netlib) BLAS is commonly ten times
# slower or more than an optimized BLAS such as OpenBLAS, MKL, or Accelerate on
# macOS. Stock R on some platforms (notably the default macOS CRAN build) ships
# the reference BLAS, so users can be paying that cost without knowing it. We
# emit a one time, suppressible note at attach so the choice is at least visible.

# TRUE if the path names a known optimized BLAS. Pure and string only so it can
# be tested directly. Tokens are the distinctive part of each vendor's real
# library filename (after symlink resolution); anything not recognized is
# treated as possibly the reference BLAS. Tokens are anchored where a bare word
# would collide: libmkl/mkl_rt not 'mkl', libgoto not 'goto', libsci_ not
# 'libsci' (vs libscience), blas_sequential/blas_openmp not bare libblas (vs
# NEC), accelerate/veclib path components not the generic libBLAS basename.
blas_looks_optimized = function(path) {
  grepl(paste0(
    'openblas|libmkl|mkl_rt|',            # OpenBLAS, Intel MKL
    'libblis|blis-|aoclblas|',            # AMD AOCL / BLIS (Flame)
    'lib[ts]?atlas|atlas/|libgoto|',      # ATLAS, GotoBLAS
    'accelerate\\.framework|veclib|',     # Apple Accelerate / vecLib
    'libsci_|libessl|armpl|',             # Cray LibSci, IBM ESSL, Arm PL
    'nvpl|nvblas|',                       # NVIDIA NVPL, NVBLAS
    'fjlapack|libfj|',                    # Fujitsu SSL2 / A64FX
    'blas_sequential|blas_openmp|',       # NEC NLC (SX-Aurora)
    'sunperf|',                           # Oracle/Sun Performance Library
    'flexiblas'),                         # FlexiBLAS dispatcher (see classify_blas)
    path, ignore.case = TRUE)
}

# a usable library path, as opposed to a bare version string or empty/NA
is_blas_path = function(x) {
  length(x) == 1L && !is.na(x) && nzchar(x) &&
    grepl('[/\\\\]|\\.(so|dylib|dll|framework)', x)
}

# Classify a single BLAS library path: TRUE reference, FALSE optimized, NA if
# it is not a usable path. Symlinks are resolved first because both Debian's
# update-alternatives and the common macOS recipe of linking libRblas at
# Accelerate point a generic name at the real library, so the optimized backend
# is only visible once the link is followed.
classify_blas = function(path) {
  if(!is_blas_path(path)) return(NA)
  path = tryCatch(normalizePath(path, mustWork = FALSE), error = function(e) path)
  if(!blas_looks_optimized(path)) return(TRUE)
  # FlexiBLAS dispatches to a backend at runtime, so a path naming flexiblas
  # does not guarantee a fast backend. Its distro default is optimized (Fedora
  # and openSUSE ship OpenBLAS as the default), but an explicit FLEXIBLAS env
  # override to the reference backend means the active BLAS really is slow.
  if(grepl('flexiblas', path, ignore.case = TRUE)) {
    fb = Sys.getenv('FLEXIBLAS')
    if(nzchar(fb) && grepl('netlib|reference', fb, ignore.case = TRUE)) return(TRUE)
  }
  FALSE
}

# TRUE if R's BLAS looks like the reference BLAS, FALSE if optimized, NA if R
# reports no usable path (e.g. some Windows builds). extSoftVersion() is light
# and reports the path on modern R; sessionInfo() is the heavier fallback.
blas_is_reference = function() {
  path = tryCatch(extSoftVersion()[['BLAS']], error = function(e) '')
  if(!is_blas_path(path)) path = tryCatch(utils::sessionInfo()$BLAS, error = function(e) '')
  classify_blas(path)
}

.onAttach = function(libname, pkgname) {
  if(isFALSE(getOption('admixtools.blas_check'))) return(invisible())
  # an advisory note must never be able to block the package from loading
  if(isTRUE(tryCatch(blas_is_reference(), error = function(e) NA)))
    packageStartupMessage(
      'admixtools does heavy dense linear algebra (qpgraph and the f-statistic ',
      'kernels) and runs much faster with an optimized BLAS such as OpenBLAS, ',
      'MKL, or Accelerate on macOS.\n',
      'R appears to be using its reference BLAS, which is commonly 10x or more ',
      'slower. See vignette("parallel", package = "admixtools") for how to ',
      'check and switch your BLAS. Set options(admixtools.blas_check = FALSE) ',
      'to silence this note.')
}
