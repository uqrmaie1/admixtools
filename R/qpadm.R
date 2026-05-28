
# Package-internal: reciprocal-condition threshold below which qpadm() flags
# f4_var as near-singular. ~sqrt(machine eps); empirically the bar where the
# pseudo-inverse fallback starts producing misleadingly tight standard
# errors. Single source of truth for qpadm()'s auto-warn gate, qpadm_sweep()'s
# verbose summary alert, and tests that pin the trigger threshold; do not
# duplicate this literal at call sites.
.rcond_concern = 1e-8


# a, b: num left - 1, num right - 1 (dimensions of f4_est)
qpwave_dof = function(a, b, r) a*b - qpadm_dof(a, b, r)
qpwave_dof = function(a, b, r) r*(a+b-r)
qpwave_dof = function(a, b, r) (a+b)*r - r^2
qpwave_dof = function(a, b, r) a*b - (a-r)*(b-r)

qpadm_dof = function(a, b, r) a*b - qpwave_dof(a, b, r)
qpadm_dof = function(a, b, r) a*b - (a+b)*r + r^2
qpadm_dof = function(a, b, r) (a-r)*(b-r)



opt_A = function(B, xmat, qinv, fudge = 0.0001) {
  # return A which minimizes covariance-weighted t(c(E)) %*% qinv %*% c(E), where E = xmat - (A %*% B)
  # A: nr * rnk
  # B: rnk * nc
  # xmat: nr * nc
  # tdim: rnk * max(nr, nc)
  # coeffs: tdim * tdim
  # rhs: rnk * nc
  # qinv: nr*nc * nr*nc
  B2 = diag(nrow(xmat)) %x% B
  coeffs = B2 %*% qinv %*% t(B2)
  rhs = c(t(xmat)) %*% qinv %*% t(B2)
  diag(coeffs) = diag(coeffs) + fudge*sum(diag(coeffs))
  A2 = solve(coeffs, rhs[1,])
  matrix(A2, nrow(xmat), byrow = TRUE)
}


opt_B = function(A, xmat, qinv, fudge = 0.0001) {
  # return B which minimizes covariance-weighted t(c(E)) %*% qinv %*% c(E), where E = xmat - (A %*% B)
  # A: nr * rnk
  # B: rnk * nc
  # xmat: nr * nc
  # tdim: rnk * max(nr, nc)
  # coeffs: tdim * tdim
  # rhs: rnk * nc
  # qinv: nr*nc * nr*nc
  A2 = A %x% diag(ncol(xmat))
  coeffs = t(A2) %*% qinv %*% A2
  rhs = c(t(xmat)) %*% qinv %*% A2
  diag(coeffs) = diag(coeffs) + fudge*sum(diag(coeffs))
  B2 = solve(coeffs, rhs[1,])
  matrix(B2, ncol = ncol(xmat), byrow = TRUE)
}


get_weights_covariance = function(f4_lo, qinv, block_lengths, fudge = 0.0001, boot = FALSE,
                                  constrained = FALSE, qpsolve = NULL) {
  rnk = dim(f4_lo)[1]-1
  numreps = dim(f4_lo)[3]
  wmat = matrix(NA, numreps, dim(f4_lo)[1])
  for(i in 1:numreps) {
    wmat[i,] = qpadm_weights(as.matrix(f4_lo[,,i]), qinv, rnk, fudge = fudge,
                             constrained = constrained, qpsolve = qpsolve)$weights
  }
  if(!boot) wmat = wmat * (numreps-1) / sqrt(numreps)
  cov(wmat)
}


qpadm_weights = function(xmat, qinv, rnk, fudge = 0.0001, iterations = 20,
                         constrained = FALSE, qpsolve = NULL) {
  if(rnk == 0) return(list(weights = 1))
  nr = nrow(xmat)
  f4svd = svd(xmat)
  B = t(f4svd$v[, seq_len(rnk), drop=FALSE])
  A = xmat %*% t(B)
  for(i in seq_len(iterations)) {
    A = opt_A(B, xmat, qinv, fudge = fudge)
    B = opt_B(A, xmat, qinv, fudge = fudge)
  }
  x = t(cbind(A, 1))
  y = c(rep(0, rnk), 1)
  rhs = crossprod(x)
  lhs = crossprod(x, y)
  if(constrained) {
    Amat = cbind(1, -1, diag(nr), -diag(nr))
    bvec = c(1, -1, rep(0, nr), -rep(1, nr))
    w = qpsolve(rhs, lhs, Amat, bvec)
  } else w = as.matrix(solve(rhs, lhs))[,1]
  weights = w/sum(w)
  namedList(weights, A, B)
}


#' Estimate admixture weights
#'
#' `qpadm` models a target population as a mixture of left (source) populations, given a set of right (outgroup) populations.
#' It can be used to estimate whether the left populations explain all genetic variation in the target population, relative to the right populations, and to estimate admixture proportions of the left populations to the target population.
#' @export
#' @param data The input data in the form of:
#' \itemize{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}}
#' \item A directory with f2 statistics
#' \item The prefix of a genotype file
#' }
#' @param left Left populations (sources)
#' @param right Right populations (outgroups)
#' @param target Target population
#' @param f4blocks Instead of f2 blocks, f4 blocks can be supplied. This is used by \code{\link{qpadm_multi}}
#' @param fudge Value added to diagonal matrix elements before inverting
#' @param fudge_twice Setting this to `TRUE` should result in p-values that better match those in the original qpAdm program
#' @param auto_only Use only chromosomes 1 to 22.
#' @param blgsize SNP block size in Morgan. Default is 0.05 (5 cM). If `blgsize` is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
#' @param poly_only Exclude sites with identical allele frequencies in all populations.
#' @param boot If `FALSE` (the default), block-jackknife resampling will be used to compute standard errors.
#' Otherwise, block-bootstrap resampling will be used to compute standard errors. If `boot` is an integer, that number
#' will specify the number of bootstrap resamplings. If `boot = TRUE`, the number of bootstrap resamplings will be
#' equal to the number of SNP blocks.
#' @param getcov Compute weights covariance. Setting `getcov = FALSE` will speed up the computation.
#' @param constrained Constrain admixture weights to be non-negative
#' @param return_f4 Return f4-statistics
#' @param cpp Use C++ functions. Setting this to `FALSE` will be slower but can help with debugging.
#' @param singular_threshold If not `NA` (the default), error out when the
#'   reciprocal condition number of the (fudged) f4 variance matrix is below
#'   this threshold. A common choice is `1e-12` -- below that the downstream
#'   `solve()` and the C++ weight optimization fall back to a pseudo-inverse,
#'   which makes the returned weight standard errors numerically tight in a
#'   way that doesn't reflect actual sampling uncertainty. The default
#'   (`NA`) preserves the historical behavior: silently fall through to the
#'   pseudo-inverse fallback. In either mode, the reciprocal condition
#'   number is reported as `f4_var_rcond` on the returned list so callers
#'   can gate downstream interpretations programmatically.
#'
#'   Low-level RcppArmadillo `solve(): system is singular; attempting approx
#'   solution` stderr lines are suppressed at compile time via
#'   `ARMA_WARN_LEVEL 1`. Those lines previously flooded stderr on batch fits
#'   at higher K without carrying any information not already exposed by the
#'   R-side diagnostics. The canonical programmatic signal is `f4_var_rcond`
#'   (always populated, regardless of `verbose`); `qpadm_multi()` and
#'   `qpadm_sweep()` force `verbose = FALSE` internally, so the interactive
#'   `warning()` only fires for direct `qpadm()` calls at `verbose = TRUE`.
#'   Batch callers should inspect `f4_var_rcond` (and `f4_var_singular_loadings`)
#'   on the returned object rather than rely on stderr / warning output.
#' @param verbose Print progress updates
#' @param ... If `data` is the prefix of genotype files, additional arguments will be passed to \code{\link{f4blockdat_from_geno}}
#' @return `qpadm` returns a list with up to four data frames describing the model fit, plus a numeric `f4_var_rcond` diagnostic:
#' \enumerate{
#' \item `f4_var_rcond`: Reciprocal condition number of the (fudged) f4 variance matrix that was inverted to compute weights and standard errors. Values near machine epsilon (~`1e-15`) indicate that the right populations are linearly dependent (e.g., a right pop is a sister to one of the sources), and the returned weight SEs may be artificially tight as a result of the pseudo-inverse fallback. See `singular_threshold` for opt-in error gating.
#' \item `f4_var_singular_loadings`: When `f4_var_rcond` is concerningly low (< `1e-8`), a tibble with one row per right population, sorted by L2 loading on the smallest right-singular vector of `f4_var`. Right pops at the top of this list are the most likely offenders -- typically a sister-clade pair (two pops with similar high loadings) or a lone right pop too close to a source. `NULL` otherwise.
#' \item `weights`: A data frame with estimated admixture proportions where each row is a left population.
#' \item `f4`: A data frame with estimated and fitted f4-statistics
#' \item `rankdrop`: A data frame describing model fits with different ranks, including *p*-values for the overall fit
#' and for nested models (comparing two models with rank difference of one). A model with `L` left populations and `R` right populations has an f4-matrix of dimensions `(L-1)*(R-1)`. If no two left population form a clade with respect to all right populations, this model will have rank `(L-1)*(R-1)`.
#'   \itemize{
#'     \item `f4rank`: Tested rank
#'     \item `dof`: Degrees of freedom of the chi-squared null distribution: `(L-1-f4rank)*(R-1-f4rank)`
#'     \item `chisq`: Chi-sqaured statistic, obtained as `E'QE`, where `E` is the difference between estimated and fitted f4-statistics, and `Q` is the f4-statistic covariance matrix.
#'     \item `p`: p-value obtained from `chisq` as `pchisq(chisq, df = dof, lower.tail = FALSE)`
#'     \item `dofdiff`: Difference in degrees of freedom between this model and the model with one less rank
#'     \item `chisqdiff`: Difference in chi-squared statistics
#'     \item `p_nested`: *p*-value testing whether the difference between two models of rank difference 1 is significant
#'   }
#' \item `popdrop`: A data frame describing model fits with different populations. Note that all models with fewer populations use the same set of SNPs as the first model.
#'   \itemize{
#'     \item `pat`: A binary code indicating which populations are present in this model. A `1` represents dropped populations. The full model is all zeros.
#'     \item `wt`: Number of populations dropped
#'     \item `dof`: Degrees of freedom of the chi-squared null distribution: `(L-1-f4rank)*(R-1-f4rank)`
#'     \item `chisq`: Chi-sqaured statistic, obtained as `E'QE`, where `E` is the difference between estimated and fitted f4-statistics, and `Q` is the f4-statistic covariance matrix.
#'     \item `p`: *p*-value obtained from `chisq` as `pchisq(chisq, df = dof, lower.tail = FALSE)`
#'     \item `f4rank`: Tested rank
#'     \item `feasible`: A model is feasible if all weights fall between 0 and 1
#'     \item `<population name>`: The weights for each population in this model
#'   }
#' }
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 10)
#' @seealso \code{\link{qpwave}}, \code{\link{lazadm}}
#' @examples
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' target = 'Denisova.DG'
#' qpadm(example_f2_blocks, left, right, target)
#' \dontrun{
#' # The original ADMIXTOOLS qpAadm program has an option called "allsnps"
#' # that selects different SNPs for each f4-statistic, which is
#' # useful when working with sparse genotype data.
#' # To get the same behavior in ADMIXTOOLS 2, supply the genotype data prefix
#' # and set `allsnps = TRUE`
#' qpadm("/my/geno/prefix", left, right, target, allsnps = TRUE)
#' }
qpadm = function(data, left, right, target, f4blocks = NULL,
                 fudge = 0.0001, fudge_twice = FALSE, auto_only = TRUE, blgsize = 0.05,
                 poly_only = FALSE, boot = FALSE, getcov = TRUE,
                 constrained = FALSE, return_f4 = FALSE, cpp = TRUE,
                 singular_threshold = NA_real_, verbose = TRUE, ...) {

  #----------------- prepare f4 stats -----------------
  f2_blocks = NULL
  args = list(...)
  if(length(args) > 0 && (!is.null(f4blocks) || !is_geno_prefix(data))) {
    stop(paste0("The following arguments are not used: '", paste0(names(args), collapse = "', '"), "'"))
  }
  if(!all(names(args) %in% names(formals(f4blockdat_from_geno)))) {
    notused = setdiff(names(args), names(formals(f4blockdat_from_geno)))
    stop(paste0("The following arguments are not recognized: '", paste0(notused, collapse = "', '"), "'"))
  }
  if(is.null(f4blocks)) {
    if(all(file.exists(left, right))) {
      left %<>% readLines
      right %<>% readLines
    }
    if(!is.null(target)) left = c(target, setdiff(left, target))
    if(is_geno_prefix(data)) {
      if(return_f4) {
        # pop1 = left[1] handles both target != NULL (line above puts target at
        # left[1]) and target = NULL (qpwave-style; the canonical anchor for
        # the f4 parameterization f4(left[1], left[-1][i]; right[1], right[-1][j])
        # used by f2blocks_to_f4blocks on the non-geno path). Previously this
        # was `pop1 = target` which silently dropped the pop1 column on the
        # target = NULL path, breaking the qpwave-with-allsnps call shape
        # documented in #69.
        popcombs = expand_grid(pop1 = left[1], pop2 = left[-1], pop3 = right, pop4 = right) %>% filter(pop3 != pop4)
        f4blockdat_all = f4blockdat_from_geno(data, popcombs = popcombs, auto_only = auto_only,
                                              blgsize = blgsize, poly_only = poly_only,  verbose = verbose, ...)
        f4blockdat = f4blockdat_all %>% filter(pop3 == right[1], pop4 != right[1])
      } else {
        f4blockdat = f4blockdat_from_geno(data, left = left, right = right, auto_only = auto_only,
                                          blgsize = blgsize, poly_only = poly_only, verbose = verbose, ...)
      }
      f4blocks = f4blockdat %>% f4blockdat_to_f4blocks()
      if(verbose) {
        snptab = f4blockdat %>% group_by(pop1, pop2, pop3, pop4) %>% summarize(n = sum(n))
        if(isTRUE(args$allsnps)) {
          alert_info('"allsnps = TRUE" uses different SNPs for each f4-statistic\n  Number of SNPs used for each f4-statistic:\n')
          snptab %>% ungroup %>% as.data.frame %>% print
        } else {
          alert_info(paste0('Number of SNPs after excluding those with missing data: ', unique(snptab$n),'\n'))
        }
      }
    } else {
      if(isTRUE(args$allsnps)) stop('"allsnps = TRUE" is only effective when reading data from genotype files!')
      if(verbose) alert_info('Computing f4 stats...\n')
      f2_blocks = get_f2(data, pops = c(left, right), auto_only = auto_only, blgsize = blgsize,
                         poly_only = poly_only, afprod = TRUE, verbose = verbose)
      f4blocks = f2blocks_to_f4blocks(f2_blocks, left, right)
    }
    f4blocks = f4blocks[left[-1], right[-1],,drop=FALSE]
  } else {
    if(!all(map_lgl(list(data, left, right), is.null)))
      stop("Can't provide 'f4blocks' and data/left/right at the same time!")
    naml = dimnames(f4blocks)[1]
    namr = dimnames(f4blocks)[2]
    left = c(names(naml), naml[[1]])
    right = c(names(namr), namr[[1]])
    if(!is.null(target)) target = left[1]
  }
  pops = c(left, right)
  if(any(duplicated(pops)))
    stop(paste0('Duplicated pops: ', paste0(unique(pops[duplicated(pops)]), collapse=', ')))

  f4stats = f4blocks %>% f4blocks_to_f4stats(boot = boot)

  f4_est = f4stats$est
  f4_var = f4stats$var
  f4_lo = f4stats$f4_lo
  block_lengths = f4stats$block_lengths
  diag(f4_var) = diag(f4_var) + fudge*sum(diag(f4_var))
  if(fudge_twice) diag(f4_var) = diag(f4_var) + fudge*sum(diag(f4_var))

  # Reciprocal condition number of (fudged) f4 variance matrix. Used as a
  # diagnostic for "are the right pops independent of the sources?".
  # Very low rcond (~1e-15 / machine epsilon) typically means one or more
  # right pops are sister to a source, which makes f4_var near-singular and
  # the downstream pseudo-inverse fallback in cpp_qpadm_weights produce
  # weight SEs that look tight but reflect numerical noise, not sampling
  # uncertainty (issue #8).
  f4_var_rcond = tryCatch(rcond(f4_var), error = function(e) NA_real_)

  # If rcond is concerningly low, compute the SVD-based attribution: the
  # smallest right-singular vector tells us which (left, right) f4 directions
  # are in the degenerate subspace. Reshape to a (R-1, L-1) matrix (f4_var's
  # vec ordering puts right-pop fast, left-pop slow -- see f4blocks_to_f4stats
  # in this file and arr3d_to_mat in utility.R) and take row norms to get
  # per-right-pop loadings on the degenerate direction. Two right pops with
  # large opposing loadings = sister-like pair; one with a dominant loading
  # alone = the lone offender.
  rcond_concern = .rcond_concern  # see definition at top of file
  # Trigger SVD attribution whenever either condition holds:
  #   (a) the auto-concern bar is hit (always-on diagnostic at very low rcond), or
  #   (b) the caller asked for fail-loud gating and the threshold actually trips
  #       (so the error message can include the attribution).
  trigger_loadings = is.finite(f4_var_rcond) &&
    (f4_var_rcond < rcond_concern ||
     (!is.na(singular_threshold) && f4_var_rcond < singular_threshold))
  f4_var_singular_loadings = NULL
  if(trigger_loadings) {
    f4_var_singular_loadings = tryCatch({
      svd_res = svd(f4_var)
      v_min = svd_res$v[, ncol(svd_res$v)]                # smallest-sigma right singular vector
      v_mat = matrix(v_min, nrow = length(right) - 1)     # rows = right[-1], cols = left[-1]
      right_norms = sqrt(rowSums(v_mat^2))
      tibble::tibble(right = right[-1], loading = right_norms) %>%
        dplyr::arrange(dplyr::desc(loading))
    }, error = function(e) NULL)
  }

  if(!is.na(singular_threshold) && is.finite(f4_var_rcond) && f4_var_rcond < singular_threshold) {
    diag_msg = ""
    if(!is.null(f4_var_singular_loadings)) {
      top = utils::head(f4_var_singular_loadings, 5)
      diag_msg = paste0(
        "\n  Right pops loading on the degenerate direction (top 5):\n",
        paste0(sprintf("    %-40s loading %.3f", top$right, top$loading), collapse = "\n"),
        "\n  Pops at the top of this list are the likely offenders -- typically a\n",
        "  sister-clade pair (two pops with similar high loadings) or a lone right\n",
        "  pop too close to a source. Drop one and refit.")
    }
    stop(sprintf(
      paste0("f4 variance matrix is near-singular (rcond = %.3g < threshold %.3g).\n",
             "  This typically means right populations aren't independent of the\n",
             "  sources (e.g. sister populations).%s\n",
             "  (Pass `singular_threshold = NA` to allow the pseudo-inverse fallback\n",
             "  anyway; weight SEs may be unreliable.)"),
      f4_var_rcond, singular_threshold, diag_msg))
  }
  if(is.finite(f4_var_rcond) && f4_var_rcond < rcond_concern && is.na(singular_threshold) && verbose) {
    # Concerning rcond, but no threshold gate set. Surface a warning so the
    # user knows the result of the silent pseudo-inverse fallback might be
    # numerically tight in a way that doesn't reflect sampling uncertainty.
    top = if(!is.null(f4_var_singular_loadings)) utils::head(f4_var_singular_loadings, 3) else NULL
    diag_msg = if(is.null(top)) "" else paste0(
      " Right pops loading on the degenerate direction: ",
      paste0(sprintf("%s (%.2f)", top$right, top$loading), collapse = ", "), ".")
    warning(sprintf(
      paste0("f4 variance matrix is near-singular (rcond = %.3g).%s ",
             "Weight SEs may be unreliable; see `?qpadm` (`singular_threshold`) for ",
             "opt-in fail-loud gating; the returned list includes the full ",
             "right-pop loading table on the `f4_var_singular_loadings` field."),
      f4_var_rcond, diag_msg))
  }

  qinv = solve(f4_var)
  out = list(f4_var_rcond = f4_var_rcond,
             f4_var_singular_loadings = f4_var_singular_loadings)

  #----------------- compute admixture weights -----------------
  if(!is.null(target)) {
    if(verbose) alert_info('Computing admixture weights...\n')
    if(cpp) {
      qpadm_weights = cpp_qpadm_weights
      get_weights_covariance = cpp_get_weights_covariance
    }
    weight = qpadm_weights(f4_est, qinv, rnk = length(left)-2, fudge = fudge, constrained = constrained,
                           qpsolve = qpsolve)$weights %>% c
    if(getcov) {
      if(verbose) alert_info('Computing standard errors...\n')
      se = sqrt(diag(get_weights_covariance(f4_lo, qinv, block_lengths, fudge = fudge, boot = boot,
                                            constrained = constrained, qpsolve = qpsolve)))
    } else se = NA
    out$weights = tibble(target, left = left[-1], weight, se) %>% mutate(z = weight/se)

    wvec = out$weights %>% select(left, weight) %>% deframe
    if(!is.null(f2_blocks)) out$f4 = fitted_f4(f2_blocks, wvec, target, left[-1], right)
    else if(return_f4) out$f4 = fitted_f4_from_f4blockdat(f4blockdat_all, wvec, target, left[-1], right)
  } else {
    if(!is.null(f2_blocks)) out$f4 = f4(f2_blocks, left[1], left[-1], right[1], right[-1], verbose = FALSE)
    else if(return_f4) out$f4 = f4blockdat %>% f4blockdat_to_f4out(boot = boot)
  }

  #----------------- compute number of admixture waves -----------------
  if(verbose) alert_info('Computing number of admixture waves...\n')
  out$rankdrop = drop_ranks(f4_est, qinv, fudge, constrained, cpp)
  if(!is.null(target)) out$popdrop = drop_pops(f4_est, qinv, fudge, constrained, cpp, left[-1])

  out
}

#' Estimate admixture waves
#'
#' `qpwave` compares two sets of populations (`left` and `right`) to each other. It estimates a lower bound on the number of admixtue waves that went from `left` into `right`, by comparing a matrix of f4-statistics to low-rank approximations. For a rank of 0 this is equivalent to testing whether `left` and `right` form clades relative to each other.
#' @export
#' @inheritParams qpadm
#' @return `qpwave` returns a list with up to two data frames describing the model fit:
#' \enumerate{
#' \item `f4` A data frame with estimated f4-statistics
#' \item `rankdrop`: A data frame describing model fits with different ranks, including *p*-values for the overall fit
#' and for nested models (comparing two models with rank difference of one). A model with `L` left populations and `R` right populations has an f4-matrix of dimensions `(L-1)*(R-1)`. If no two left population form a clade with respect to all right populations, this model will have rank `(L-1)*(R-1)`.
#'   \itemize{
#'     \item `f4rank`: Tested rank
#'     \item `dof`: Degrees of freedom of the chi-squared null distribution: `(L-1-f4rank)*(R-1-f4rank)`
#'     \item `chisq`: Chi-sqaured statistic, obtained as `E'QE`, where `E` is the difference between estimated and fitted f4-statistics, and `Q` is the f4-statistic covariance matrix.
#'     \item `p`: p-value obtained from `chisq` as `pchisq(chisq, df = dof, lower.tail = FALSE)`
#'     \item `dofdiff`: Difference in degrees of freedom between this model and the model with one less rank
#'     \item `chisqdiff`: Difference in chi-squared statistics
#'     \item `p_nested`: *p*-value testing whether the difference between two models of rank difference 1 is significant
#'   }
#' }
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 10)
#' @seealso \code{\link{qpadm}}
#' @examples
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' qpwave(example_f2_blocks, left, right)
qpwave = function(data, left, right,
                  fudge = 0.0001, auto_only = TRUE, blgsize = 0.05, poly_only = FALSE, boot = FALSE,
                  constrained = FALSE, cpp = TRUE, singular_threshold = NA_real_, verbose = TRUE)
  qpadm(data = data, left = left, right = right, target = NULL,
        fudge = fudge, auto_only = auto_only, blgsize = blgsize, poly_only = poly_only,  boot = boot,
        constrained = constrained, cpp = cpp,
        singular_threshold = singular_threshold, verbose = verbose)



f2blocks_to_f4stats = function(f2_blocks, left, right, boot = FALSE) {

  f2blocks_to_f4blocks(f2_blocks, left, right) %>%
    f4blocks_to_f4stats(boot = boot)
}

f2blocks_to_f4blocks = function(f2_blocks, left, right) {

  if(length(right) < 2) stop('Not enough right populations!')
  if(length(left) < 2) stop('Not enough left populations!')

  nr = length(right) - 1
  nl = length(left) - 1
  f4_blocks = (f2_blocks[left[-1], rep(right[1], nr), , drop = FALSE] +
                 f2_blocks[rep(left[1], nl), right[-1], , drop = FALSE] -
                 f2_blocks[rep(left[1], nl), rep(right[1], nr), , drop = FALSE] -
                 f2_blocks[left[-1], right[-1], , drop = FALSE])/2
  dimnames(f4_blocks)[[2]] = right[-1]
  names(dimnames(f4_blocks))[1:2] = c(left[1], right[1])
  f4_blocks
}


f4blocks_to_f4stats = function(f4_blocks, boot = FALSE) {

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo_nafix)
  statfun = ifelse(boot, boot_pairarr_stats, jack_pairarr_stats)

  f4_lo = f4_blocks %>% samplefun
  block_lengths = parse_number(dimnames(f4_lo)[[3]])

  out = f4_lo %>% statfun(block_lengths)
  #out = f4_lo %>% jack_pairarr_stats(block_lengths)
  out$f4_lo = f4_lo
  out$block_lengths = block_lengths
  out
}

#' Turn f4 block data to 3d array
#'
#' @param f4blockdat f4 block data frame generated by \code{\link{f4blockdat_from_geno}}
#' @param remove_na Remove blocks with missing values
#' @export
f4blockdat_to_f4blocks = function(f4blockdat, remove_na = TRUE) {
  # assumes pop1 = left[1], pop2 = left[-1], pop3 = right[1], pop4 = right[-1]
  stopifnot(length(unique(f4blockdat$pop1)) == 1 && length(unique(f4blockdat$pop3)) == 1)
  stopifnot(all(c('block', 'n', paste0('pop', 1:4)) %in% names(f4blockdat)))

  f4blockdat %<>% arrange(block, pop4, pop2)
  p1 = f4blockdat$pop1[1]
  p3 = f4blockdat$pop3[1]
  p2 = unique(f4blockdat$pop2)
  p4 = unique(f4blockdat$pop4)
  bl = f4blockdat %>% select(block, n) %>%
    group_by(block) %>% slice_max(n, with_ties = F) %$% n %>% paste0('l', .)

  stopifnot(nrow(f4blockdat) == length(p2)*length(p4)*length(bl))

  arr = array(f4blockdat$est, c(length(p2), length(p4), length(bl)), list(p2, p4, bl) %>% purrr::set_names(c(p1, p3, 'bl')))
  if(remove_na) {
    keep = apply(arr, 3, function(x) sum(is.na(x)) == 0)
    arr = arr[,,keep, drop = FALSE]
    if(mean(keep) < 0.5) warning(paste0('Discarding ', sum(!keep), ' block(s) due to missing values!'))
  }
  arr
}

f4blockdat_to_f4out = function(f4blockdat, boot) {

  samplefun = ifelse(boot, function(...) est_to_boo_dat(...), est_to_loo_dat)
  datstatfun = ifelse(boot, boot_dat_stats, jack_dat_stats)
  popcomb = f4blockdat %>% filter(block == 1) %>% select(pop1:pop4)
  totn = f4blockdat %>%
    group_by(pop1, pop2, pop3, pop4) %>%
    summarize(n = sum(n, na.rm=T), .groups = 'drop')

  if(!'den' %in% names(f4blockdat)) {
    f4out = f4blockdat %>%
      group_by(pop1, pop2, pop3, pop4) %>%
      samplefun()
  } else {
    f4out = f4blockdat %>%
      rename(num = est) %>%
      pivot_longer(c(num, den), names_to = 'type', values_to = 'est') %>%
      group_by(type, pop1, pop2, pop3, pop4) %>%
      samplefun() %>%
      pivot_wider(names_from = type, values_from = c(est, loo)) %>%
      mutate(loo = loo_num/loo_den)
  }
  out = f4out %>%
    datstatfun %>%
    ungroup %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    transmute(pop1, pop2, pop3, pop4, est, se, z, p) %>%
    left_join(totn, by = paste0('pop', 1:4))
  popcomb %>% left_join(out, by = paste0('pop', 1:4))
}

f3blockdat_to_f3out = function(f3blockdat, boot) {

  if('numer' %in% names(f3blockdat)) {
    if(boot) warning('boot argument will being ignored')
    return(f3blockdat_to_f3out_numden(f3blockdat))
  }
  samplefun = ifelse(boot, function(...) est_to_boo_dat(...), est_to_loo_dat)
  datstatfun = ifelse(boot, boot_dat_stats, jack_dat_stats)
  totn = f3blockdat %>%
    group_by(pop1, pop2, pop3) %>%
    summarize(n = sum(n, na.rm=T))

  f3out = f3blockdat %>%
      group_by(pop1, pop2, pop3) %>%
      samplefun()

  f3out %>%
    datstatfun %>%
    ungroup %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    transmute(pop1, pop2, pop3, est, se, z, p) %>%
    left_join(totn, by = c('pop1', 'pop2', 'pop3'))
}

f3blockdat_to_f3out_numden = function(f3blockdat) {

  totn = f3blockdat %>%
    group_by(pop1, pop2, pop3) %>%
    summarize(n = sum(n, na.rm=T))

  loo = f3blockdat %>% select(-est) %>%
    pivot_longer(c(numer, denom), names_to = 'type', values_to = 'est') %>%
    group_by(pop1, pop2, pop3, type) %>%
    est_to_loo_dat() %>% select(-est) %>%
    pivot_wider(id_cols = c(pop1:pop3, block, n, length), names_from = 'type', values_from = 'loo') %>%
    mutate(loo = numer/denom)

  loo %>%
    jack_dat_stats %>%
    ungroup %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    transmute(pop1, pop2, pop3, est, se, z, p) %>%
    left_join(totn, by = c('pop1', 'pop2', 'pop3'))
}



lazadm_old = function(f2_data, left, right, target = NULL,
                      boot = FALSE, constrained = TRUE) {

  #----------------- prepare f4 stats -----------------
  if(is.null(target)) {
    left %<>% readLines
    right %<>% readLines
    target = left[1]
    left = left[-1]
  }
  pops = c(target, left, right)
  stopifnot(!any(duplicated(pops)))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo_nafix)
  f2_blocks = get_f2(f2_data, pops, afprod = TRUE) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  f2_mat = apply(f2_blocks, 1:2, weighted.mean, block_lengths)

  r = 1:length(right)
  og_indices = expand.grid(r, r, r) %>%
    filter(Var1 != Var2, Var1 != Var3, Var2 < Var3)

  pos1 = target
  pos2 = right[og_indices[,1]]
  pos3 = right[og_indices[,2]]
  pos4 = right[og_indices[,3]]

  y = f2_mat[cbind(pos1, pos4)] +
      f2_mat[cbind(pos2, pos3)] -
      f2_mat[cbind(pos1, pos3)] -
      f2_mat[cbind(pos2, pos4)]

  x = f2_mat[pos4, left] +
      f2_mat[cbind(pos2, pos3)] -
      f2_mat[pos3, left] -
      f2_mat[cbind(pos2, pos4)]

  #----------------- compute admixture weights -----------------
  lhs = crossprod(x, y)
  rhs = crossprod(x)
  nc = length(left)

  if(constrained) weight = qpsolve(rhs, lhs, diag(nc), rep(0, nc))
  else weight = solve(rhs, lhs)[,1]
  weight = weight/sum(weight)

  tibble(target, left, weight)
}


#' Estimate admixture weights
#'
#' Models target as a mixture of left populations, and outgroup right populations. Uses Lazaridis method
#' based non-negative least squares of f4 matrix.
#' @export
#' @inheritParams qpadm
#' @return `lazadm` returns a data frame with weights and standard errors for each left population
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European
#' languages in Europe.} Nature (SI 9)
#' @seealso \code{\link{qpadm}}
#' @examples
#' target = 'Denisova.DG'
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' lazadm(example_f2_blocks, left, right, target)
#' lazadm(example_f2_blocks, left, right, target, constrained = FALSE)
lazadm = function(data, left, right, target,
                  boot = FALSE, constrained = TRUE) {

  #----------------- prepare f4 stats -----------------
  if(is.null(target)) {
    left %<>% readLines
    right %<>% readLines
    target = left[1]
    left = left[-1]
  }
  pops = c(target, left, right)
  stopifnot(!any(duplicated(pops)))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo_nafix)
  statfun = ifelse(boot, boot_mat_stats, jack_mat_stats)
  f2_blocks = get_f2(data, pops, afprod = TRUE)
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])
  loo_blocks = f2_blocks %>% samplefun
  loo_blocks %<>% abind::abind(apply(f2_blocks, 1:2, weighted.mean, block_lengths))

  numblocks = length(block_lengths) + 1

  r = match(right, pops)
  og_indices = expand.grid(r, r, r) %>%
    filter(Var1 != Var2, Var1 != Var3, Var2 < Var3)
  ncomb = nrow(og_indices)
  nleft = length(left)
  blocknums = rep(1:numblocks, each = ncomb)

  pos1 = match(target, pops)
  pos2 = og_indices[,1]
  pos3 = og_indices[,2]
  pos4 = og_indices[,3]

  ymat = matrix(loo_blocks[cbind(pos1, pos4, blocknums)] +
                loo_blocks[cbind(pos2, pos3, blocknums)] -
                loo_blocks[cbind(pos1, pos3, blocknums)] -
                loo_blocks[cbind(pos2, pos4, blocknums)], ncomb)

  xarr = loo_blocks[pos4, left,] +
    array(loo_blocks[cbind(pos2, pos3, rep(blocknums, each = nleft))], c(ncomb, nleft, numblocks)) -
    loo_blocks[pos3, left,] -
    array(loo_blocks[cbind(pos2, pos4, rep(blocknums, each = nleft))], c(ncomb, nleft, numblocks))

  #----------------- compute admixture weights -----------------
  lhs = matrix(NA, nleft, numblocks)
  rhs = array(NA, c(nleft, nleft, numblocks))
  for(i in 1:numblocks) {
    lhs[,i] = crossprod(xarr[,,i], ymat[,i])
    rhs[,,i] = crossprod(xarr[,,i])
  }

  if(constrained) {
    Amat = cbind(1, -1, diag(nleft), -diag(nleft))
    bvec = c(1, -1, rep(0, nleft), -rep(1, nleft))
    fun = function(rhs, lhs) qpsolve(rhs, lhs, Amat, bvec)
  } else fun = function(rhs, lhs) solve(rhs, lhs)[,1]

  wmat = sapply(1:numblocks, function(i) {
    w = fun(rhs[,,i], lhs[,i,drop = FALSE])
    if(sum(w) > 1e-10) w = w/sum(w)
    w
  })
  stats = wmat[,-numblocks] %>% statfun(block_lengths, tot = wmat[,numblocks])
  weight = stats$est
  se = sqrt(diag(stats$var))

  tibble(target, left, weight, se) %>% mutate(z = weight/se)
}



drop_ranks = function(f4_est, qinv, fudge, constrained, cpp) {
  # drops rank and fits qpadm model

  rnk = nrow(f4_est) - 1
  fitrank = function(x) qpadm_fit(f4_est, qinv, x, fudge = fudge,
                                  constrained = constrained, cpp = cpp, addweights = FALSE)
  #rankdrop = map_dfr(rev(seq_len(max(1, rnk))-(rnk==0)), fitrank) %>%
  rankdrop = map_dfr(rnk:0, fitrank) %>%
    mutate(dofdiff = lead(dof)-dof, chisqdiff = lead(chisq)-chisq,
           p_nested = pchisq(chisqdiff, dofdiff, lower.tail = FALSE))
  rankdrop
}


drop_pops = function(f4_est, qinv, fudge, constrained, cpp, left) {
  # drops each subset of left populations and fits qpadm model
  #if(nrow(f4_est) == 1) return()
  fitpop = function(x, y) qpadm_fit(x, y, nrow(x)-1, fudge = fudge, constrained = constrained,
                                    cpp = cpp, addweights = TRUE)
  nc = ncol(f4_est)
  popdrop = left %>% power_set %>% enframe(name = 'i', value = 'pop') %>%
    #filter(map_int(pop, length) > 1) %>%
    mutate(f4mat = map(pop, ~f4_est[., , drop=F]),
           ind = map(pop, ~((rep(match(., left), each=nc)-1)*nc+(1:nc))),
           qinvs = map(ind, ~qinv[., ., drop=F]),
           fd = map2(f4mat, qinvs, fitpop)) %$% fd %>% bind_rows

  wmat = popdrop %>% select(-1:-4) %>% as.matrix
  pat = apply(is.na(wmat)+0, 1, paste0, collapse='')
  feasible = apply(wmat, 1, function(x) all(between(na.omit(x), 0, 1)))
  popdrop %<>% add_column(pat, feasible)

  child_patterns = function(pat) {
    if(length(pat) == 0) return(list())
    ones = str_locate_all(pat, '0')[[1]][,1]
    map(ones, ~{pat2 = pat; str_sub(pat2, ., .) = 1; pat2})
  }

  rnk = nrow(f4_est) - 1
  nested = popdrop %>% filter(f4rank == rnk-1) %>% select(pat, dof, chisq, feasible)
  for(i in rev(seq_len(max(0,rnk-1)))) {
    children = child_patterns(nested$pat)
    thisrnk = popdrop %>% filter(f4rank == i, pat %in% children, feasible) %>%
      top_n(1, -jitter(chisq)) %>% select(pat, dof, chisq, feasible)
    nested %<>% bind_rows(thisrnk)
  }
  nested %<>% arrange(-dof) %>%
    mutate(dofdiff = dof-lead(dof), chisqdiff = chisq-lead(chisq),
           p_nested = pchisq(chisqdiff, dofdiff, lower.tail = FALSE)) %>%
    transmute(best = TRUE, pat, dofdiff, chisqdiff, p_nested)

  popdrop %<>% left_join(nested, by = 'pat') %>% mutate(wt = rnk-f4rank) %>%
    select(pat, wt, dof, chisq, p, everything(),
           feasible, best, dofdiff, chisqdiff, p_nested) %>%
    arrange(dof, pat)
  popdrop
}



qpadm_evaluate_fit = function(xmat, qinv, A, B, f4rank) {

  res = t(xmat - A %*% B)
  chisq = (t(c(res)) %*% qinv %*% c(res))[,1]
  dof = qpadm_dof(nrow(A), ncol(B), f4rank)
  p = if(dof == 0) NA else pchisq(chisq, df = dof, lower.tail = FALSE)
  tibble(f4rank, dof, chisq, p)
}


qpadm_fit = function(xmat, qinv, rnk, fudge = 0.0001, iterations = 20,
                     constrained = FALSE, cpp = TRUE, addweights = FALSE) {
  # returns one-row data frame with fit of one qpadm model

  if(cpp) qpadm_weights = cpp_qpadm_weights
  if(rnk == 0) {
    fit = list(A = matrix(0, nrow(xmat)), B = t(matrix(0, ncol(xmat))), weights = 1)
  } else {
    fit = qpadm_weights(xmat, qinv, rnk, fudge = fudge, constrained = constrained, qpsolve = qpsolve)
  }
  out = qpadm_evaluate_fit(xmat, qinv, fit$A, fit$B, rnk)
  if(addweights && length(fit$weights) == nrow(xmat)) {
    out %<>% bind_cols(fit$weights %>% t %>% as_tibble(.name_repair = ~rownames(xmat)))
  }
  out
}

add_weighted_f2 = function(f2_blocks, weights) {

  nam = dimnames(f2_blocks)[[1]]
  npop = dim(f2_blocks)[1]
  nblocks = dim(f2_blocks)[3]
  stopifnot(all(names(weights) %in% nam))
  matchedweights = rep(0, npop)
  matchedweights[match(names(weights), nam)] = weights
  add = apply(f2_blocks * matchedweights, 2:3, sum, na.rm=T)
  f2_blocks = abind::abind(f2_blocks, add, along = 1)
  f2_blocks = abind::abind(f2_blocks, rbind(add, 0), along = 2)
  dimnames(f2_blocks)[[1]] = dimnames(f2_blocks)[[2]] = c(nam, 'fit')
  f2_blocks
}


fitted_f4 = function(f2_blocks, weights, target, left, right) {

  weights = weights/sum(weights)
  pops = names(weights)
  f2_blocks_plus = add_weighted_f2(f2_blocks, weights)
  fitf4 = f4(f2_blocks_plus, target, c(left, 'fit'), right, right, verbose = FALSE) %>% filter(pop3 != pop4)
  fitf4 %>% left_join(enframe(weights, name = 'pop2', value = 'weight'), by = 'pop2') %>%
    arrange(pop1, pop3, pop4, pop2)
}

fitted_f4_from_f4blockdat = function(f4blockdat, weights, target, left, right) {

  weights = weights[sort(names(weights))]
  weights = weights/sum(weights)
  fitf4 = f4blockdat %>% arrange(block, pop3, pop4, pop2) %>%
    summarize(est = sum(est*weights), pop1 = pop1[1], pop2 = 'fit', n = n[1], length = length[1], .by=c(pop3,pop4,block))
  bind_rows(f4blockdat, fitf4) %>%
    f4blockdat_to_f4out(boot = FALSE) %>%
    arrange(pop1, pop3, pop4, pop2)
}

#' Faster version of \code{\link{qpadm}} with reduced output
#'
#' Models target as a mixture of left populations, given a set of outgroup right populations.
#' Can be used to estimate admixture proportions, and to estimate the number of independent
#' admixture events.
#' @export
#' @inheritParams qpadm
#' @param f2_data Blocked f2-statistics (3d array), a directory path, or a genotype file prefix.
#' @param rnk Rank of f4-matrix. Defaults to one less than full rank.
#' @param weights Return weights (default = `FALSE`)
#' @return Data frame with `f4rank`, `dof`, `chisq`, `p`, `feasible`
#' @seealso \code{\link{qpadm}}
#' @examples
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' target = 'Denisova.DG'
#' qpadm_p(example_f2_blocks, left, right, target)
qpadm_p = function(f2_data, left, right, target = NULL, auto_only = TRUE, fudge = 0.0001, boot = FALSE,
                   constrained = FALSE, rnk = length(setdiff(left, target)) - 1, cpp = TRUE,
                   weights = FALSE, f4blocks = NULL) {

  if(is.null(f4blocks)) {
    if(!is.null(target)) left = c(target, setdiff(left, target))
    f2_blocks = get_f2(f2_data, pops = left, pops2 = right, auto_only = auto_only, afprod = TRUE)
    f4dat = f2blocks_to_f4stats(f2_blocks, left, right, boot = boot)
  } else {
    f4dat = f4blocks_to_f4stats(f4blocks)
    rnk = dim(f4blocks)[1] - 1
  }
  f4_est = f4dat$est
  f4_var = f4dat$var
  diag(f4_var) = diag(f4_var) + fudge*sum(diag(f4_var))
  qinv = solve(f4_var)
  out = qpadm_fit(f4_est, qinv, rnk, fudge = fudge,
                  constrained = constrained, cpp = cpp, addweights = TRUE)
  w = out %>% select(-1:-4) %>% as.matrix %>% c
  if(!weights) out = out %>% select(1:4)
  out %>%
    mutate(feasible = all(between(w, 0, 1)))
}


#' Test if two sets of populations form two clades
#'
#' A thin wrapper around \code{\link{qpadm_p}} with `rnk` set to zero
#' @export
#' @inheritParams qpadm
#' @param f2_data Blocked f2-statistics (3d array), a directory path, or a genotype file prefix.
test_cladality = function(f2_data, left, right, fudge = 0.0001, boot = FALSE, cpp = TRUE) {
  qpadm_p(f2_data, left, right, fudge = fudge, boot = boot, rnk = 0, cpp = cpp)
}


find_right = function(f2_blocks, target, pops) {
  f4p = f4(f2_blocks, verbose = FALSE) %>%
    select(pop1:pop4, p)
}


all_comb = function(pops, target = NULL, left = NULL, right = NULL) {
  stopifnot(length(union(pops, target)) > 5)
  stopifnot(is.null(left) || is.null(right) || is.null(target))
  if(is.null(target)) target = setdiff(pops, union(left, right))
  if(!is.null(left)) return(bind_rows(map_dfr(target, ~all_comb_right(., left, right = setdiff(pops, c(., left))))))
  if(!is.null(right)) return(bind_rows(map_dfr(target, ~all_comb_left(., left = setdiff(pops, c(., right)), right))))
  lr = map(target, ~setdiff(pops, .))
  lrcomb = map(lr, all_lr)
  len = map_dbl(lrcomb, ~length(.$left))
  t = rep(target, len)
  l = flatten(map(lrcomb, 1))
  r = flatten(map(lrcomb, 2))
  tibble(target = t, left = l, right = r)
}

all_comb_right = function(target, left, right) {
  nleft = length(left)
  r = power_set(right) %>% discard(~length(.) <= nleft)
  tibble(target, left = list(left), right = r)
}
all_comb_left = function(target, left, right) {
  nright = length(right)
  l = power_set(left) %>% keep(~between(length(.), 2, nright-1))
  tibble(target, left = l, right = list(right))
}

all_lr = function(pops) {
  # splits pops into all possible combinations of left and right, and unused
  npops = length(pops)
  left = power_set(pops) %>% discard(~length(.) >= npops/2 | length(.) < 2)
  right = map(left, ~power_set(setdiff(pops, .)))
  right2 = map2(left, right, ~discard(.y, function(x) length(x) <= length(.x)))
  l = rep(left, map_dbl(right2, length))
  r = flatten(right2)
  list(left = l, right = r)
}

all_lr2 = function(pops, rightfix = 0) {
  # splits pops into all possible combinations of left and right
  maxleft = min(length(pops), floor((length(pops) + rightfix - 1)/2))
  left = map(1:maxleft, ~combn(pops, ., simplify = F)) %>% flatten
  right = map(left, ~setdiff(pops, .))
  namedList(left, right)
}


#' Compute all pairwise qpwave p-values
#'
#' For all pairs of left populations qpwave rank 0 Chi-squared statistics and p-values will be computed. This tests for each pair of left populations whether it forms a clade with respect to the right populations.
#' @export
#' @param f2_data A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}. Alternatively, a directory with f2 statistics.
#' @param left Left populations
#' @param right Right populations
#' @examples
#' \dontrun{
#' left = pops[5:7]
#' right = pops[1:4]
#' f2_blocks = f2_from_precomp('/my/f2/dir/', pops = left, pops2 = right, afprod = TRUE)
#' qpwave_pairs(f2_blocks, left, right)
#' # If f2-stats are too big to load them into memory,
#' # the following will read the data for each pair from disk:
#' qpwave_pairs('/my/f2/dir/', left, right)
#' }
qpwave_pairs = function(f2_data, left, right) {
  expand_grid(pop1 = left, pop2 = left) %>%
    filter(pop1 < pop2) %>%
    mutate(out = furrr::future_map2(pop1, pop2, ~qpadm_p(f2_data, .y, right, .x, rnk = 0))) %>%
    unnest_wider(out) %>%
    select(pop1, pop2, chisq, p) %>%
    bind_rows(rename(., pop1 = pop2, pop2 = pop1)) %>%
    arrange(pop1, pop2)
}



#' Compute p-values for many qpadm models
#'
#' This functions evaluates many qpadm models simultaneously by keeping the target population
#' and the `rightfix` populations fixed, and distributing the `leftright` populations by keeping some
#' in the set of left population and adding the remaining populations to the right populations.
#' (See details for an example of how models are generated)
#' @export
#' @param f2_blocks 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' @param leftright Populations which will be distributed between left and right
#' @param rightfix Populations which will be on the right side in all models
#' @param target Target population
#' @param full_results Return all output items which are returned by \code{\link{qpadm}}. By default (`full_results = FALSE`), weights and several other statistics will not be computed for each model, making it faster and the output more readable. If `full_results = TRUE`, the output will be a nested data frame where each row is one `qpadm` model, and each column has one data frame item from the regular qpadm output (`weights`, `f4`, `rankdrop`, `popdrop`).
#' @param verbose Print progress updates
#' @details If `leftright` consists of the populations L1, L2, L3, L4; `rightfix` is the population R; and `target` is T,
#' the following models will be genrated: \cr\cr

#' (left), (right), (target) \cr
#' (L1), (L2, L3, L4, R), (T) \cr
#' (L2), (L1, L3, L4, R), (T) \cr
#' (L3), (L1, L2, L4, R), (T) \cr
#' (L4), (L1, L2, L3, R), (T) \cr
#' (L1, L2), (L3, L4, R), (T) \cr
#' (L1, L3), (L2, L4, R), (T) \cr
#' (L1, L4), (L2, L3, R), (T) \cr
#' (L2, L3), (L1, L4, R), (T) \cr
#' (L2, L4), (L1, L3, R), (T) \cr
#' (L3, L4), (L1, L2, R), (T) \cr
#'
#' @return A data frame with Chi-squared statistics and p-values for each population combination
#' @examples
#' \dontrun{
#' pops = dimnames(example_f2_blocks)[[1]]
#' qpadm_rotate(example_f2_blocks, leftright = pops[1:4],
#'              target = pops[5], rightfix = pops[6:7])
#' }
qpadm_rotate = function(f2_blocks, leftright, target, rightfix = NULL, full_results = FALSE, verbose = TRUE) {

  lr = all_lr2(leftright, length(rightfix))
  if(verbose) alert_info(paste0('Evaluating ', length(lr[[1]]), ' models...\n'))
  qpadm_eval_rotate(f2_blocks, target, lr, rightfix, full_results = full_results, verbose = verbose)
}

#' Rotate populations between left and right
#'
#' This functions creates a data frame with population combinations which can be used as the input for \code{\link{qpadm_multi}}
#' @export
#' @param leftright Populations which will be distributed between left and right
#' @param rightfix Populations which will be on the right side in all models
#' @param target Target population
#' @return A data frame with Chi-squared statistics and p-values for each population combination
#' @examples
#' \dontrun{
#' pops = dimnames(example_f2_blocks)[[1]]
#' rotate_models(leftright = pops[1:4],
#'              target = pops[5], rightfix = pops[6:7])
#' }
rotate_models = function(leftright, target, rightfix = NULL) {

  all_lr2(leftright, length(rightfix)) %>% as_tibble %>% rowwise %>%
    mutate(right = list(c(right, rightfix)), target = target) %>% ungroup
}


qpadm_eval_rotate = function(f2_blocks, target, leftright_dat, rightfix, full_results = FALSE, verbose = TRUE) {
  if(full_results) fun = function(...) qpadm(..., verbose = FALSE)
  else fun = qpadm_p
  leftright_dat %>%
    as_tibble %>%
    rowwise %>% mutate(right = list(c(rightfix, right))) %>% ungroup %>%
    mutate(res = furrr::future_map2(left, right, ~fun(f2_blocks, .x, .y, target),
                                    .progress = verbose, .options = furrr::furrr_options(seed = TRUE))) %>%
    unnest_wider(res) #%>%
  #mutate(chisq = map(rankdrop, 'chisq') %>% map_dbl(1)) %>%
  #arrange(chisq)
}




qpadmmodels_to_popcombs = function(models) {
  # takes data frame of left, right, target populations
  # returns all f4 popcombs numbered by model

  if('target' %in% names(models)) models$left = map2(models$left, models$target, ~union(.y, .x))
  models %>% as_tibble %>% mutate(model = 1:n()) %>% rowwise %>%
    transmute(model, pop1 = left[1], pop2 = list(left[-1]), pop3 = right[1], pop4 = list(right[-1])) %>%
    unnest_longer(pop2) %>% unnest_longer(pop4)
}


#' Run multiple qpadm models
#'
#' This function runs multiple qpadm models, re-using f4-statistics where possible. Supports parallel evaluation of models, which can be turned on with `future::plan('multisession')` or similar, and turned off with `future::plan('sequential')`.
#' @export
#' @param data The input data in the form of:
#' \itemize{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}}
#' \item A directory with f2 statistics
#' \item The prefix of a genotype file
#' }
#' @param models A nested list (or data frame) with qpadm models. It should consist of two or three other named lists (or columns) containing `left`, `right`, (and `target`) populations.
#' @param allsnps Use all SNPs with allele frequency estimates in every population of any given population quadruple. If `FALSE` (the default) only SNPs which are present in all populations in `popcombs` (or any given model in it) will be used. When there are populations with lots of (non-randomly) missing data, `allsnps = TRUE` can lead to false positive results. This option only has an effect if `data` is the prefix of a genotype file. If `data` are f2-statistics, the behavior will be determined by the options that were used in computing the f2-statistics.
#' @param full_results Return the full qpadm output list for each model; if `FALSE`, return only summary statistics (default `TRUE`). Note: the lean (`FALSE`) path dispatches to [qpadm_p()], which omits `f4_var_rcond` and `f4_var_singular_loadings` from each per-fit result.
#' @param verbose Print progress updates
#' @param ... Further arguments forwarded to downstream functions. When `data`
#'   is a genotype-file prefix, kwargs matching [f4blockdat_from_geno()]'s
#'   formals are routed to the preflight f4 computation. Kwargs matching the
#'   per-fit function's formals (\code{\link{qpadm}} when `full_results=TRUE`,
#'   [qpadm_p()] when `FALSE`) are routed to each per-model fit.
#'   Unrecognised names raise a single error at entry naming both surfaces.
#'   `singular_threshold` (a [qpadm()] formal) is forwarded but be aware that
#'   if any model trips its threshold, [qpadm()] raises an error, the entire
#'   `qpadm_multi` call aborts via `furrr::future_map`, and no results are
#'   returned. For "mark and continue" semantics over a batch, omit
#'   `singular_threshold` and filter on `f4_var_rcond` after the call.
#' @return A list where each element is the output of one qpadm model.
#' @examples
#' \dontrun{
#' # the following specifies two models: one with 2/3/1 and one with
#' # 1/2/1 left/right/target populations
#' models = tibble(
#'            left = list(c('pop1', 'pop2'), c('pop3')),
#'            right = list(c('pop5', 'pop6', 'pop7'), c('pop7', 'pop8')),
#'            target = c('pop10', 'pop10'))
#' results = qpadm_multi('/my/geno/prefix', models)
#' }
qpadm_multi = function(data, models, allsnps = FALSE, full_results = TRUE, verbose = TRUE, ...) {

  if(!'left' %in% names(models) || !'right' %in% names(models))
    stop("'models' should have elements 'left' and 'right'!")
  if(length(unique(map_dbl(models, length))) != 1)
    stop("'left', 'right' (and 'target') should be of equal length!")

  if('target' %in% names(models)) models$left = map2(models$left, models$target, ~union(.y, .x))
  if(!'tibble' %in% class(models)) models %<>% as_tibble

  # Single-pass kwarg validation + dispatch. Selects the fit-function up
  # front so qpadm-only kwargs (when full_results=TRUE) and qpadm_p-only
  # kwargs (when FALSE) are correctly routed and validated. The dispatch
  # also validates that no reserved (qpadm_multi-internal) names appear in
  # `...` and that no geno-only kwargs are supplied when `data` is not a
  # genotype-file prefix (would otherwise be silently dropped).
  fit_fun     = if(full_results) qpadm else qpadm_p
  on_f2_path  = !is_geno_prefix(data)
  dispatch    = .qpadm_multi_dispatch_dots(list(...), fit_fun, on_f2_path)

  if(is_geno_prefix(data)) {
    popcombs = qpadmmodels_to_popcombs(models)
    f4blockdat = do.call(f4blockdat_from_geno,
                         c(list(data, popcombs, allsnps = allsnps, verbose = verbose),
                           dispatch$geno))
    # split(.$model) returns a named list. Names are kept here as per-fit
    # labels (used by the tryCatch wrapper below to identify which combo
    # failed). qpadm_sweep strips these names at its vapply layer so they
    # don't leak into output column attributes.
    f4blocks = f4blockdat %>% split(.$model) %>%
      furrr::future_map(quietly(f4blockdat_to_f4blocks),
                        .options = furrr::furrr_options(seed = TRUE)) %>%
      map('result')
  } else {
    if(verbose && allsnps) alert_warning('allsnps = TRUE is not effective when using precomputed f2 statistics\n')
    pops = models %>% select(any_of(c('target', 'left'))) %>% unlist %>% unique
    pops2 = models %>% select(right) %>% unlist %>% unique
    f2blocks = data %>% get_f2(pops, pops2, afprod = TRUE)
    f4blocks = models %>% rowwise %>% mutate(f4blocks = list(f2blocks_to_f4blocks(f2blocks, left, right))) %>% pull
  }
  if('target' %in% names(models)) models %<>% rowwise %>% mutate(left = list(setdiff(left, target))) %>% ungroup

  if(verbose) alert_info('Running models...\n')

  # Inner closure: only fit-fun formals are forwarded via dispatch$fit. The
  # tryCatch wrapper catches per-fit errors (e.g., singular_threshold trips
  # raising stop() inside qpadm()) and rewraps them with the model label
  # via rlang::abort(parent = e), which preserves the original condition
  # chain (class, call) so caller-side class-specific tryCatch handlers
  # still fire. `data` slot is named explicitly (not positional) so that
  # the names-strip cannot mistakenly displace it through do.call's
  # positional-fill behavior.
  fit_call = function(.x, .label) {
    base_args = list(data = NULL, left = NULL, right = NULL,
                     target = TRUE, f4blocks = .x)
    if(full_results) base_args$verbose = FALSE
    tryCatch(
      do.call(fit_fun, c(base_args, dispatch$fit)),
      error = function(e)
        rlang::abort(paste0("qpadm fit failed for model '", .label, "': ",
                            conditionMessage(e)),
                     parent = e))
  }

  # f4blocks always has names on the geno path (from split(.$model)) and
  # never on the f2 path. Explicit length-aware check rather than %||%,
  # which only triggers on strict NULL — `names()` of an unnamed list is
  # NULL, but a future change that returns a names-attributed-but-empty
  # `character(0)` would fall through %||% with character(0), producing a
  # zero-length label vector that future_map2 silently accepts.
  labels = names(f4blocks)
  if(is.null(labels) || length(labels) != length(f4blocks))
    labels = as.character(seq_along(f4blocks))
  out = furrr::future_map2(f4blocks, labels, fit_call,
                           .progress = verbose,
                           .options = furrr::furrr_options(seed = TRUE))

  if(!full_results) out = bind_cols(models, out %<>% bind_rows)
  out
}


# Names that qpadm_multi consumes internally rather than via `...`.
# Passing any of these via `...` is a user error (they collide with the
# positional / named args qpadm_multi already supplies at do.call time —
# the round-4 dispatch helper allowed routing them silently, which round-5
# review surfaced as the source of confusing duplicate-name errors AND
# silent positional NULL-displacement bugs):
#
#   data / pref / f2_data  : data argument across the 3 downstream callees
#   left / right / target  : filled from `models`
#   f4blocks / popcombs    : computed internally
#   verbose                : qpadm_multi forces FALSE on the inner qpadm
#
# Each must be reachable through its proper argument (`data`, `models`, or
# qpadm_multi's own formals) rather than through `...`.
.QPADM_MULTI_RESERVED = c("data", "pref", "f2_data",
                          "left", "right", "target",
                          "f4blocks", "popcombs",
                          "verbose")


# Splits qpadm_multi's `...` into per-callee buckets and validates upfront.
#
# qpadm_multi forwards user kwargs to two downstream callees with disjoint
# and overlapping formal sets:
#
#   * `f4blockdat_from_geno()`  (the geno-prefix preflight). 15 formals.
#   * `fit_fun` (the inner per-fit call). Either `qpadm()` (18 formals) when
#     `full_results=TRUE`, or `qpadm_p()` (12 formals) when FALSE.
#
# Three validation layers, each raising a single clear error at
# qpadm_multi's entry point rather than allowing the bad kwarg to surface
# as a cryptic downstream do.call error:
#
#   1. Reserved names (see `.QPADM_MULTI_RESERVED`) are rejected — they
#      would collide with the positional / named args qpadm_multi binds
#      internally at do.call time.
#   2. Unknown names (not in the union of non-reserved geno + fit formals)
#      are rejected.
#   3. Geno-only names (in geno_known but not fit_known) are rejected when
#      `on_f2_path = TRUE` — they would otherwise be silently dropped
#      because the f2-data path never calls f4blockdat_from_geno.
#
# Returns a list with two buckets (`$geno`, `$fit`), each holding only the
# kwargs whose names match that callee's non-reserved formals. Each
# callsite uses `do.call(callee, c(positional_args, dispatch$bucket))`.
#
# Empty names (positional kwargs supplied via `...`) are ignored — qpadm_multi
# has no semantic for them and the documented contract is named-kwargs-only.
.qpadm_multi_dispatch_dots = function(dots, fit_fun, on_f2_path) {
  geno_known = setdiff(names(formals(f4blockdat_from_geno)), .QPADM_MULTI_RESERVED)
  fit_known  = setdiff(names(formals(fit_fun)),              .QPADM_MULTI_RESERVED)
  all_known  = union(geno_known, fit_known)

  named = if(is.null(names(dots))) rep("", length(dots)) else names(dots)

  # Layer 1: reserved names.
  reserved_in_dots = intersect(named, .QPADM_MULTI_RESERVED)
  if(length(reserved_in_dots) > 0)
    stop("Reserved argument(s) to qpadm_multi(): '",
         paste(reserved_in_dots, collapse = "', '"),
         "'. These are filled internally from `data`, `models`, or qpadm_multi's ",
         "own formals; pass them via the proper argument instead of via `...`.",
         call. = FALSE)

  # Layer 2: unknown names.
  unknown = setdiff(named[nzchar(named)], all_known)
  if(length(unknown) > 0)
    stop("Unknown argument(s) to qpadm_multi(): '",
         paste(unknown, collapse = "', '"),
         "'. Valid kwargs match the non-reserved formals of f4blockdat_from_geno() ",
         "(preflight, geno-path only) or the per-fit function (qpadm() when ",
         "full_results=TRUE, qpadm_p() when FALSE). See ?qpadm_multi.",
         call. = FALSE)

  # Layer 3: geno-only on f2 path (silent-drop prevention).
  if(isTRUE(on_f2_path)) {
    geno_only      = setdiff(geno_known, fit_known)
    geno_only_dots = intersect(named, geno_only)
    if(length(geno_only_dots) > 0)
      stop("Argument(s) '", paste(geno_only_dots, collapse = "', '"),
           "' apply only when `data` is a genotype file prefix. ",
           "On the f2-data path they would be silently ignored. ",
           "Remove them or pass a genotype prefix.", call. = FALSE)
  }

  list(geno = dots[named %in% geno_known],
       fit  = dots[named %in% fit_known])
}


#' Sweep qpadm over a Cartesian product of targets, source-sets, and right-sets
#'
#' Patterson-style sweeps fit qpadm for every combination of
#' (target, source-set, right-set). Each invocation through [qpadm()] would re-load
#' the f2 cache from disk; `qpadm_sweep()` loads it once via [qpadm_multi()] and
#' returns a flat tibble with one row per combination, suitable for
#' filtering / ranking model fits across a sweep without unnesting nested lists.
#'
#' This is a convenience wrapper around [qpadm_multi()] that adds:
#' * named source-sets and right-sets so each combination is labelled in the output
#' * implicit Cartesian product over `(targets x source_sets x right_sets)`
#' * a flat tibble result with top-level columns extracted from each model fit
#'
#' For per-model parallel evaluation, set `future::plan('multisession')` before calling.
#'
#' @export
#' @inheritParams qpadm_multi
#' @param targets Character vector of target populations.
#' @param source_sets Named list of character vectors; each element is a candidate
#'   set of source ("left") populations for one model. If unnamed, names default
#'   to `S1`, `S2`, .... Empty names are an error.
#' @param right_sets Named list of character vectors; each element is a candidate
#'   set of "right" / outgroup populations. If unnamed, names default to
#'   `R1`, `R2`, .... Empty names are an error.
#' @param full_results If `TRUE` (the default), the returned tibble includes
#'   list-columns `weights`, `rankdrop`, `popdrop`, and
#'   `f4_var_singular_loadings` with the full per-model output of [qpadm()].
#'   If `FALSE`, only the flat summary columns are returned (including the
#'   always-on `f4_var_rcond` diagnostic).
#' @return A tibble with one row per `(target, source_set, right_set)` combination
#'   and columns:
#'   \itemize{
#'     \item `target`, `source_set`, `right_set`: identifiers for the combination
#'     \item `left`, `right`: list-columns with the source / right pops for this model
#'     \item `f4rank`: tested rank in the top row of [qpadm()]'s `rankdrop` (= `length(left) - 1`)
#'     \item `p`, `chisq`, `dof`: top-row of `rankdrop` (the "auto" model fit)
#'     \item `feasible`: `TRUE` if all weights are between 0 and 1
#'     \item `f4_var_rcond`: scalar reciprocal-condition diagnostic of the f4
#'       variance matrix used in the GLS solve. Values near machine epsilon
#'       (`~1e-15`) signal near-singular `f4_var` and warrant inspection of
#'       `f4_var_singular_loadings`. Always returned (not gated on
#'       `full_results`).
#'     \item `weights`: list-column with the per-source `weight` / `se` / `z` tibble (`full_results = TRUE`)
#'     \item `rankdrop`: list-column with the full rankdrop table (`full_results = TRUE`)
#'     \item `popdrop`: list-column with [qpadm()]'s leave-one-out per-source-pop
#'       fit table (`full_results = TRUE`). Always populated since `qpadm_sweep`
#'       always supplies a target.
#'     \item `f4_var_singular_loadings`: list-column (`full_results = TRUE`).
#'       Cells are tibbles identifying which right population is
#'       loading-heaviest on the singular direction; populated when the
#'       auto-bar fires (`f4_var_rcond < 1e-8`). Cells are `NULL` when the
#'       auto-bar does not fire (the common clean-data case) and also `NULL`
#'       if [qpadm()]'s SVD inside the loadings computation fails on a
#'       pathological matrix. Useful for SVD-guided right-pop pruning.
#'   }
#'
#' Kwargs in `...` are forwarded to [qpadm_multi()], which validates them
#' upfront against the union of formals for [f4blockdat_from_geno()] (the
#' geno-path preflight) and [qpadm()] (the per-fit call), then routes each
#' to its respective downstream callee. Unrecognised names raise a single
#' clear error at entry — see [qpadm_multi()]'s `...` docs for the full
#' contract.
#'
#' `singular_threshold` (a [qpadm()] formal, forwarded via `...`) does not
#' produce a populated `f4_var_singular_loadings` cell in the sweep:
#' [qpadm()] raises an error when the threshold trips, propagates through
#' `furrr::future_map`, and aborts the entire sweep. For "mark and continue"
#' semantics, omit `singular_threshold` and filter on `f4_var_rcond` after
#' the sweep returns.
#' @seealso [qpadm()], [qpadm_multi()]
#' @examples
#' \dontrun{
#' # Run 3 targets * 2 source-sets * 2 right-sets = 12 qpadm models from one f2 dir.
#' targets    = c("Patterson_England_IA", "Patterson_England_BA",
#'                "Patterson_England_C_EBA")
#' sources    = list(canonical_3way = c("WHGA", "Balkan_N", "OldSteppe"),
#'                   with_ehg       = c("WHGA", "Balkan_N", "OldSteppe", "EHG_Karelia"))
#' rights     = list(distal_4pop    = c("OldAfrica", "WHGB", "Turkey_N", "Russia_Afanasievo"),
#'                   distal_refined = c("OldAfrica", "WHGB", "Turkey_N", "Russia_Afanasievo",
#'                                      "Iran_GanjDareh_N"))
#' res = qpadm_sweep("f2_dir/", targets, sources, rights)
#' res %>% arrange(target, p)
#' }
qpadm_sweep = function(data, targets, source_sets, right_sets,
                       allsnps = FALSE, full_results = TRUE, verbose = TRUE, ...) {

  if(length(targets) < 1) stop("'targets' must be a non-empty character vector")
  if(!is.list(source_sets) || length(source_sets) < 1)
    stop("'source_sets' must be a non-empty list of character vectors")
  if(!is.list(right_sets) || length(right_sets) < 1)
    stop("'right_sets' must be a non-empty list of character vectors")
  if(is.null(names(source_sets))) names(source_sets) = paste0("S", seq_along(source_sets))
  if(is.null(names(right_sets)))  names(right_sets)  = paste0("R", seq_along(right_sets))
  if(any(!nzchar(names(source_sets)))) stop("All 'source_sets' entries must be named")
  if(any(!nzchar(names(right_sets))))  stop("All 'right_sets' entries must be named")
  if(any(duplicated(names(source_sets)))) stop("'source_sets' names must be unique")
  if(any(duplicated(names(right_sets))))  stop("'right_sets' names must be unique")
  if(!all(vapply(source_sets, is.character, logical(1))))
    stop("All 'source_sets' entries must be character vectors")
  if(!all(vapply(right_sets, is.character, logical(1))))
    stop("All 'right_sets' entries must be character vectors")
  if(any(duplicated(targets)))
    warning("'targets' contains duplicates; the same model will be fit multiple times")

  combos = expand.grid(target     = as.character(targets),
                       source_set = names(source_sets),
                       right_set  = names(right_sets),
                       stringsAsFactors = FALSE,
                       KEEP.OUT.ATTRS = FALSE)
  combos$left  = unname(source_sets[combos$source_set])
  combos$right = unname(right_sets[combos$right_set])

  if(verbose) alert_info(paste0("Running ", nrow(combos),
                                " qpadm models (", length(targets), " target * ",
                                length(source_sets), " source-set * ",
                                length(right_sets), " right-set)...\n"))

  models = tibble::tibble(left   = combos$left,
                          right  = combos$right,
                          target = combos$target)
  # full_results = TRUE is required here so we can extract the summary columns
  # (p, chisq, feasible via weights) below. The outer full_results param controls
  # whether the list-columns survive into the returned tibble, not what qpadm_multi
  # computes internally.
  fits = qpadm_multi(data, models, allsnps = allsnps,
                     full_results = TRUE, verbose = verbose, ...)

  # Flatten each qpadm() result into the top row of its rankdrop + the weights tibble.
  top_row = function(x) if(is.null(x) || nrow(x) == 0) NULL else x[1, ]
  out = tibble::tibble(
    target     = combos$target,
    source_set = combos$source_set,
    right_set  = combos$right_set,
    left       = combos$left,
    right      = combos$right,
    f4rank     = vapply(fits, function(f) {
                   r = top_row(f$rankdrop); if(is.null(r)) NA_integer_ else as.integer(r$f4rank) },
                   integer(1)),
    p          = vapply(fits, function(f) {
                   r = top_row(f$rankdrop); if(is.null(r)) NA_real_ else as.numeric(r$p) },
                   numeric(1)),
    chisq      = vapply(fits, function(f) {
                   r = top_row(f$rankdrop); if(is.null(r)) NA_real_ else as.numeric(r$chisq) },
                   numeric(1)),
    dof        = vapply(fits, function(f) {
                   r = top_row(f$rankdrop); if(is.null(r)) NA_integer_ else as.integer(r$dof) },
                   integer(1)),
    feasible   = vapply(fits, function(f) {
                   w = f$weights; if(is.null(w)) NA else all(w$weight >= 0 & w$weight <= 1) },
                   logical(1)),
    # f4_var_rcond is a scalar diagnostic. It is surfaced on the always-on
    # flat-summary surface (alongside p / chisq / dof / feasible) rather than
    # gated on full_results because callers using full_results = FALSE for
    # smaller output tibbles still want to filter sweep rows by rank
    # deficiency. Note: full_results = FALSE only suppresses output columns;
    # the per-model qpadm() (including its f4_var_rcond computation) always
    # runs because qpadm_multi is called with full_results = TRUE above. Use
    # `[[` rather than `$` so a future longer slot name (e.g.
    # `f4_var_rcond_threshold`) can't partial-match silently. `r[[1]]`
    # rather than `as.numeric(r)` preserves any future names/attributes on
    # the upstream scalar (avoids attr-mismatch diagnostics in
    # expect_identical-based round-trip tests).
    f4_var_rcond = vapply(fits, function(f) {
                   r = f[['f4_var_rcond']]
                   if(is.null(r) || length(r) != 1L) NA_real_
                   else r[[1]] },
                   numeric(1)))

  # Strip names that may have leaked through the wrapper chain. qpadm_multi
  # keeps names on f4blocks to label tryCatch errors with the failing combo,
  # so vapply preserves them into the output tibble's atomic columns. Strip
  # them here so res$<col>[1] is identical to a direct qpadm()'s scalar
  # under attribute-strict comparison (expect_identical / identical()).
  # Programmatic loop over all atomic columns rather than a hardcoded list,
  # so future-added vapply columns get the same treatment automatically.
  for(col in names(out))
    if(is.atomic(out[[col]]) && !is.null(names(out[[col]])))
      names(out[[col]]) = NULL

  if(full_results) {
    # `unname(lapply(...))` strips the per-fit list names (same source as the
    # vapply-column names above: qpadm_multi keeps fit names for error
    # reporting). List-columns of name-bearing lists print confusingly and
    # break expect_identical against direct qpadm() returns.
    out$weights  = unname(lapply(fits, `[[`, 'weights'))
    out$rankdrop = unname(lapply(fits, `[[`, 'rankdrop'))
    # popdrop is qpadm()'s table of "what happens when each left pop is
    # dropped" — populated whenever target is non-null (always TRUE for
    # qpadm_sweep, which always sets target = combos$target). Surfaced as a
    # list-column for symmetry with the other diagnostic tables.
    out$popdrop = unname(lapply(fits, `[[`, 'popdrop'))
    # f4_var_singular_loadings is a tibble (or NULL) — gated on full_results
    # because of variable per-row payload. Each cell can be NULL for two
    # reasons: the rcond gate did not fire (the common case on clean data),
    # OR the gate fired but svd() inside qpadm errored on a pathological
    # f4_var matrix (qpadm wraps the SVD in tryCatch with NULL fallback).
    out$f4_var_singular_loadings = unname(lapply(fits, `[[`, 'f4_var_singular_loadings'))
  }

  # Sweep-level summary of rank deficiency. qpadm_multi forces verbose=FALSE
  # on each inner qpadm() call (so the per-model near-singular warning at
  # qpadm.R:317-322 never reaches the user). Surface a single sweep-level
  # alert here so callers running verbose=TRUE see how many rows tripped
  # the auto-bar without having to filter the f4_var_rcond column manually.
  # When full_results=FALSE the loadings list-column is absent — adapt the
  # alert text so it doesn't point users at a non-existent column.
  #
  # Routed through cli::cli_warn (stderr via the R condition system) rather
  # than alert_warning (cat to stdout). Stdout-mixing would corrupt the
  # captured data stream when callers redirect (Rscript ... > results.txt
  # 2> progress.txt), and warnings semantically belong on stderr.
  if(verbose) {
    n_concerning = sum(!is.na(out$f4_var_rcond) & out$f4_var_rcond < .rcond_concern)
    if(n_concerning > 0) {
      # cli >= 3.4 interprets `{.foo}` inside cli_warn/cli_inform glue as a
      # style reference (e.g., `{.code x}` styles `x` as code), so a dot-
      # prefixed variable name like `.rcond_concern` errors with "Invalid
      # cli literal: ... starts with a dot." Bind to a non-dot local before
      # interpolating. The {.code $f4_var_singular_loadings} substitution is
      # legitimate styling ({.code} is a real cli style).
      rcond_concern_local = .rcond_concern
      n_rows = nrow(out)
      loadings_hint = if(full_results)
        " Inspect {.code $f4_var_singular_loadings} for offending right pops."
      else
        " Re-run with {.code full_results = TRUE} to see {.code $f4_var_singular_loadings}."
      cli::cli_warn(c("!" = paste0("{n_concerning} of {n_rows} sweep row(s) ",
                                   "have f4_var_rcond < {rcond_concern_local} ",
                                   "(near-singular f4_var)."),
                      "i" = loadings_hint))
    }
  }

  out
}


# ---------------------------------------------------------------------------
# qpadm_with_pruning: SVD-guided automated right-pop pruning
# ---------------------------------------------------------------------------

#' Iteratively prune sister-clade right pops from a qpadm fit
#'
#' Wraps [qpadm()] in a loop that drops the highest-loading right population
#' from the f4-variance singular vector until the fit's reciprocal condition
#' number clears `singular_threshold`. Codifies the hand-curation defense
#' against near-singular f4 covariance matrices (the diagnostic surface from
#' [qpadm()]'s `f4_var_rcond` + `f4_var_singular_loadings` outputs) as a
#' first-class function with a load-bearing audit trail.
#'
#' Two pruning strategies are supported:
#'
#' \describe{
#'   \item{`"greedy"` (default)}{At each iteration, drop the single right pop
#'     with the largest absolute loading on the smallest right-singular vector
#'     of `f4_var`. Cost: one [qpadm()] fit per iteration. Almost always
#'     correct.}
#'   \item{`"lookahead"`}{At each iteration, evaluate the top
#'     `lookahead_top_j` candidate drops by actually performing the post-drop
#'     [qpadm()] fit, and pick the candidate that maximizes post-drop
#'     `f4_var_rcond`. Cost: `lookahead_top_j + 1` fits per iteration (e.g.
#'     4x greedy at the default `lookahead_top_j = 3`). Protects against the
#'     "sister-pair masking" edge case: when two right pops have near-equal
#'     loadings but very different post-drop rconds (e.g. one is the true
#'     sister-of-source while the other was sister-of-the-first), greedy
#'     `which.max` can pick the wrong one and chain into over-pruning;
#'     lookahead sees the post-drop landscape and picks the genuinely
#'     informative drop.}
#' }
#'
#' The pruner commits to a drop at each step and does not backtrack across
#' iterations. For tightly-constrained panels where a globally-optimal subset
#' matters more than runtime, a full branch-and-bound enumeration would be the
#' next escalation; that is out of scope for v1.
#'
#' **`right[1]` is the reference anchor and cannot be pruned.** [qpadm()]'s f4
#' parameterization is `f4(left[1], left[-1]; right[1], right[-1])`, and the
#' SVD attribution covers only `right[-1]`. A near-collinear `right[1]` would
#' not appear in the loadings table and the pruner cannot detect or remove it.
#' Convention is to place a phylogenetically-distant, unambiguous outgroup
#' (e.g. `Chimp.REF`) at `right[1]`. If you suspect `right[1]` is the
#' offender, try rotating it with a different pop and re-running.
#'
#' **Performance note for genotype-prefix data.** When `data` is a file-system
#' prefix rather than a pre-computed `f2_blocks` array, each inner [qpadm()]
#' call reads from disk. With `lookahead_top_j = 3` and ten iterations, that
#' is up to 40 disk reads. Pre-compute with [extract_f2()] and pass the
#' resulting array for batch work.
#'
#' Defensive guards return early with `converged = FALSE` and a populated
#' `reason` slot (rather than guessing) if [qpadm()]'s diagnostic surface
#' returns something the pruner cannot interpret: an empty / `NULL`
#' loadings table, a named drop pop that's not in the active right set
#' (stale loading from a prior fit), the floor `min_right_pops` being hit
#' before convergence, or the iteration cap being reached.
#'
#' @export
#' @inheritParams qpadm
#' @param target Target population, or `NULL` for a qpwave-style (no-target)
#'   fit. Passed directly to [qpadm()]. Default `NULL`.
#' @param right Character vector of initial right ("outgroup") populations.
#'   `right[1]` is the reference anchor and cannot be pruned (see Details).
#'   Pruning drops from `right[-1]`; the final retained set is reported as
#'   `right_final` on the returned object.
#' @param singular_threshold Reciprocal condition number floor; the loop
#'   terminates with `converged = TRUE` once `f4_var_rcond` clears this bar.
#'   Unlike [qpadm()]'s `singular_threshold` (which raises an error on trip),
#'   here the threshold is the convergence criterion, not a fail-loud gate.
#'   Default `1e-12`. Values above `1e-8` trigger a warning because inner
#'   fits may not populate loadings in the `(1e-8, singular_threshold)` band.
#' @param min_right_pops Minimum size of the right set; the loop terminates
#'   with `converged = FALSE` and `reason = "floor"` if pruning would
#'   reduce `right` below this floor without converging. Defaults to
#'   `length(left) + (!is.null(target))`, the structural minimum for qpadm
#'   identification. When `reason = "floor"`, the returned fit corresponds to
#'   the last (near-singular) right set; inspect `$f4_var_rcond` before using
#'   the weights.
#' @param max_iterations Hard cap on iterations; defaults to
#'   `length(right) - min_right_pops + 1`. The default is the largest
#'   number of drops the floor allows, so it should never trip; this is a
#'   belt-and-suspenders bound against a logic bug.
#' @param strategy Pruning strategy: `"greedy"` (default) or `"lookahead"`.
#'   See Details.
#' @param lookahead_top_j Only used when `strategy = "lookahead"`. Number
#'   of top candidates (by absolute loading) to evaluate per iteration.
#'   Default `3`. Higher values trade runtime for thoroughness; values above
#'   ~5 rarely help on real panels.
#' @param verbose If `TRUE`, emits a `cli` info line per iteration naming
#'   the dropped pop and rcond. Default `FALSE` (batch-friendly).
#' @param ... Additional arguments forwarded to each inner [qpadm()] call.
#' @return A list with the slots of the final [qpadm()] fit (so `f4`,
#'   `rankdrop`, `f4_var_rcond`, `f4_var_singular_loadings` are present, plus
#'   `weights` and `popdrop` when `target` is non-`NULL`; on the qpwave-style
#'   `target = NULL` path [qpadm()] omits `weights` and `popdrop`) plus four
#'   pruner-specific slots. As a special case, when `reason = "inner_error"`
#'   the qpadm slots are absent entirely (`fit` was `NULL`); only the four
#'   pruner-specific slots below are returned. The four slots are:
#' \itemize{
#'   \item `pruning_trail`: tibble with one row per drop, columns `iteration`
#'     (integer), `dropped` (character, the pop dropped), `loading` (double,
#'     absolute loading of the dropped pop on the smallest singular vector
#'     at drop time), `rcond_before` (double, `f4_var_rcond` of the fit
#'     just before the drop), `rcond_after` (double, `f4_var_rcond` of the
#'     fit after the drop is committed). Empty tibble if the initial fit
#'     was already well-conditioned.
#'   \item `right_final`: character vector of right pops in the final fit
#'     (i.e. `setdiff(right, dropped pops in trail)` in original order).
#'   \item `converged`: `TRUE` if the loop exited with `f4_var_rcond >=
#'     singular_threshold`, `FALSE` otherwise.
#'   \item `reason`: character scalar explaining the terminal state. One of
#'     `"converged"` (the success case), `"floor"` (would have dropped below
#'     `min_right_pops`; returned fit is near-singular), `"iter_cap"`
#'     (defensive; max_iterations reached), `"null_loadings"` (loadings table
#'     missing or empty; can happen when `singular_threshold > 1e-8`),
#'     `"stale_pop"` (the named drop pop was not in the active right set;
#'     defensive guard), `"inner_error"` (an inner [qpadm()] call errored;
#'     a warning is emitted and the partial trail is preserved).
#' }
#'
#'   For a converged fit, calling `qpadm(data, left, right_final, target, ...)`
#'   directly produces identical `weights`, `f4`, `rankdrop`, `popdrop`,
#'   and `f4_var_rcond` to the returned object (the pruner emits the same
#'   call shape for the terminal fit).
#' @seealso [qpadm()], [qpadm_multi()], [qpadm_sweep()]
#' @examples
#' \dontrun{
#' # Pruning loop converges in one step on a sister-pair fixture:
#' fit = qpadm_with_pruning(
#'   example_f2_blocks,
#'   left   = c("Russia_Ust_Ishim.DG", "Vindija.DG"),
#'   right  = c("Chimp.REF", "Altai_Neanderthal.DG", "Mbuti.DG",
#'              "Switzerland_Bichon.SG"),
#'   target = "Denisova.DG")
#' fit$converged          # TRUE
#' fit$pruning_trail      # tibble: which pop was dropped, when, why
#' fit$right_final        # the converged right set
#' }
qpadm_with_pruning = function(data, left, right, target = NULL,
                              singular_threshold = 1e-12,
                              min_right_pops = NULL,
                              max_iterations = NULL,
                              strategy = c("greedy", "lookahead"),
                              lookahead_top_j = 3,
                              verbose = FALSE,
                              ...) {

  strategy = match.arg(strategy)

  # ---- argument validation ----
  if(!is.numeric(singular_threshold) || length(singular_threshold) != 1 ||
     is.na(singular_threshold) || singular_threshold <= 0)
    stop("'singular_threshold' must be a single positive numeric. ",
         "Use a fail-loud value like 1e-12; see ?qpadm's singular_threshold.")

  # Warn when singular_threshold exceeds the auto-concern bar (1e-8). In that
  # band, inner qpadm() calls won't populate f4_var_singular_loadings (the
  # trigger_loadings condition inside qpadm() fires only at rcond < 1e-8 OR
  # when the inner singular_threshold trips -- and we force the inner threshold
  # to NA to avoid mid-loop errors). The pruner would exit with
  # reason = "null_loadings" rather than pruning, which is confusing. The
  # default 1e-12 is always below the concern bar, so this fires only for
  # unusual caller-specified thresholds.
  if(singular_threshold > .rcond_concern)
    cli::cli_warn(
      c("!" = paste0("'singular_threshold' = ", singular_threshold,
                     " is above the auto-concern bar (", .rcond_concern, "). ",
                     "Inner qpadm() fits may not populate loadings for rcond in ",
                     "(", .rcond_concern, ", ", singular_threshold, "), causing ",
                     "premature reason = 'null_loadings' exits."),
        "i" = "The default singular_threshold = 1e-12 is always below the concern bar."))

  # min_right_pops: default to the structural minimum for qpadm identification.
  # qpadm's f4 parameterization is f4(left[1], left[-1][i]; right[1], right[-1][j]),
  # so the f4_var matrix has shape (R-1) x (R-1). With L = length(left) + 1
  # (target prepended), the system is identified when R >= L, i.e. length(right)
  # >= length(left) + (!is.null(target)). We use this as the default floor.
  if(is.null(min_right_pops)) {
    min_right_pops = as.integer(length(left) + (!is.null(target)))
    # Never below 2: qpadm requires length(right) >= 2 to build f4_var at all
    # (right[-1] must be non-empty), so a floor of 1 would let the loop drop to
    # a 1-pop set that errors inside qpadm rather than tripping the floor guard.
    min_right_pops = max(min_right_pops, 2L)
  } else {
    if(!is.numeric(min_right_pops) || length(min_right_pops) != 1 ||
       is.na(min_right_pops) || min_right_pops < 1)
      stop("'min_right_pops' must be a single positive integer.")
    # Catch non-integer values (e.g. 3.7) that as.integer() would silently truncate.
    if(min_right_pops != as.integer(min_right_pops))
      stop("'min_right_pops' must be a whole number; ", min_right_pops, " is not.")
    min_right_pops = as.integer(min_right_pops)
  }

  if(!is.character(right) || length(right) < min_right_pops)
    stop("'right' must be a character vector of at least ", min_right_pops,
         " populations (min_right_pops). Got ", length(right), ".")

  if(is.null(max_iterations)) {
    max_iterations = length(right) - min_right_pops + 1L
  } else {
    if(!is.numeric(max_iterations) || length(max_iterations) != 1 ||
       is.na(max_iterations) || max_iterations < 1)
      stop("'max_iterations' must be a single positive integer or NULL.")
    if(max_iterations != as.integer(max_iterations))
      stop("'max_iterations' must be a whole number; ", max_iterations, " is not.")
    max_iterations = as.integer(max_iterations)
  }
  if(strategy == "lookahead") {
    if(!is.numeric(lookahead_top_j) || length(lookahead_top_j) != 1 ||
       is.na(lookahead_top_j) || lookahead_top_j < 1)
      stop("'lookahead_top_j' must be a single positive integer.")
    if(lookahead_top_j != as.integer(lookahead_top_j))
      stop("'lookahead_top_j' must be a whole number; ", lookahead_top_j, " is not.")
    lookahead_top_j = as.integer(lookahead_top_j)
  }

  # ---- single-shot qpadm wrapper -------------------------------------
  # Each inner call runs with verbose = FALSE (the loop's cli output owns
  # the user-facing surface). qpadm's `singular_threshold` is also forced
  # off because the OUTER convergence test treats it as a non-fatal floor,
  # while qpadm()'s singular_threshold raises an error on trip; we don't
  # want that to derail the loop. Both names are formals of this function,
  # so R's argument matching ensures users can't sneak them into `...`.
  user_dots = list(...)
  fit_one = function(right_now) {
    do.call(qpadm,
            c(list(data = data, left = left, right = right_now,
                   target = target, verbose = FALSE,
                   singular_threshold = NA_real_),
              user_dots))
  }

  # ---- pruning loop ----
  trail_rows = list()  # accumulated row-list; cheaper than incremental rbind
  right_now = right
  reason    = NA_character_
  fit       = NULL

  for(iter in seq_len(max_iterations)) {
    # Iteration 1: let an inner qpadm() error propagate. A failure on the
    # initial right set means the caller's inputs are bad (misspelled pop, too
    # few right pops, ...) and a loud error is far more useful than a shell
    # return with reason = "inner_error" and no qpadm fields.
    # Iteration 2+: a failure means a *drop* produced an unfittable set, so
    # catch it, warn, and return the partial trail rather than losing progress.
    if(iter == 1L) {
      fit = fit_one(right_now)
    } else {
      fit = tryCatch(fit_one(right_now), error = function(e) {
        cli::cli_warn(
          c("!" = paste0("[qpadm_with_pruning] inner qpadm() errored at iter ",
                         iter, "; returning partial trail."),
            "i" = conditionMessage(e)))
        NULL
      })
      if(is.null(fit)) {
        reason = "inner_error"
        break
      }
    }

    # Success?
    if(isTRUE(is.finite(fit$f4_var_rcond)) &&
       fit$f4_var_rcond >= singular_threshold) {
      reason = "converged"
      break
    }

    # Floor?
    if(length(right_now) <= min_right_pops) {
      reason = "floor"
      break
    }

    # Loadings table available?
    loadings = fit$f4_var_singular_loadings
    if(is.null(loadings) || nrow(loadings) == 0) {
      reason = "null_loadings"
      break
    }

    # Pick which pop to drop.
    if(strategy == "greedy") {
      # Guard against a degenerate loadings table where all entries are NA
      # (defensive; qpadm's SVD tryCatch returns NULL on error, so this should
      # be unreachable in practice, but a future qpadm change could alter that).
      if(all(is.na(loadings$loading))) {
        reason = "null_loadings"
        break
      }
      best_idx = which.max(abs(loadings$loading))
      pop_to_drop = loadings$right[best_idx]
      drop_loading = abs(loadings$loading[best_idx])
      rcond_after = NA_real_  # not measured for greedy at attribution time
    } else {
      # "lookahead": evaluate top-J candidates by absolute loading; pick the
      # one whose post-drop fit has the largest f4_var_rcond.
      ordered = loadings[order(abs(loadings$loading), decreasing = TRUE), ]
      n_cand = min(lookahead_top_j, nrow(ordered))
      candidates = ordered$right[seq_len(n_cand)]
      cand_loadings = abs(ordered$loading[seq_len(n_cand)])
      cand_rconds = vapply(candidates, function(p) {
        cand_right = setdiff(right_now, p)
        # suppressWarnings: probe fits on a near-singular candidate may emit
        # qpadm's near-singular R-level warning; these are expected during
        # lookahead scoring and would flood the user's console.
        cand_fit = tryCatch(
          suppressWarnings(fit_one(cand_right)),
          error = function(e) NULL)
        if(is.null(cand_fit) || !is.finite(cand_fit$f4_var_rcond %||% NA_real_))
          NA_real_
        else
          cand_fit$f4_var_rcond
      }, numeric(1))
      # vapply names its result by the input candidates vector; strip so
      # rcond_after on the trail is a plain unnamed double per row.
      names(cand_rconds) = NULL
      if(all(is.na(cand_rconds))) {
        # All probes errored or returned non-finite rcond. Emit a breadcrumb
        # so the user knows lookahead degenerated, then fall back to the
        # loading-ranked (greedy) choice.
        cli::cli_warn(
          c("!" = paste0("[qpadm_with_pruning] lookahead iter ", iter,
                         ": all ", n_cand, " probe fit(s) failed; ",
                         "falling back to loading-ranked greedy choice."),
            "i" = paste0("Candidates: ",
                          paste(candidates, collapse = ", "))),
          class = "qpadm_with_pruning_lookahead_all_na")
        best_idx    = 1L
        rcond_after = NA_real_
      } else {
        best_idx    = which.max(cand_rconds)
        rcond_after = cand_rconds[best_idx]
      }
      pop_to_drop  = candidates[best_idx]
      drop_loading = cand_loadings[best_idx]
    }

    # Stale pop guard (defensive: the named drop must be in the active set).
    if(!pop_to_drop %in% right_now) {
      reason = "stale_pop"
      break
    }

    if(verbose) {
      cli::cli_inform(c("i" = paste0("[qpadm_with_pruning] iter ", iter,
                                     ": rcond = ",
                                     format(fit$f4_var_rcond, digits = 3),
                                     " < ", singular_threshold,
                                     ", dropping '", pop_to_drop, "' ",
                                     "(loading = ",
                                     format(drop_loading, digits = 3), ")")))
    }

    trail_rows[[iter]] = list(iteration = iter, dropped = pop_to_drop,
                              loading = drop_loading,
                              rcond_before = fit$f4_var_rcond,
                              rcond_after = rcond_after)
    right_now = setdiff(right_now, pop_to_drop)
  }

  # If we fell off the end of the loop without setting reason, the iteration
  # cap was reached. Unreachable with the default max_iterations (the floor
  # fires on the last allowed iteration), but reachable when the caller passes
  # an explicit max_iterations smaller than that default.
  #
  # The loop commits its drop at the *end* of each iteration body and only
  # fits at the *top*, so on a cap exit the last drop was applied to right_now
  # but never fitted: `fit` lags `right_final` by one drop. Re-fit on the final
  # right_now to re-synchronise, backfill the last trail row's rcond_after, and
  # reclassify as converged if that final drop happened to clear the threshold.
  if(is.na(reason)) {
    final_fit = tryCatch(fit_one(right_now), error = function(e) {
      cli::cli_warn(
        c("!" = paste0("[qpadm_with_pruning] reconciliation fit on the final ",
                       "right set errored after the iteration cap."),
          "i" = conditionMessage(e)))
      NULL
    })
    if(is.null(final_fit)) {
      # The final drop produced an unfittable set. Report inner_error; the
      # returned object carries no qpadm fields (see @return).
      fit = NULL
      reason = "inner_error"
    } else {
      fit = final_fit
      if(length(trail_rows) > 0)
        trail_rows[[length(trail_rows)]]$rcond_after = fit$f4_var_rcond
      reason = if(isTRUE(is.finite(fit$f4_var_rcond)) &&
                  fit$f4_var_rcond >= singular_threshold)
        "converged" else "iter_cap"
    }
  }

  # ---- assemble trail tibble ----
  if(length(trail_rows) == 0) {
    trail = tibble(iteration = integer(0),
                   dropped = character(0),
                   loading = double(0),
                   rcond_before = double(0),
                   rcond_after = double(0))
  } else {
    trail = bind_rows(trail_rows)
  }

  # ---- return ----
  # Floor exit note: when reason = "floor", the returned `fit` is the last
  # fit_one() result on the current right_now (the set that triggered the
  # floor check). It is near-singular by construction. Weights are present but
  # their SEs may be unreliable; inspect $f4_var_rcond and consider reducing
  # min_right_pops or expanding the right set.
  c(fit, list(pruning_trail = trail,
              right_final = right_now,
              converged = (reason == "converged"),
              reason = reason))
}

