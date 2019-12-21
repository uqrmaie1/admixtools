

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
  diag(coeffs) = diag(coeffs) + fudge
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
  diag(coeffs) = diag(coeffs) + fudge
  B2 = solve(coeffs, rhs[1,])
  #t(matrix(B2, ncol(xmat)))
  matrix(B2, ncol = ncol(xmat), byrow = TRUE)
}


get_weights_covariance = function(f4_blocks, qinv, cpp = FALSE, fudge = 0.0001) {

  rnk = dim(f4_blocks)[1]-1
  numblocks = dim(f4_blocks)[3]
  wmat = matrix(NA, numblocks, dim(f4_blocks)[1])
  if(cpp) qpadm_weights = cpp_qpadm_weights
  for(i in 1:numblocks) {
    wmat[i,] = qpadm_weights(as.matrix(f4_blocks[,,i]), qinv, rnk, fudge = fudge)
  }
  jackmeans = colMeans(wmat)
  mnc = jackmeans - t(wmat)
  (numblocks-1)/numblocks * (mnc %*% t(mnc))
}


qpadm_weights = function(xmat, qinv, rnk, fudge = 0.0001) {
  f4svd = svd(xmat)
  B = t(f4svd$v[, 1:rnk, drop=FALSE])
  A = xmat %*% t(B)
  for(i in 1:20) {
    A = opt_A(B, xmat, qinv, fudge = fudge)
    B = opt_B(A, xmat, qinv, fudge = fudge)
  }
  w = solve(t(cbind(A, 1)), c(rep(0, rnk), 1))
  w/sum(w)
}


#' Estimate admixture weights
#'
#' Models target as a mixture of left populations, and outgroup right populations.
#' @export
#' @param target the target population
#' @param left the source population
#' @param right outgroup populations
#' @param f2_blocks 3d array of block-jackknife leave-one-block-out estimates of f2 statistics. output of \code{\link{afs_to_f2_blocks}}. they are weighted by inverse of outgroup heterozygosity, if outgroup was specified.
#' @param block_lengths the jackknife block lengths used in computing the f2 statistics. see \code{\link{get_block_lengths}}.
#' @param f2_dir a directory with f2 statistics for each population pair in the graph. must contain 'block_lengths.RData'.
#' @param fudge value added to diagonal matrix elements before inverting
#' @param getcov should standard errors be returned? Setting this to FALSE makes this function much faster.
#' @param cpp should optimization be done using C++ or R function? cpp = TRUE is much faster.
#' @return a data frame with weights and standard errors for each left population
#' @examples
#' target = 'Denisova.DG'
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' qpadm(target, left, right, f2_blocks, block_lengths)
qpadm = function(target, left, right, f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL, fudge = 0.0001, getcov = TRUE, cpp = TRUE) {

  allpops = c(target, left, right)

  #----------------- read f-stats -----------------
  if(is.null(f2_dir) & (is.null(f2_blocks) | is.null(block_lengths))) stop('You have to provide an f2_dir argument, or f2_blocks and block_lengths!')
  if(!is.null(f2_dir) & (is.null(f2_blocks) | is.null(block_lengths))) {
    f2_blocks = read_f2(f2_dir, pops = allpops)
    load(paste0(f2_dir, '/block_lengths.RData'))
  }

  #----------------- process f-stats -----------------
  f2nam = dimnames(f2_blocks)[[1]]
  stopifnot(length(unique(c(target, left, right))) == length(allpops))
  stopifnot(all(allpops %in% f2nam))
  stopifnot(length(right) > length(left))

  nam = intersect(f2nam, allpops)
  f2_blocks = rray::rray(f2_blocks[nam,nam,], dim_names = list(nam, nam, NULL))

  f4_blocks = (f2_blocks[target, right[-1], ] +
                 f2_blocks[left, right[1], ] -
                 f2_blocks[target, right[1], ] -
                 f2_blocks[left, right[-1], ])/2

  sts = bj_pairarr_stats(f4_blocks, block_lengths)
  f4_jest = sts[[1]]
  f4_jvar = sts[[2]]

  #----------------- compute admixture weights -----------------
  diag(f4_jvar) = diag(f4_jvar) + fudge
  qinv = solve(f4_jvar)
  rnk = length(left)-1
  if(cpp) {
    qpadm_weights = cpp_qpadm_weights
    get_weights_covariance = cpp_get_weights_covariance
  }
  weight = qpadm_weights(f4_jest, qinv, rnk, fudge = fudge) %>% c
  if(getcov) se = sqrt(diag(get_weights_covariance(f4_blocks, qinv)))
  else se = rep(NA, length(weight))

  tibble(target, left, weight, se)
}


#' Wrapper function around the original qpAdm program
#'
#' This requires a working installation of qpAdm, which will be called using \code{\link{system}}
#'
#' @param target the target population
#' @param left the source population
#' @param right outgroup populations
#' @param bin path to the qpAdm binary file
#' @param pref path to and prefix of the packedancestrymap genotype files
#' @param outdir the output directory. files \code{out}, \code{parfile}, \code{leftlist}, \code{rightlist} will be overwritten
#' @param printonly should the command be printed or executed?
#' @return if not printonly, a data frame with parsed qpAdm output
#' @export
#' @examples
#' \dontrun{
#' target = 'Denisova.DG'
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' qpadm_wrapper(target, left, right,
#'   bin = 'path/to/qpAdm', pref = 'path/to/packedancestrymap_prefix',
#'   env = 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/blas/')
#' }
qpadm_wrapper = function(target, left, right, bin, pref, outdir='.', printonly=FALSE, env='') {

  parfile = paste0('genotypename: ', pref, '.geno\n',
                   'snpname: ', pref, '.snp\n',
                   'indivname: ', pref, '.ind\n',
                   'popleft: ', outdir, '/leftlist\n',
                   'popright: ', outdir, '/rightlist\n',
                   'details: YES\n',
                   'hashcheck: NO')

  write(parfile, paste0(outdir, '/parfile'))
  write(c(target, left), paste0(outdir, '/leftlist'))
  write(right, paste0(outdir, '/rightlist'))

  cmd = paste0(env,' ', bin, ' -p ', outdir, '/parfile > ', outdir, '/out')

  if(printonly) {
    print(cmd)
  } else {
    system(cmd)
    return(parse_qpadm_output(paste0(outdir, '/out')))
  }
}


