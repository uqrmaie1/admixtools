

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

# todo: this need to incorporate block_lengths; cpp version too
get_weights_covariance = function(f4_blocks, qinv, cpp = FALSE, fudge = 0.0001,
                                  constrained = FALSE, qpsolve = NULL) {

  rnk = dim(f4_blocks)[1]-1
  numblocks = dim(f4_blocks)[3]
  wmat = matrix(NA, numblocks, dim(f4_blocks)[1])
  if(cpp) qpadm_weights = cpp_qpadm_weights
  for(i in 1:numblocks) {
    wmat[i,] = qpadm_weights(as.matrix(f4_blocks[,,i]), qinv, rnk,
                             fudge = fudge, constrained = constrained, qpsolve = qpsolve)
  }
  jackmeans = colMeans(wmat)
  mnc = jackmeans - t(wmat)
  (numblocks-1)/numblocks * (mnc %*% t(mnc))
}


qpadm_weights = function(xmat, qinv, rnk, fudge = 0.0001,
                         iterations = 20, constrained = FALSE, qpsolve = NULL) {
  f4svd = svd(xmat)
  B = t(f4svd$v[, 1:rnk, drop=FALSE])
  A = xmat %*% t(B)
  for(i in 1:iterations) {
    A = opt_A(B, xmat, qinv, fudge = fudge)
    B = opt_B(A, xmat, qinv, fudge = fudge)
  }
  rhs = t(cbind(A, 1))
  lhs = c(rep(0, rnk), 1)
  #w = solve(t(cbind(A, 1)), c(rep(0, rnk), 1))
  nc = ncol(rhs)
  if(constrained) w = -qpsolve(rhs, lhs, -diag(nc), rep(0, nc))
  else w = solve(rhs, lhs)[,1]
  w/sum(w)
}

qpadm_f2 = function(target, left, right, f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL, f2_denom = 1) {

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
  f2_blocks = f2_blocks / f2_denom
  namedList(f2_blocks, block_lengths)
}

#' Estimate admixture weights
#'
#' Models target as a mixture of left populations, and outgroup right populations.
#' @export
#' @param target target population
#' @param left source populations
#' @param right outgroup populations
#' @param f2_blocks 3d array of block-jackknife leave-one-block-out estimates of f2 statistics. output of \code{\link{afs_to_f2_blocks}}
#' @param block_lengths the jackknife block lengths used in computing the f2 statistics. see \code{\link{get_block_lengths}}.
#' @param f2_dir a directory with f2 statistics for each population pair in the graph. must contain \code{block_lengths.RData}.
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param fudge value added to diagonal matrix elements before inverting
#' @param getcov should standard errors be returned? Setting this to \code{FALSE} makes this function much faster.
#' @param constrained if \code{FALSE} (default), admixture weights can be negative. if \code{TRUE}, they will all be non-negative, as in \code{\link{lazadm}}
#' @param cpp should optimization be done using C++ or R function? \code{cpp = TRUE} is much faster.
#' @return a data frame with weights and standard errors for each left population
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European languages in Europe.} Nature (SI 10)
#' @seealso \code{\link{lazadm}}
#' @examples
#' target = 'Denisova.DG'
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' qpadm(target, left, right, example_f2_blocks, example_block_lengths)
qpadm = function(target, left, right, f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL, f2_denom = 1, fudge = 0.0001, getcov = TRUE, constrained = FALSE, cpp = TRUE) {

  #----------------- prepare f4 stats -----------------
  f2dat = qpadm_f2(target, left, right, f2_blocks, block_lengths, f2_dir, f2_denom)
  f2_blocks = f2dat$f2_blocks
  block_lengths = f2dat$block_lengths

  f4_blocks = (f2_blocks[target, right[-1], ] +
               f2_blocks[left, right[1], ] -
               f2_blocks[target, right[1], ] -
               f2_blocks[left, right[-1], ])/2

  f4dat = bj_pairarr_stats(f4_blocks, block_lengths)
  f4_jest = f4dat$jest
  f4_jvar = f4dat$jvar

  #----------------- compute admixture weights -----------------
  diag(f4_jvar) = diag(f4_jvar) + fudge
  qinv = solve(f4_jvar)
  rnk = length(left)-1
  if(cpp) {
    qpadm_weights = cpp_qpadm_weights
    get_weights_covariance = cpp_get_weights_covariance
  }
  qpsolve = function(...) quadprog::solve.QP(...)$solution
  weight = qpadm_weights(f4_jest, qinv, rnk, fudge = fudge,
                         constrained = constrained, qpsolve = qpsolve) %>% c
  if(getcov) se = sqrt(diag(get_weights_covariance(f4_blocks, qinv, constrained = constrained,
                                                   qpsolve = qpsolve)))
  else se = rep(NA, length(weight))

  tibble(target, left, weight, se)
}


#' Wrapper function around the original qpAdm program
#'
#' This requires a working installation of qpAdm, which will be called using \code{\link{system}}
#'
#' @param target target population
#' @param left source populations (or leftlist file)
#' @param right outgroup populations (or rightlist file)
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
qpadm_wrapper = function(target = NULL, left = NULL, right = NULL, bin, pref = NULL,
                         outdir = './', parfile = NULL, printonly=FALSE, env='') {

  stopifnot(!is.null(parfile) & is.null(c(target, left, right)) |
            !is.null(pref) & !is.null(left) & !is.null(right))
  stopifnot(file.exists(str_replace(bin, '.+ ', '')))
  stopifnot(!is.null(target) | all(file.exists(c(left, right))))

  if(is.null(parfile)) {

    if(!is.null(target)) {
      leftfile = paste0(outdir, '/leftlist')
      rightfile = paste0(outdir, '/rightlist')
      write(c(target, left), leftfile)
      write(right, rightfile)
    } else {
      leftfile = left
      rightfile = right
    }

    pref = normalizePath(pref, mustWork = FALSE)
    parfile = paste0('genotypename: ', pref, '.geno\n',
                     'snpname: ', pref, '.snp\n',
                     'indivname: ', pref, '.ind\n',
                     'popleft: ', leftfile, '\n',
                     'popright: ', rightfile, '\n',
                     'details: YES\n',
                     'hashcheck: NO')

    parfilename = paste0(outdir, '/parfile')
    write(parfile, parfilename)

  } else {
    parfilename = parfile
  }

  cmd = paste0(env,' ', bin, ' -p ', parfilename, ' > ', outdir, '/out')

  if(printonly) {
    print(cmd)
  } else {
    system(cmd)
    return(parse_qpadm_output(paste0(outdir, '/out')))
  }
}


#' Estimate admixture weights
#'
#' Models target as a mixture of left populations, and outgroup right populations. Uses Lazaridis method based non-negative least squares of f4 matrix.
#' @export
#' @param target target population
#' @param left source populations
#' @param right outgroup populations
#' @param f2_blocks 3d array of block-jackknife leave-one-block-out estimates of f2 statistics. output of \code{\link{afs_to_f2_blocks}}
#' @param block_lengths the jackknife block lengths used in computing the f2 statistics. see \code{\link{get_block_lengths}}.
#' @param f2_dir a directory with f2 statistics for each population pair in the graph. must contain \code{block_lengths.RData}.
#' @param f2_denom scales f2-statistics. A value of around 0.278 converts F2 to Fst.
#' @param getcov should standard errors be returned? Currently not implemented.
#' @param constrained if \code{TRUE} (default), admixture weights will all be non-negative. if \code{FALSE}, they can be negative, as in \code{\link{qpadm}}
#' @return a data frame with weights and standard errors for each left population
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @references Haak, W. et al. (2015) \emph{Massive migration from the steppe was a source for Indo-European languages in Europe.} Nature (SI 9)
#' @seealso \code{\link{qpadm}}
#' @examples
#' target = 'Denisova.DG'
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' lazadm(target, left, right, example_f2_blocks, example_block_lengths)
#' lazadm(target, left, right, example_f2_blocks, example_block_lengths, constrained = FALSE)
lazadm = function(target, left, right, f2_blocks = NULL, block_lengths = NULL,
                  f2_dir = NULL, f2_denom = 1, getcov = FALSE, constrained = TRUE) {

  #----------------- prepare f4 stats -----------------
  f2dat = qpadm_f2(target, left, right, f2_blocks, block_lengths, f2_dir, f2_denom)
  f2_blocks = f2dat$f2_blocks
  block_lengths = f2dat$block_lengths

  f2_mat = apply(f2_blocks, 1:2, weighted.mean, block_lengths)

  ri = 1:length(right)
  og_indices = expand.grid(ri, ri, ri) %>%
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

  if(constrained) weight = -quadprog::solve.QP(rhs, lhs, -diag(nc), rep(0, nc))$solution
  else weight = solve(rhs, lhs)[,1]
  weight = weight/sum(weight)

  # todo: implement this
  if(getcov) se = rep(NA, length(weight))
  else se = rep(NA, length(weight))

  tibble(target, left, weight, se)
}
