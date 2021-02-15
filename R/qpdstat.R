

f2_f3 = function(f2_12, f2_13, f2_23) (f2_12 + f2_13 - f2_23) / 2
f2_f4 = function(f2_14, f2_23, f2_13, f2_24) (f2_14 + f2_23 - f2_13 - f2_24) / 2


#' Estimate f2 statistics
#'
#' Computes f2 statistics from f2 blocks of the form \eqn{f2(A, B)}
#' @export
#' @param data Input data in one of three forms: \enumerate{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}} (fastest option)
#' \item A directory which contains pre-computed f2-statistics
#' \item The prefix of genotype files (slowest option)
#' }
#' @param pop1 One of the following four:
#' \enumerate{
#' \item `NULL`: all possible population combinations will be returned
#' \item A vector of population labels. All combinations with the other `pop` arguments will be returned
#' \item A matrix with population combinations to be tested, with one population per column and one
#' combination per row. Other `pop` arguments will be ignored.
#' \item the location of a file (`poplistname` or `popfilename`) which specifies the populations or
#' population combinations to be tested. Other `pop` arguments will be ignored.
#' }
#' @param pop2 A vector of population labels
#' @param boot If `FALSE` (the default), block-jackknife resampling will be used to compute standard errors.
#' Otherwise, block-bootstrap resampling will be used to compute standard errors. If `boot` is an integer, that number
#' will specify the number of bootstrap resamplings. If `boot = TRUE`, the number of bootstrap resamplings will be
#' equal to the number of SNP blocks.
#' @param sure The number of population combinations can get very large. This is a safety option that stops you
#' from accidently computing all combinations if that number is large.
#' @param unique_only If `TRUE` (the default), redundant combinations will be excluded
#' @param verbose Print progress updates
#' @return `f2` returns a data frame with f2 statistics
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history} Genetics
#' @references Peter, B. (2016) \emph{Admixture, Population Structure, and F-Statistics} Genetics
#' @examples
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' f2(example_f2_blocks, pop1, pop2)
#' \dontrun{
#' f2(f2_dir, pop1, pop2)
#' }
f2 = function(data, pop1 = NULL, pop2 = NULL,
              boot = FALSE, sure = FALSE, unique_only = TRUE, verbose = FALSE) {

  out = fstat_get_popcombs(data, pop1 = pop1, pop2 = pop2,
                           sure = sure, unique_only = unique_only, fnum = 2)
  pops = unique(c(out$pop1, out$pop2))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, cpp_boot_vec_stats, cpp_jack_vec_stats)
  f2_blocks = get_f2(data, pops) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  #----------------- compute f2 -----------------
  if(verbose) alert_info('Computing f2-statistics\n')

  out %<>% group_by(pop1, pop2) %>%
    summarize(f2dat = list(f2_blocks[pop1, pop2, ])) %>% ungroup %>%
    mutate(sts = map(f2dat, ~statfun(., block_lengths)), est = map_dbl(sts, 'est'), var = map_dbl(sts, 'var')) %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    select(-f2dat, -var, -sts)

  out
}

#' Compute Fst
#'
#' This function reads Fst from a directory with precomputed f2-statistics, and turns per-block data
#' into estimates and standard errors for each population pair. See `details` for how Fst is computed.
#' @export
#' @inheritParams f2
#' @param data A directory which contains pre-computed f2- and fst-statistics
#' @details The Hudson Fst estimator used here is described in the two publications below.
#' For two populations with estimated allele frequency vectors `p1` and `p2`,
#' and allele count vectors `n1` and `n2`, it is calculated as follows:\cr\cr
#' `num = (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)`\cr
#' `denom = p1 + p2 - 2*p1*p2`\cr
#' `fst = mean(num)/mean(denom)`\cr\cr
#' This is done independently for each SNP block, and is stored on disk for each population pair.
#' Jackknifing or bootstrapping across these per-block estimates yields the overall estimates and standard errors.
#' @references Reich, D. (2009) \emph{Reconstructing Indian population history} Nature
#' @references Bhatia, G. (2013) \emph{Estimating and interpreting Fst: the impact of rare variants} Genome Research
#' @examples
#' \dontrun{
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' fst(f2_dir, pop1, pop2)
#' }
fst = function(data, pop1 = NULL, pop2 = NULL,
               boot = FALSE, verbose = FALSE) {

  out = fstat_get_popcombs(data, pop1 = pop1, pop2 = pop2,
                           sure = TRUE, unique_only = TRUE, fnum = 2)

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, cpp_boot_vec_stats, cpp_jack_vec_stats)
  f2_blocks = f2_from_precomp(data, pops = pop1, pops2 = pop2, fst = TRUE) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  out %>% group_by(pop1, pop2) %>%
    summarize(f2dat = list(f2_blocks[pop1, pop2, ])) %>% ungroup %>%
    mutate(sts = map(f2dat, ~statfun(., block_lengths)), est = map_dbl(sts, 'est'), var = map_dbl(sts, 'var')) %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    select(-f2dat, -var, -sts)
}



#' Estimate f3 statistics
#'
#' Computes f3 statistics from f2 blocks of the form \eqn{f3(A; B, C)}. Equivalent to
#' \eqn{(f2(A, B) + f2(A, C) - f2(B, C)) / 2} and to \eqn{f4(A, B; A, C)}
#' @export
#' @param pop3 A vector of population labels
#' @inheritParams f2
#' @return `qp3pop` returns a data frame with f3 statistics
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history} Genetics
#' @references Peter, B. (2016) \emph{Admixture, Population Structure, and F-Statistics} Genetics
#' @aliases f3
#' @section Alias:
#' `f3`
#' @examples
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' pop3 = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' qp3pop(example_f2_blocks, pop1, pop2, pop3)
#' \dontrun{
#' qp3pop(f2_dir, pop1, pop2, pop3)
#' }
qp3pop = function(data, pop1 = NULL, pop2 = NULL, pop3 = NULL,
                  boot = FALSE, sure = FALSE, unique_only = TRUE, verbose = FALSE) {

  stopifnot(is.null(pop2) & is.null(pop3) | !is.null(pop2) & !is.null(pop3))
  stopifnot(!is_geno_prefix(data) || !is.null(pop1))

  out = fstat_get_popcombs(data, pop1 = pop1, pop2 = pop2, pop3 = pop3,
                           sure = sure, unique_only = unique_only, fnum = 3)
  pops = unique(c(out$pop1, out$pop2, out$pop3))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, cpp_boot_vec_stats, cpp_jack_vec_stats)
  f2_blocks = get_f2(data, pops) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  #----------------- compute f3 -----------------
  if(verbose) alert_info('Computing f3-statistics\n')

  out %<>% group_by(pop1, pop2, pop3) %>%
    summarize(f3dat = list(f2_f3(f2_blocks[pop1, pop2, ],
                                 f2_blocks[pop1, pop3, ],
                                 f2_blocks[pop2, pop3, ]))) %>% ungroup %>%
    mutate(sts = map(f3dat, ~statfun(., block_lengths)), est = map_dbl(sts, 'est'), var = map_dbl(sts, 'var')) %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    select(-f3dat, -var, -sts)

  out
}
#' @export
f3 = qp3pop


#' Estimate f4 statistics
#'
#' Computes f4 statistics from f2 blocks of the form \eqn{f4(A, B; C, D)}. Equivalent to
#' \eqn{(f2(A, D) + f2(B, C) - f2(A, C) - f2(B, D)) / 2}
#' @export
#' @inheritParams f2
#' @param pop3 A vector of population labels
#' @param pop4 A vector of population labels
#' @param comb Generate all combinations of `pop1`, `pop2`, `pop3`, `pop4`. If `FALSE`, `pop1`, `pop2`, `pop3`, `pop4` should all be vectors of the same length.
#' @param blgsize SNP block size in Morgan. Default is 0.05 (50 cM). Only used when `data` is the prefix of genotype files
#' @param block_lengths Vector with lengths of each jackknife block. \code{sum(block_lengths)} has to
#' match the number of SNPs. only used when `data` is the prefix of genotype files
#' @param f4mode If `TRUE`: f4 is computed from allele frequencies `a`, `b`, `c`, and `d` as `(a-b)*(c-d)`. if `FALSE`, D-statistics are computed instead, defined as `(a-b)*(c-d) / ((a + b - 2*a*b) * (c + d - 2*c*d))`, which is the same as `(P(ABBA) - P(BABA)) / (P(ABBA) + P(BABA))`. `f4mode = FALSE` is only available when `data` is the prefix of genotype files
#' @param afprod Compute f4 from allele frequency products instead of f2. Only used if `data` is a directory with precomputed data.
#' @param cpp Use C++ functions. Setting this to `FALSE` will be slower but can help with debugging.
#' @param ... Additional arguments passed to \code{\link{f4blockdat_from_geno}} if `data` is a genotype file prefix or \code{\link{f2_from_precomp}} if `data` is a directory with f2-statistics
#' @return `qpdstat` returns a data frame with f4 statistics
#' @aliases f4
#' @section Alias:
#' `f4`
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history} Genetics
#' @references Peter, B. (2016) \emph{Admixture, Population Structure, and F-Statistics} Genetics
#' @examples
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' pop3 = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' pop4 = 'Switzerland_Bichon.SG'
#' qpdstat(example_f2_blocks, pop1, pop2, pop3, pop4)
#' \dontrun{
#' qpdstat(f2_dir, pop1, pop2, pop3, pop4)
#' }
qpdstat = function(data, pop1 = NULL, pop2 = NULL, pop3 = NULL, pop4 = NULL,
                   boot = FALSE, sure = FALSE, unique_only = TRUE,
                   comb = TRUE, blgsize = NULL, block_lengths = NULL, f4mode = TRUE,
                   afprod = TRUE, cpp = TRUE, verbose = is.character(data), ...) {

  stopifnot(is.null(pop2) & is.null(pop3) & is.null(pop4) |
            !is.null(pop2) & !is.null(pop3) & !is.null(pop4))

  if(!comb) {
    stopifnot(!is.null(pop2))
    stopifnot(length(unique(length(pop1), length(pop2), length(pop3), length(pop4))) == 1)
    out = tibble(pop1, pop2, pop3, pop4)
  } else {
    if(verbose) alert_info('Getting population combinations...\n')
    out = fstat_get_popcombs(data, pop1, pop2, pop3, pop4,
                             sure = sure, unique_only = unique_only, fnum = 4)
    if(verbose) alert_info(paste0(nrow(out), ' population combinations found\n'))
  }
  pops = unique(c(out$pop1, out$pop2, out$pop3, out$pop4))
  pops1 = unique(c(out$pop1, out$pop2))
  pops2 = unique(c(out$pop3, out$pop4))

  if(is_geno_prefix(data)) {
    if(verbose) alert_info('Computing from f4 from genotype data...\n')
    return(qpdstat_geno(data, out, blgsize = ifelse(is.null(blgsize), 0.05, blgsize),
                        f4mode = f4mode, block_lengths = block_lengths, boot = boot,
                        allsnps = TRUE, verbose = verbose, ...))
  }

  if(cpp) {
    boot_vec_stats = cpp_boot_vec_stats
    jack_vec_stats = cpp_jack_vec_stats
  }
  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, boot_vec_stats, jack_vec_stats)

  if(verbose) alert_info(paste0('Loading f2 data for ', length(pops1)*length(pops2), ' population pairs...\n'))
  f2_blocks = get_f2(data, pops1, pops2 = pops2, afprod = afprod, ...) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  #----------------- compute f4 -----------------
  if(verbose) alert_info('Computing f4-statistics\n')

  out %<>% rowwise %>%
    mutate(f4dat = list(f2_f4(f2_blocks[pop1, pop4, ],
                              f2_blocks[pop2, pop3, ],
                              f2_blocks[pop1, pop3, ],
                              f2_blocks[pop2, pop4, ]))) %>% ungroup %>%
    mutate(sts = map(f4dat, ~statfun(., block_lengths)), est = map_dbl(sts, 'est'), var = map_dbl(sts, 'var')) %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    select(-f4dat, -var, -sts)
  out
}

#' @export
f4 = qpdstat


fstat_get_popcombs = function(f2_data = NULL, pop1 = NULL, pop2 = NULL, pop3 = NULL, pop4 = NULL,
                              sure = FALSE, unique_only = TRUE, fnum = NULL) {
  # used by f2, f3, and f4 function
  # returns data frame 'out' with pop combs
  stopifnot(!is.null(pop1) | !is.null(f2_data))

  #----------------- make combinations -----------------
  out = NULL
  nam = c('pop1', 'pop2', 'pop3', 'pop4')[1:fnum]
  maxcomb = 1e6
  if(is_geno_prefix(f2_data) && is.null(pop1)) {
    if(is_ancestrymap_prefix(f2_data)) {
      indend = '.ind'
      popcol = 3
    } else {
      indend = '.fam'
      popcol = 1
    }
    pop1 = read_table2(paste0(f2_data, indend), col_types = cols(), col_names = FALSE, progress = FALSE)[[popcol]]
  }
  if(!is.null(pop2)) {
    ncomb = length(pop1) * length(pop2) * max(1, length(pop3)) * max(1, length(pop4))
    if(ncomb > maxcomb & !sure) {
      stop(paste0('If you really want to compute ', ncomb,
                  ' f-statistics, run this again with "sure = TRUE".'))
    }
    out = expand_grid(pop1, pop2, pop3, pop4)
  } else if(is.null(pop1)) {
    if(is.character(f2_data)) pop1 = list.dirs(f2_data, full.names=FALSE, recursive=FALSE)
    else pop1 = dimnames(f2_data)[[1]]
  } else if(is.character(pop1) && file.exists(pop1)) {
    pop1 = read_table2(pop1, col_names = FALSE)
    if(ncol(pop1) == 1) {
      pop1 = pop1[[1]]
    } else {
      out = pop1 %>% set_colnames(nam)
    }
  } else if('data.frame' %in% class(pop1) || is.matrix(pop1)) {
    if(ncol(pop1) != fnum) stop(paste0("Wrong number of columns in 'pop1'! (Is ",ncol(pop1)," should be ",fnum,")"))
    out = pop1 %>% set_colnames(nam) %>% as_tibble
  }

  if(is.null(out)) {
    # pop1 is a character vector at this point
    if(unique_only) {
      ncomb = choose(length(pop1), fnum)*(fnum-1)
      if(ncomb > maxcomb & !sure) {
        stop(paste0('If you really want to compute ', ncomb,
                    ' f-statistics, run this again with "sure = TRUE", or select your populations or combinations of interest.'))
      }
      if(length(pop1) < fnum) stop('Not enough populations!')
      outmat = t(combn(pop1, fnum))
      nr = nrow(outmat)
      if(fnum == 2) {
        out = outmat %>% set_colnames(nam) %>% as_tibble
      } else if(fnum == 3) {
        out = rbind(outmat[,1:3], outmat[,c(2,3,1)], outmat[,c(3,1,2)]) %>%
          set_colnames(nam) %>% as_tibble
      } else if(fnum == 4) out = rbind(outmat[,1:4], outmat[,c(1,3,2,4)], outmat[,c(1,4,2,3)]) %>%
        set_colnames(nam) %>% as_tibble() %>% slice(rep(1:nr, each=3) + (0:2)*nr)
      else {
        stop('fnum should be 2, 3, or 4!')
      }
    } else {
      # if(fnum^length(pop1) > maxcomb & !sure) {
      #   stop(paste0('If you really want to compute close to ', fnum^length(pop1),
      #               ' f-statistics, run this again with "sure = TRUE". Or specify more than just pop1.'))
      # }
      if(fnum == 2) out = expand_grid(pop1 = pop1, pop2 = pop1)
      else if(fnum == 3) out = expand_grid(pop1 = pop1, pop2 = pop1, pop3 = pop1)
      else if(fnum == 4) out = as.data.frame(t(combn(pop1, 2)), stringsAsFactors = FALSE) %>%
          expand_grid(x1=., x2=.) %>%
          {quietly(flatten_dfc)(.)$result} %>%
          set_colnames(paste0('pop', 1:4))
      else stop('fnum should be 2, 3, or 4!')
    }
  }
  out %>% distinct
}


qpfstats = function(f2_blocks) {

  pops = dimnames(f2_blocks)[[1]]
  pairs = t(combn(sort(pops), 2))
  pp = paste(pairs[,1], pairs[,2])

  a = f4(f2_blocks, verbose = FALSE) %>% select(pop1:pop4, est, se, z) %>% mutate(est = ifelse((pop1 > pop2) == (pop3 > pop4), est, -est), f13 = paste(pmin(pop1, pop3), pmax(pop1, pop3)), f24 = paste(pmin(pop2, pop4), pmax(pop2, pop4)), f14 = paste(pmin(pop1, pop4), pmax(pop1, pop4)), f23 = paste(pmin(pop2, pop3), pmax(pop2, pop3))) %>% expand_grid(pp) %>% mutate(coef = ifelse(pp == f14 | pp == f23, 1, ifelse(pp  == f13 | pp == f24, -1, 0))) %>% pivot_wider(pop1:z, names_from = pp, values_from = coef) %>%
    #mutate(est = -est*2) %>%
    select(-pop1:-pop4) %>% lm(as.formula(paste0('est ~ 0 + ', paste('`', colnames(.)[-1:-3], '`', sep = '', collapse = ' + '))), data = .)


  summary(a)$coefficients %>% as_tibble(rownames = 'pp') %>% mutate(pp = str_replace_all(pp, '`', '')) %>% separate(pp, c('pop1', 'pop2'), sep = ' ') %>% left_join(f2(f2_blocks[sort(pops), sort(pops),]), by = c('pop1', 'pop2')) %>% ggplot(aes(Estimate, est)) + geom_point() + geom_abline()


  ff2 = f2(f2_blocks[sort(pops), sort(pops),]) %>% bind_rows(rename(., pop1=pop2, pop2=pop1)) %>% rename(p1 = pop1, p2 = pop2)

}



gmat_to_aftable = function(gmat, popvec) {
  # raw genotype matrix, not corrected for ploidy, nind x nsnp
  rowsum(gmat, popvec, na.rm = TRUE) / rowsum((!is.na(gmat))+0, popvec) / 2
}


qpdstat_geno = function(pref, popcombs, blgsize = 0.05, block_lengths = NULL,
                        f4mode = TRUE, boot = FALSE, allsnps = FALSE, verbose = TRUE, ...) {

  pref = normalizePath(pref, mustWork = FALSE)
  f4blockdat = f4blockdat_from_geno(pref, popcombs, blgsize = blgsize, block_lengths = block_lengths,
                                    f4mode = f4mode, allsnps = allsnps, verbose = verbose, ...)

  if(verbose) alert_info('Summarize across blocks...\n')
  out = f4blockdat %>% f4blockdat_to_f4out(boot)
  popcombs %>% left_join(out, by = c('pop1', 'pop2', 'pop3', 'pop4'))
}


#' Get per-block f4-statistics
#'
#' This function turns per-block f2-statistics into per-block f4-statistics of the form `f4(pop1, pop2; pop3, pop4)`
#' @export
#' @param f2_data A 3d array with blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{f2_from_geno}}
#' Alternatively, a directory with precomputed data. See \code{\link{extract_f2}} and \code{\link{extract_counts}}.
#' @param pop1 Either the name(s) of the first population(s), or a four column matrix with the names of all four populations.
#' @param pop2 Population 2 (same length as `pop1`)
#' @param pop3 Population 3 (same length as `pop1`)
#' @param pop4 Population 4 (same length as `pop1`)
#' @return A matrix of per-block f4-statistics (`popcomb x block`)
f4_from_f2 = function(f2_data, pop1, pop2 = NULL, pop3 = NULL, pop4 = NULL) {
  if(is.matrix(pop1)) {
    stopifnot(is.null(c(pop2, pop3, pop4)))
    pop2 = pop1[,2]
    pop3 = pop1[,3]
    pop4 = pop1[,4]
    pop1 = pop1[,1]
  } else stopifnot(length(unique(c(length(pop1), length(pop2), length(pop3), length(pop4)))) == 1)
  f2_blocks = get_f2(f2_data, union(pop1, pop2), union(pop3, pop4), afprod = TRUE)
  map(seq_along(pop1), ~{
    f2_f4(f2_blocks[pop1[.], pop4[.],],
          f2_blocks[pop2[.], pop3[.],],
          f2_blocks[pop1[.], pop3[.],],
          f2_blocks[pop2[.], pop4[.],]) %>% unname
  }) %>% do.call(rbind, .)
}


#' Estimate admixture proportions via f4 ratios
#'
#' @export
#' @param data Input data in one of three forms:
#' \enumerate{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{extract_f2}} (fastest option)
#' \item A directory which contains pre-computed f2-statistics
#' \item The prefix of genotype files (slowest option)
#' }
#' @param pops A vector of 5 populations or a five column population matrix.
#' The following ratios will be computed: `f4(1, 2; 3, 4)/f4(1, 2; 5, 4)`
#' @param boot If `FALSE` (the default), block-jackknife resampling will be used to compute standard errors.
#' Otherwise, block-bootstrap resampling will be used to compute standard errors. If `boot` is an integer, that number
#' will specify the number of bootstrap resamplings. If `boot = TRUE`, the number of bootstrap resamplings will be
#' equal to the number of SNP blocks.
#' @param verbose Print progress updates
#' @return `qpf4ratio` returns a data frame with f4 ratios
qpf4ratio = function(data, pops, boot = FALSE, verbose = FALSE) {

  if(!is.matrix(pops)) pops %<>% t
  if(ncol(pops) != 5) stop("'pops' should be a vector of length 5, or a matrix with 5 columns.")

  f2_blocks = get_f2(data, pops, afprod = TRUE, verbose = verbose)

  samplefun = ifelse(boot, function(x, ...) est_to_boo(x, boot, ...), est_to_loo)
  statfun = ifelse(boot, boot_mat_stats, jack_mat_stats)

  block_lengths = parse_number(dimnames(f2_blocks)[[3]])
  f4_num = f4_from_f2(f2_blocks, pops[,1], pops[,2], pops[,3], pops[,4])
  f4_den = f4_from_f2(f2_blocks, pops[,1], pops[,2], pops[,5], pops[,4])

  thresh = 1e-6
  setmiss = abs(f4_den) < thresh
  f4_den[setmiss] = f4_num[setmiss] = NA
  totnum = weighted_row_means(f4_num, block_lengths, na.rm = TRUE)
  totden = weighted_row_means(f4_den, block_lengths, na.rm = TRUE)
  tot = totnum/totden
  f4_num_loo = f4_num %>% samplefun(block_lengths)
  f4_den_loo = f4_den %>% samplefun(block_lengths)
  if(boot) block_lengths = parse_number(dimnames(f4_num_loo)[[2]])
  stats = (f4_num_loo/f4_den_loo) %>% statfun(block_lengths, tot = tot)

  pops %>%
    as_tibble(.name_repair = ~paste0('pop', 1:5)) %>%
    mutate(alpha = stats$est, se = sqrt(diag(stats$var)), z = alpha/se)
}


