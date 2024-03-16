

f2_f3 = function(f2_12, f2_13, f2_23) (f2_12 + f2_13 - f2_23) / 2
f2_f4 = function(f2_14, f2_23, f2_13, f2_24) (f2_14 + f2_23 - f2_13 - f2_24) / 2


#' Estimate f2 statistics
#'
#' Computes f2 statistics from f2 blocks of the form \eqn{f2(A, B)}
#' @export
#' @param data Input data in one of three forms: \enumerate{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{f2_from_geno}} (fastest option)
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
#' @param ... Additional arguments passed to \code{\link{f2_from_geno}} when `data` is a genotype prefix
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
              boot = FALSE, sure = FALSE, unique_only = TRUE, verbose = FALSE, ...) {

  out = fstat_get_popcombs(data, pop1 = pop1, pop2 = pop2,
                           sure = sure, unique_only = unique_only, fnum = 2)
  pops = unique(c(out$pop1, out$pop2))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, cpp_boot_vec_stats, cpp_jack_vec_stats)
  f2_blocks = get_f2(data, pops, verbose = verbose, ...) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  #----------------- compute f2 -----------------
  if(verbose) alert_info('Computing f2-statistics\n')

  out %<>% group_by(pop1, pop2) %>%
    summarize(f2dat = list(f2_blocks[pop1, pop2, ])) %>% ungroup %>%
    mutate(sts = map(f2dat, ~statfun(., block_lengths)), est = map_dbl(sts, 'est'), var = map_dbl(sts, 'var')) %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    select(pop1, pop2, est, se)

  out
}

#' Compute Fst
#'
#' This function reads Fst from a directory with precomputed f2-statistics, and turns per-block data
#' into estimates and standard errors for each population pair. See `details` for how Fst is computed.
#' @export
#' @param data Input data in one of three forms:
#' \enumerate{
#' \item A 3d array of blocked Fst, output of \code{\link{f2_from_precomp}} with `fst = TRUE`
#' \item A directory which contains pre-computed Fst
#' \item The prefix of genotype files
#' }
#' @inheritParams f2
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
               boot = FALSE, verbose = FALSE, ...) {

  out = fstat_get_popcombs(data, pop1 = pop1, pop2 = pop2,
                           sure = TRUE, unique_only = TRUE, fnum = 2)

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, cpp_boot_vec_stats, cpp_jack_vec_stats)
  #f2fun = if(is_geno_prefix(data)) f2_from_geno else f2_from_precomp
  ell = list(...)
  ell$fst = TRUE
  ell$f2_data = data
  ell$pops = pop1
  ell$pops2 = pop2
  f2_blocks = do.call(get_f2, ell) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  out %>% group_by(pop1, pop2) %>%
    summarize(f2dat = list(f2_blocks[pop1, pop2, ])) %>% ungroup %>%
    mutate(sts = map(f2dat, ~statfun(., block_lengths)), est = map_dbl(sts, 'est'), var = map_dbl(sts, 'var')) %>%
    mutate(se = sqrt(var), z = est/se, p = ztop(z)) %>%
    select(pop1, pop2, est, se)
}



#' Estimate f3 statistics
#'
#' Computes f3 statistics of the form \eqn{f3(A; B, C)}. When using the same SNPs for all populations, this is equivalent to
#' \eqn{(f2(A, B) + f2(A, C) - f2(B, C)) / 2} and to \eqn{f4(A, B; A, C)}
#' @export
#' @inheritParams f2
#' @param pop3 A vector of population labels
#' @param ... Additional arguments passed to \code{\link{f3blockdat_from_geno}} if `data` is a genotype prefix, or to \code{\link{get_f2}} otherwise
#' @details
#' There are several arguments that can be passed via ... which affect the estimated f3-statistics.
#' The default options are the same as in the original qp3pop program,
#' but some options are not effective when using precomputed f2-statistics. See `examples` for more information.
#'
#'
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
#'
#' pops = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' qp3pop(example_f2_blocks, pops)
#' qp3pop(example_f2_blocks, pops, unique_only = FALSE)
#' \dontrun{
#' qp3pop(f2_dir, pop1, pop2, pop3)
#'
#' Below are three scenarios, and in each one `qp3pop()` and `qp3pop_wrapper()`
#' should give the same or very similar estimates. Note that to compute `f3(A; B, C)`,
#' `qp3pop_wrapper()` expects the populations to be in the order `B`, `C`, `A`.
#'
#' prefix = '/path/to/geno/prefix'
#' qp3popbin = '/path/to/AdmixTools/bin/qp3Pop'
#' pops = dimnames(example_f2_blocks)[[1]]

# target diploid
# outgroupmode NO (this is the default when passing a geno file prefix)
#' qp3pop_wrapper(prefix, pops[2], pops[3], pops[1], bin = qp3popbin, outgroupmode = FALSE)
#' qp3pop(prefix, pops[1], pops[2], pops[3])
#' qp3pop(prefix, pops[1], pops[2], pops[3], poly_only = TRUE)

# outgroupmode YES (this is the only option with precomputed f2-stats)
#' qp3pop_wrapper(prefix, pops[2], pops[3], pops[1], bin = qp3popbin, outgroupmode = TRUE)
#' qp3pop(prefix, pops[1], pops[2], pops[3], outgroupmode = TRUE)
#' f2b = f2_from_geno(prefix, pops = pops[1:3], poly_only = FALSE)
#' qp3pop(f2b, pops[1], pops[2], pops[3])

# target pseudodiploid (no heterozygotes means heterozygosity rate correction is not possible)
#' qp3pop_wrapper(prefix, pops[1], pops[3], pops[2], bin = qp3popbin, outgroupmode = TRUE)
#' qp3pop(prefix, pops[2], pops[1], pops[3], outgroupmode = TRUE, apply_corr = FALSE)
#' }
qp3pop = function(data, pop1 = NULL, pop2 = NULL, pop3 = NULL,
                  boot = FALSE, sure = FALSE, unique_only = TRUE,
                  blgsize = NULL, block_lengths = NULL, verbose = FALSE, ...) {

  stopifnot(is.null(pop2) & is.null(pop3) | !is.null(pop2) & !is.null(pop3))
  stopifnot(!is_geno_prefix(data) || !is.null(pop1))

  out = fstat_get_popcombs(data, pop1 = pop1, pop2 = pop2, pop3 = pop3,
                           sure = sure, unique_only = unique_only, fnum = 3)

  if(is_geno_prefix(data)) {
    if(verbose) alert_info('Computing from f3 from genotype data...\n')
    return(qp3pop_geno(data, out, blgsize = ifelse(is.null(blgsize), 0.05, blgsize),
                       block_lengths = block_lengths, boot = boot,
                       verbose = verbose, ...))
  }

  pops = unique(c(out$pop1, out$pop2, out$pop3))
  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, cpp_boot_vec_stats, cpp_jack_vec_stats)
  f2_blocks = get_f2(data, pops, verbose = verbose, ...) %>% samplefun
  if(is.null(block_lengths)) block_lengths = parse_number(dimnames(f2_blocks)[[3]])

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
#' Computes f4-statistics of the form \eqn{f4(A, B; C, D)}. For allele frequencies `a`, `b`, `c`, `d`, `f4(A, B; C, D)` is computed as the average of `(a-b)*(c-d)` over all SNPs in each SNP block. This is equivalent to \eqn{(f2(A, D) + f2(B, C) - f2(A, C) - f2(B, D)) / 2} (assuming no missing data). The input of this function can either be a 3d array of f2-statistics generated by \code{\link{f2_from_precomp}} or \code{\link{f2_from_geno}}, a directory with f2-statistics, or the prefix of genotype files. Computing f4 from genotype files directly is slower, but provides more flexibility in dealing with missing data (see details).
#' @export
#' @inheritParams f2
#' @param pop3 A vector of population labels
#' @param pop4 A vector of population labels
#' @param comb Generate all combinations of `pop1`, `pop2`, `pop3`, `pop4`. If `FALSE`, `pop1`, `pop2`, `pop3`, `pop4` should all be vectors of the same length.
#' @param blgsize SNP block size in Morgan. Default is 0.05 (5 cM). Only used when `data` is the prefix of genotype files
#' @param block_lengths Vector with lengths of each jackknife block. \code{sum(block_lengths)} has to
#' match the number of SNPs. only used when `data` is the prefix of genotype files
#' @param f4mode Set this to `FALSE` to compute D-statistics instead of f4. This only has an effect if the first argument is a genotype prefix. D-statistics are computed as `(a-b)*(c-d) / ((a + b - 2*a*b) * (c + d - 2*c*d))`, which is the same as `(P(ABBA) - P(BABA)) / (P(ABBA) + P(BABA))`
#' @param afprod Compute f4 from allele frequency products instead of f2 (default `TRUE`). Only used if `data` is a directory with precomputed data. For populations with lots of missing data, this option reduces bias that can result from setting `maxmiss` to values greater than 0. In all other cases it should not make a difference.
#' @param cpp Use C++ functions. Setting this to `FALSE` will be slower but can help with debugging.
#' @param ... Additional arguments passed to \code{\link{f4blockdat_from_geno}} if `data` is a genotype file prefix or \code{\link{f2_from_precomp}} if `data` is a directory with f2-statistics
#' @details f4- and D-statistics are informative about how four populations are related to one another. Estimates of f4 are unbiased as long as assumptions about SNP ascertainment and mutation rates are met. Missing data can violate the ascertainment assumptions: an f4-statistic may be significantly different when it is calculated from all non-missing SNPs, compared to what it would be if it were calculated from all SNPs in the genome. However, because this difference is often small, f4 is often calculated using samples or populations with missing data, on a particular subset of all SNPs. There are different strategies for choosing the SNPs in this case, and these strategies differ in how many SNPs they use, how likely they lead to bias, and whether pre-computed f2-statistics can be used.
#' * Use the same SNPs for every f4-statistic\cr
#' This is the most conservative option, but also the option which will use the smallest number of SNPs, which may result in a lack of power. It is the default option when pre-computing f2-statistics (`maxmiss = 0` in \code{\link{extract_f2}} or \code{\link{f2_from_geno}}).
#' * Use different SNPs for each f4-statistic\cr
#' This option strikes a balance between avoiding bias and using a larger number of SNPs. For each f4-statistic it selects all SNPs which are present in all four populations. This option only works when the first argument is a genotype prefix (it doesn't work with pre-computed f2-statistics). In that case, this option is used by default. To turn it off and instead use the same SNPs for every f4-statistic, set `allsnps = FALSE`. This option is the default option in the original qpDstat program (it is in fact the only option for selecting SNPs in the original qpDstat; however in the original qpAdm and qpGraph programs, this mode of selecting SNPs can be enabled with the `allsnps` option, whereas the default mode in qpAdm and qpGraph is equivalent to `maxmiss = 0`)
#' * Use different SNPs for each f2-statistic\cr
#' This is the least conservative option, but also the option which uses most of the available information. It makes it possible to use pre-computed f2-statistics for a large number of populations without losing a large number of SNPs. To use this option, set the `maxmiss` parameter to a value greater than 0 (and not larger than 1) in \code{\link{extract_f2}} or \code{\link{f2_from_geno}}. When using this option, be aware that bias is possible, in particular for f4-statistics where some populations have large amounts of missing data. To reduce the bias that can result from using this option, you may want to combine it with using the option `afprod = TRUE` in \code{\link{f2_from_precomp}}.
#' In summary, whenever you work with populations with missing data, there is no guarantee that f4- or D-statistics involving these populations are not skewed in some way. If you choose to analyze these populations anyway, and you decide which SNPs to use, there is a trade-off between maximizing power and minimizing the risk of bias. One strategy might be to first use the least conservative option (setting `maxmiss = 1` in \code{\link{extract_f2}}) to get an overview, and then spot-check individual results using more conservative options.
#'
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
#' \dontrun{
#' # compute D-statistics instead
#' qpdstat("/geno/prefix", pop1, pop2, pop3, pop4)
#' }
#' \dontrun{
#' # make a data frame with the population combinations for which f4 should be computed
#' combinations = tibble(pop1 = pop1, pop2 = pop2[1], pop3 = pop3, pop4 = pop4)
#' qpdstat(example_f2_blocks, combinations)
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
                        verbose = verbose, ...))
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
  # if(is_geno_prefix(f2_data) && is.null(pop1)) {
  #   if(is_ancestrymap_prefix(f2_data)) {
  #     indend = '.ind'
  #     popcol = 3
  #   } else {
  #     indend = '.fam'
  #     popcol = 1
  #   }
  #   pop1 = read_table2(paste0(f2_data, indend), col_types = cols(), col_names = FALSE, progress = FALSE)[[popcol]]
  # }
  if(!is.null(pop2)) {
    ncomb = length(pop1) * length(pop2) * max(1, length(pop3)) * max(1, length(pop4))
    if(ncomb > maxcomb & !sure) {
      stop(paste0('If you really want to compute ', ncomb,
                  ' f-statistics, run this again with "sure = TRUE".'))
    }
    out = expand_grid(pop1, pop2, pop3, pop4)
  } else if(is.null(pop1)) {
    if(is_precomp_dir(f2_data)) pop1 = list.dirs(f2_data, full.names=FALSE, recursive=FALSE)
    else if(is_plink_prefix(f2_data)) pop1 = unique(read_table2(paste0(f2_data,'.fam'), col_names=F, col_types = cols())$X1)
    else if(is_ancestrymap_prefix(f2_data)) pop1 = unique(read_table2(paste0(f2_data,'.ind'), col_names=F, col_types = cols())$X3)
    else pop1 = dimnames(f2_data)[[1]]
  } else if(is.character(pop1)[1] && file.exists(pop1[1])) {
    pop1 = read_table2(pop1, col_names = FALSE)
    if(ncol(pop1) == 1) {
      pop1 = pop1[[1]]
    } else {
      out = pop1 %>% set_colnames(nam)
    }
  } else if('data.frame' %in% class(pop1) || is.matrix(pop1)) {
    if(ncol(pop1) != fnum) stop(paste0("Wrong number of columns in 'pop1'! (Is ",ncol(pop1)," should be ",fnum,")"))
    if(is.matrix(pop1)) out = pop1 %>% set_colnames(nam) %>% as_tibble
    else out = pop1 %>% select(paste0('pop', 1:fnum))
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



gmat_to_aftable = function(gmat, popvec) {
  # raw genotype matrix, not corrected for ploidy, nind x nsnp
  rowsum(gmat, popvec, na.rm = TRUE) / rowsum((!is.na(gmat))+0, popvec) / 2
}


qpdstat_geno = function(pref, popcombs, blgsize = 0.05, block_lengths = NULL,
                        f4mode = TRUE, boot = FALSE, allsnps = TRUE, poly_only = FALSE, verbose = TRUE, ...) {

  pref = normalizePath(pref, mustWork = FALSE)
  f4blockdat = f4blockdat_from_geno(pref, popcombs, blgsize = blgsize, block_lengths = block_lengths,
                                    f4mode = f4mode, allsnps = allsnps, poly_only = poly_only, verbose = verbose, ...)

  if(verbose) alert_info('Summarize across blocks...\n')
  out = f4blockdat %>% f4blockdat_to_f4out(boot)
  popcombs %>% left_join(out, by = c('pop1', 'pop2', 'pop3', 'pop4'))
}

qp3pop_geno = function(pref, popcombs, blgsize = 0.05, block_lengths = NULL,
                       boot = FALSE, allsnps = TRUE, poly_only = FALSE, verbose = TRUE, ...) {

  if(!all(...names() %in% names(formals(f3blockdat_from_geno)))) {
    notused = setdiff(...names(), names(formals(f3blockdat_from_geno)))
    stop(paste0("The following arguments are not recognized: '", paste0(notused, collapse = "', '"), "'"))
  }
  pref = normalizePath(pref, mustWork = FALSE)
  f3blockdat = f3blockdat_from_geno(pref, popcombs, blgsize = blgsize, block_lengths = block_lengths,
                                    allsnps = allsnps, poly_only = poly_only, verbose = verbose, ...)

  if(verbose) alert_info('Summarize across blocks...\n')
  out = f3blockdat %>% f3blockdat_to_f3out(boot)
  popcombs %>% left_join(out, by = c('pop1', 'pop2', 'pop3'))
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
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{f2_from_geno}} (fastest option)
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

#' Estimate f4 differences
#'
#' @export
#' @param data Input data in one of three forms:
#' \enumerate{
#' \item A 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}} or \code{\link{f2_from_geno}} (fastest option)
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
#' @return `qpf4diff` returns a data frame with f4 ratios
qpf4diff = function(data, pops, boot = FALSE, verbose = FALSE) {

  if(!is.matrix(pops)) pops %<>% t
  if(ncol(pops) != 5) stop("'pops' should be a vector of length 5, or a matrix with 5 columns.")

  f2_blocks = get_f2(data, pops, apply_corr = FALSE, poly_only=FALSE, verbose = verbose)

  samplefun = ifelse(boot, function(x, ...) est_to_boo(x, boot, ...), est_to_loo)
  statfun = ifelse(boot, boot_mat_stats, jack_mat_stats)

  block_lengths = parse_number(dimnames(f2_blocks)[[3]])
  f4_1 = f4_from_f2(f2_blocks, pops[,1], pops[,2], pops[,3], pops[,4])
  f4_2 = f4_from_f2(f2_blocks, pops[,1], pops[,2], pops[,5], pops[,4])

  # tot1 = weighted_row_means(f4_1, block_lengths, na.rm = TRUE)
  # tot2 = weighted_row_means(f4_2, block_lengths, na.rm = TRUE)
  # tot = tot1-tot2
  f4diff_loo = (f4_1-f4_2) %>% samplefun(block_lengths)
  if(boot) block_lengths = parse_number(dimnames(f4_1)[[2]])
  stats = f4diff_loo %>% statfun(block_lengths)

  pops %>%
    as_tibble(.name_repair = ~paste0('pop', 1:5)) %>%
    mutate(f4diff = stats$est, se = sqrt(diag(stats$var)), z = f4diff/se, n = count_snps(f2_blocks))
}



#' Compute f4 from allele frequencies
#'
#' @export
#' @param afdat A data frame with allele frequencies and SNP metadata. Can be grouped.
#' @param popcombs A data frame with population combinations. Columns `pop1` to `pop4`
#' @examples
#' \dontrun{
#' # Compute f4 for all mutatation classes separately
#' afs = plink_to_afs('/my/geno/prefix', pops = c('p1', 'p2', 'p3', 'p4', 'p5'))
#' afdat = bind_cols(afs$snpfile, afs$afs %>% as_tibble()) %>%
#'         mutate(gr = paste0(pmin(A1, A2), pmax(A1, A2))) %>%
#'         group_by(gr)
#' popcombs = tibble(pop1 = c('p1', 'p5'), pop2 = 'p2', pop3 = 'p3', pop4 = 'p4')
#' out = f4_from_afdat(afdat, popcombs)
#' out %>% ggplot(aes(gr, est)) + geom_point() +
#'           geom_errorbar(aes(ymin = est - se, ymax = est + se)) +
#'           facet_wrap(~paste(pop1, pop2, pop3, pop4), scales = 'free')
#' }
f4_from_afdat = function(afdat, popcombs) {

  for(i in 1:nrow(popcombs)) {
    p1 = popcombs$pop1[i]
    p2 = popcombs$pop2[i]
    p3 = popcombs$pop3[i]
    p4 = popcombs$pop4[i]
    afdat %<>% mutate(!!paste0('f4_', i) := (!!sym(p1)-!!sym(p2))*(!!sym(p3)-!!sym(p4)))
  }
  gr = groups(afdat)
  afdat %>% select(!!!gr, any_of(c('CHR', 'cm', 'POS', 'block')), starts_with('f4_')) %>%
    snpdat_to_jackest %>% arrange(.col, !!!gr) %>% left_join(popcombs %>% mutate(.col = paste0('f4_', 1:n()))) %>%
    transmute(!!!gr, pop1, pop2, pop3, pop4, est, se = sqrt(var), z = est/se, cnt)
}

