

f3_from_f2 = function(f2_12, f2_13, f2_23) (f2_12 + f2_13 - f2_23) / 2
f4_from_f2 = function(f2_14, f2_23, f2_13, f2_24) (f2_14 + f2_23 - f2_13 - f2_24) / 2




#' Estimate f2 statistics
#'
#' Computes f2 statistics from f2 blocks of the form \eqn{f2(A, B)}
#' @export
#' @param f2_data a 3d array of blocked f2 statistics, output of \code{\link{f2_from_precomp}}.
#' alternatively, a directory with precomputed data. see \code{\link{extract_f2}} and \code{\link{extract_counts}}.
#' @param pop1 one of the following four:
#' \enumerate{
#' \item `NULL`: all possible pairs of the populations in `f2_blocks` will be returned
#' \item a vector of population labels
#' \item a data frame with population combinations to be tested, with one population per column and one
#' combination per row. Other `pop` arguments will be ignored.
#' \item the location of a file (`poplistname` or `popfilename`) which specifies the populations or
#' population combinations to be tested. Other `pop` arguments will be ignored.
#' }
#' @param pop2 a vector of population labels
#' @param f2_denom scales f2-statistics. 1 correspondes to `f4mode: YES`. 1/4.75 is similar to `f4mode: NO`.
#' @param boot If `FALSE` (the default), each block will be left out at a time and the covariance matrix
#' of the f statistics will be computed using block-jackknife. Otherwise bootstrap resampling is performed `n` times,
#' where `n` is either equal to `boot` if it is an integer, or equal to the number of blocks if `boot` is `TRUE`.
#' The the covariance matrix of the f statistics will be computed using bootstrapping.
#' @param sure The number of population combinations can get very large. This is a safety option that stops you
#' from accidently computing all combinations if that number is large.
#' @param unique_only If `TRUE` (the default), redundant combinations will be returned as well.
#' @param verbose print progress updates
#' @return a data frame with f2 statistics
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @examples
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' f2(example_f2_blocks, pop1, pop2)
#' \dontrun{
#' f2(f2_dir, pop1, pop2)
#' }
f2 = function(f2_data, pop1 = NULL, pop2 = NULL,
              f2_denom = 1, boot = FALSE, sure = FALSE, unique_only = TRUE, verbose = FALSE) {

  out = fstat_get_popcombs(f2_data = f2_data, pop1 = pop1, pop2 = pop2,
                           sure = sure, unique_only = unique_only, fnum = 2)
  pops = unique(c(out$pop1, out$pop2))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, boot_mat_stats, jack_mat_stats)
  f2_blocks = get_f2(f2_data, pops, f2_denom) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  #----------------- compute f2 -----------------
  if(verbose) alert_info('Computing f2-statistics\n')

  out %<>% group_by(pop1, pop2) %>%
    summarize(f2dat = list(f2_blocks[pop1, pop2, ])) %>% ungroup %>%
    mutate(sts = map(f2dat, ~statfun(t(.), block_lengths))) %>%
    unnest_wider(sts) %>%
    mutate(se = map_dbl(var[,1], sqrt), z = est/se, p = ztop(z)) %>%
    select(-f2dat, -var)

  out
}



#' Estimate f3 statistics
#'
#' Computes f3 statistics from f2 blocks of the form \eqn{f3(A; B, C)}. Equivalent to
#' \eqn{(f2(A, B) + f2(A, C) - f2(B, C)) / 2} and to \eqn{f4(A, B; A, C)}
#' @export
#' @param pop3 a vector of population labels
#' @inheritParams f2
#' @return a data frame with f3 statistics
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
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
qp3pop = function(f2_data, pop1 = NULL, pop2 = NULL, pop3 = NULL,
                  f2_denom = 1, boot = FALSE, sure = FALSE, unique_only = TRUE, verbose = FALSE) {

  stopifnot(is.null(pop2) & is.null(pop3) |
              !is.null(pop2) & !is.null(pop3))

  out = fstat_get_popcombs(f2_data = f2_data, pop1 = pop1, pop2 = pop2, pop3 = pop3,
                           sure = sure, unique_only = unique_only, fnum = 3)
  pops = unique(c(out$pop1, out$pop2, out$pop3))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo)
  statfun = ifelse(boot, boot_mat_stats, jack_mat_stats)
  f2_blocks = get_f2(f2_data, pops, f2_denom) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  #----------------- compute f3 -----------------
  if(verbose) alert_info('Computing f3-statistics\n')

  out %<>% group_by(pop1, pop2, pop3) %>%
    summarize(f3dat = list(f3_from_f2(f2_blocks[pop1, pop2, ],
                                     f2_blocks[pop1, pop3, ],
                                     f2_blocks[pop2, pop3, ]))) %>% ungroup %>%
    mutate(sts = map(f3dat, ~statfun(t(.), block_lengths))) %>%
    unnest_wider(sts) %>%
    mutate(se = map_dbl(var[,1], sqrt), z = est/se, p = ztop(z)) %>%
    select(-f3dat, -var)

  out
}
#' @export
f3 = qp3pop


#' Estimate f4 statistics
#'
#' Computes f4 statistics from f2 blocks of the form \eqn{f4(A, B; C, D)}. Equivalent to
#' \eqn{(f2(A, D) + f2(B, C) - f2(A, C) - f2(B, D)) / 2}
#' @export
#' @param f2_data f2 data in one of the following formats
#' \enumerate{
#' \item A 3d array of block-jackknife leave-one-block-out estimates of f2 statistics,
#' output of \code{\link{extract_f2}} and \code{\link{extract_counts}}
#' \item A directory with f2 statistics
#' \item Prefix of a packedancestrymap file. This is the slowest option, but allows to compute f4-statistics based on all non-missing SNPs in each population quadruple. This can be more precise in the presence of large amounts of missing data.
#' }
#' @param pop3 a vector of population labels
#' @param pop4 a vector of population labels
#' @param dist genetic distance in Morgan. Default is 0.05 (50 cM). only used when `f2_data` is the prefix of packedancestrymap files
#' @param block_lengths vector with lengths of each jackknife block. \code{sum(block_lengths)} has to
#' match the number of SNPs. only used when `f2_data` is the prefix of packedancestrymap files
#' @param f4mode if `TRUE`: f4 is computed from allele frequencies `a`, `b`, `c`, and `d` as `(a-b)*(c-d)`. if `FALSE`, D-statistics are computed instead, defined as `(a-b)*(c-d) / ((a + b - 2*a*b) * (c + d - 2*c*d))`. `f4mode = FALSE` is only available when `f2_data` is the prefix of packedancestrymap files
#' @inheritParams f2
#' @return a data frame with f4 statistics
#' @aliases f4
#' @section Alias:
#' `f4`
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @examples
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' pop3 = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' pop4 = 'Switzerland_Bichon.SG'
#' qpdstat(example_f2_blocks, pop1, pop2, pop3, pop4)
#' \dontrun{
#' qpdstat(f2_dir, pop1, pop2, pop3, pop4)
#' }
qpdstat = function(f2_data, pop1 = NULL, pop2 = NULL, pop3 = NULL, pop4 = NULL,
                   f2_denom = 1, boot = FALSE, sure = FALSE, unique_only = TRUE,
                   comb = TRUE, dist = NULL, block_lengths = NULL, f4mode = TRUE, verbose = FALSE) {

  stopifnot(is.null(pop2) & is.null(pop3) & is.null(pop4) |
            !is.null(pop2) & !is.null(pop3) & !is.null(pop4))

  if(is_packedancestrymap_prefix(f2_data)) {
    if(verbose) alert_info('Computing from f4 from genotype data...\n')
    return(f4_from_geno(f2_data, pop1, pop2, pop3, pop4, dist = ifelse(is.null(dist), 0.05, dist),
                        f4mode = f4mode, block_lengths = block_lengths, verbose = TRUE))
  }

  if(!comb) {
    stopifnot(!is.null(pop2))
    stopifnot(length(unique(length(pop1), length(pop2), length(pop3), length(pop4))) == 1)
    out = tibble(pop1, pop2, pop3, pop4)
  } else {
    out = fstat_get_popcombs(f2_data = f2_data, pop1 = pop1, pop2 = pop2, pop3 = pop3, pop4 = pop4,
                             sure = sure, unique_only = unique_only, fnum = 4)
  }
  pops = unique(c(out$pop1, out$pop2, out$pop3, out$pop4))
  pops1 = unique(c(out$pop1, out$pop2))
  pops2 = unique(c(out$pop3, out$pop4))

  samplefun = ifelse(boot, function(x) est_to_boo(x, boot), est_to_loo_nafix)
  statfun = ifelse(boot, boot_mat_stats, jack_mat_stats)
  #f2_blocks = get_f2(f2_data, pops, f2_denom) %>% samplefun
  f2_blocks = get_f2(f2_data, pops1, f2_denom, pops2 = pops2) %>% samplefun
  block_lengths = parse_number(dimnames(f2_blocks)[[3]])

  #----------------- compute f4 -----------------
  if(verbose) alert_info('Computing f4-statistics\n')

  out %<>% rowwise %>%
    mutate(f4dat = list(f4_from_f2(f2_blocks[pop1, pop4, ],
                                   f2_blocks[pop2, pop3, ],
                                   f2_blocks[pop1, pop3, ],
                                   f2_blocks[pop2, pop4, ]))) %>% ungroup %>%
    mutate(sts = map(f4dat, ~statfun(t(.), block_lengths))) %>%
    unnest_wider(sts) %>%
    mutate(se = map_dbl(var[,1], sqrt), z = est/se, p = ztop(z)) %>%
    select(-f4dat, -var)
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
  maxcomb = 1e5
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
  } else if(!'data.frame' %in% class(pop1) && file.exists(pop1)) {
    pop1 = read_table2(pop1, col_names = FALSE)
    if(ncol(pop1) == 1) {
      pop1 = pop1[[1]]
    } else {
      out = pop1 %>% set_colnames(nam)
    }
  } else if('data.frame' %in% class(pop1)) {
    out = pop1 %>% set_colnames(nam)
  }

  if(is.null(out)) {
    # pop1 is a character vector at this point
    if(unique_only) {
      ncomb = choose(length(pop1), fnum)*(fnum-1)
      if(ncomb > maxcomb & !sure) {
        stop(paste0('If you really want to compute ', ncomb,
                    ' f-statistics, run this again with "sure = TRUE". Or specify more than just pop1.'))
      }
      if(length(pop1) < fnum) stop('Not enough populations!')
      outmat = t(combn(pop1, fnum))
      nr = nrow(outmat)
      if(fnum == 2) {
        out = outmat %>% as_tibble %>% set_colnames(nam)
      } else if(fnum == 3) {
        out = rbind(outmat[,1:3], outmat[,c(2,3,1)], outmat[,c(3,1,2)]) %>%
          as_tibble %>% set_colnames(nam)
      } else if(fnum == 4) out = rbind(outmat[,1:4], outmat[,c(1,3,2,4)], outmat[,c(1,4,2,3)]) %>%
        as_tibble() %>% set_colnames(nam) %>% slice(rep(1:nr, each=3) + (0:2)*nr)
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
      else if(fnum == 4) out = as_tibble(t(combn(pop1, 2))) %>%
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

  a = f4(f2_blocks) %>% select(pop1:pop4, est, se, z) %>% mutate(est = ifelse((pop1 > pop2) == (pop3 > pop4), est, -est), f13 = paste(pmin(pop1, pop3), pmax(pop1, pop3)), f24 = paste(pmin(pop2, pop4), pmax(pop2, pop4)), f14 = paste(pmin(pop1, pop4), pmax(pop1, pop4)), f23 = paste(pmin(pop2, pop3), pmax(pop2, pop3))) %>% expand_grid(pp) %>% mutate(coef = ifelse(pp == f14 | pp == f23, 1, ifelse(pp  == f13 | pp == f24, -1, 0))) %>% pivot_wider(pop1:z, names_from = pp, values_from = coef) %>% mutate(est = -est*2) %>% select(-pop1:-pop4) %>% lm(as.formula(paste0('est ~ 0 + ', paste('`', colnames(.)[-1:-3], '`', sep = '', collapse = ' + '))), weights = 1/se, data = .)


  summary(a)$coefficients %>% as_tibble(rownames = 'pp') %>% mutate(pp = str_replace_all(pp, '`', '')) %>% separate(pp, c('pop1', 'pop2'), sep = ' ') %>% left_join(f2(f2_blocks[sort(pops), sort(pops),]), by = c('pop1', 'pop2')) %>% ggplot(aes(Estimate, est)) + geom_point() + geom_abline()


  ff2 = f2(f2_blocks[sort(pops), sort(pops),]) %>% bind_rows(rename(., pop1=pop2, pop2=pop1)) %>% rename(p1 = pop1, p2 = pop2)

}




