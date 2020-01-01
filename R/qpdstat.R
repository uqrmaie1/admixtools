


# den = (w + x - 2*w*x) * (y + z -2*y*z)
#
# w = y
# x = z
#
# den = (y + z - 2*y*z) * (y + z -2*y*z)
# den = (y + z - 2*y*z)^2
# den = y^2 + z^2 - 4*y^2*z^2 + 2*y*z - 4*y^2*z - 4*y*z^2
# den = (y + z)^2 - 4*y*z * (y*z + y + z)


f3_from_f2 = function(f2_12, f2_13, f2_23) (f2_12 + f2_13 - f2_23) / 2
f4_from_f2 = function(f2_14, f2_23, f2_13, f2_24) (f2_14 + f2_23 - f2_13 - f2_24) / 2


#' Wrapper function around the original qp3Pop program
#'
#' Computes f3 statistics of the form \eqn{F_3(P_1; P_2, P_3)}. Equivalent to \eqn{(F_2(P_1, P_2) + F_2(P_1, P_3) - F_2(P_2, P_3)) / 2} and to \eqn{F_4(P_1, P_2; P_1, P_3)}. Requires a working installation of qp3Pop, which will be called using \code{\link{system}}
#' @export
#' @param source1 one of the following four:
#' \enumerate{
#' \item \code{NULL}: populations will be read from \code{poplistname} or \code{popfilename} specified in \code{parfile}
#' \item a vector of population labels
#' \item a data frame in which the first four columns specify the population quadruples to be tested. \code{source2}, \code{target} will be ignored.
#' \item the location of a file (\code{poplistname} or \code{popfilename}) which specifies the populations or population quadruples to be tested. \code{source2} and \code{target} will be ignored.
#' }
#' @param source2 a vector of population labels
#' @param target a vector of population labels
#' @param bin path to the qp3Pop binary file
#' @param pref path to and prefix of the packedancestrymap genotype files
#' @param outdir the output directory. files \code{out}, \code{parfile}, \code{poplistname}, \code{popfilename} may be overwritten
#' @param parfile qp3Pop parameter file. If this is specified, \code{source1}, \code{source2}, \code{target} will be ignored.
#' @param printonly should the command be printed or executed?
#' @return If \code{printonly}, the \code{qp3Pop} command, otherwise a data frame with parsed \code{qp3Pop} output
#' @examples
#' \dontrun{
#' target = 'Denisova.DG'
#' source1 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' source2 = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' qp3pop_wrapper(source1, source2, target,
#'   bin = 'path/to/qp3Pop', pref = 'path/to/packedancestrymap_prefix',
#'   env = 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/blas/')
#'
#' qp3pop_wrapper(bin = 'path/to/qp3Pop', parfile = 'path/to/parfile')
#' }
qp3pop_wrapper = function(source1, source2 = NULL, target = NULL, bin = './qp3Pop',
                          pref = NULL, outdir='.', parfile = NULL,
                          inbreed = 'NO', outgroupmode = 'YES',
                          printonly=FALSE, env='') {

  stopifnot(!is.null(parfile) | !is.null(pref) & !is.null(source1))
  stopifnot(file.exists(str_replace(bin, '.+ ', '')))
  stopifnot(is.null(source1) | is.null(source2) & is.null(target) |
              !is.null(source2) & !is.null(target))

  if(is.null(parfile)) {

    popfilename = paste0(outdir, '/popfilename')
    if(!is.null(source2)) {
      expand_grid(source1, source2, target) %>% write_tsv(popfilename, col_names = FALSE)
    } else if('data.frame' %in% class(source1)) {
      source1 %>% select(1:3) %>% write_tsv(popfilename, col_names = FALSE)
    } else {
      stopifnot(file.exists(source1))
      popfilename = source1
    }
    pref = normalizePath(pref, mustWork = FALSE)
    parfile = paste0('genotypename: ', pref, '.geno\n',
                     'snpname: ', pref, '.snp\n',
                     'indivname: ', pref, '.ind\n',
                     'popfilename: ', popfilename, '\n',
                   # 'f4mode: ', f4mode, '\n',
                     'inbreed: ', inbreed, '\n',
                     'outgroupmode: ', outgroupmode, '\n',
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
    return(parse_qp3pop_output(paste0(outdir, '/out')))
  }
}



#' Wrapper function around the original qpDstat program
#'
#' This requires a working installation of qpDstat, which will be called using \code{\link{system}}
#' @export
#' @param pop1 one of the following four:
#' \enumerate{
#' \item \code{NULL}: populations will be read from \code{poplistname} or \code{popfilename} specified in \code{parfile}
#' \item a vector of population labels
#' \item a data frame in which the first four columns specify the population quadruples to be tested. \code{pop2}, \code{pop3}, \code{pop4} will be ignored.
#' \item the location of a file (\code{poplistname} or \code{popfilename}) which specifies the populations or population quadruples to be tested. \code{pop2}, \code{pop3}, \code{pop4} will be ignored.
#' }
#' @param pop2 a vector of population labels
#' @param pop3 a vector of population labels
#' @param pop4 a vector of population labels
#' @param bin path to the qpDstat binary file
#' @param pref path to and prefix of the packedancestrymap genotype files
#' @param outdir the output directory. files \code{out}, \code{parfile}, \code{poplistname}, \code{popfilename} may be overwritten
#' @param parfile qpDstat parameter file. If this is specified, \code{pop}, \code{pop2}, \code{pop3}, and \code{pop4} will be ignored.
#' @param printonly should the command be printed or executed?
#' @return If \code{printonly}, the \code{qpDstat} command, otherwise a data frame with parsed \code{qpDstat} output
#' @examples
#' \dontrun{
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' pop3 = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' pop4 = 'Switzerland_Bichon.SG'
#' qpdstat_wrapper(pop1, pop2, pop3, pop4,
#'   bin = 'path/to/qpDstat', pref = 'path/to/packedancestrymap_prefix',
#'   env = 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/blas/')
#'
#' qpdstat_wrapper(bin = 'path/to/qpDstat', parfile = 'path/to/parfile')
#' }
qpdstat_wrapper = function(pop1 = NULL, pop2 = NULL, pop3 = NULL, pop4 = NULL, bin = './qpDstat',
                           pref = NULL, outdir='.', parfile = NULL, f4mode = 'YES', printonly=FALSE, env='') {

  stopifnot(!is.null(parfile) | !is.null(pref) & !is.null(pop1))
  stopifnot(file.exists(str_replace(bin, '.+ ', '')))
  stopifnot(is.null(pop2) & is.null(pop3) & is.null(pop4) |
            !is.null(pop2) & !is.null(pop3) & !is.null(pop4))

  if(is.null(parfile)) {

    popfiletype = 'popfilename'
    popfilename = paste0(outdir, '/popfilename')
    if(!is.null(pop2)) {
      expand_grid(pop1, pop2, pop3, pop4) %>% write_tsv(popfilename, col_names = FALSE)
    } else if('data.frame' %in% class(pop1)) {
      pop1 %>% select(1:4) %>% write_tsv(popfilename, col_names = FALSE)
    } else if(length(pop1) == 1) {
      stopifnot(file.exists(pop1))
      dat = read_table2(pop1, col_names = FALSE)
      if(ncol(dat) < 4) popfiletype = 'poplistname'
      popfilename = pop1
    } else {
      popfiletype = 'poplistname'
      popfilename = paste0(outdir, '/poplistname')
      write(pop1, popfilename)
    }
    pref = normalizePath(pref, mustWork = FALSE)
    parfile = paste0('genotypename: ', pref, '.geno\n',
                     'snpname: ', pref, '.snp\n',
                     'indivname: ', pref, '.ind\n',
                     popfiletype, ': ', popfilename, '\n',
                     'f4mode: ', f4mode, '\n',
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
    return(parse_qpdstat_output(paste0(outdir, '/out')))
  }
}


#' Estimate f2 statistics
#'
#' Computes f2 statistics from f2 jackknife blocks of the form \eqn{F_2(P_1, P_2)}
#' @export
#' @param pop1 one of the following four:
#' \enumerate{
#' \item \code{NULL}: all possible pairs of the populations in \code{f2_blocks} will be returned
#' \item a vector of population labels
#' \item a data frame in which the first four columns specify the population quadruples to be tested. \code{pop2} will be ignored.
#' \item the location of a file (\code{poplistname} or \code{popfilename}) which specifies the populations or population quadruples to be tested. \code{pop2} will be ignored.
#' }
#' @param pop2 a vector of population labels
#' @param pop3 a vector of population labels
#' @param f2_blocks 3d array of block-jackknife leave-one-block-out estimates of f2 statistics. output of \code{\link{afs_to_f2_blocks}}. they are weighted by inverse of outgroup heterozygosity, if outgroup was specified.
#' @param block_lengths the jackknife block lengths used in computing the f2 statistics. see \code{\link{get_block_lengths}}.
#' @param f2_dir a directory with f2 statistics for each population pair in the graph. must contain 'block_lengths.RData'.
#' @param f2_denom scales f2-statistics. 1 correspondes to \code{f4mode: YES}. 1/4.75 is similar to \code{f4mode: NO}.
#' @return If \code{printonly}, the \code{qpDstat} command, otherwise a data frame with parsed \code{qpDstat} output
#' @examples
#' \dontrun{
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' f2(pop1, pop2, f2_blocks = f2_blocks, block_lengths = block_lengths)
#' f2(pop1, pop2, f2_dir = f2_dir)
#' }
f2 = function(pop1 = NULL, pop2 = NULL,
              f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL, f2_denom = 1,
              sure = FALSE, unique_only = TRUE, verbose = FALSE) {

  out = fstat_get_popcombs(pop1 = pop1, pop2 = pop2,
                           f2_blocks = f2_blocks, sure = sure, unique_only = unique_only,
                           fnum = 2)
  pops = unique(c(out$pop1, out$pop2))
  f2s = fstat_get_f2(pops = pops, f2_blocks = f2_blocks, block_lengths = block_lengths,
                     f2_dir = f2_dir, f2_denom = f2_denom)
  f2_blocks_scaled = f2s[[1]]
  block_lengths = f2s[[2]]
  f2nam = dimnames(f2_blocks_scaled)[[1]]

  #----------------- compute f2 -----------------
  if(verbose) alert_info('Computing f2-statistics\n')

  out %<>% group_by(pop1, pop2) %>%
    summarize(jack = list(f2_blocks_scaled[pop1, pop2, ])) %>% ungroup %>%
    mutate(sts = map(jack, ~bj_mat_stats(t(.), block_lengths))) %>%
    unnest_wider(sts) %>%
    mutate(estimate = jest, se = map_dbl(jvar[,1], sqrt), z = estimate/se, p.value = ztop(z)) %>%
    select(-jack, -jest, -jvar)

  out

}



#' Estimate f3 statistics
#'
#' Computes f3 statistics from f2 jackknife blocks of the form \eqn{F_3(P_1; P_2, P_3)}. Equivalent to  \eqn{(F_2(P_1, P_2) + F_2(P_1, P_3) - F_2(P_2, P_3)) / 2} and to \eqn{F_4(P_1, P_2; P_1, P_3)}
#' @export
#' @param pop1 one of the following four:
#' \enumerate{
#' \item \code{NULL}: all possible triples of the populations in \code{f2_blocks} will be returned
#' \item a vector of population labels
#' \item a data frame in which the first four columns specify the population quadruples to be tested. \code{pop2} and \code{pop3} will be ignored.
#' \item the location of a file (\code{poplistname} or \code{popfilename}) which specifies the populations or population quadruples to be tested. \code{pop2} and \code{pop3} will be ignored.
#' }
#' @param pop2 a vector of population labels
#' @param pop3 a vector of population labels
#' @param f2_blocks 3d array of block-jackknife leave-one-block-out estimates of f2 statistics. output of \code{\link{afs_to_f2_blocks}}. they are weighted by inverse of outgroup heterozygosity, if outgroup was specified.
#' @param block_lengths the jackknife block lengths used in computing the f2 statistics. see \code{\link{get_block_lengths}}.
#' @param f2_dir a directory with f2 statistics for each population pair in the graph. must contain 'block_lengths.RData'.
#' @param f2_denom scales f2-statistics. 1 correspondes to \code{f4mode: YES}. 1/4.75 is similar to \code{f4mode: NO}
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @return If \code{printonly}, the \code{qp3Pop} command, otherwise a data frame with parsed \code{qp3Pop} output
#' @aliases f3
#' @section Alias:
#' \code{f3}
#' @examples
#' \dontrun{
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' pop3 = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' qp3pop(pop1, pop2, pop3,
#'   f2_blocks = f2_blocks, block_lengths = block_lengths)
#' qp3pop(pop1, pop2, pop3, f2_dir = f2_dir)
#' }
qp3pop = function(pop1 = NULL, pop2 = NULL, pop3 = NULL,
              f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL, f2_denom = 1,
              sure = FALSE, unique_only = TRUE, verbose = FALSE) {

  stopifnot(is.null(pop2) & is.null(pop3) |
              !is.null(pop2) & !is.null(pop3))

  out = fstat_get_popcombs(pop1 = pop1, pop2 = pop2, pop3 = pop3,
                           f2_blocks = f2_blocks, sure = sure, unique_only = unique_only,
                           fnum = 3)
  pops = unique(c(out$pop1, out$pop2, out$pop3))
  f2s = fstat_get_f2(pops = pops, f2_blocks = f2_blocks, block_lengths = block_lengths,
                     f2_dir = f2_dir, f2_denom = f2_denom)
  f2_blocks_scaled = f2s[[1]]
  block_lengths = f2s[[2]]
  f2nam = dimnames(f2_blocks_scaled)[[1]]

  #----------------- compute f3 -----------------
  if(verbose) alert_info('Computing f3-statistics\n')

  out %<>% group_by(pop1, pop2, pop3) %>%
    summarize(jack = list(f3_from_f2(f2_blocks_scaled[pop1, pop2, ],
                                     f2_blocks_scaled[pop1, pop3, ],
                                     f2_blocks_scaled[pop2, pop3, ]))) %>% ungroup %>%
    mutate(sts = map(jack, ~bj_mat_stats(t(.), block_lengths))) %>%
    unnest_wider(sts) %>%
    mutate(estimate = jest, se = map_dbl(jvar[,1], sqrt), z = estimate/se, p.value = ztop(z)) %>%
    select(-jack, -jest, -jvar)

  out

}
#' @export
f3 = qp3pop


#' Estimate f4 statistics
#'
#' Computes f4 statistics from f2 jackknife blocks of the form \eqn{F_4(P_1, P_2; P_3, P_4)}. Equivalent to \eqn{(F_2(P_1, P_4) + F_2(P_2, P_3) - F_2(P_1, P_3) - F_2(P_2, P_4)) / 2}
#' @export
#' @param pop1 one of the following four:
#' \enumerate{
#' \item \code{NULL}: all possible quadruples of the populations in \code{f2_blocks} will be returned
#' \item a vector of population labels
#' \item a data frame in which the first four columns specify the population quadruples to be tested. \code{pop2}, \code{pop3}, \code{pop4} will be ignored.
#' \item the location of a file (\code{poplistname} or \code{popfilename}) which specifies the populations or population quadruples to be tested. \code{pop2}, \code{pop3}, \code{pop4} will be ignored.
#' }
#' @param pop2 a vector of population labels
#' @param pop3 a vector of population labels
#' @param pop4 a vector of population labels
#' @param f2_blocks 3d array of block-jackknife leave-one-block-out estimates of f2 statistics. output of \code{\link{afs_to_f2_blocks}}. they are weighted by inverse of outgroup heterozygosity, if outgroup was specified.
#' @param block_lengths the jackknife block lengths used in computing the f2 statistics. see \code{\link{get_block_lengths}}.
#' @param f2_dir a directory with f2 statistics for each population pair in the graph. must contain 'block_lengths.RData'.
#' @param f2_denom scales f2-statistics. 1 correspondes to \code{f4mode: YES}. 1/4.75 is similar to \code{f4mode: NO}
#' @return If \code{printonly}, the \code{qpDstat} command, otherwise a data frame with parsed \code{qpDstat} output
#' @aliases f4
#' @section Alias:
#' \code{f4}
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#' @examples
#' \dontrun{
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' pop3 = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' pop4 = 'Switzerland_Bichon.SG'
#' qpdstat(pop1, pop2, pop3, pop4,
#'   f2_blocks = f2_blocks, block_lengths = block_lengths)
#' qpdstat(pop1, pop2, pop3, pop4, f2_dir = f2_dir)
#' }
qpdstat = function(pop1 = NULL, pop2 = NULL, pop3 = NULL, pop4 = NULL,
                   f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL, f2_denom = 1,
                   sure = FALSE, unique_only = TRUE, verbose = FALSE) {

  stopifnot(is.null(pop2) & is.null(pop3) & is.null(pop4) |
            !is.null(pop2) & !is.null(pop3) & !is.null(pop4))

  out = fstat_get_popcombs(pop1 = pop1, pop2 = pop2, pop3 = pop3, pop4 = pop4,
                           f2_blocks = f2_blocks, sure = sure, unique_only = unique_only,
                           fnum = 4)
  pops = unique(c(out$pop1, out$pop2, out$pop3, out$pop4))
  f2s = fstat_get_f2(pops = pops, f2_blocks = f2_blocks, block_lengths = block_lengths,
                     f2_dir = f2_dir, f2_denom = f2_denom)
  f2_blocks_scaled = f2s[[1]]
  block_lengths = f2s[[2]]
  f2nam = dimnames(f2_blocks_scaled)[[1]]

  #----------------- compute f4 -----------------
  if(verbose) alert_info('Computing f4-statistics\n')

  out %<>% group_by(pop1, pop2, pop3, pop4) %>%
    summarize(jack = list(f4_from_f2(f2_blocks_scaled[pop1, pop4, ],
                                     f2_blocks_scaled[pop2, pop3, ],
                                     f2_blocks_scaled[pop1, pop3, ],
                                     f2_blocks_scaled[pop2, pop4, ]))) %>% ungroup %>%
    mutate(sts = map(jack, ~bj_mat_stats(t(.), block_lengths))) %>%
    unnest_wider(sts) %>%
    mutate(estimate = jest, se = map_dbl(jvar[,1], sqrt), z = estimate/se, p.value = ztop(z)) %>%
    select(-jack, -jest, -jvar)

  #----- get output in same order as qpDstat ----
  # out %>% mutate_at(vars(starts_with('pop')), ~(ordered(., levels = f2nam))) %>%
  #   group_by(pop1, pop2, pop3, pop4) %>%
  #   mutate(grp = paste(sort(c(pop2, pop3, pop4)), collapse = ' ')) %>%
  #   group_by(grp) %>% arrange(pop1, grp, pop2, .by_group=FALSE) %>% ungroup %>%
  #   mutate_at(vars(starts_with('pop')), ~(as.character(.))) %>% select(-grp)
  out
}
#' @export
f4 = qpdstat


fstat_get_popcombs = function(pop1 = NULL, pop2 = NULL, pop3 = NULL, pop4 = NULL,
              f2_blocks = NULL, sure = FALSE, unique_only = TRUE, fnum = NULL) {
  # used by f2, f3, and f4 function
  # returns data frame 'out' with pop combs
  stopifnot(!is.null(pop1) | !is.null(f2_blocks))

  #----------------- make combinations -----------------
  out = NULL
  nam = c('pop1', 'pop2', 'pop3', 'pop4')[1:fnum]
  maxcomb = 1e5
  if(!is.null(pop2)) {
    ncomb = length(pop1) * length(pop2) * max(1, length(pop3)) * max(1, length(pop4))
    if(ncomb > maxcomb & !sure) stop(paste0('If you really want to compute ', ncomb, ' f-statistics, run this again with "sure = TRUE".'))
    out = expand_grid(pop1, pop2, pop3, pop4)
  } else if(is.null(pop1)) {
    pop1 = dimnames(f2_blocks)[[1]]
  } else if(!'data.frame' %in% class(pop1)) {
    pop1 = read_table2(pop1, col_names = FALSE)
    if(ncol(pop1) == 1) {
      pop1 = pop1[[1]]
    } else {
      out = pop1 %>% set_colnames(nam)
    }
  } else {
    out = pop1 %>% set_colnames(nam)
  }

  if(is.null(out)) {
    # pop1 is a character vector at this point
    if(unique_only) {
      ncomb = choose(length(pop1), fnum)*(fnum-1)
      if(ncomb > maxcomb & !sure) stop(paste0('If you really want to compute ', ncomb, ' f-statistics, run this again with "sure = TRUE". Or specify more than just pop1.'))
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
      if(fnum^length(pop1) > maxcomb & !sure) stop(paste0('If you really want to compute close to ', fnum^length(pop1), ' f-statistics, run this again with "sure = TRUE". Or specify more than just pop1.'))
      if(fnum == 2) out = expand_grid(pop1 = pop1, pop2 = pop1)
      else if(fnum == 3) out = expand_grid(pop1 = pop1, pop2 = pop1, pop3 = pop1)
      else if(fnum == 4) out = expand_grid(pop1 = pop1, pop2 = pop1, pop3 = pop1, pop4 = pop1) %>% filter(pop1 != pop3, pop2 != pop4, pop1 != pop2 | pop3 != pop4)
      else stop('fnum should be 2, 3, or 4!')
    }
  }
  out
}


fstat_get_f2 = function(pops, f2_blocks = NULL, block_lengths = NULL, f2_dir = NULL, f2_denom = 1) {
  # used by f2, f3, and f4 function
  # returns list with 'f2_blocks_scaled' and 'block_lengths'

  #----------------- read f-stats -----------------
  if(is.null(f2_dir) & (is.null(f2_blocks) | is.null(block_lengths))) stop('You have to provide an f2_dir argument, or f2_blocks and block_lengths!')
  if(!is.null(f2_dir) & (is.null(f2_blocks) | is.null(block_lengths))) {
    f2_blocks = read_f2(f2_dir, pops = pops)
    load(paste0(f2_dir, '/block_lengths.RData'))
  }

  #----------------- process f-stats -----------------
  f2nam = dimnames(f2_blocks)[[1]]
  stopifnot(all(pops %in% f2nam))
  f2_blocks_scaled = f2_blocks / f2_denom

  namedList(f2_blocks_scaled, block_lengths)

}

