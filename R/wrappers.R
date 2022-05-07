yesno = function(arg) {
  if(isTRUE(arg) || arg == 'YES') return('YES')
  else if(isFALSE(arg) || arg == 'NO') return('NO')
  stop("Argument has to be 'YES', 'NO', TRUE, or FALSE!")
}

#' Wrapper function around the original qp3Pop program
#'
#' Computes f3 statistics of the form \eqn{f3(A; B, C)}. Equivalent to \eqn{(f2(A, B) + f2(A, C) - f2(B, C)) / 2}
#' and to \eqn{f4(A, B; A, C)}. Requires a working installation of qp3Pop, which will be called
#' using \code{\link{system}}
#' @export
#' @param pref Path to and prefix of the packedancestrymap genotype files
#' @param source1 One of the following four:
#' \enumerate{
#' \item \code{NULL}: Populations will be read from \code{poplistname} or \code{popfilename} specified in \code{parfile}
#' \item A vector of population labels
#' \item A data frame in which the first four columns specify the population triples to be tested.
#' \code{source2}, \code{target} will be ignored.
#' \item The location of a file (\code{poplistname} or \code{popfilename}) which specifies the populations or
#' population combinations to be tested. \code{source2} and \code{target} will be ignored.
#' }
#' @param source2 A vector of population labels
#' @param target A vector of population labels
#' @param bin Path to the qp3Pop binary file
#' @param outdir Output directory. files \code{out}, \code{parfile}, \code{poplistname},
#' \code{popfilename} may be overwritten
#' @param parfile qp3Pop parameter file. If this is specified, \code{source1}, \code{source2},
#' \code{target} will be ignored.
#' @param inbreed inbreed
#' @param outgroupmode outgroupmode
#' @param f4mode f4mode
#' @param printonly Should the command be printed or executed?
#' @param env Export environmental variables. See examples.
#' @param verbose Print progress updates
#' @return If \code{printonly}, the \code{qp3Pop} command, otherwise a data frame with parsed \code{qp3Pop} output
#' @examples
#' \dontrun{
#' source1 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' source2 = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' target = 'Denisova.DG'
#' qp3pop_wrapper('genotype_prefix', source1, source2, target,
#'   bin = 'path/to/qp3Pop')
#' }
qp3pop_wrapper = function(pref, source1, source2 = NULL, target = NULL, bin = '~np29/o2bin/qp3Pop',
                          outdir = '.', parfile = NULL,
                          inbreed = 'NO', outgroupmode = 'YES', f4mode = 'YES',
                          printonly = FALSE, env = '', verbose = TRUE) {

  stopifnot(!is.null(parfile) | !is.null(pref) & !is.null(source1))
  stopifnot(file.exists(str_replace(bin, '.+ ', '')))
  stopifnot(is.null(source1) | is.null(source2) & is.null(target) |
              !is.null(source2) & !is.null(target))

  outdir = normalizePath(outdir, mustWork = FALSE)
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
                     'f4mode: ', f4mode %>% yesno, '\n',
                     'inbreed: ', inbreed %>% yesno, '\n',
                     'outgroupmode: ', outgroupmode %>% yesno, '\n',
                     'details: YES\n')
    parfilename = paste0(outdir, '/parfile')
    write(parfile, parfilename)

  } else {
    parfilename = parfile
  }

  cmd = paste0(env,' ', bin, ' -p ', parfilename, ' > ', outdir, '/out')
  outfile = paste0(outdir, '/out')
  run_admixtools(cmd, parse_qp3pop_output, outfile, printonly, verbose)
}



#' Wrapper function around the original qpDstat program
#'
#' This requires a working installation of qpDstat, which will be called using \code{\link{system}}
#' @export
#' @param pref Path to and prefix of the packedancestrymap genotype files
#' @param pop1 One of the following four:
#' \enumerate{
#' \item \code{NULL}: populations will be read from \code{poplistname} or \code{popfilename} specified in \code{parfile}
#' \item A vector of population labels
#' \item A data frame in which the first four columns specify the population combinations to be tested.
#' \code{pop2}, \code{pop3}, \code{pop4} will be ignored.
#' \item the location of a file (\code{poplistname} or \code{popfilename}) which specifies
#' the populations or population combinations to be tested. \code{pop2}, \code{pop3}, \code{pop4} will be ignored.
#' }
#' @param pop2 A vector of population labels
#' @param pop3 A vector of population labels
#' @param pop4 A vector of population labels
#' @param bin Path to the qpDstat binary file
#' @param outdir Output directory. files \code{out}, \code{parfile}, \code{poplistname},
#' \code{popfilename} may be overwritten
#' @param parfile qpDstat parameter file. If this is specified, \code{pop}, \code{pop2}, \code{pop3},
#' and \code{pop4} will be ignored.
#' @param f4mode f4mode
#' @param inbreed inbreed
#' @param printonly Should the command be printed or executed?
#' @param env Export environmental variables. See examples.
#' @param verbose Print progress updates
#' @return If \code{printonly}, the \code{qpDstat} command, otherwise a data frame with parsed \code{qpDstat} output
#' @examples
#' \dontrun{
#' pop1 = 'Denisova.DG'
#' pop2 = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' pop3 = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG')
#' pop4 = 'Switzerland_Bichon.SG'
#' qpdstat_wrapper('genotype_prefix', pop1, pop2, pop3, pop4,
#'   bin = 'path/to/qpDstat', pref = 'path/to/packedancestrymap_prefix')
#' }
qpdstat_wrapper = function(pref, pop1, pop2 = NULL, pop3 = NULL, pop4 = NULL,
                           bin = '~np29/o2bin/qpDstat', outdir='.', parfile = NULL,
                           f4mode = 'YES', inbreed = 'NO',
                           printonly=FALSE, env='', verbose = TRUE) {

  stopifnot(!is.null(parfile) | !is.null(pref) & !is.null(pop1))
  stopifnot(file.exists(str_replace(bin, '.+ ', '')))
  stopifnot(is.null(pop2) & is.null(pop3) & is.null(pop4) |
              !is.null(pop2) & !is.null(pop3) & !is.null(pop4))

  outdir = normalizePath(outdir, mustWork = FALSE)
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
                     'f4mode: ', f4mode %>% yesno, '\n',
                     'inbreed: ', inbreed %>% yesno, '\n',
                     'details: YES\n')
    parfilename = paste0(outdir, '/parfile')
    write(parfile, parfilename)

  } else {
    parfilename = parfile
  }

  cmd = paste0(env,' ', bin, ' -p ', parfilename, ' > ', outdir, '/out')
  outfile = paste0(outdir, '/out')
  run_admixtools(cmd, parse_qpdstat_output, outfile, printonly, verbose)
}


#' Wrapper function around the original qpF4ratio program
#'
#' This requires a working installation of qpF4ratio, which will be called using \code{\link{system}}
#' @export
#' @param pref Path to and prefix of the packedancestrymap genotype files
#' @param pops A vector of five populations, or a 5 x n matrix with population names. For each line `alpha` will be computed as `f4(1,2; 3,4)/f4(1,2; 5,4)`
#' @param bin Path to the qpF4ratio binary file
#' @param outdir Output directory. files \code{out}, \code{parfile}, \code{poplistname},
#' \code{popfilename} may be overwritten
#' @param parfile qpF4ratio parameter file. If this is specified, `pops` will be ignored.
#' @param blgsize blgsize
#' @param fancyf4 fancyf4
#' @param printonly Should the command be printed or executed?
#' @param env Export environmental variables. See examples.
#' @param verbose Print progress updates
#' @return If \code{printonly}, the \code{qpF4ratio} command, otherwise a data frame with parsed \code{qpF4ratio} output
#' @examples
#' \dontrun{
#' pops = c('Denisova.DG', 'Altai_Neanderthal.DG', 'Vindija.DG', 'Chimp.REF', 'Mbuti.DG')
#' qpf4ratio_wrapper('genotype_prefix', pops, bin = 'path/to/qpDstat')
#' }
qpf4ratio_wrapper = function(pref, pops, bin = '~np29/o2bin/qpF4ratio', outdir='.', parfile = NULL,
                             blgsize = 0.05, fancyf4 = 'YES', printonly = FALSE, env='', verbose = TRUE) {

  stopifnot(!is.null(parfile) || !is.null(pref))
  stopifnot(file.exists(str_replace(bin, '.+ ', '')))
  outdir = normalizePath(outdir, mustWork = FALSE)

  if(is.null(parfile)) {

    if(!is.matrix(pops)) pops %<>% t

    popfilename = paste0(outdir, '/popfilename')
    write_delim(data.frame(pops[,1:4,drop=F], '::', pops[,c(1,2,5,4),drop=F]), popfilename, delim = ' ', col_names = FALSE)

    pref = normalizePath(pref, mustWork = FALSE)

    parfile = paste0('genotypename: ', pref, '.geno\n',
                     'snpname: ', pref, '.snp\n',
                     'indivname: ', pref, '.ind\n',
                     'blgsize: ', blgsize, '\n',
                     'fancyf4: ', fancyf4 %>% yesno, '\n',
                     'popfilename: ', popfilename, '\n')
    parfilename = paste0(outdir, '/parfile')
    write(parfile, parfilename)

  } else {
    parfilename = parfile
  }

  cmd = paste0(env,' ', bin, ' -p ', parfilename, ' > ', outdir, '/out')
  outfile = paste0(outdir, '/out')
  run_admixtools(cmd, parse_qpf4ratio_output, outfile, printonly, verbose)

}


#' Wrapper function around the original qpAdm program
#'
#' This requires a working installation of qpAdm, which will be called using \code{\link{system}}
#'
#' @export
#' @param pref Path to and prefix of the packedancestrymap genotype files
#' @param target Target population
#' @param left Left populations (or leftlist file)
#' @param right Right populations (or rightlist file)
#' @param bin Path to the qpAdm binary file
#' @param outdir Output directory. files \code{out}, \code{parfile}, \code{leftlist},
#' \code{rightlist} will be overwritten
#' @param parfile qpAdm parameter file
#' @param allsnps allsnps
#' @param blgsize blgsize
#' @param fancyf4 fancyf4
#' @param f4mode f4mode
#' @param inbreed inbreed
#' @param printonly Should the command be printed or executed?
#' @param env Export environmental variables. See examples.
#' @param verbose Print progress updates
#' @return If not printonly, a data frame with parsed qpAdm output
#' @examples
#' \dontrun{
#' left = c('Altai_Neanderthal.DG', 'Vindija.DG')
#' right = c('Chimp.REF', 'Mbuti.DG', 'Russia_Ust_Ishim.DG', 'Switzerland_Bichon.SG')
#' target = 'Denisova.DG'
#' qpadm_wrapper('genotype_prefix', left, right, target,
#'   bin = 'path/to/qpAdm')
#' }
qpadm_wrapper = function(pref, left, right, target = NULL, bin = '~np29/o2bin/qpAdm',
                         outdir = './', parfile = NULL, allsnps = 'NO', blgsize = 0.05, fancyf4 = 'NO',
                         f4mode = 'YES', inbreed = 'NO', printonly = FALSE, env = '', verbose = TRUE) {

  stopifnot(!is.null(parfile) & is.null(c(target, left, right)) |
            !is.null(pref) & !is.null(left) & !is.null(right))
  stopifnot(file.exists(str_replace(bin, '.+ ', '')))
  #stopifnot(!is.null(target) | all(file.exists(c(left, right))))

  if(!is.null(target)) left = c(target, setdiff(left, target))
  outdir = normalizePath(outdir, mustWork = FALSE)
  if(is.null(parfile)) {

    if(all(file.exists(left, right))) {
      leftfile = left
      rightfile = right
    } else {
      leftfile = paste0(outdir, '/leftlist')
      rightfile = paste0(outdir, '/rightlist')
      write(left, leftfile)
      write(right, rightfile)
    }

    pref = normalizePath(pref, mustWork = FALSE)
    parfile = paste0('genotypename: ', pref, '.geno\n',
                     'snpname: ', pref, '.snp\n',
                     'indivname: ', pref, '.ind\n',
                     'popleft: ', leftfile, '\n',
                     'popright: ', rightfile, '\n',
                     'allsnps: ', allsnps %>% yesno, '\n',
                     'blgsize: ', blgsize, '\n',
                     'fancyf4: ', fancyf4 %>% yesno, '\n',
                     'f4mode: ', f4mode %>% yesno, '\n',
                     'inbreed: ', inbreed %>% yesno, '\n',
                     'details: YES\n',
                     'fstdetails: YES\n')

    parfilename = paste0(outdir, '/parfile')
    write(parfile, parfilename)

  } else {
    parfilename = parfile
  }

  cmd = paste0(env,' ', bin, ' -p ', parfilename, ' > ', outdir, '/out')
  outfile = paste0(outdir, '/out')
  run_admixtools(cmd, parse_qpadm_output, outfile, printonly, verbose)
}

#' Wrapper function around the original qpWave program
#'
#'
#' @export
#' @inheritParams qpadm
qpwave_wrapper = qpadm_wrapper


#' Wrapper function around the original qpGraph program
#' @export
#' @param pref Prefix of the packedancestrymap format genotype files.
#' @param graph An admixture graph or qpGraph graph file
#' @param bin Location of the qpGraph binary
#' @param parfile qpGraph parameter file
#' @param outdir Output directory
#' @param printonly Should output be executed or the command just be printed?
#' @param badsnps badsnps
#' @param lambdascale lambdascale
#' @param inbreed inbreed
#' @param diag diag
#' @param blgsize blgsize
#' @param outpop outgroup population
#' @param loadf3 loadf3
#' @param lsqmode least-squares mode. sets the offdiagonal elements of the block-jackknife covariance matrix to zero.
#' @param fstdmode fstdmode
#' @param hires hires
#' @param forcezmode forcezmode
#' @param zthresh zthresh
#' @param allsnps allsnps
#' @param oldallsnps oldallsnps
#' @param doanalysis doanalysis
#' @param bigiter bigiter
#' @param initmix initmix
#' @param env Export environmental variables. See examples.
#' @param verbose Print progress updates
#' @return A list with parsed qpGraph output
#' \enumerate{
#' \item \code{edges}: data frame
#' \item \code{score}: scalar
#' \item \code{f2}: data frame
#' }
#' @examples
#' \dontrun{
#' qpgraph_wrapper('genotype_prefix', example_graph,
#'                  bin = 'path/to/qpGraph')
#' }
qpgraph_wrapper = function(pref, graph, bin = '~np29/o2bin/qpGraph', parfile = NULL, outdir = '.',
                           printonly = FALSE, badsnps = NULL, lambdascale = -1, inbreed = 'NO',
                           diag = 0.0001, blgsize = 0.05, outpop = 'NULL', loadf3 = NULL,
                           lsqmode = 'NO', fstdmode = 'NO', hires = 'NO', forcezmode = 'NO', zthresh = 0,
                           allsnps = 'NO', oldallsnps = 'NO', doanalysis = 'YES',
                           bigiter = 100, initmix = 0, env = '', verbose = TRUE) {
  # wrapper around AdmixTools qpGraph
  # makes parfile and graphfile
  stopifnot(!is.null(parfile) | !is.null(pref))

  outdir = normalizePath(outdir, mustWork = FALSE)
  if(is.null(parfile)) {
    pref = normalizePath(pref, mustWork = FALSE) %>% str_replace('(?<!\\\\) ', '\\\\ ')
    badsnpfile = paste0(outdir, '/badsnps')
    fstatsfile = paste0(outdir, '/fstats.out')
    parfile = paste0('indivname:       ', pref, '.ind\n',
                     'snpname:         ', pref, '.snp\n',
                     'genotypename:    ', pref, '.geno\n',
                     'outpop:         ', outpop, '\n',
                     'blgsize: ', blgsize, '\n',
                     'details: YES\n',
                     'fstdetails: YES\n',
                     'diag: ', diag, '\n',
                     'lsqmode: ', lsqmode %>% yesno, '\n',
                     'fstdmode: ', fstdmode %>% yesno, '\n',
                     'hires: ', hires %>% yesno, '\n',
                     'forcezmode: ', forcezmode %>% yesno, '\n',
                     'zthresh: ', zthresh, '\n',
                     'allsnps: ', allsnps %>% yesno, '\n',
                     'oldallsnps: ', oldallsnps %>% yesno, '\n',
                     'lambdascale: ', lambdascale, '\n',
                     'fstatsoutname: ', fstatsfile, '\n',
                     'badsnpname: ', badsnpfile, '\n',
                     'inbreed: ', inbreed %>% yesno, '\n',
                     'bigiter: ', bigiter, '\n',
                     'initmix: ', initmix, '\n',
                     'doanalysis: ', doanalysis %>% yesno, '\n')
    if(!is.null(loadf3) && loadf3 == 'YES') parfile %<>% paste0('fstatsname: ', fstatsfile, '\n')
    if(allsnps == 'YES') parfile %<>% paste0('allsnps: YES\nloadf3: YES\n')
    if(oldallsnps == 'YES') parfile %<>% paste0('oldallsnps: YES\nallsnps: YES\nloadf3: NO\n')

    pf = paste0(outdir, '/parfile')
    write(parfile, pf)
    write(badsnps, badsnpfile)
  } else {
    pf = parfile
  }

  if(!'character' %in% class(graph)) {

    if(class(graph)[1] == 'igraph') {
      igraph = graph
      graph = igraph::as_edgelist(graph)
    } else {
      igraph = igraph::graph_from_edgelist(as.matrix(graph)[,1:2])
    }
    edg = as_tibble(graph, .name_repair = ~c('from', 'to'))
    edg %<>% group_by(to) %>% mutate(type = ifelse(n()==1, 'edge', 'admix')) %>% ungroup
    e1 = (edg %>% filter(type == 'edge'))$from
    e2 = (edg %>% filter(type == 'edge'))$to
    a1 = (edg %>% filter(type == 'admix'))$from
    a2 = (edg %>% filter(type == 'admix'))$to
    leaves = get_leafnames(igraph)
    root = setdiff(edg$from, edg$to)
    admix = tibble()
    for(m in unique(a2)) {
      admix %<>% bind_rows(tibble(v1='admix', v2=m, v3=(edg %>% filter(to == m))$from[1],
                                  v4=(edg %>% filter(to == m))$from[2]))
    }

    simfile = tibble(v1 = c('root'), v2 = root, v3='', v4='') %>%
      bind_rows(tibble(v1 = 'label', v2=leaves, v3=leaves, v4='')) %>%
      bind_rows(tibble(v1='edge', v2=paste0('e', 1:length(e1)), v3=e1, v4=e2)) %>%
      bind_rows(admix)

    gf = paste0(outdir, '/graphfile')
    simfile %>% write_tsv(gf, col_names=F)
  } else {
    stopifnot(file.exists(graph))
    gf = normalizePath(graph)
  }

  qpgraph_wrapper2(bin = paste0(env,' ', bin), parfile = pf, graphfile = gf,
                   outfile = paste0(outdir, '/out'), printonly = printonly, verbose = verbose)
}


qpgraph_wrapper2 = function(bin='./qpGraph', parfile='./parfile', graphfile='./graphfile',
                            outfile='./out', printonly=FALSE, verbose = TRUE) {
  # wrapper around AdmixTools qpGraph
  # input is locations of parfile and graphfile
  # output is parsed output

  cmd = paste0(bin, ' -p ', parfile, ' -g ', graphfile, ' > ', outfile)
  run_admixtools(cmd, parse_qpgraph_output, outfile, printonly, verbose)
}



qpfstats_wrapper = function(pref, pops, bin = '~np29/o2bin/qpfstats', parfile = NULL,
                            outdir = './', allsnps = 'NO', inbreed = 'NO', scale = 'NO', details = 'YES', hires = 'YES',
                            printonly = FALSE, env = '', verbose = TRUE) {

  stopifnot(is.null(parfile) || is.null(pops))
  stopifnot(file.exists(str_replace(bin, '.+ ', '')))

  outdir = normalizePath(outdir, mustWork = FALSE)

  if(is.null(parfile)) {

    poplistfile = paste0(outdir, '/poplist')
    fstatsfile = paste0(outdir, '/fstats')
    writeLines(pops, poplistfile)

    pref = normalizePath(pref, mustWork = FALSE)
    parfile = paste0('genotypename: ', pref, '.geno\n',
                     'snpname: ', pref, '.snp\n',
                     'indivname: ', pref, '.ind\n',
                     'poplistname: ', poplistfile, '\n',
                     'fstatsoutname: ', fstatsfile, '\n',
                     'allsnps: ', allsnps %>% yesno, '\n',
                     'inbreed: ', inbreed %>% yesno, '\n',
                     'details: ', details %>% yesno, '\n',
                     'hires: ', hires %>% yesno, '\n',
                     'scale: ', scale %>% yesno, '\n')

    parfilename = paste0(outdir, '/parfile')
    write(parfile, parfilename)

  } else {
    parfilename = parfile
  }

  cmd = paste0(env,' ', bin, ' -p ', parfilename, ' > ', outdir, '/out')
  outfile = paste0(outdir, '/out')
  run_admixtools(cmd, parse_qpfstats_output, outfile, printonly, verbose)
}


run_admixtools = function(cmd, parsefun, outfile, printonly, verbose) {
  if(printonly) {
    print(cmd)
  } else{
    if(verbose) {
      alert_info('Running admixtools:\n')
      cat(paste0(cmd, '\n'))
    }
    system(cmd)
    if(verbose) alert_info('Parsing output\n')
    return(parsefun(outfile))
  }
}

#' @export
parse_fstats = function(outfile, denom1 = 1e3, denom2 = 1e7) {
  # parse Nick's qpGraph fstats file

  dat = read_lines(outfile) %>% str_squish
  basepop = dat[1] %>% str_replace('.+basepop: ', '') %>% str_replace(' .+', '')
  dat = dat[-1]

  f3 = dat %>% {str_count(., ' ') == 2} %>% `[`(dat, .) %>% enframe %>%
    separate(value, c('pop2', 'pop3', 'est'), ' ', T, T)
  f3var = dat %>% {str_count(., ' ') > 2} %>% `[`(dat, .) %>% enframe %>%
    separate(value, c('pop1', 'pop2', 'pop3', 'pop4', 'se2'), ' ', T, T)
  pops = c(basepop, unique(c(f3$pop2, f3$pop3)))
  npop = length(pops)
  npair = choose(npop, 2)
  npair2 = choose(choose(npop,2)+1, 2)
  p1 = split(1:npair2, rep(1:npair, npair:1))
  p2 = sapply(2:npair, function(i) rep(i, i-1)-npair + cumsum(npair:(npair-i+2)))
  indices = unname(unlist(interleave(p1, p2)))

  f3_est = (f3$est/denom1) %>% structure(pops = pops)
  f3varmat = matrix(f3var$se2[indices]/denom2, npair, npair)*10
  ppinv = solve(f3varmat) %>% structure(pops = pops)
  f3out = f3 %>% mutate(pop1 = basepop, est=est/denom1, se = sqrt(diag(f3varmat))) %>%
    select(pop1, pop2, pop3, est, se)

  namedList(f3_est, ppinv, pops, f3out)
}

# for Nick
write_qpgraph_output = function(fit, outfile = './out', decimals = 3, labels = NULL, counts = NULL) {

  sep = '\t'
  paste('final_score:', fit$score, sep = sep) %>% write(outfile)
  if('f2' %in% names(fit)) {
    write('f2:', outfile, append = TRUE)
    fit$f2 %>% write_delim(outfile, delim = sep, append = TRUE, col_names = TRUE)
  }
  if('f3' %in% names(fit)) {
    write('f3:', outfile, append = TRUE)
    fit$f3 %>% write_delim(outfile, delim = sep, append = TRUE, col_names = TRUE)
  }

  c('graph:', fit_to_qpgraph_format(fit$edges, decimals = decimals, sep = sep)) %>%
    write(outfile, append = TRUE)
}

fit_to_qpgraph_format = function(edges, decimals = 3, sep = '\t') {
  # edges is data frame with columns type, from, to, weight
  leaves = setdiff(edges$to, edges$from)
  vertex = paste('vertex', unique(c(t(as.matrix(cbind(edges$from, edges$to))))), '0', sep = sep)
  label = paste('label', leaves, leaves, sep = sep)
  norm = edges %>% filter(type == 'edge') %$%
    paste('edge', paste0(from, '_', to), from, to, round(weight, decimals), sep = sep)
  admix = edges %>% filter(type == 'admix') %>% group_by(to) %>%
    summarize(from = paste(from, collapse = sep),
              weight = paste(round(weight, decimals), collapse = sep)) %$%
    paste(rep('admix', length(to)), to, from, weight, sep = sep)
  c(vertex, label, norm, admix)
}

igraph_to_qpgraph = function(graph, outfile, sep = '\t') {
  edges = igraph::as_edgelist(graph) %>%
    as_tibble(.name_repair = ~c('from', 'to')) %>%
    add_count(to) %>%
    mutate(type = ifelse(n == 2, 'admix', 'edge'))
  leaves = setdiff(edges$to, edges$from)
  root = paste('root', setdiff(edges$from, edges$to), sep = sep)
  label = paste('label', leaves, leaves, sep = sep)
  norm = edges %>% filter(type == 'edge') %$%
    paste('edge', paste0(from, '_', to), from, to, sep = sep)
  admix = edges %>% filter(type == 'admix') %>% group_by(to) %>%
    summarize(from = paste(from, collapse = sep)) %$%
    paste(rep('admix', length(to)), to, from, sep = sep)
  c(root, label, norm, admix) %>%
  write(outfile)
}

#' Read qpGraph output file
#' @export
#' @param outfile output file generated by qpGraph.
#' @return list of output data.
parse_qpgraph_output = function(outfile, logfile = outfile) {
  # reads qpGraph output file
  # returns list of three objects:
  # 'edges': data.frame of branch lengths and admixture weights
  # 'score': best fit score
  # 'f2': data.frame of estimated and fitted f2 values

  # addition of outliers currently not safe when two population have the same three-letter-prefix

  edges = parse_qpgraph_output_edges(outfile)

  dat = readLines(logfile) %>% str_squish() %>% enframe(value = 'X1') %>% select(X1)

  score = (dat %>% filter(grepl('^final score', X1)) %>%
             separate('X1', c('a', 'b', 'score'), sep=' +', convert = T, extra='drop', fill='right'))$score

  pops = dat$X1 %>% str_subset('^population:') %>% str_squish %>% word(3)
  if(length(pops) == 1 && is.na(pops)) {
    pshort = dat$X1 %>% str_subset(' f2: ') %>% str_squish
    pshort = c(word(pshort[1], 1), unique(word(pshort, 2)))
    pfull = edges %>% edges_to_igraph() %>% get_leafnames()
    pfull = pfull[-1]
    pre3 = str_sub(pfull, 1,3)
    if(any(duplicated(pre3))) stop('Prefixes duplicated! Need to fix parse_qpgraph_output')
    stopifnot(isTRUE(all.equal(sort(pre3), sort(pshort))))
    pops = pfull[match(pshort, pre3)]
  }
  numpop = length(pops)
  numpair = choose(numpop, 2)

  f2 = dat %>% filter(grepl(' f2: ', X1)) %>%
    separate('X1', c('pop1', 'pop2', 'fst','fit','est','diff','se','z'), sep=' +', convert = TRUE) %>%
    select(-fst) %>%
    mutate(pop1 = rep(head(pops, -1), (numpop-1):1), pop2 = unlist(map(2:numpop, ~pops[.:numpop])))

  f3dat = dat %>% filter(grepl(' ff3fit: ', X1))
  if(nrow(f3dat) > 0) {
    f3 = f3dat %>%
      separate('X1', c('pop2', 'pop3', 'ff3fit','fit','est'), sep=' +', convert = TRUE) %>%
      select(-ff3fit) %>%
      mutate(pop2 = rep(pops[-1], each = numpop-1), pop3 = rep(pops[-1], numpop-1))
  } else {
    mul = 1000
    f3eststart = str_which(dat$X1, '^ff3:')[1]+2
    f3estend = str_which(dat$X1, '^ff3sig\\*10:')[1]-2
    f3sigstart = f3estend+4
    f3sigend = f3sigstart + (f3estend-f3eststart)
    f3fitstart = str_which(dat$X1, '^ff3fit:')[1]+2
    f3fitend = f3fitstart + (f3estend-f3eststart)
    fn = function(start, end) dat %>% slice(start:end) %>% pull(X1) %>% str_squish() %>%
      str_sub(3) %>% str_split(' ') %>% unlist %>% as.numeric
    f3 = expand_grid(pop2 = pops, pop3 = pops) %>% mutate(pop1 = pops[1], .before = 1) %>%
      mutate(est = fn(f3eststart, f3estend)/mul, fit = fn(f3fitstart, f3fitend)/mul, se = fn(f3sigstart, f3sigend)/mul/10) %>%
      mutate(diff = est - fit, z = diff/se)

    #f3 = NULL
  }

  outlierstart = str_which(dat$X1, '^outliers:')[1]+2
  outlierend = str_which(dat$X1, '^worst f-stat:')[1]-3

  numalloutliers = choose(choose(numpop, 2)+1, 2)
  outliersmissing = outlierend - outlierstart + 1 != numalloutliers
  if(outliersmissing) warning('Some outliers are missing. No f4-stats returned.')

  #if(is.na(pops)[1]) pops = dat$X1 %>% str_subset('^zzaddw') %>% str_split(' ') %>% map(2) %>% unlist %>% unique
  denom = 1000
  amb = names(which(table(str_sub(pops, 1, 3)) > 1))
  if(length(amb) > 0) warning(paste('Ambiguous populations ommited from outliers: ', paste0(amb, collapse = ', ')))

  if(!outliersmissing) {

  poppairs = t(combn(pops, 2))
  popquads = cbind(poppairs[rep(1:numpair, numpair:1),], poppairs[unlist(map(1:numpair, ~(.:numpair))),])

  outliers = dat$X1[outlierstart:outlierend] %>%
    str_split(' +') %>%
    do.call(rbind, .) %>%
    `[`(,-1:-4) %>%
    cbind(popquads, .) %>%
    as_tibble(.name_repair = ~c(paste0('pop', 1:4), c('fit', 'est', 'diff', 'se', 'z'))) %>%
    type_convert(col_types = cols()) %>%
    mutate(diff = -diff, z = -z) #%>%
    #filter(length(intersect(amb, c(pop1, pop2, pop3, pop4))) == 0)
  #f3 %<>% left_join(outliers %>% filter(pop1 == pop3, pop1 == str_sub(pops[1], 1, 3)) %>%
  #                    select(pop2, pop4, diff, se, z), by = c('pop2', 'pop3'='pop4'))

  # if(length(unique(str_sub(pops, 1, 3))) == length(unique(pops))) {
  #   outliers %<>% mutate(across(pop1:pop4, ~pops[match(., str_sub(pops, 1, 3))]))
  # }
  } else {
    outliers = NULL
  }

  #f2 %<>% mutate(pop1 = rep(head(pops, -1), (numpop-1):1), pop2 = unlist(map(2:numpop, ~pops[.:numpop])))
  #f3 %<>% mutate(pop2 = rep(pops[-1], each = numpop-1), pop3 = rep(pops[-1], numpop-1))

  namedList(edges, score, f2, f3, f4 = outliers)
}

parse_qpgraph_output_edges = function(outfile, uselabels = FALSE) {

  dat = readLines(outfile) %>% str_squish() %>% enframe(value = 'X1') %>% select(X1)

  nammap = dat %>%
    filter(grepl('^label', X1)) %>%
    separate('X1', c('x', 'label', 'name'), sep=' +', extra = 'drop') %>%
    select(label, name) %>%
    deframe

  edges = dat %>%
    filter(grepl('^edge|^ledge|^redge|^admix', X1)) %>%
    separate('X1', c('type', 'name', 'from', 'to', 'weight', 'w2'),
             sep=' +', convert = T, extra='drop', fill='right')
  admix1 = edges %>% filter(type=='admix') %>%
    mutate(type='aedge', to=name, name='', w2=NA)
  admix2 = edges %>% filter(type=='admix') %>%
    mutate(type='aedge', from=to, to=name, name='', weight=w2, w2=NA)
  edges %<>%
    bind_rows(admix1) %>%
    bind_rows(admix2) %>%
    filter(!type == 'admix') %>%
    mutate(type = ifelse(type=='aedge', 'admix', 'edge')) %>%
    select(from, to, type, weight)
  if(!uselabels) edges %<>% mutate(from = recode(from, !!!nammap), to = recode(to, !!!nammap))
  edges
}



#' Read qpGraph graph file
#' @export
#' @param graphfile File with admixture graph in qpGraph format.
#' @param split_multi Split multifurcations
#' @param igraph Convert to igraph format
#' @param uselabels Should labels or full names be used? Defaults to full names.
#' @return Graph represented as two column edge matrix. Can have four columns if edges are locked
parse_qpgraph_graphfile = function(graphfile, split_multi = TRUE, igraph = FALSE, uselabels = FALSE) {
  # reads graph in qpGraph format
  # returns edge matrix (adjacency list)
  lines = read_lines(graphfile) %>% enframe %>% transmute(V1 = value)
  nammap = lines %>% filter(grepl('^label', V1)) %>%
    separate('V1', c('type', 'name', 'label'), sep = '\\s+', extra = 'drop') %>%
    select(name, label) %>% deframe

  dat = lines %>%
    filter(grepl('^edge|^redge|^ledge|^admix|^lock', V1)) %>%
    separate('V1', c('type', 'name', 'from', 'to', 'l1', 'l2'),
             sep = '\\s+', extra = 'drop', fill = 'right') %>%
    mutate(type = recode(type, ledge = 'edge', redge = 'edge'))
  locks = dat %>% filter(type=='lock') %$% name
  admix = dat %>% filter(type == 'admix')
  admix1 = admix %>% mutate(type = 'edge', to = name,
                            lower = ifelse(name %in% locks, as.numeric(l1)/100, NA),
                            upper = lower, name = '')
  admix2 = admix %>% mutate(type = 'edge', from = to, to = name,
                            lower = ifelse(name %in% locks, as.numeric(l2)/100, NA),
                            upper = lower, name = '')
  out = dat %>%
    filter(type != 'lock') %>%
    bind_rows(admix1) %>%
    bind_rows(admix2) %>%
    filter(type != 'admix') %>%
    select(from, to, lower, upper)
  if(!uselabels) out %<>% mutate(from = recode(from, !!!nammap), to = recode(to, !!!nammap))
  if(all(is.na(out$lower)) && all(is.na(out$upper))) out %<>% select(-lower, -upper)
  #out %<>% as.matrix
  if(split_multi) out %<>% split_multifurcations()
  if(igraph) out %<>% edges_to_igraph()
  out
}


# Read qpGraph parameter file
# @export
# @param parfile parameter file for qpGraph.
# @return named list of qpGraph parameters.
parse_qpgraph_parfile = function(parfile) {
  # reads qpGraph parfile
  # returns named list of parameters
  # all genotype files have to have same prefix

  dat = read_table2(parfile, comment = '#', col_names = c('par', 'value'), col_types = cols())
  dat %<>% mutate(par = str_replace_all(par, c(':$'='', 'genotypename'='pref')),
                  value = str_replace_all(value,
                                          c('\\.geno$'='',
                                            'S1'=filter(dat, par=='S1:')$value[1],
                                            'DIR'=filter(dat, par=='DIR:')$value[1])),
                  value=ifelse(value=='YES', TRUE, ifelse(value=='NO', FALSE, value))) %>%
    filter(!par %in% c('DIR', 'S1', 'indivname', 'snpname')) %>%
    t %>%
    as_tibble()
  dat %>%
    set_colnames(slice(dat, 1)) %>%
    slice(-1) %>%
    type_convert(col_types = cols()) %>%
    as.list
}


# Read qpDstat parameter file
# @export
# @param parfile parameter file for qpDstat.
# @return named list of qpGraph parameters.
parse_qpdstat_parfile = function(parfile) {
  # reads qpGraph parfile
  # returns named list of parameters
  # all genotype files have to have same prefix

  dat = read_table2(parfile, comment = '#', col_names = c('par', 'value'), col_types = cols())
  dat %>% mutate(par = str_replace_all(par, c(':$'='', 'genotypename'='pref')),
                 value = str_replace_all(value,
                                         c('\\.geno$'='',
                                           'SSS'=filter(dat, par=='SSS:')$value[1],
                                           'DIR'=filter(dat, par=='DIR:')$value[1])),
                 value=ifelse(value=='YES', TRUE, ifelse(value=='NO', FALSE, value))) %>%
    filter(!par %in% c('DIR', 'SSS', 'indivname', 'snpname')) %>%
    deframe %>%
    as.list %>%
    as_tibble()
}

#' Read qpF4ratio output file
#' @export
#' @param outfile Output file generated by qpF4ratio
#' @return Data frame with output data
parse_qpf4ratio_output = function(parfile) {
  parfile %>%
    readLines %>%
    str_subset(' result: ') %>%
    str_split(' +') %>%
    map(`[`, c(3:6,10,12:15)) %>%
    do.call('rbind', .) %>%
    set_colnames(c(paste0('pop', 1:5), 'alpha', 'se', 'z', 'n')) %>%
    as_tibble %>%
    type_convert(col_types = cols())
}

#' Read qpAdm output file
#' @export
#' @param outfile Output file generated by qpAdm.
#' @return Data frame with output data.
parse_qpadm_output = function(outfile) {
  # reads qpAdm output file

  dat = read_lines(outfile)

  lstart = str_which(dat, 'left pops:')[1]+1
  rstart = str_which(dat, 'right pops:')[1]+1
  lend = str_which(dat, '^$') %>% magrittr::extract(. > lstart) %>% head(1)-1
  rend = str_which(dat, '^$') %>% magrittr::extract(. > rstart) %>% head(1)-1
  coefstart = str_which(dat, '^best coefficients:')[1]
  sigstart = str_which(dat, 'fixed pat')[1]+1
  sigend = str_which(dat, '^best pat:')[1]-1

  target = dat[lstart]
  left = dat[lstart:lend][-1]
  right = dat[rstart:rend]
  coefs = dat[coefstart:(coefstart+2)]
  coefs %<>% str_split(' +') %>%
    map(~tail(., length(left)+1) %>% head(-1) %>% as.numeric %>% set_names(left)) %>%
    set_names(c('weight', 'mean', 'se')) %>% as_tibble %>% mutate(z = mean/se)
  weights = tibble(target, left) %>% bind_cols(coefs)

  popdrop = do.call(rbind, str_split(dat[sigstart:sigend], ' +')) %>%
    as.data.frame(stringsAsFactors=F) %>% select(-1) %>%
    set_colnames(c('pat', 'wt', 'dof', 'chisq', 'p', left, 'feasible')) %>%
    mutate(feasible = feasible != 'infeasible', across(!c('pat', 'feasible'), as.numeric)) %>% as_tibble

  namedList(weights, popdrop)
}


#' Read qpDstat output file
#' @export
#' @param outfile output file generated by qpDstat.
#' @return tibble with output data.
parse_qpdstat_output = function(outfile) {

  dat = read_lines(outfile) %>% str_subset('^result: ') %>% str_squish
  stopifnot(length(dat) > 0)
  nc = dat[[1]] %>% str_split(' +', simplify = T) %>% ncol
  nam1 = c('result', 'pop1', 'pop2', 'pop3', 'pop4', 'f4', 'Z')
  nam2 = c('BABA', 'ABBA', 'numsnps')
  if(nc == 10) nam = c(nam1, nam2)
  else if(nc == 11) nam = c(nam1, 'best', nam2)
  else stop('Unexpected number of columns in output!')

  dat %>% as_tibble %>% separate(value, nam, sep = ' +', convert = TRUE) %>%
    select(-result) %>% mutate(se = f4/Z, p.value = ztop(Z))
}


#' Read qp3Pop output file
#' @export
#' @param outfile output file generated by qp3Pop.
#' @return tibble with output data.
parse_qp3pop_output = function(outfile) {

  dat = read_lines(outfile) %>% str_replace(' no data', 'result:') %>% str_subset('^ result: ') %>% str_squish
  stopifnot(length(dat) > 0)
  nam = c('result', 'source1', 'source2', 'target', 'f3', 'se', 'Z', 'numsnps')

  dat %>% as_tibble %>% separate(value, nam, sep = ' +', convert = TRUE) %>%
    select(-result) %>% mutate(se = f3/Z, p.value = ztop(Z))
}


parse_qpfstats_output = function(outfile) {

  dat = read_lines(outfile)
  pops = dat %>% str_subset('^population: ') %>% str_squish %>% str_split(' ') %>% map(3) %>% unlist

  nam1 = c('p11', 'p12', 'p21', 'p22', 'f21', 'f22', 'se')
  nam2 = c('p1', 'p2', 'p3', 'p4', 'f2', 'se', 'z')
  nam3 = c('p1', 'p2', 'f21', 'se1', 'f22', 'se2')

  djest = dat %>% str_subset('^jest')
  if(length(djest) > 0) jest = djest %>% str_squish %>% str_split(' ') %>% do.call(rbind, .) %>% `[`(,-1:-3) %>% apply(2, as.numeric) %>% as_tibble(.name_repair = ~nam1) %>% mutate(across(p11:p22, ~pops[.+1]))
  else jest = NULL

  dfstat = dat %>% str_subset('^fstat:')
  if(length(dfstat) > 0) fstat = dfstat %>% str_squish %>% str_split(' ') %>% do.call(rbind, .) %>% `[`(,-1) %>% apply(2, as.numeric) %>% as_tibble(.name_repair = ~nam2) %>% mutate(across(p1:p4, ~pops[.+1])) %>%
    rename(pop1 = p1, pop2 = p2, pop3 = p3, pop4 = p4)
  else fstat = NULL

  basis = dat %>% str_subset('^basis:') %>% str_squish %>% str_split(' ') %>% do.call(rbind, .) %>% `[`(,-c(1,6)) %>% apply(2, as.numeric) %>% as_tibble(.name_repair = ~nam3) %>% mutate(across(p1:p2, ~pops[.+1])) %>%
    rename(pop1 = p1, pop2 = p2)

  namedList(jest, fstat, basis)
}

parse_qpfstats = function(outfile) {

  dat = read_lines(outfile) %>% str_squish()
  len = str_split(dat, '\\s+') %>% map_dbl(length)
  basepop = str_split(dat[1], '\\s+')[[1]][3]

  f3est = dat[len == 3] %>% str_split('\\s+') %>% do.call(rbind, .) %>%
    as_tibble(.name_repair = ~c('pop2', 'pop3', 'est')) %>%
    add_column(pop1 = basepop, .before = 1) %>% mutate(est = as.numeric(est)/1e3)

  f3cov = dat[len == 5] %>% str_split('\\s+') %>% do.call(rbind, .) %>%
    as_tibble(.name_repair = ~c('pop11', 'pop12', 'pop21', 'pop22', 'cov')) %>%
    mutate(cov = as.numeric(cov)/1e6)
  f3se = f3cov %>% filter(pop11 == pop21, pop12 == pop22) %>%
    transmute(pop2 = pop11, pop3 = pop12, se = sqrt(cov))

  f3est %<>% left_join(f3se, by = c('pop2', 'pop3'))

  namedList(f3est, f3cov)
}

parse_qpfmv_output = function(outfile) {

  outfile %>%
    read_lines %>%
    str_subset('^result:') %>%
    str_split('\\s+') %>%
    do.call(rbind, .) %>%
    as_tibble(.name_repair = ~c('X1', paste0('pop', 1:4), 'est', 'z', 'X2')) %>%
    select(-X1, -X2) %>%
    mutate(est = as.numeric(est), z = as.numeric(z))
}

# for Nick
parse_qpff3base_output = function(outfile, denom = 1000) {

  dat = read_lines(outfile) %>% str_squish
  stopifnot(length(dat) > 0)

  popstart = str_which(dat, '^end of inpack')[1]+1
  popend = str_which(dat, '^outpop:.+basepop:')[1]-1
  fststart = str_which(dat, '^fst:')[1]+2
  fstend = str_which(dat, '^f2:')[1]-2
  f2start = fstend+4
  f2end = str_which(dat, '^ff3 \\(unscaled\\):')[1]-2

  pops = dat[popstart:popend] %>% str_replace('.+ ', '')
  fst = dat[fststart:fstend] %>% enframe %>% separate(value, c('pop1', pops), ' ', T, T) %>%
    select(-name) %>% mutate(pop1 = pops) %>% pivot_longer(-pop1, 'pop2', values_to = 'fst') %>% mutate(fst = fst/denom)
  f2 = dat[f2start:f2end] %>% enframe %>% separate(value, c('pop1', pops), ' ', T, T) %>%
    select(-name) %>% mutate(pop1 = pops) %>% pivot_longer(-pop1, 'pop2', values_to = 'f2') %>% mutate(f2 = f2/denom)

  fst %>% left_join(f2, by = c('pop1', 'pop2'))
}

#' Get pseudo dates for graph nodes
#'
#' This function assigns a date to each node in an admixture graph and is used in \code{\link{msprime_sim}}. 
#' The date will correspond to the y-coordinate of each node used for plotting in \code{\link{plotly_graph}} 
#' unless `fix` option is set \code{TRUE}, in which case all leaf nodes returned by \code{\link{get_leafnames()}} will be set to 0.
#' @export
#' @param graph An admixture graph
#' @param time Scalar by which y-coordinate values will be multiplied to get dates
#' @param fix Boolean specifying if leaf nodes will be fixed to 0. If `TRUE` (default), all samples will be drawn at the end of the simulation (i.e., from today).
#' @return A named vector with pseudo dates for each graph node
pseudo_dates = function(graph, time = 1000, fix=TRUE) {
  edges = graph %>% as_edgelist()
  pdat = graph_to_plotdat(edges)$eg
  out = bind_rows(transmute(pdat, name, y), transmute(pdat, name = to, y = yend)) %>%
    mutate(y = y - min(y)) %>% distinct %>% deframe %>% multiply_by(time)
  if (isTRUE(fix))  out[get_leafnames(graph)] = 0
  out
}

#' Simulate an admixture graph in msprime
#'
#' This function generates an msprime simulation script, and optionally executes it in python
#' @export
#' @param graph Graph as `igraph` object or edge list with columns 'from' and 'to'. If it's an edge list with column 'weight' (possibly from a fitted graph), the admixture weights will be used. Otherwise, all admixture edges will have weight 0.5.
#' @param outpref Prefix of output files
#' @param nsnps The number of SNPs to simulate. All SNPs will be simulated independently of each other.
#' @param neff Effective population size. If a scalar value, it will be constant across all populations. Alternatively, it can be a named vector with a different value for each population.
#' @param ind_per_pop The number of individuals to simulate for each population. If a scalar value, it will be constant across all populations. Alternatively, it can be a named vector with a different value for each population.
#' @param mutation_rate Mutation rate per site per generation. The default is set to a high value (0.001) to obtain more polymorphic SNPs in order to speed up the simulation.
#' @param time Either a scalar value (1000 by default) which the dates generated by \code{\link{pseudo_dates}}, or a named vector with dates for each graph node.
#' @param admix_default The default admixture proportions for all admixture nodes
#' @param run If `FALSE`, the function will terminate after writing the msprime script. If `TRUE`, it will try and execute the function with the default python installation. If you want to use some other python installation, you can set `run = /my/python`.
#' @param numcores The number of cores to use when simulating data.
#' @return The file name and path of the simulation script
#' @examples
#' \dontrun{
#' results = qpgraph(example_f2_blocks, example_graph)
#' msprime_sim(results$edges)
#' }
msprime_sim = function(graph, outpref = 'msprime_sim', nsnps = 1000, neff = 1000, ind_per_pop = 1,
                       mutation_rate = 1e-3, time = 1000, admix_default = 0.5, run = FALSE, numcores = NULL, shorten_admixed_leaves = FALSE) {

  outpref %<>% normalizePath(mustWork = FALSE)
  if('igraph' %in% class(graph)) edges = as_edgelist(graph) %>% as_tibble(.name_repair = ~c('from', 'to'))
  else if(is.matrix(graph) && ncol(graph) == 2) edges = as_tibble(graph, .name_repair = ~c('from', 'to'))
  else edges = graph

  edges %<>% add_count(to) %>%
    mutate(type = ifelse(n > 1, 'admix', 'normal')) %>% select(-n)
  adm = edges %>% filter(type == 'admix') %>% pull(to)
  graph = igraph::graph_from_edgelist(as.matrix(edges[,1:2]))
  nodes = names(V(graph))
  leaves = get_leafnames(graph)
  leaves = intersect(nodes, leaves)

  if(length(time) == 1) {

    dates = pseudo_dates(graph, time)
    #dates[leaves] = dates[leaves] * 0.5
    if(shorten_admixed_leaves) {
      # simulates the behavior of treemix where drift after admixture is not allowed
      admleaves = edges %>% filter(to %in% leaves, from %in% adm) %>% select(from, to) %>% deframe
      dates[admleaves] = 0.999 * dates[names(admleaves)]
    }
    #edges %<>% mutate(date = dates[edges$to]*time)
    edges %<>% mutate(date = dates[edges$from])

  } else {

    edges %<>% left_join(enframe(time, 'from', 'date'), by = 'from')
    dates = time

  }
  edges %<>% mutate(source = match(to, nodes)-1, dest = match(from, nodes)-1) %>%
    arrange(date)
  #dates[leaves[1]] = 0

  # continue here: popsize controls neff and adm weights?
  if(!'weight' %in% names(edges)) edges %<>% group_by(to) %>%
    mutate(weight = ifelse(type == 'admix', ifelse(from == min(from), admix_default, 1-admix_default), 1)) %>% ungroup
  popsize = edges %>% transmute(to, w = ifelse(type == 'admix', 1, 1/weight)) %>% distinct %>% deframe
  popsize[setdiff(edges$from, edges$to)] = 1
  edges %<>% group_by(to) %>% mutate(i = 1:n(), weight = ifelse(type == 'admix' & i == 1, weight, 1)) %>% ungroup
  #date = tibble(to = nodes) %>% left_join(edges %>% select(to, date) %>% distinct, by = 'to') %>% deframe %>% replace_na(0)
  #date = date - time

  out = "import math\nimport numpy\nimport msprime\nimport multiprocessing\n"
  out = "import math\nimport numpy\nimport msprime\nimport multiprocess\n"

  out = paste0(out, '\nnsnps = int(', nsnps, ')')
  if(length(neff) == 1) initsize = round(neff*replace_na(popsize[nodes], neff*100), 2)
  else if(all(nodes %in% names(neff))) initsize = neff[nodes]
  else stop("'neff' has to be a single number or a named vector with a number for each node!")
  out = paste0(out, '\npops = [\n', paste0('  msprime.PopulationConfiguration(initial_size = ', initsize,
                                        ')', c(rep(',', length(nodes)-1), ''),
                                        ' #',(1:length(nodes))-1, ' ', nodes, '\n', collapse = '') ,']\n')
  out = paste0(out, '\ndate = [', paste0(dates[nodes], collapse = ', '),']\n')

  lnum = match(leaves, nodes) - 1

  if(length(ind_per_pop) > 1) {
    if(!isTRUE(all.equal(sort(names(ind_per_pop)), sort(leaves)))) stop("'ind_per_pop' has to be a single number or a named vector with a number for each leaf node!")

  out = paste0(out,
  '\nind_per_pop = [',paste0(ind_per_pop[leaves], collapse = ', '),']',
  '\nss = [j for l in ind_per_pop for j in range(l)]',
  '\npp = list(numpy.repeat([', paste0(lnum, collapse = ', '), '], ind_per_pop))',
  '\npp2 = list(numpy.repeat(["', paste0(leaves, collapse = '", "'), '"], ind_per_pop))',
  '\nsamples = [msprime.Sample(pp[i], 0) for i in range(len(pp)) for j in range(2)]',
  '\nindnam = [str(pp2[i]) + "_" + str(ss[i]+1) for i in range(len(pp2))]')

  } else {

    indnam = paste(rep(leaves, each = ind_per_pop), seq_len(ind_per_pop), sep = '_')
    out = paste0(out, '\nindnam = ["', paste0(indnam, collapse = '", "'), '"]')
    out = paste0(out, '\n\nsamples = [msprime.Sample(j, max(0, date[j])) for j in [', paste(lnum, collapse = ', '),'] for i in range(',2*ind_per_pop,')]')
  }

  out = paste0(out, '\n\nevents = [\n',
               paste0('  msprime.MassMigration(time = ',
                      edges$date, ', source = ', edges$source,', destination = ', edges$dest,
                      ', proportion = ', edges$weight,')',
                      collapse = ',\n'), '\n]\n')

  out = paste0(out, "\ngt = numpy.zeros((int(nsnps),len(samples)//2), 'int')")

  out = paste0(out, "\n\ndef f(i):", "\n  tree_sequence = msprime.simulate(population_configurations = pops, samples = samples, demographic_events = events, mutation_rate = ", mutation_rate,")",
               "\n  if tree_sequence.genotype_matrix().shape[0] > 0:",
               '\n    return (tree_sequence.genotype_matrix()[0,range(0,len(samples),2)] + tree_sequence.genotype_matrix()[0,range(1,len(samples),2)])',
               "\n  else:",
               "\n    return numpy.zeros(len(samples)//2, 'int')")

  if(is.null(numcores)) numcores = 100
  out = paste0(out, "\n\np = multiprocess.Pool(int(min(multiprocess.cpu_count(),", numcores,")))")
  out = paste0(out, "\ngt = numpy.array(p.map(f, range(nsnps)))")

  # out = paste0(out, "\nfor i in range(int(",nsnps,")):")
  # out = paste0(out, paste0("\n  tree_sequence = msprime.simulate(population_configurations = pops, samples = samples, demographic_events = events, mutation_rate = ", mutation_rate,")"))
  # out = paste0(out, "\n  if tree_sequence.genotype_matrix().shape[0] > 0:")
  # out = paste0(out, '\n    gt[i,:] = (tree_sequence.genotype_matrix()[0,range(0,len(samples),2)] + tree_sequence.genotype_matrix()[0,range(1,len(samples),2)])')
  # out = paste0(out, "\n  else:")
  # out = paste0(out, "\n    gt[i,:] = numpy.zeros((1, len(samples)//2))")

  out = paste0(out, "\n\nnumpy.savetxt('",outpref,".geno', gt, '%d', '')")

  out = paste0(out, "\n\nwith open('",outpref,".snp', 'w') as f:")
  out = paste0(out, "\n  for i in range(nsnps):")
  out = paste0(out, "\n    bytes = f.write('rs'+str(i+1)+'\\t1\\t' + str(i*50/float(nsnps-1)) + '\\t' + str(i*100) + '\\tA\\tC\\n')")
  out = paste0(out, "\n\nwith open('",outpref,".ind', 'w') as f:")
  out = paste0(out, "\n  for i in range(len(indnam)):")
  out = paste0(out, "\n    bytes = f.write(indnam[i] + '\\tU\\t' + indnam[i].rstrip(\"1234567890\").rstrip(\"_\") + '\\n')\n")

  outfilename = paste0(outpref, '.py')
  writeLines(out, outfilename)


  if(run != FALSE) {
    if(isTRUE(run)) run = 'python'
    system(paste(run, outfilename))
  }
  tools::file_path_as_absolute(outfilename)
}

#' Simulate an admixture graph in msprime
#'
#' This function generates an msprime simulation script, executes it in python,
#' and turns the resulting genotype data into f2-statistics
#' @inheritParams msprime_sim
#' @seealso \code{\link{msprime_sim}}
#' @export
f2_from_msprime = function(..., blgsize = 0.05, cleanup = TRUE, verbose = TRUE) {

  ell = list(...)
  if(!'outpref' %in% names(ell)) ell$outpref = tempfile()
  if(verbose) alert_info('Simulating data...\n')

  do.call(msprime_sim, ell)

  if(verbose) alert_info('Reading data...\n')
  f2_blocks = f2_from_geno(ell$outpref, verbose = verbose, auto_only = FALSE, blgsize = blgsize)
  if(cleanup) file.remove(paste0(ell$outpref, c('.geno', '.ind', '.snp', '.py')))
  f2_blocks
}

# newick_to_edges = function(newwick) {
#
#   tree = ape::read.tree(text = newwick)
#   edges = tree$edge
#   edges = apply(edges, 2, function(x) paste0('l', x))
#   edges[tree$edge[,2] <= length(tree$tip.label), 2] = tree$tip.label
#   edges
# }

plot_ellout = function(evecfile, ellfile) {

  evec = read_table2(evecfile, col_names = FALSE, skip = 1) %>%
    transmute(id = X1, pop = X4, x = X2, y = X3)
  dd = read_lines(ellfile)
  samples = dd %>% str_subset('^sample:') %>% str_squish %>% word(2)
  num = dd %>% str_subset('^ell') %>% str_squish %>% str_split(' ') %>% do.call(rbind, .) %>%
    `[`(,-c(1,7)) %>% apply(2,as.numeric)
  elldat = as_tibble(num, .name_repair = ~c('x0', 'y0', 'a', 'b', 'angle', 'ci')) %>%
    mutate(id = samples, .before = 1)
  evec %>% left_join(elldat %>% select(-x0, -y0), by = 'id') %>%
    mutate(col = ifelse(is.na(a), pop, id)) %>%
    ggplot() + geom_point(aes(x, y, col = col, shape = pop), size=1) +
    ggforce::geom_ellipse(aes(x0=x, y0=y, a=a, b=b, angle=angle, fill = col, col = col), size=0.05, alpha = 0.05) +
    theme(panel.background = element_blank(), axis.line = element_line(), legend.title = element_blank()) + xlab('Eigenvector 1') + ylab('Eigenvector 2') + guides(col='none')

}

parse_treemix_treeout = function(treeout) {

  dat = read_lines(treeout)
  edges = newick_to_edges(dat[1])
  graph = graph_from_edgelist(edges)
  leaves = get_leafnames(graph)

  nammap = graph %>% igraph::distances(setdiff(names(V(graph)), leaves), leaves, mode = 'out') %>%
    apply(1, function(x) paste0(sort(names(which(is.finite(x)))), collapse=',')) %>%
    enframe %>% select(2:1) %>% deframe %>% c(purrr::set_names(leaves))

  newedges = str_split(dat[-1], ' ') %>%
    do.call(rbind, .) %>%
    `[`(,5:6) %>%
    as_tibble(.name_repair = ~c('from', 'to')) %>%
    mutate(i = 1:n()) %>%
    pivot_longer(c(from, to), names_to = 'k', values_to = 'v') %>%
    mutate(v = str_split(v, ','),
           v = map(v, ~str_replace_all(., '\\(|:.+', '')),
           v = map_chr(v, ~paste(sort(.), collapse = ','))) %>%
    pivot_wider(names_from = k, values_from = v) %>%
    transmute(from = nammap[from], to = nammap[to]) %>%
    mutate(inboth = to %in% .$from) %>%
    arrange(inboth, from, to)

  graph %>% insert_admix(source_to = newedges$from, dest_to = newedges$to, substitute = FALSE)
}



parse_treemix = function(stem, split = FALSE) {

  vert = read_table2(paste0(stem, '.vertices.gz'), col_names = FALSE, col_types = cols(.default = 'c')) %>%
    transmute(num = X1, nam = ifelse(is.na(X2), paste0('n',num), X2)) %>% deframe
  edges = read_table2(paste0(stem, '.edges.gz'), col_names = FALSE, col_types = cols(.default = 'c')) %>%
    transmute(from = vert[X1], to = vert[X2])
  if(split) edges %<>% as.matrix %>% split_multifurcations %>% select(1:2)
  edges %>% as.matrix %>% graph_from_edgelist()
}


f2_from_fastsimcoal = function(graph, nblocks = 100, verbose = TRUE, ...) {

  #sfs = map(seq_len(nblocks), ~fastsimcoal_sim_dsfs(graph, ...), ...)
  sfs = list()
  run = list(...)$run
  if(is.null(run)) stop("Need to pass fastsimcoal path to 'run'!")
  if(!file.exists(run)) stop(paste0(run, ' not found!'))
  for(i in seq_len(nblocks)) {
    sfs[[i]] = fastsimcoal_sim_dsfs(graph, ..., verbose = FALSE)
  }
  if(verbose) alert_info('fastsimcoal sfs done. Computing f2 from sfs...\n')
  sfs %>% map(sfs_to_f2) %>% bind_rows(.id = 'block') %>% f2dat_to_f2blocks()
}

f2dat_to_f2blocks = function(f2dat, fill_diag = TRUE) {
  # f2dat has columns pop1, pop2, f2, block

  pops = sort(union(f2dat$pop1, f2dat$pop2))
  nblocks = length(unique(f2dat$block))
  if('length' %in% names(f2dat)) {
    bl = f2dat %>% filter(pop1 == f2dat$pop1[1], pop2 == f2dat$pop2[1]) %>% pull(length) %>% paste0('l', .)
  } else {
    rep('l1', nblocks)
  }
  out = f2dat %>% select(pop1, pop2, f2, block) %>% mutate(block = as.numeric(block)) %>%
    bind_rows(rename(., pop1 = pop2, pop2 = pop1)) %>% distinct
  if(fill_diag) {
    out %<>% bind_rows(expand_grid(block = seq_len(nblocks), pop1 = pops, f2 = 0) %>% mutate(pop2 = pop1))
  }
  out %>%
    arrange(block, pop1, pop2) %$%
    array(f2, c(length(pops), length(pops), nblocks), list(pops, pops, bl))
}

f2dat_to_f2blocks2 = function(f2dat, nblocks = 1000, cv = 0.1) {

  pops = sort(union(f2dat$pop1, f2dat$pop2))
  f2dat %>% select(pop1, pop2, f2) %>% filter(pop1 < pop2) %>% expand_grid(block = seq_len(nblocks)) %>%
    rowwise %>% mutate(f2 = rnorm(1, f2, f2*cv)) %>% ungroup %>%
    bind_rows(rename(., pop1 = pop2, pop2 = pop1)) %>%
    bind_rows(expand_grid(block = seq_len(nblocks), pop1 = pops, f2 = 0) %>% mutate(pop2 = pop1)) %>%
    arrange(block, pop1, pop2) %$%
    array(f2, c(length(pops), length(pops), nblocks), list(pops, pops, rep('l1', nblocks)))
}

write_fastsimcoal_obs = function(afs, outfile = stdout()) {

  popcounts = apply(afs$counts, 2, max)
  out = "1 observations. No. of demes and sample sizes are on next line\n"
  out %<>% paste0(length(popcounts), '\t', paste(popcounts, collapse = '\t'), '\n')
  out %<>% paste0(paste(admixtools:::joint_sfs(afs) %>% pull(n), collapse = '\t'))

  writeLines(out, outfile)
}


write_fastsimcoal_files = function(graph, outpref, tpl = FALSE, nsnps = 10000000,
                                   recombination_rate = 0.00000001, mutation_rate = 0.00000002, time = 1000) {

  pops = get_leafnames(graph)
  nodes = names(V(graph))
  npop = length(pops)
  nnode = length(nodes)
  dates = if(length(time) > 1) time else pseudo_dates(graph, time) %>% enframe('node', 'date')
  edges = graph %>% as_edgelist() %>% as_tibble(.name_repair = ~c('from', 'to')) %>% add_count(to) %>% mutate(type = ifelse(n > 1, 'admix', 'normal')) %>% select(-n) %>% left_join(dates %>% transmute(from = node, fromdate = date)) %>% left_join(dates %>% transmute(to = node, todate = date)) %>% suppressMessages() %>%
    mutate(weight = ifelse(type == 'admix' & !duplicated(to), 0.5, 1)) %>%
    mutate(time = ifelse(type == 'admix', fromdate, fromdate), source = match(to, nodes)-1, sink = match(from, nodes)-1, migrants = weight, size = 1, gr = 0, migmat = 0,
           gr = ifelse(type == 'admix', 0, paste0('G_', 1:nrow(.), '$')),
           epar = ifelse(type == 'admix', paste0('T_', from, '$'), paste0('T_', from, '$')),
           epar2 = paste0('T_', to, '$'),
           rule = ifelse(to %in% pops, '', paste0(epar2, ' <= ', epar))) %>% arrange(time)
  nedge = nrow(edges)


  if(tpl) {
    nadm = edges %>% filter(migrants < 1) %>% nrow
    apar = paste0('A', seq_len(nadm), '$')[seq_len(nadm)]
    edges$migrants[edges$migrants < 1] = apar
    #edges$time = edges$epar
    edges$time = edges$todate

    est = "// Priors and rules file
  // *********************

  [PARAMETERS]
  //#isInt? #name   #dist.#min  #max
  //all Ns are in number of haploid individuals\n"

    apars = paste0('0\t', apar, '\tunif\t0\t1\toutput', collapse = '\n')
    #epars = paste0('1\t', unique(edges$epar), '\tlogunif\t10\t100000\toutput', collapse = '\n')
    gpars = paste0('0\t', unique(edges$gr)[1], '\tunif\t-1\t1\toutput', collapse = '\n')

    est %<>% paste0(gpars, '\n')
    if(nadm > 0) est %<>% paste0(apars, '\n')
    est %<>% paste0("\n[RULES]\n")
    #est %<>% paste0(paste0(unique(edges$rule), collapse = '\n'))
    est %<>% paste0("\n\n[COMPLEX PARAMETERS]\n")

    writeLines(est, paste0(outpref, '.est'))

  }

  events = edges %$% paste(time, source, sink, migrants, size, 0, migmat) %>% paste(collapse = '\n')

  out = "//Number of population samples (demes)\n"
  out %<>% paste0(nnode)
  out %<>% paste0('\n//Population effective sizes (number of genes)\n')
  out %<>% paste0(paste0(rep(10000, nnode), collapse = '\n'))
  out %<>% paste0('\n//Sample sizes\n')
  out %<>% paste0(paste0(ifelse(nodes %in% pops, 2, 0), collapse = '\n'))
  out %<>% paste0('\n//Growth rates: negative growth implies population expansion\n')
  out %<>% paste0(paste0(rep(0, nnode), collapse = '\n'))
  out %<>% paste0('\n//Number of migration matrices : 0 implies no migration between demes\n')
  out %<>% paste0(0)
  out %<>% paste0('\n//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 4 historical event\n')
  out %<>% paste0(nedge, ' historical event\n')
  out %<>% paste0(events)
  out %<>% paste0('\n//Number of independent loci [chromosome]\n')
  out %<>% paste0('1 0')
  out %<>% paste0('\n//Per chromosome: Number of linkage blocks\n')
  out %<>% paste0('1')
  out %<>% paste0('\n//per Block: data type, num loci, rec. rate and mut rate + optional parameters\n')
  if(tpl) out %<>% paste0('FREQ ', nsnps,' ', recombination_rate,' ', mutation_rate,'\n')
  else out %<>% paste0('SNP ', nsnps,' ', recombination_rate,' ', mutation_rate,'\n')
  #out %<>% paste0('DNA 10000000 0.00000001 0.00000002 0.33\n')
  #out %<>% paste0('STANDARD 10000000 0.00000001 0.00000002\n')

  writeLines(out, paste0(outpref, if(tpl) '.tpl' else '.par'))
}



parse_fastsimcoal_dsfs = function(obs_file, popcounts) {
  # reads fastsimcoal DSFS.obs file
  # popcounts is vector of haplotype counts for each population

  cnt = obs_file %>%
    read_lines %>%
    pluck(3) %>%
    str_split('\t') %>%
    pluck(1) %>%
    as.numeric
  expand_grid(!!!popcounts %>% map(~0:.)) %>% mutate(n = cnt)
}

fastsimcoal_sim_dsfs = function(graph, outpref, run = './fsc26', popcounts = 1, nsnps = 10000000,
                                recombination_rate = 0.00000001, mutation_rate = 0.00000002, verbose = TRUE) {

  write_fastsimcoal_files(graph, outpref, tpl = FALSE, nsnps = nsnps, recombination_rate = recombination_rate,
                          mutation_rate = mutation_rate)
  if(file.exists(run)) {
    obs_file = paste0(outpref, '/', basename(outpref), '_DSFS.obs')
    if(file.exists(obs_file)) file.remove(obs_file)
    else dir.create(outpref, showWarnings = FALSE)
    cmd = paste0(run, ' -i ', outpref,'.par -n1 -d -u -s0 -k 100000000 > ', outpref, '.txt')
    if(verbose) alert_info(paste0('Running "', cmd, '"...\n'))
    olddir = getwd()
    system(paste0('cd ', dirname(outpref), '; ', cmd, '; cd ', olddir))
    pops = get_leafnames(graph)
    if(length(popcounts) == 1) popcounts = rep(popcounts*2, length(pops)) %>% set_names(pops)
    dsfs = parse_fastsimcoal_dsfs(obs_file, popcounts)
    return(dsfs)
  } else {
    if(verbose) alert_info(paste0('file ', run, ' not found. Files written to ', outpref, '\n'))
  }
}




