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
    if(!is.null(lambdascale)) parfile %<>% paste0('lambdascale: ', lambdascale, '\n')
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
    edg = as_tibble(graph) %>% set_colnames(c('from', 'to'))
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
    set_colnames(c('from', 'to')) %>%
    as_tibble %>%
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
parse_qpgraph_output = function(outfile) {
  # reads qpGraph output file
  # returns list of three objects:
  # 'edges': data.frame of branch lengths and admixture weights
  # 'score': best fit score
  # 'f2': data.frame of estimated and fitted f2 values

  # addition of outliers currently not safe when two population have the same three-letter-prefix

  edges = parse_qpgraph_output_edges(outfile)

  dat = readLines(outfile) %>% str_squish() %>% enframe(value = 'X1') %>% select(X1)

  score = (dat %>% filter(grepl('^final score', X1)) %>%
             separate('X1', c('a', 'b', 'score'), sep=' +', convert = T, extra='drop', fill='right'))$score

  pops = dat$X1 %>% str_subset('^population:') %>% str_squish %>% word(3)
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
    f3 = NULL
  }

  outlierstart = str_which(dat$X1, '^outliers:')[1]+2
  outlierend = str_which(dat$X1, '^worst f-stat:')[1]-3

  numalloutliers = choose(choose(numpop, 2)+1, 2)
  if(outlierend - outlierstart + 1 != numalloutliers) warning('Outliers are missing')

  if(is.na(pops)[1]) pops = dat$X1 %>% str_subset('^zzaddw') %>% str_split(' ') %>% map(2) %>% unlist %>% unique
  denom = 1000
  amb = names(which(table(str_sub(pops, 1, 3)) > 1))
  if(length(amb) > 0) warning(paste('Ambiguous populations ommited from outliers: ', paste0(amb, collapse = ', ')))

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

  #f2 %<>% mutate(pop1 = rep(head(pops, -1), (numpop-1):1), pop2 = unlist(map(2:numpop, ~pops[.:numpop])))
  #f3 %<>% mutate(pop2 = rep(pops[-1], each = numpop-1), pop3 = rep(pops[-1], numpop-1))

  namedList(edges, score, f2, f3, f4 = outliers)
}

parse_qpgraph_output_edges = function(outfile) {

  dat = readLines(outfile) %>% str_squish() %>% enframe(value = 'X1') %>% select(X1)

  edges = dat %>%
    filter(grepl('^ledge|^redge|^admix', X1)) %>%
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

}



#' Read qpGraph graph file
#' @export
#' @param graphfile File with admixture graph in qpGraph format.
#' @param split_multi Split multifurcations
#' @param igraph Convert to igraph format
#' @return Graph represented as two column edge matrix. Can have four columns if edges are locked
parse_qpgraph_graphfile = function(graphfile, split_multi = TRUE, igraph = FALSE) {
  # reads graph in qpGraph format
  # returns edge matrix (adjacency list)
  lines = read_lines(graphfile) %>% enframe %>% transmute(V1 = value)
  namemap = lines %>% filter(grepl('^label', V1)) %>%
    separate('V1', c('type', 'label', 'name'), sep = '\\s+', extra = 'drop') %>%
    select(-type) %>% deframe

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
    mutate(from = recode(from, !!!namemap),
           to = recode(to, !!!namemap)) %>%
    select(from, to, lower, upper)
  if(all(is.na(out$lower)) && all(is.na(out$upper))) out %<>% select(-lower, -upper)
  out %<>% as.matrix
  if(split_multi) out %<>% split_multifurcations
  if(igraph) out = graph_from_edgelist(as.matrix(out[,1:2]))
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


graph_to_lineages = function(graph, node = V(graph)[1], num = 1) {
  # adds vertex attribute 'lineage' to graph, for 'msprime_sim'

  lin = vertex_attr(graph, 'lineage', node)
  if(!is.null(lin) && !is.na(lin)) return(graph)
  graph = set_vertex_attr(graph, 'lineage', node, num)
  children = neighbors(graph, node, mode = 'out')
  if(length(children) == 0) {
    return(graph)
  } else if(length(children) == 1) {
    if(length(neighbors(graph, children, mode = 'in')) == 1) {
      graph = graph_to_lineages(graph, children[1], num)
    } else {
      maxlin = max(vertex_attr(graph, 'lineage'), na.rm = TRUE)
      graph = graph_to_lineages(graph, children[1], maxlin + 1)
    }
  } else if(length(children) == 2) {
    graph = graph_to_lineages(graph, children[1], num)
    maxlin = max(vertex_attr(graph, 'lineage'), na.rm = TRUE)
    graph = graph_to_lineages(graph, children[2], maxlin + 1)
  }
  graph
}

#' Generate msprime simulation code from admixture graph
#' @export
#' @param edges Graph as edge list with columns 'from' and 'to' (optionally 'weight')
#' @return msprime python code backbone
#' @examples
#' \dontrun{
#' results = qpgraph(example_f2_blocks, example_graph)
#' msprime_sim(results$edges)
#' }
msprime_sim = function(edges, outfilename = 'msprime_sim.py',
                       vcffilename = 'msprimesim.vcf', nind = 1, popsize = 1000, len = 1e9,
                       recombination_rate = 1e-9, mutation_rate = 1e-9, time = 100) {

  if(is.matrix(edges) && ncol(edges) == 2) edges %<>% as_tibble(.name_repair = ~c('from', 'to'))
  pdat = graph_to_plotdat(edges)$eg
  dates = bind_rows(transmute(pdat, name, y), transmute(pdat, name = to, y = yend)) %>% distinct %>% deframe

  edges %<>% add_count(to) %>%
    mutate(from = str_replace_all(from, '\\.', ''), to = str_replace_all(to, '\\.', ''),
           type = ifelse(n > 1, 'admix', 'normal')) %>% select(-n)
  if(!'weight' %in% names(edges)) edges %<>% mutate(weight = 1)
  graph = igraph::graph_from_edgelist(as.matrix(edges[,1:2]))
  leaves = get_leafnames(graph)
  nodes = names(V(graph))[-1]
  lineages = vertex_attr(graph_to_lineages(graph), 'lineage')
  edges2 = edges %>%
    mutate(nfrom = match(from, names(V(graph))),
           nto = match(to, names(V(graph))),
           lfrom = lineages[nfrom],
           lto = lineages[nto],
           time = nfrom) %>%
    filter(lfrom != lto) %>%
    group_by(to) %>%
    mutate(weight2 = ifelse(n() == 2 & time == max(time), weight, 1)) %>%
    ungroup

  out = "import math
import msprime
import argparse

parser = argparse.ArgumentParser(description = 'Simulating pseudo sequences by chromosomes')
parser.add_argument('-c', '--CHR', type = int, default = 1, help = 'Chromosomes number')
parser.add_argument('-m', '--migration_percent', type = int, default = 5, help = 'Gene flow percent')
parser.add_argument('-i', '--iteration', type = int, default = 1, help ='Iteration and seed')
args = parser.parse_args()
"

  out = paste0(out, '\n#eff. pop. sizes\n', paste0(nodes, ' = ',popsize,'\n', collapse = ''))
  out = paste0(out, '\ngen = 29\n')
  out = paste0(out, 'pops = [\n', paste0('\tmsprime.PopulationConfiguration(initial_size = ', nodes,
                                        ')', c(rep(',', length(nodes)-1), ''),
                                        ' #',(1:length(nodes))-1,'\n', collapse = '') ,']\n')
  out = paste0(out, '\n#ind. dates\nsamples = [\n', paste0('\tmsprime.Sample(', (1:length(nodes))-1,
                                                           ', ',dates[nodes]*time,'/gen)',
                                                           c(rep(',', length(nodes)-1), ' '),
                                                           ' #', nodes, '\n', collapse = ''), ']\n')
  out = paste0(out, '\n#pop. split dates\n', paste0('T_', edges2$from, ' = ',dates[edges2$from]*time,'\n', collapse = ''))
  out = paste0(out, '\nevents = [\n',
               paste0('\tmsprime.MassMigration(time = T_',
                      edges2$from, '/gen, source = ', edges2$lto-1,', destination = ', edges2$lfrom-1,
                      ', proportion = ', ifelse(edges2$type == 'admix', edges2$weight2, 1),')',
                      collapse = ',\n'), '\n]\n')

  out = paste0(out, paste0("\ntree_sequence = msprime.simulate(population_configurations = pops,
                                 length = ",len,",
                                 samples = samples,
                                 demographic_events = events,
                                 recombination_rate = ", recombination_rate,",
                                 mutation_rate = ", mutation_rate,",
                                 random_seed = args.iteration)

with open('",vcffilename,"', 'w') as f:
  tree_sequence.write_vcf(f, ploidy = 2)"))

  writeLines(out, outfilename)
}


