#' @import dplyr
#' @import readr
#' @import purrr
#' @import tidyr
#' @import stringr
#' @import ggplot2
#' @import utils
#' @importFrom magrittr set_colnames "%$%" "%<>%"
#' @importFrom rray rray rray_transpose
#' @importFrom lobstr obj_size
#' @importFrom abind abind
#' @importFrom crayon blue red green bold italic
#' @importFrom tibble as_tibble deframe enframe add_column
#' @importFrom stats na.omit setNames runif
#' @importFrom grDevices hcl
#' @importFrom rlang .data
#' @importFrom plotly plot_ly add_trace add_markers layout
#' @importFrom Rcpp cppFunction
#' @importFrom igraph V E neighbors subcomponent get.edge.ids degree incident_edges all_simple_paths graph_from_edgelist as_edgelist
#' @useDynLib admixtools

TRUE

#' Read allele frequencies from packedancestrymap genotype file
#' @export
#' @param pref prefix of packedancestrymap files (has to end in .geno, .ind, .snp).
#' @param pops vector of populations from which to compute allele frequencies.
#' @param inds vector of samples from which to compute allele frequencies. \code{pops} or \code{inds} have to be specified.
#' @param blocksize number of SNPs read in each block.
#' @param na.action what to do with missing SNPs. \code{none} output table will have missing data. \code{impute} allele frequencies will be imputed from other populations. \code{remove} SNPs with missing data will be removed.
#' @return a tibble with allele frequency data. first six columns are snpfile, remaining columns are allele frequencies for reach population.
packedancestrymap_to_aftable = function(pref, pops=NULL, inds=NULL, blocksize=1000, na.action='none', return_matrix=FALSE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # pops is vector of populations for which to calculate AFs
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population
  # na.action should be 'none' (default), 'mean impute', or 'remove SNP'

  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = cols())
  snpfile = read_table2(paste0(pref, '.snp'), col_names = FALSE, col_types = cols()) %>% set_colnames(c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2'))
  if(!is.null(inds)) {
    indfile$X3 = indfile$X1
    indfile$X3[!indfile$X3 %in% inds] = NA
    pops = inds
  }
  if(is.null(pops) & is.null(inds)) pops = unique(indfile$X3)
  pops = unique(na.omit(pops))
  stopifnot(all(pops %in% indfile$X3))
  fl = paste0(pref, '.geno')
  conn = file(fl, 'rb')
  on.exit(close(conn))
  hd = strsplit(readBin(conn, 'character', n = 1), ' +')[[1]]
  close(conn)
  nind = as.numeric(hd[2])
  nsnp = as.numeric(hd[3])
  numpop = length(pops)
  popind = which(indfile$X3 %in% pops)
  popind2 = na.omit(popind)

  alert_info(paste0(basename(pref), '.geno has ', nind, ' samples and ', nsnp, ' SNPs.\n'))
  alert_info(paste0('Calculating allele frequencies from ', length(popind2), ' samples in ', numpop, ' populations.\n'))
  alert_info(paste0('Expected size of allele frequency data: ', round((nsnp*numpop*8+nsnp*112)/1e6), ' MB\n'))
  # 8, 112: estimated scaling factors for AF columns and annotation columns

  rlen = file.info(fl)$size/(nsnp+1)
  conn = file(fl, 'rb')
  invisible(readBin(conn, 'raw', n = rlen))
  afmatrix = matrix(NA, nsnp, numpop)
  colnames(afmatrix) = pops
  popind3 = c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`))
  popindmat = matrix(sapply(1:numpop, function(i) indfile$X3[popind2] == pops[i]), ncol=numpop)
  cnt = 1
  sumna = 0
  while(cnt <= nsnp) {
    if(cnt+blocksize > nsnp) {
      blocksize = nsnp-cnt+1
      popind3 = sort(c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`)))
    }
    bitmat = matrix(as.integer(rawToBits(readBin(conn, 'raw', n = rlen*blocksize))), ncol=8, byrow = TRUE)
    gmat = matrix(c(t(bitmat[,c(8,6,4,2)]*2+bitmat[,c(7,5,3,1)]))[popind3], ncol=blocksize)
    gmat[gmat==3]=NA # assuming non-missing genotypes are 0, 1, 2 missing is 3
    popfreqs = sapply(1:numpop, function(i) colMeans(gmat[popindmat[,i],, drop=FALSE], na.rm=TRUE)/2)
    popfreqs[is.nan(popfreqs)] = NA
    sumna = sumna + sum(is.na(popfreqs))
    afmatrix[cnt:(cnt+blocksize-1),] = popfreqs
    cat(paste0('\r', (cnt-1)/1e3, 'k SNPs read...'))
    cnt = cnt+blocksize
  }
  cat('\n')
  alert_success(paste0(cnt-1, ' SNPs read in total\n'))
  alert_warning(paste0(sumna, ' allele frequencies are missing (on average ', round(sumna/numpop), ' per population)\n'))

  if(na.action == 'impute') {
    afmatrix = mean_impute(afmatrix, by=1)
    isnan = is.nan(afmatrix)
    afmatrix[isnan] = NA
    allna = apply(afmatrix, 1, function(x) all(is.na(x)))
    afmatrix = afmatrix[!allna,]
    snpfile = snpfile[!allna,]
    alert_danger(paste0(sumna-sum(isnan), ' allele frequencies were imputed, (on average ', round((sumna-sum(isnan))/numpop),' per population) ', sum(allna), ' SNPs were removed\n'))
  }
  if(na.action == 'remove') {
    anyna = apply(afmatrix, 1, function(x) any(is.na(x)))
    afmatrix = afmatrix[!anyna,]
    snpfile = snpfile[!anyna,]
    alert_danger(paste0(sum(anyna), ' SNPs were removed, ', sum(!anyna), ' SNPs remain\n'))
  }
  if(return_matrix) return(afmatrix)
  snpfile %>% bind_cols(as_tibble(afmatrix))
}

#' Reads qpGraph output file
#' @export
#' @param outfile output file generated by qpGraph.
#' @return list of output data.
parse_qpgraph_output = function(outfile) {
  # reads qpGraph output file
  # returns list of three objects:
  # 'edges': data.frame of branch lengths and admixture weights
  # 'score': best fit score
  # 'f2': data.frame of estimated and fitted f2 values

  dat = read_table(outfile, col_names=F, col_types = cols(), guess_max = 1e6)

  edges = dat %>%
    filter(grepl('^ledge|^redge|^admix', .data$X1)) %>%
    separate('X1', c('type', 'name', 'from', 'to', 'weight', 'w2'),
             sep=' +', convert = T, extra='drop', fill='right')
  admix1 = edges %>% filter(.data$type=='admix') %>%
    mutate(type='aedge', to=.data$name, name='', w2=NA)
  admix2 = edges %>% filter(.data$type=='admix') %>%
    mutate(type='aedge', from=.data$to, to=.data$name, name='', weight=.data$w2, w2=NA)
  edges %<>%
    bind_rows(admix1) %>%
    bind_rows(admix2) %>%
    filter(!.data$type == 'admix') %>%
    mutate(type = ifelse(.data$type=='aedge', 'admix', 'edge')) %>%
    select(-'w2', -'name')

  score = (dat %>% filter(grepl('^final score', .data$X1)) %>%
             separate('X1', c('a', 'b', 'score'), sep=' +', convert = T, extra='drop', fill='right'))$score

  f2 = dat %>% filter(grepl(' f2: ', .data$X1)) %>%
    separate('X1', c('pop1', 'pop2', 'fst','f2fit','f2est','diff','se','z'), sep=' +', convert = TRUE) %>%
    select(-.data$fst)

  namedList(edges, score, f2)
}

#' Reads graph in qpGraph format
#' @export
#' @param graphfile file with admixture graph in qpGraph format.
#' @return graph represented as two column edge matrix.
parse_qpgraph_graphfile = function(graphfile) {
  # reads graph in qpGraph format
  # returns edge matrix (adjacency list)

  dat = read_lines(graphfile) %>%
    tibble %>% set_colnames('V1') %>%
    filter(grepl('^edge|admix', .data$V1)) %>%
    separate('V1', c('type', 'name', 'from', 'to'), sep = '\\s+')
  admix1 = dat %>% filter(.data$type=='admix') %>% mutate(type='edge', to=.data$name, name='')
  admix2 = dat %>% filter(.data$type=='admix') %>% mutate(type='edge', from=.data$to, to=.data$name, name='')
  dat %>%
    bind_rows(admix1) %>%
    bind_rows(admix2) %>%
    filter(.data$type != 'admix') %>%
    select(.data$from, .data$to) %>% as.matrix
}

#' Reads qpGraph parameter file
#' @export
#' @param parfile parameter file for qpGraph.
#' @return named list of qpGraph parameters.
parse_qpgraph_parfile = function(parfile) {
  # reads qpGraph parfile
  # returns named list of parameters
  # all genotype files have to have same prefix

  dat = read_table2(parfile, comment = '#', col_names = c('par', 'value'), col_types = cols())
  dat %<>% mutate(par = str_replace_all(.data$par, c(':$'='', 'genotypename'='pref')),
           value = str_replace_all(.data$value,
                                   c('\\.geno$'='',
                                     'S1'=filter(dat, .data$par=='S1:')$value[1],
                                     'DIR'=filter(dat, .data$par=='DIR:')$value[1])),
           value=ifelse(.data$value=='YES', TRUE, ifelse(.data$value=='NO', FALSE, .data$value))) %>%
    filter(!.data$par %in% c('DIR', 'S1', 'indivname', 'snpname')) %>%
    t %>%
    as_tibble()
  dat %>%
    set_colnames(slice(dat, 1)) %>%
    slice(-1) %>%
    type_convert(col_types = cols()) %>%
    as.list
}



#' Reads qpAdm output file
#' @export
#' @param outfile output file generated by qpAdm.
#' @return tibble with output data.
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
  coefs %<>% str_split(' +') %>% map(~tail(., length(left)+1) %>% head(-1) %>% as.numeric %>% set_names(left)) %>%
    set_names(c('weight', 'jmean', 'jse')) %>% as_tibble
  weights = tibble(target, left) %>% bind_cols(coefs)

  sig = do.call(rbind, str_split(dat[sigstart:sigend], ' +')) %>% as.data.frame(stringsAsFactors=F) %>% select(-1) %>% set_colnames(c('pattern', 'wt', 'dof', 'chisq', 'tail prob', left, 'feasible')) %>% mutate(feasible = feasible != 'infeasible')

  namedList(weights, sig)
}

#' Write f2 block jackknife estimates to disk
#'
#' This function takes a 3d array of f2 block jackknife estimates, splits it by population pair, and writes each pair to a separate \code{.RData} file.
#' @export
#' @param f2_block 3d array with f2 block jackknife estimates. The first two dimensions of \code{f2_block} have to have population names. Each file will be stored under \code{{path}/{pop1}/{pop2}.RData}.
#' @param path directory into which to write the files.
#' @param overwrite should existing \code{.RData} files be overwritten?
#' @return NULL
write_f2 = function(f2_block, path, overwrite=FALSE) {

  if(!dir.exists(path)) dir.create(path)
  d1 = dim(f2_block)[1]
  d2 = dim(f2_block)[2]
  nam1 = dimnames(f2_block)[[1]]
  nam2 = dimnames(f2_block)[[2]]
  for(i in seq_len(d1)) {
    for(j in seq_len(d2)) {
      pop1 = min(nam1[i], nam2[j])
      pop2 = max(nam1[i], nam2[j])
      f2 = as.vector(f2_block[i, j, ])
      dir = paste0(path, '/', pop1, '/')
      fl = paste0(dir, pop2, '.RData')
      if(!dir.exists(dir)) dir.create(dir)
      if(!file.exists(fl) | overwrite) save(f2, pop1, pop2, file=fl)
    }
  }
}

#' @export
read_f2 = function(path, pops = NULL) {
  # reads f2 block jackknife .RData files and returns 3d array
  # pops is vector of populations which should be read. defaults to all populations.
  if(is.null(pops)) pops = list.dirs(path, full.names=FALSE, recursive=FALSE)
  load(paste0(path, '/', pops[1], '/', pops[1], '.RData'))
  numblocks = length(f2)
  numpops = length(pops)
  f2_blocks = array(0, c(numpops, numpops, numblocks))
  dimnames(f2_blocks)[[1]] = dimnames(f2_blocks)[[2]] = pops
  for(pop1 in pops) {
    for(pop2 in pops) {
      if(pop1 <= pop2) {
        fl = paste0(path, '/', pop1, '/', pop2, '.RData')
        load(fl)
        if(any(is.na(f2))) warn(paste0('missing values in ', pop1, ' - ', pop2, '!'))
        f2_blocks[pop1, pop2, ] = f2_blocks[pop2, pop1, ] = f2
      }
    }
  }
  f2_blocks
}


split_afmat = function(afmat, pops_per_block, outprefix) {
  # splits afmat into numparts parts, and saves them to {outprefix}_{part}.RData
  npops = ncol(afmat)
  starts = seq(1, npops, pops_per_block)
  numparts = length(starts)
  ends = c(lead(starts)[-numparts]-1, npops)
  for(i in seq_len(numparts)) {
    cat(paste0('\rpart ', i, ' of ', numparts))
    afs = afmat[, starts[i]:ends[i]]
    save(afs, file = paste0(outprefix, '_', i, '.RData'))
  }
  cat('\n')
}

