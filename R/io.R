
#' Read allele frequencies from packedancestrymap files
#' @export
#' @param pref prefix of packedancestrymap files (files have to end in \code{.geno}, \code{.ind}, \code{.snp}).
#' @param pops vector of populations from which to compute allele frequencies.
#' @param inds vector of samples from which to compute allele frequencies. \code{pops} and \code{inds} cannot be specified at the same time. If none are specified, populations will be extracted from the third column in the \code{.ind} file
#' @param blocksize number of SNPs read in each block.
#' @param na.action what to do with SNPs with missing data.
#' \itemize{
#' \item \code{none} (default) output table will have missing data
#' \item \code{impute} allele frequencies will be imputed from other populations
#' \item \code{remove} SNPs with missing data will be removed
#' }
#' @param return_matrix return data as list of matrices
#' @return a list with two items: allele frequency data and allele counts. Unless \code{return_matrix} is specified, first six columns are bimfile, remaining columns are data for each population.
#' @examples
#' \dontrun{
#' afdat = packedancestrymap_to_aftable(prefix, pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
packedancestrymap_to_aftable = function(pref, pops=NULL, inds=NULL, blocksize=1000, na.action='none', return_matrix = FALSE, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # pops is vector of populations for which to calculate AFs
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population
  # na.action should be 'none' (default), 'mean impute', or 'remove SNP'

  stopifnot(is.null(pops) | is.null(inds))
  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = cols(), progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = cols(), progress = FALSE)
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

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nind, ' samples and ', nsnp, ' SNPs.\n'))
    alert_info(paste0('Calculating allele frequencies from ', length(popind2), ' samples in ', numpop, ' populations.\n'))
    alert_info(paste0('Expected size of allele frequency data: ', round((nsnp*numpop*8+nsnp*112)/1e6), ' MB\n'))
    # 8, 112: estimated scaling factors for AF columns and annotation columns
  }

  rlen = file.info(fl)$size/(nsnp+1)
  conn = file(fl, 'rb')
  invisible(readBin(conn, 'raw', n = rlen))
  afmatrix = countmatrix = matrix(NA, nsnp, numpop)
  #countarr = array(NA, dim = c(3, nsnps, numpop))
  colnames(afmatrix) = colnames(countmatrix) = pops
  rownames(afmatrix) = rownames(countmatrix) = snpfile$SNP
  popind3 = c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`))
  popindmat = matrix(sapply(1:numpop, function(i) indfile$X3[popind2] == pops[i]), ncol=numpop)
  cnt = 1
  sumna = 0
  #count_nonmissing = function(x) sum(!is.na(x))
  count_gt = function(x) {x=na.omit(x); c(sum(x==0), sum(x==1), sum(x==2))}
  while(cnt <= nsnp) {
    if(cnt+blocksize > nsnp) {
      blocksize = nsnp-cnt+1
      popind3 = sort(c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`)))
    }
    bitmat = matrix(as.integer(rawToBits(readBin(conn, 'raw', n = rlen*blocksize))), ncol=8, byrow = TRUE)
    gmat = matrix(c(t(bitmat[,c(8,6,4,2)]*2+bitmat[,c(7,5,3,1)]))[popind3], ncol=blocksize)
    gmat[gmat==3]=NA # assuming non-missing genotypes are 0, 1, 2 missing is 3
    if(cnt == 1) {
      ploidy = apply(gmat, 1, function(x) length(unique(na.omit(x)))-1)
      alert_info(paste0('Detected ', sum(ploidy == 2), ' diploid samples and ',
                        sum(ploidy == 1), ' pseudohaploid samples.\n'))
    }
    popfreqs = popcounts = matrix(NA, blocksize, numpop)
    for(i in seq_len(numpop)) {
      gm = gmat[popindmat[,i],, drop=FALSE]
      popfreqs[,i] = colMeans(gm, na.rm=TRUE)/2
      popcounts[,i] = colSums((!is.na(gm))*ploidy[popindmat[,i]])
      #popcounts[,i] = colSums((!is.na(gm)))
    }
    popfreqs[is.nan(popfreqs)] = NA
    sumna = sumna + sum(is.na(popfreqs))
    afmatrix[cnt:(cnt+blocksize-1),] = popfreqs
    countmatrix[cnt:(cnt+blocksize-1),] = popcounts
    if(verbose) cat(paste0('\r', (cnt-1)/1e3, 'k SNPs read...'))
    cnt = cnt+blocksize
  }
  if(verbose) {
    cat('\n')
    alert_success(paste0(cnt-1, ' SNPs read in total\n'))
    alert_warning(paste0(sumna, ' allele frequencies are missing (on average ', round(sumna/numpop), ' per population)\n'))
  }
  keepsnps = snpfile %>% mutate(i = 1:n()) %>% filter(as.numeric(gsub('[a-zA-Z]+', '', CHR)) %in% 1:22) %$% i
  outdat = treat_missing(afmatrix[keepsnps,], countmatrix[keepsnps,], snpfile[keepsnps,], na.action = na.action, verbose = verbose)

  outlist = list(afs = outdat$afmatrix, counts = outdat$countmatrix)
  if(!return_matrix) outlist = map(outlist, ~bind_cols(outdat$snpfile, as_tibble(.)))
  outlist$snpfile = outdat$snpfile
  outlist
}

#' Read genotype data from packedancestrymap files
#' @export
#' @param pref prefix of packedancestrymap files (files have to end in \code{.geno}, \code{.ind}, \code{.snp}).
#' @param inds optional vector of samples to read in
#' @param blocksize number of SNPs read in each block.
#' @param na.action what to do with SNPs with missing data.
#' \itemize{
#' \item \code{none} (default) output table will have missing data
#' \item \code{impute} allele frequencies will be imputed from other populations
#' \item \code{remove} SNPs with missing data will be removed
#' }
#' @return a list with the genotype data matrix, the \code{.ind} file, and the \code{.snp} file
#' @examples
#' \dontrun{
#' samples = c('Ind1', 'Ind2', 'Ind3')
#' geno = read_packedancestrymap(prefix, samples)
#' }
read_packedancestrymap = function(pref, inds=NULL, blocksize = 1000, na.action = 'none',
                                  auto_only = TRUE, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population
  # na.action should be 'none' (default), 'mean impute', or 'remove SNP'

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = cols(), progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = cols(), progress = FALSE)
  indfile$X3 = indfile$X1
  if(!is.null(inds)) {
    indfile$X3[!indfile$X3 %in% inds] = NA
  } else {
    inds = indfile$X1
  }
  fl = paste0(pref, '.geno')
  conn = file(fl, 'rb')
  on.exit(close(conn))
  hd = strsplit(readBin(conn, 'character', n = 1), ' +')[[1]]
  close(conn)
  nind = as.numeric(hd[2])
  nsnp = as.numeric(hd[3])

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nind, ' samples and ', nsnp, ' SNPs.\n'))
    alert_info(paste0('Expected size of genotype data: ', round((nsnp*nind*8+nsnp*112)/1e6), ' MB\n'))
    # 8, 112: estimated scaling factors for AF columns and annotation columns
  }

  rlen = file.info(fl)$size/(nsnp+1)
  conn = file(fl, 'rb')
  invisible(readBin(conn, 'raw', n = rlen))
  afmatrix = countmatrix = matrix(NA, nsnp, nind)
  colnames(afmatrix) = colnames(countmatrix) = inds
  rownames(afmatrix) = rownames(countmatrix) = snpfile$SNP
  popind2 = which(!is.na(indfile$X3))
  popind3 = c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`))
  popindmat = matrix(sapply(1:nind, function(i) indfile$X3[popind2] == inds[i]), ncol=nind)
  cnt = 1
  sumna = 0
  count_nonmissing = function(x) sum(!is.na(x))
  while(cnt <= nsnp) {
    if(cnt+blocksize > nsnp) {
      blocksize = nsnp-cnt+1
      popind3 = sort(c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`)))
    }
    bitmat = matrix(as.integer(rawToBits(readBin(conn, 'raw', n = rlen*blocksize))), ncol=8, byrow = TRUE)
    gmat = matrix(c(t(bitmat[,c(8,6,4,2)]*2+bitmat[,c(7,5,3,1)]))[popind3], ncol=blocksize)
    gmat[gmat==3]=NA # assuming non-missing genotypes are 0, 1, 2 missing is 3
    popfreqs = sapply(1:nind, function(i) colMeans(gmat[popindmat[,i],, drop=FALSE], na.rm=TRUE))
    popfreqs[is.nan(popfreqs)] = NA
    sumna = sumna + sum(is.na(popfreqs))
    afmatrix[cnt:(cnt+blocksize-1),] = popfreqs
    if(verbose) cat(paste0('\r', (cnt-1)/1e3, 'k SNPs read...'))
    cnt = cnt+blocksize
  }
  if(verbose) {
    cat('\n')
    alert_success(paste0(cnt-1, ' SNPs read in total\n'))
    alert_warning(paste0(sumna, ' genotypes are missing (on average ', round(sumna/nind), ' per sample)\n'))
  }
  if(auto_only) keepsnps = snpfile %>% mutate(i = 1:n()) %>% filter(as.numeric(gsub('[a-zA-Z]+', '', CHR)) %in% 1:22) %$% i
  else keepsnps = seq_len(nrow(snpfile))
  outdat = treat_missing(afmatrix[keepsnps,], NULL, snpfile[keepsnps,], na.action = na.action, verbose = verbose)
  outlist = list(geno = outdat$afmatrix, ind = indfile %>% filter(!is.na(X3)), snp = snpfile %>% slice(keepsnps))
  outlist
}



#' Read allele frequencies from PLINK files
#'
#' @export
#' @param pref prefix of PLINK files (files have to end in \code{.bed}, \code{.bim}, \code{.fam}).
#' @param pops vector of populations from which to compute allele frequencies.
#' @param inds vector of samples from which to compute allele frequencies. \code{pops} and \code{inds} cannot be specified at the same time. If none are specified, populations will be extracted from the first column in the \code{.fam} file
#' @param na.action what to do with SNPs with missing data.
#' \itemize{
#' \item \code{none} (default) output table will have missing data
#' \item \code{impute} allele frequencies will be imputed from other populations
#' \item \code{remove} SNPs with missing data will be removed
#' }
#' @param return_matrix return data as list of matrices
#' @return a list with two items: allele frequency data and individual counts. Unless \code{return_matrix} is specified, first six columns are bimfile, remaining columns are data for each population.
#' @examples
#' \dontrun{
#' afdat = plink_to_aftable(prefix, pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
plink_to_aftable = function(pref, pops = NULL, inds = NULL, na.action = 'none',
                            return_matrix = FALSE, verbose = FALSE) {
  # This is based on Gad Abraham's "plink2R" package
  # Modified to return per-group allele frequencies rather than raw genotypes.

  stopifnot(is.null(pops) | is.null(inds))
  if(verbose) alert_info('Reading allele frequencies from PLINK file...\n')

  bedfile = paste0(pref, '.bed')
  famfile = paste0(pref, '.fam')
  bimfile = paste0(pref, '.bim')
  nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
  bim = read_table2(bimfile, col_names = nam, progress = FALSE, col_types = cols())
  fam = read_table2(famfile, col_names = FALSE, progress = FALSE, col_types = cols())
  if(is.null(pops)) pops = fam$X1
  indvec = pop_indices(fam, pops, inds)
  indvec2 = which(indvec > 0)
  keepinds = unique(fam[[is.null(pops)+1]][indvec > 0])
  keepsnps = bim %>% mutate(i = 1:n()) %>% filter(as.numeric(gsub('[a-zA-Z]+', '', CHR)) %in% 1:22) %$% i
  afs = read_plink_afs_cpp(normalizePath(bedfile), indvec, indvec2, verbose = verbose)
  afmatrix = afs[[1]][keepsnps,]
  countmatrix = afs[[2]][keepsnps,]
  rownames(afmatrix) = rownames(countmatrix) = bim$SNP[keepsnps]
  colnames(afmatrix) = colnames(countmatrix) = keepinds

  outdat = treat_missing(afmatrix, countmatrix, bim[keepsnps,], na.action = na.action, verbose = verbose)
  outlist = list(afs = outdat$afmatrix, counts = outdat$countmatrix)
  if(!return_matrix) outlist = map(outlist, ~bind_cols(outdat$snpfile, as_tibble(.)))
  outlist$snpfile = outdat$snpfile
  outlist
}

pop_indices = function(famdat, pops=NULL, inds=NULL) {
  # returns vector with length nindiv where each entry corresponds to a population (1 to npop).
  # 0 means do not use
  # famdat has two columns: group ID, IID
  famdat %<>% select(1:2) %>% set_colnames(c('FID', 'IID'))
  stopifnot(all(pops %in% famdat$FID))
  stopifnot(all(inds %in% famdat$IID))
  if(is.null(pops)) {
    pops = famdat$IID
    famdat$FID = pops
  }
  if(is.null(inds)) inds = famdat$IID
  famdat %>% mutate(index = ifelse(FID %in% pops & IID %in% inds, FID, NA)) %$%
    replace_na(as.numeric(factor(index, levels = unique(index))), 0)
}


#' @export
read_plink = function(pref, inds = NULL, na.action = 'none', auto_only = TRUE, verbose = FALSE) {
  # This is based on Gad Abraham's "plink2R" package, but genotypes are m x n, not n x m
  # Modified to return per-group allele frequencies rather than raw genotypes.
  stopifnot(!is.null(auto_only))
  if(verbose) alert_info('Reading PLINK files...\n')

  bedfile = paste0(pref, '.bed')
  famfile = paste0(pref, '.fam')
  bimfile = paste0(pref, '.bim')
  nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
  bim = read_table2(bimfile, col_names = nam, col_types = cols(), progress = FALSE)
  fam = read_table2(famfile, col_names = FALSE, col_types = cols(), progress = FALSE)
  if(is.null(inds)) inds = fam[[2]]
  indvec2 = which(fam[[2]] %in% inds)
  keepinds = fam[[2]][indvec2]

  if(auto_only) keepsnps = bim %>% mutate(i = 1:n()) %>% filter(as.numeric(gsub('[a-zA-Z]+', '', CHR)) %in% 1:22) %$% i
  else keepsnps = 1:nrow(bim)
  bim %<>% slice(keepsnps)
  g = read_plink_afs_cpp(normalizePath(bedfile), indvec=seq_len(nrow(fam)), indvec2, verbose = verbose)[[1]]
  g = g[keepsnps,] * 2
  rownames(g) = bim$SNP
  colnames(g) = keepinds

  outdat = treat_missing(g, NULL, bim, na.action = na.action, verbose = verbose)
  outlist = list(bed = outdat$afmatrix, fam = fam[indvec2,], bim = outdat$snpfile)
  outlist
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
    separate('X1', c('pop1', 'pop2', 'fst','fit','est2','diff','se','z'), sep=' +', convert = TRUE) %>%
    select(-.data$fst)

  f3 = dat %>% filter(grepl(' ff3fit: ', .data$X1)) %>%
    separate('X1', c('pop2', 'pop3', 'ff3fit','fit','est'), sep=' +', convert = TRUE) %>%
    select(-.data$ff3fit)

  fststart = str_which(dat$X1, '^fst:')[1]+2
  fstend = str_which(dat$X1, '^f2:')[1]-1
  f2start = fstend+3
  f2end = str_which(dat$X1, '^ff3:')[1]-1
  pops = dat$X1 %>% str_subset('^population:') %>% str_squish %>% word(3)
  denom = 1000

  f21 = dat %>% slice(f2start:f2end) %>% separate(X1, c('pop1', pops), ' +', T, T) %>% mutate(pop1 = pops) %>% pivot_longer(-pop1, 'pop2', values_to = 'est') %>% mutate(est = est/denom)
  fst = dat %>% slice(fststart:fstend) %>% separate(X1, c('pop1', pops), ' +', T, T) %>% mutate(pop1 = pops) %>% pivot_longer(-pop1, 'pop2', values_to = 'fst') %>% mutate(fst = fst/denom)
  f2before = f21 %>% left_join(fst, by = c('pop1', 'pop2')) %>%
    mutate(pop1 = str_sub(pop1, 1, 3), pop2 = str_sub(pop2, 1, 3))

  f2 %<>% left_join(f2before, by = c('pop1', 'pop2'))

  namedList(edges, score, f2, f3)
}

#' Read qpGraph graph file
#' @export
#' @param graphfile file with admixture graph in qpGraph format.
#' @return graph represented as two column edge matrix.
parse_qpgraph_graphfile = function(graphfile) {
  # reads graph in qpGraph format
  # returns edge matrix (adjacency list)
  lines = read_lines(graphfile) %>%
    tibble %>% set_colnames('V1')
  namemap = lines %>% filter(grepl('^label', .data$V1)) %>%
    separate('V1', c('type', 'label', 'name'), sep = '\\s+', extra = 'drop') %>%
    select(-type) %>% deframe


  dat = lines %>%
    filter(grepl('^edge|redge|ledge|admix', .data$V1)) %>%
    separate('V1', c('type', 'name', 'from', 'to'), sep = '\\s+', extra = 'drop') %>%
    mutate(type = recode(type, ledge = 'edge', redge = 'edge'))
  admix1 = dat %>% filter(.data$type=='admix') %>% mutate(type='edge', to=.data$name, name='')
  admix2 = dat %>% filter(.data$type=='admix') %>% mutate(type='edge', from=.data$to, to=.data$name, name='')
  dat %>%
    bind_rows(admix1) %>%
    bind_rows(admix2) %>%
    filter(.data$type != 'admix') %>%
    mutate(from = recode(from, !!!namemap),
           to = recode(to, !!!namemap)) %>%
    select(.data$from, .data$to) %>% as.matrix
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

# Read qpDstat parameter file
# @export
# @param parfile parameter file for qpDstat.
# @return named list of qpGraph parameters.
parse_qpdstat_parfile = function(parfile) {
  # reads qpGraph parfile
  # returns named list of parameters
  # all genotype files have to have same prefix

  dat = read_table2(parfile, comment = '#', col_names = c('par', 'value'), col_types = cols())
  dat %>% mutate(par = str_replace_all(.data$par, c(':$'='', 'genotypename'='pref')),
                  value = str_replace_all(.data$value,
                                          c('\\.geno$'='',
                                            'SSS'=filter(dat, .data$par=='SSS:')$value[1],
                                            'DIR'=filter(dat, .data$par=='DIR:')$value[1])),
                  value=ifelse(.data$value=='YES', TRUE, ifelse(.data$value=='NO', FALSE, .data$value))) %>%
    filter(!.data$par %in% c('DIR', 'SSS', 'indivname', 'snpname')) %>%
    deframe %>%
    as.list %>%
    as_tibble()
}



#' Read qpAdm output file
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
    set_names(c('weight', 'jmean', 'jse')) %>% as_tibble %>% mutate(z = jmean/jse)
  weights = tibble(target, left) %>% bind_cols(coefs)

  sig = do.call(rbind, str_split(dat[sigstart:sigend], ' +')) %>% as.data.frame(stringsAsFactors=F) %>% select(-1) %>% set_colnames(c('pattern', 'wt', 'dof', 'chisq', 'tail prob', left, 'feasible')) %>% mutate(feasible = feasible != 'infeasible')

  namedList(weights, sig)
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
  fst = dat[fststart:fstend] %>% enframe %>% separate(value, c('pop1', pops), ' ', T, T) %>% select(-name) %>% mutate(pop1 = pops) %>% pivot_longer(-pop1, 'pop2', values_to = 'fst') %>% mutate(fst = fst/denom)
  f2 = dat[f2start:f2end] %>% enframe %>% separate(value, c('pop1', pops), ' ', T, T) %>% select(-name) %>% mutate(pop1 = pops) %>% pivot_longer(-pop1, 'pop2', values_to = 'f2') %>% mutate(f2 = f2/denom)

  fst %>% left_join(f2, by = c('pop1', 'pop2'))

}



#' Write f2 block jackknife estimates to disk
#'
#' This function takes a 3d array of f2 block jackknife estimates, splits it by population pair, and writes each pair to a separate \code{.RData} file.
#' @export
#' @param f2_block 3d array with f2 block jackknife estimates. The first two dimensions of \code{f2_block} have to have population names. Each file will be stored under \code{{outdir}/{pop1}/{pop2}.RData}.
#' @param outdir directory into which to write the files.
#' @param overwrite should existing \code{.RData} files be overwritten?
#' @return NULL
#' @seealso \code{\link{read_f2}}
#' @examples
#' \dontrun{
#' write_f2(f2_block, outdir = 'path/to/f2stats/')
#' }
write_f2 = function(f2_block, outdir, overwrite=FALSE) {

  if(!dir.exists(outdir)) dir.create(outdir)
  d1 = dim(f2_block)[1]
  d2 = dim(f2_block)[2]
  nam1 = dimnames(f2_block)[[1]]
  nam2 = dimnames(f2_block)[[2]]
  for(i in seq_len(d1)) {
    for(j in seq_len(d2)) {
      pop1 = min(nam1[i], nam2[j])
      pop2 = max(nam1[i], nam2[j])
      f2 = as.vector(f2_block[i, j, ])
      dir = paste0(outdir, '/', pop1, '/')
      fl = paste0(dir, pop2, '.rds')
      if(!dir.exists(dir)) dir.create(dir)
      if(!file.exists(fl) | overwrite) saveRDS(f2, file=fl)
    }
  }
}

#' Read f2 block jackknife estimates from disk
#'
#' This function reads f2 block jackknife estimates which were writtend to disk by \code{\link{write_f2}} and returns a 3d array of f2 block jackknife estimates.
#' @export
#' @param f2_dir directory from which to read files
#' @param pops the populations for which f2 statistics should be read. Defaults to all populations, which may require a lot of memory.
#' @return a 3d array of block jackknife estimates
#' @seealso \code{\link{write_f2}}
#' @examples
#' \dontrun{
#' read_f2(f2_dir, pops = c('Zerg', 'Terran', 'Protoss'))
#' }
read_f2 = function(f2_dir, pops = NULL) {
  # reads f2 block jackknife .RData files and returns 3d array
  # pops is vector of populations which should be read. defaults to all populations.
  remove_na = is.null(pops)
  if(is.null(pops)) {
    pops = list.dirs(f2_dir, full.names = FALSE, recursive = FALSE)
  }
  old = file.exists(paste0(f2_dir, '/block_lengths.RData'))
  if(old) {
    load(paste0(f2_dir, '/block_lengths.RData'))
  } else {
    block_lengths = readRDS(paste0(f2_dir, '/block_lengths.rds'))
  }
  numblocks = length(block_lengths)
  numpops = length(pops)
  f2_blocks = array(NA, c(numpops, numpops, numblocks))
  f2_blocks = structure(f2_blocks, block_lengths = block_lengths)
  dimnames(f2_blocks)[[1]] = dimnames(f2_blocks)[[2]] = pops
  for(pop1 in pops) {
    for(pop2 in pops) {
      if(pop1 <= pop2) {
        if(old) {
          fl = paste0(f2_dir, '/', pop1, '/', pop2, '.RData')
          if(file.exists(fl)) {
            load(fl)
            f2_blocks[pop1, pop2, ] = f2_blocks[pop2, pop1, ] = f2
            if(any(is.na(f2))) warning(paste0('missing values in ', pop1, ' - ', pop2, '!'))
          } else warning(paste0('file ', fl, ' not found!'))
        } else {
          fl = paste0(f2_dir, '/', pop1, '/', pop2, '.rds')
          if(file.exists(fl)) {
            f2 = readRDS(fl)
            f2_blocks[pop1, pop2, ] = f2_blocks[pop2, pop1, ] = f2
            if(any(is.na(f2))) warning(paste0('missing values in ', pop1, ' - ', pop2, '!'))
          } else warning(paste0('file ', fl, ' not found!'))
        }
      }
    }
  }
  if(remove_na) {
    keep = apply(f2_blocks, 1, function(x) sum(is.na(x)) == 0)
    f2_blocks = f2_blocks[keep, keep, ]
  }
  f2_blocks
}


write_pairdat = function(aa_arr, nn_arr, outdir, overwrite = FALSE) {

  if(!dir.exists(outdir)) dir.create(outdir)
  d1 = dim(aa_arr)[1]
  d2 = dim(aa_arr)[2]
  nam1 = dimnames(aa_arr)[[1]]
  nam2 = dimnames(aa_arr)[[2]]
  for(i in seq_len(d1)) {
    for(j in seq_len(d2)) {
      ind1 = min(nam1[i], nam2[j])
      ind2 = max(nam1[i], nam2[j])
      dir = paste0(outdir, '/', ind1, '/')
      fl = paste0(dir, ind2, '.rds')
      if(!dir.exists(dir)) dir.create(dir)
      prods = cbind(aa = aa_arr[i, j, ], nn = nn_arr[i, j, ])
      if(!file.exists(fl) | overwrite) saveRDS(prods, file=fl)
    }
  }

}

# #' Read block_lengths from disk
# #'
# #' @export
# #' @param f2_dir directory from which to read files
# #' @return a vector with block_lengths
# #' @examples
# #' \dontrun{
# #' read_bl(f2_dir)
# #' }
# read_bl = function(f2_dir) {
#   load(paste0(f2_dir, '/block_lengths.RData'))
#   block_lengths
# }


#' Write allele frequency estimates to disk
#'
#' This function takes a 3d array of f2 block jackknife estimates, splits it by population pair, and writes each pair to a separate \code{.RData} file.
#' @export
#' @param afmat matrix with allele frequency estimates. Has to have column names. Data for each population will be stored in a variable called \code{allele_frequencies} under \code{{path}/{pop}.RData}
#' @param path directory into which to write the files.
#' @param countmat optional matrix with allele counts for each SNP and population. Dimensions have to match \code{afmat}. Data will be stored in same file (\code{{path}/{pop}.RData}) in the variable \code{num_individuals}
#' @param overwrite should existing \code{.RData} files be overwritten?
#' @seealso \code{\link{read_f2}}
#' @examples
#' \dontrun{
#' write_afs(afmat, path = 'path/to/afs/')
#' }
write_afs = function(afmat, path, countmat = NULL, overwrite = FALSE) {
  # needs to be updated

  stopifnot(is.null(countmat) || all(dim(afmat) == dim(countmat)))
  if(!dir.exists(path)) dir.create(path)
  pops = colnames(afmat)

  for(i in seq_len(pops)) {
    allele_frequencies = afmat[,i]
    num_alleles = countmat[,i]
    population = pops[i]
    f2 = as.vector(f2_block[i, j, ])
    fl = paste0(path, '/', pop, '.RData')
    if(!file.exists(fl) | overwrite) save(allele_frequencies, num_alleles, population, file=fl)
  }
}


#' Split an allele frequency matrix into blocks
#'
#' This function splits a large allele frequency matrix into smaller blocks with \code{pops_per_block} populations per block, and saves them as \code{.RData} files in \code{outprefix}
#' @export
#' @param afmat the allele frequency matrix to be split
#' @param pops_per_block the number of populations per block
#' @param outprefix the prefix of the output files. Directories are not created on the fly.
#' @seealso \code{\link{packedancestrymap_to_aftable}}, \code{\link{write_split_f2_block}}
#' @examples
#' \dontrun{
#' afmatall = packedancestrymap_to_aftable('path/to/packedancestrymap_prefix', allpopulations,
#'                                          na.action = 'none', return_matrix = TRUE)
#' split_afmat(afmatall, pops_per_block = 20, outprefix = 'afmatall_split_v42.1/afmatall_')
#' }
split_afmat = function(afmat, pops_per_block, outprefix, verbose = TRUE) {
  # splits afmat into numparts parts, and saves them to {outprefix}_{part}.RData
  npops = ncol(afmat)
  starts = seq(1, npops, pops_per_block)
  numparts = length(starts)
  ends = c(lead(starts)[-numparts]-1, npops)
  for(i in seq_len(numparts)) {
    if(verbose) cat(paste0('\rpart ', i, ' of ', numparts))
    afs = afmat[, starts[i]:ends[i]]
    save(afs, file = paste0(outprefix, i, '.RData'))
  }
  if(verbose) cat('\n')
}


#' Write f2 jackknife blocks to disk
#'
#' Prepare data for various admixtools functions. Reads data from packedancestrymap or PLINK files, and computes allele frequencies and f2 block jackknife statistics for selected populations. Data is either written to disk in files for each population and population pair, or returned as list. This function calls \code{\link{packedancestrymap_to_aftable}} / \code{\link{plink_to_aftable}} and \code{\link{afs_to_f2_blocks}}. The default treatment of missing data is \code{na.action = 'none'} which is similar to \code{useallsnps: YES}. To get results that are closer to \code{useallsnps: NO}, use \code{na.action = 'remove'}.
#' @export
#' @param pref prefix of packedancestrymap or PLINK files. packedancestrymap has to end in \code{.geno}, \code{.snp}, \code{.ind}, PLINK has to end in \code{.bed}, \code{.bim}, \code{.fam}
#' @param outdir the directory to which to write data
#' @param pops the populations for which data should be extracted. Defaults to all populations, which may require a lot of memory. If \code{NULL}, each individual will be its own population.
#' @param inds the individuals for which data should be extracted. \code{pops} and \code{inds} cannot be specified at the same time. If none are specified, all populations will be extracted
#' @param na.action what to do with SNPs with missing data.
#' \itemize{
#' \item \code{none} (default) output table will have missing data
#' \item \code{impute} allele frequencies will be imputed from other populations
#' \item \code{remove} SNPs with missing data will be removed
#' }
#' @param maxmem split up allele frequency data into blocks, if memory requirements exceed \code{maxmem} MB.
#' @param overwrite should existing files be overwritten?
#' @param loo if \code{TRUE} (the default), return the total estimates minus the estimates from each block. if \code{FALSE}, return the estimate from each block.
#' @param verbose print progress updates
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' extract_f2(pref, pops = c('Protoss', 'Zerg', 'Terran'))
#' }
extract_f2 = function(pref, outdir, pops = NULL, na.action = 'none', maxmem = 8000,
                      overwrite = FALSE, loo = TRUE, verbose = TRUE) {

  outdir = normalizePath(outdir, mustWork = FALSE)
  afdat = anygeno_to_aftable(pref, pops, na.action = na.action, verbose = verbose)
  block_lengths = get_block_lengths(afdat$snpfile)
  f2_blocks = afs_to_f2_blocks(afdat$afs, afdat$counts, block_lengths,
                               outdir = outdir, overwrite = overwrite,
                               maxmem = maxmem, loo = loo, verbose = verbose)
  allele_counts = colMeans(afdat$counts, na.rm = TRUE)
  saveRDS(block_lengths, file = paste0(outdir, '/block_lengths.rds'))
  saveRDS(allele_counts, file = paste0(outdir, '/allele_counts.rds'))
  if(verbose) alert_info(paste0('Data written to ', outdir, '/\n'))
}

anygeno_to_aftable = function(pref, pops = NULL, na.action = 'none', verbose = TRUE) {

  if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
    geno_to_aftable = packedancestrymap_to_aftable
  } else if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) {
    geno_to_aftable = plink_to_aftable
  } else stop('Genotype files not found!')
  afdat = geno_to_aftable(pref, pops, na.action = na.action,
                          return_matrix = TRUE, verbose = verbose)
  afdat
}

read_anygeno = function(pref, inds = NULL, na.action = 'none', verbose = TRUE) {

  if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
    read_geno = function(...) {g = read_packedancestrymap(...); list(bed = g$geno, fam = g$ind, bim = g$snp)}
  } else if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) {
    read_geno = read_plink
  } else stop('Genotype files not found!')
  read_geno(pref, inds, na.action = na.action, verbose = verbose)

}

#' Write data for individual pairs to disk
#'
#' @export
extract_indpairs = function(pref, outdir, inds = NULL, na.action = 'none',
                           maxmem = 8000, overwrite = FALSE, verbose = TRUE) {
  # continue here / not tested

  g = read_anygeno(pref, inds, na.action = na.action, verbose = verbose)

  block_lengths = get_block_lengths(g$bim)
  xmat_to_pairdat(g$bed, block_lengths,
                  outdir = outdir, overwrite = overwrite,
                  maxmem = maxmem, verbose = verbose)
}



#' Read f2 jackknife blocks from genotype data
#'
#' Prepare data for various admixtools functions. Reads data from packedancestrymap or PLINK files, and computes allele frequencies and f2 block jackknife statistics for selected populations. Data is either written to disk in files for each population and population pair, or returned as list. This function calls \code{\link{packedancestrymap_to_aftable}} / \code{\link{plink_to_aftable}} and \code{\link{afs_to_f2_blocks}}. The default treatment of missing data is \code{na.action = 'none'} which is similar to \code{useallsnps: YES}. To get results that are closer to \code{useallsnps: NO}, use \code{na.action = 'remove'}.
#' @export
#' @param pref prefix of packedancestrymap or PLINK files. packedancestrymap has to end in \code{.geno}, \code{.snp}, \code{.ind}, PLINK has to end in \code{.bed}, \code{.bim}, \code{.fam}
#' @param pops the populations for which data should be extracted. Defaults to all populations, which may require a lot of memory. If \code{NULL}, each individual will be its own population.
#' @param inds the individuals for which data should be extracted. \code{pops} and \code{inds} cannot be specified at the same time. If none are specified, all populations will be extracted
#' @param na.action what to do with SNPs with missing data.
#' \itemize{
#' \item \code{none} (default) output table will have missing data
#' \item \code{impute} allele frequencies will be imputed from other populations
#' \item \code{remove} SNPs with missing data will be removed
#' }
#' @param maxmem split up allele frequency data into blocks, if memory requirements exceed \code{maxmem} MB.
#' @param verbose print progress updates
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' f2_blocks = f2_from_geno(pref, pops = c('Protoss', 'Zerg', 'Terran'))
#' }
f2_from_geno = function(pref, poplist = NULL, pops = NULL, inds = NULL,
                        na.action = 'none', maxmem = 8000, verbose = TRUE) {
  # this should become a replacement for extract_data, with bypassing writing to disk and reading genotypes straight to afs
  # put back inds
  # think of how what to do with poplist. maybe only use inds first, then after reading geno group by pop?

  afdat = anygeno_to_aftable(pref, pops, na.action = na.action, verbose = verbose)
  afs = afdat$afs
  counts = afdat$counts
  #hets = which(apply(afs, 1, function(x) sum(na.omit(x) > 0 & na.omit(x) < 1)) > 0)
  #afs = afs[hets,]
  #counts = counts[hets,]
  #block_lengths = get_block_lengths(afdat$snpfile[hets,])
  block_lengths = get_block_lengths(afdat$snpfile)
  f2_blocks = afs_to_f2_blocks(afs, counts, block_lengths,
                               maxmem = maxmem, loo = TRUE, verbose = verbose)
  f2_blocks
}

# this should not really be needed
f2_from_geno_indivs = function(pref, poplist = NULL, pops = NULL, inds = NULL, na.action = 'none', maxmem = 8000, verbose = TRUE) {

  # think of how what to do with poplist. maybe only use inds first, then after reading geno group by pop?

  if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
    read_geno = function(...) {g = read_packedancestrymap(...); list(bed = g$geno, bim = g$snp, fam = g$ind)}
  } else if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) {
    read_geno = read_plink
  } else stop('Genotype files not found!')
  g = read_geno(pref, inds, na.action = na.action, verbose = verbose)

  block_lengths = get_block_lengths(g$bim)
  pairdat = xmat_to_pairdat(g$bed, block_lengths, maxmem = maxmem, verbose = verbose)

  f2_blocks = indpairs_to_f2blocks(pairdat$indivs, pairdat$pairs, poplist, block_lengths)
  f2_blocks
}



#' Read f2 jackknife blocks from disk
#'
#' @export
#' @param dir directory with precomputed f2 jackknife blocks
#' @param poplist the populations for which data should be extracted. Defaults to all populations, which may require a lot of memory. If \code{NULL}, each individual will be its own population.
#' @param verbose print progress updates
#' @examples
#' \dontrun{
#' dir = 'my/f2/dir/'
#' f2_blocks = f2_from_precomp(dir, pops = c('Protoss', 'Zerg', 'Terran'))
#' }
f2_from_precomp = function(dir, poplist = NULL, verbose = TRUE) {

  if(file.exists(paste0(dir, '/indivs.rds'))) {
    if(verbose) alert_info(paste0('Reading precomputed data for ', nrow(poplist), ' individuals...\n'))
    f2_blocks = f2_from_precomp_indivs(dir, poplist = poplist)$f2_blocks
  } else {
    if(!is.character(poplist)) poplist = poplist$pop
    if(verbose) alert_info(paste0('Reading precomputed data for ', length(poplist), ' populations...\n'))
    f2_blocks = read_f2(dir, poplist)
  }
  f2_blocks
}


#' @export
f2_from_precomp_indivs = function(dir, poplist = NULL) {

  remove_na = FALSE
  if(is.null(poplist)) {
    ind = pop = list.dirs(dir, full.names = FALSE, recursive = FALSE)
    poplist = tibble(pop, ind)
    remove_na = TRUE
  }
  inds = poplist$ind
  indivs = readRDS(paste0(dir, '/indivs.rds'))
  block_lengths = readRDS(paste0(dir, '/block_lengths.rds'))
  bl = seq_along(block_lengths)
  pairs = tibble()
  for(ind1 in inds) {
    for(ind2 in inds) {
      if(ind1 <= ind2) {
        fl = paste0(dir, '/', ind1, '/', ind2, '.rds')
        if(file.exists(fl)) {
          prods = readRDS(fl)
          if(any(is.na(prods))) warning(paste0('missing values in ', ind1, ' - ', ind2, '!'))
          pairs %<>% bind_rows(tibble(ind1, ind2, bl) %>% bind_cols(as_tibble(prods)))
        } else warning(paste0('file ', fl, ' not found!'))
      }
    }
  }
  f2_blocks = indpairs_to_f2blocks(indivs, pairs, poplist, block_lengths)
  if(remove_na) {
    keep = apply(f2_blocks, 1, function(x) sum(is.na(x)) == 0)
    f2_blocks = f2_blocks[keep, keep, ]
  }
  #f2_blocks
  namedList(f2_blocks, indivs, pairs, poplist)
}

#' @export
read_indivs = function(pairs, dir, add) {
  # reads data for selected individuals (and pairs)
  # updates (returns) pairs (indivs should contain all sampels on disk)

  stopifnot(length(intersect(pairs$ind1, add)) == 0)
  block_lengths = readRDS(paste0(dir, '/block_lengths.rds'))
  bl = seq_along(block_lengths)
  inds = union(pairs$ind1, add)

  for(ind1 in inds) {
    for(ind2 in inds) {
      if(ind1 <= ind2 && (ind1 %in% add || ind2 %in% add)) {
        fl = paste0(dir, '/', ind1, '/', ind2, '.rds')
        if(file.exists(fl)) {
          prods = readRDS(fl)
          if(any(is.na(prods))) warning(paste0('missing values in ', ind1, ' - ', ind2, '!'))
          newdat = tibble(ind1, ind2, bl) %>% bind_cols(as_tibble(prods))
          #newdat %<>% bind_rows(newdat %>% rename(ind1 = ind2, ind2 = ind1))
          pairs %<>% bind_rows(newdat)
        } else warning(paste0('file ', fl, ' not found!'))
      }
    }
  }
  pairs
}


xmat_to_pairdat = function(xmat, block_lengths, f2_denom = 1, maxmem = 8000,
                           outdir = NULL, overwrite = FALSE, verbose = TRUE) {
  # input is genotype matrix
  # either writes data to outdir, or returns it
  # output is either a list with a, n, aa, nn information for all indivs and pairs, or nothing

  nr = nrow(xmat)
  nc = ncol(xmat)
  mem1 = lobstr::obj_size(xmat)
  mem2 = mem1*nc*2
  numsplits = ceiling(mem2/1e6/maxmem)
  width = ceiling(nc/numsplits)
  starts = seq(1, nc, width)
  numsplits2 = length(starts)
  ends = c(lead(starts)[-numsplits2]-1, nc)

  if(verbose) {
    alert_info(paste0('genotype matrix for ', nr, ' SNPs and ', nc, ' samples is ', round(mem1/1e6), ' MB\n'))
    alert_warning(paste0('matrices of pairwise stats for all SNPs and sample pairs require ', round(mem2/1e6), ' MB\n'))
    if(numsplits2 > 1) alert_info(paste0('splitting into ', numsplits2, ' blocks of ', width, ' samples and up to ', maxmem, ' MB (', choose(numsplits2+1,2), ' block pairs)\n'))
  }

  xmat %<>% fix_ploidy
  ploidy = attr(xmat, 'ploidy')
  bids = rep(seq_along(block_lengths), block_lengths)

  indivs = as_tibble(xmat) %>%
    rownames_to_column(var='SNP') %>%
    pivot_longer(-SNP, names_to = 'ind', values_to = 'a') %>%
    mutate(n = 1-is.na(a), n = n*ploidy[ind], a = replace_na(a, 0)) %>%
    add_column(bl = rep(bids, each=nc), .before = 'ind') %>%
    group_by(bl, ind) %>% summarize(a = mean(a), n = mean(n)) %>% ungroup

  cmb = combn(0:numsplits2, 2)+(1:0)
  aa_list = nn_list = replicate(numsplits2, list())
  prodrray = function(m1, m2) rray(t(m1), dim=c(ncol(m1), 1, nrow(m1))) *
    rray(t(m2), dim=c(1, ncol(m2), nrow(m1)))

  for(i in 1:ncol(cmb)) {
    if(numsplits2 > 1 & verbose) cat(paste0('\rsample pair block ', i, ' out of ', ncol(cmb)))
    c1 = cmb[1,i]
    c2 = cmb[2,i]
    s1 = starts[c1]:ends[c1]
    s2 = starts[c2]:ends[c2]

    a1 = xmat[,s1, drop=F]
    a2 = xmat[,s2, drop=F]
    n1 = (!is.na(a1)) * rep(ploidy[s1], each = nr)
    n2 = (!is.na(a2)) * rep(ploidy[s2], each = nr)
    a1 %<>% replace_na(0)
    a2 %<>% replace_na(0)
    aa_arr = prodrray(a1, a2)
    nn_arr = prodrray(n1, n2)
    dimnames(aa_arr)[1:2] = dimnames(nn_arr)[1:2] = list(colnames(a1), colnames(a2))
    aa_subblock = block_arr_mean(aa_arr, block_lengths)
    nn_subblock = block_arr_mean(nn_arr, block_lengths)
    if(!is.null(outdir)) {
      write_pairdat(aa_subblock, nn_subblock, outdir = outdir, overwrite = overwrite)
    } else {
      aa_list[[c1]][[c2]] = aa_subblock
      aa_list[[c2]][[c1]] = aperm(aa_list[[c1]][[c2]], c(2,1,3))
      nn_list[[c1]][[c2]] = nn_subblock
      nn_list[[c2]][[c1]] = aperm(nn_list[[c1]][[c2]], c(2,1,3))
    }
  }
  if(numsplits2 > 1 & verbose) cat('\n')
  if(!is.null(outdir)) {
    ifile = paste0(outdir, '/indivs.rds')
    bfile = paste0(outdir, '/block_lengths.rds')
    if(overwrite || !file.exists(ifile)) saveRDS(indivs, file = ifile)
    if(overwrite || !file.exists(bfile)) saveRDS(block_lengths, file = bfile)

  } else {

    assemble_arrays = function(l) do.call(abind, list(lapply(l, function(x) do.call(abind, list(x, along=2))), along=1))
    aa_arr_full = assemble_arrays(aa_list) %>% structure(block_lengths = block_lengths)
    nn_arr_full = assemble_arrays(nn_list) %>% structure(block_lengths = block_lengths)
    dimnames(aa_arr_full) = dimnames(nn_arr_full) = list(ind1 = colnames(xmat), ind2 = colnames(xmat), bl = 1:length(block_lengths))

    pairs = aa_arr_full %>% as.array %>% as.tbl_cube(met_name = 'aa') %>% as_tibble %>% mutate(nn = nn_arr_full %>% as.array %>% as.tbl_cube(met_name = 'nn') %>% as_tibble %$% nn)

    return(namedList(indivs, pairs))

  }

}


