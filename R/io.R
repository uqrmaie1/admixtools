
#' Read allele frequencies from packedancestrymap genotype file
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
#' @return a list with two items: allele frequency data and individual counts. Unless \code{return_matrix} is specified, first six columns are bimfile, remaining columns are data for each population.
#' @examples
#' \dontrun{
#' afdat = packedancestrymap_to_aftable(prefix, pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
packedancestrymap_to_aftable = function(pref, pops=NULL, inds=NULL, blocksize=1000, na.action='none', return_matrix=FALSE, verbose = TRUE) {
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
  colnames(afmatrix) = colnames(countmatrix) = pops
  rownames(afmatrix) = rownames(countmatrix) = snpfile$SNP
  popind3 = c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`))
  popindmat = matrix(sapply(1:numpop, function(i) indfile$X3[popind2] == pops[i]), ncol=numpop)
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
    popfreqs = sapply(1:numpop, function(i) colMeans(gmat[popindmat[,i],, drop=FALSE], na.rm=TRUE)/2)
    popcounts = sapply(1:numpop, function(i) apply(gmat[popindmat[,i],, drop=FALSE], 2, count_nonmissing))
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

  outdat = treat_missing(afmatrix, countmatrix, snpfile, na.action = na.action, verbose = verbose)
  outlist = list(afs = outdat$afmatrix, counts = outdat$countmatrix)
  if(!return_matrix) outlist = map(outlist, ~bind_cols(outdat$snpfile, as_tibble(.)))
  outlist$snpfile = outdat$snpfile
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

  bim = read_table2(bimfile, col_names = FALSE, col_types = cols(), progress = FALSE)
  fam = read_table2(famfile, col_names = FALSE, col_types = cols(), progress = FALSE)
  indvec = pop_indices(fam, pops, inds)
  indvec2 = which(indvec > 0)
  keep = unique(fam[[is.null(pops)+1]][indvec > 0])
  afs = read_plink_afs_cpp(normalizePath(bedfile), indvec, indvec2, verbose = verbose)
  afmatrix = afs[[1]]
  countmatrix = afs[[2]]
  rownames(afmatrix) = rownames(countmatrix) = bim[[2]]
  colnames(afmatrix) = colnames(countmatrix) = keep

  outdat = treat_missing(afmatrix, countmatrix, bim, na.action = na.action, verbose = verbose)
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


treat_missing = function(afmatrix, countmatrix = NULL, snpfile = NULL,
                         na.action = 'none', verbose = TRUE) {
  discard = FALSE
  if(na.action == 'impute') {
    afmatrix = mean_impute(afmatrix, by=1)
    isnan = is.nan(afmatrix)
    afmatrix[isnan] = NA
    discard = apply(afmatrix, 1, function(x) all(is.na(x)))
    if(verbose) alert_danger(paste0(sumna-sum(isnan), ' allele frequencies were imputed, (on average ', round((sumna-sum(isnan))/numpop),' per population) ', sum(discard), ' SNPs were removed\n'))
  }
  if(na.action == 'remove') {
    discard = apply(afmatrix, 1, function(x) any(is.na(x)))
    if(verbose) alert_danger(paste0(sum(discard), ' SNPs were removed, ', sum(!discard), ' SNPs remain\n'))
  }
  afmatrix = afmatrix[!discard,]
  countmatrix = countmatrix[!discard,]
  snpfile = snpfile[!discard,]

  namedList(afmatrix, countmatrix, snpfile)
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
    separate('X1', c('pop1', 'pop2', 'fst','f2fit','f2est','diff','se','z'), sep=' +', convert = TRUE) %>%
    select(-.data$fst)

  namedList(edges, score, f2)
}

#' Read qpGraph graph file
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
    set_names(c('weight', 'jmean', 'jse')) %>% as_tibble
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
      fl = paste0(dir, pop2, '.RData')
      if(!dir.exists(dir)) dir.create(dir)
      if(!file.exists(fl) | overwrite) save(f2, pop1, pop2, file=fl)
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
  if(is.null(pops)) pops = list.dirs(f2_dir, full.names=FALSE, recursive=FALSE)
  load(paste0(f2_dir, '/', pops[1], '/', pops[1], '.RData'))
  numblocks = length(f2)
  numpops = length(pops)
  f2_blocks = array(0, c(numpops, numpops, numblocks))
  dimnames(f2_blocks)[[1]] = dimnames(f2_blocks)[[2]] = pops
  for(pop1 in pops) {
    for(pop2 in pops) {
      if(pop1 <= pop2) {
        fl = paste0(f2_dir, '/', pop1, '/', pop2, '.RData')
        load(fl)
        if(any(is.na(f2))) warn(paste0('missing values in ', pop1, ' - ', pop2, '!'))
        f2_blocks[pop1, pop2, ] = f2_blocks[pop2, pop1, ] = f2
      }
    }
  }
  f2_blocks
}


#' Write allele frequency estimates to disk
#'
#' This function takes a 3d array of f2 block jackknife estimates, splits it by population pair, and writes each pair to a separate \code{.RData} file.
#' @export
#' @param afmat matrix with allele frequency estimates. Has to have column names. Data for each population will be stored in a variable called \code{allele_frequencies} under \code{{path}/{pop}.RData}
#' @param path directory into which to write the files.
#' @param countmat optional matrix with individual counts for each SNP and population. Dimensions have to match \code{afmat}. Data will be stored in same file (\code{{path}/{pop}.RData}) in the variable \code{num_individuals}
#' @param overwrite should existing \code{.RData} files be overwritten?
#' @seealso \code{\link{read_f2}}
#' @examples
#' \dontrun{
#' write_afs(afmat, path = 'path/to/afs/')
#' }
write_afs = function(afmat, path, countmat = NULL, overwrite=FALSE) {

  stopifnot(is.null(countmat) || all(dim(afmat) == dim(countmat)))
  if(!dir.exists(path)) dir.create(path)
  pops = colnames(afmat)

  for(i in seq_len(pops)) {
    allele_frequencies = afmat[,i]
    num_individuals = countmat[,i]
    population = pops[i]
    f2 = as.vector(f2_block[i, j, ])
    fl = paste0(path, '/', pop, '.RData')
    if(!file.exists(fl) | overwrite) save(allele_frequencies, num_individuals, population, file=fl)
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


#' Extract subset of genotype data
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
#' @param verbose print progress updates
#' @return Nothing, if \code{outdir} is specified. Otherwise a list with \code{f2_blocks} and \code{block_lengths}
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' data = extract_data(pref, pops = c('Protoss', 'Zerg', 'Terran'))
#' f2_blocks = data$f2_blocks
#' block_lengths = data$block_lengths
#' }
extract_data = function(pref, outdir = NULL, pops = NULL, inds = NULL,
                        na.action = 'none', maxmem = 1000, overwrite = FALSE,
                        verbose = TRUE) {

  stopifnot(is.null(pops) | is.null(inds))
  if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
    geno_to_aftable = packedancestrymap_to_aftable
    snpend = '.snp'
    nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  } else if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) {
    geno_to_aftable = plink_to_aftable
    snpend = '.bim'
    nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
  } else stop('Genotype files not found!')
  snpfile = read_table2(paste0(pref, snpend), col_names = nam,
                        col_types = cols(), progress = FALSE)
  afdat = geno_to_aftable(pref, pops, inds, na.action = na.action,
                          return_matrix = TRUE, verbose = verbose)
  afs = afdat$afs
  pop_counts = ceiling(colMeans(afdat$counts))
  block_lengths = get_block_lengths(afdat$snpfile)
  f2_blocks = afs_to_f2_blocks(afs, pop_counts, block_lengths,
                               outdir = outdir, overwrite = overwrite,
                               maxmem = maxmem, verbose = verbose)
  if(is.null(outdir)) return(namedList(f2_blocks, block_lengths, pop_counts))
  save(block_lengths, file = paste0(outdir, '/block_lengths.RData'))
  save(pop_counts, file = paste0(outdir, '/pop_counts.RData'))
}

