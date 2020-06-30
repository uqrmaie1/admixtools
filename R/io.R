
#' Read allele frequencies from packedancestrymap files
#'
#' This is currently slower than `plink_to_aftable` because it is implemented only in `R`, not in `C++`.
#' @export
#' @param pref prefix of packedancestrymap files (files have to end in `.geno`, `.ind`, `.snp`).
#' @param inds vector of samples from which to compute allele frequencies.
#' @param pops vector of populations from which to compute allele frequencies. If `NULL` (default), populations will be extracted from the third column in the `.ind` file
#' @param blocksize number of SNPs read in each block.
#' @param verbose print progress updates
#' @return a list with two items: allele frequency data and allele counts.
#' @examples
#' \dontrun{
#' afdat = packedancestrymap_to_aftable(prefix, pops = pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
packedancestrymap_to_aftable = function(pref, inds = NULL, pops = NULL, blocksize = 1000, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # pops is vector of populations for which to calculate AFs
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population

  stopifnot(is.null(pops) || is.null(inds) || length(inds) == length(pops))

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = cols(), progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = cols(), progress = FALSE)
  if(is.null(inds)) {
    inds = indfile$X1
  } else {
    indfile$X1[!indfile$X1 %in% inds] = NA
  }
  if(is.null(pops)) {
    pops = indfile$X3
    indfile$X3[!indfile$X1 %in% inds] = NA
  } else {
    indfile$X3[!is.na(indfile$X1)] = pops
  }
  upops = unique(na.omit(pops))

  fl = paste0(pref, '.geno')
  conn = file(fl, 'rb')
  on.exit(close(conn))
  hd = strsplit(readBin(conn, 'character', n = 1), ' +')[[1]]
  close(conn)
  nind = as.numeric(hd[2])
  nsnp = as.numeric(hd[3])
  numpop = length(upops)
  popind = which(indfile$X3 %in% upops)
  popind2 = na.omit(popind)

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nind, ' samples and ', nsnp, ' SNPs\n'))
    alert_info(paste0('Calculating allele frequencies from ', length(popind2), ' samples in ', numpop, ' populations\n'))
    alert_info(paste0('Expected size of allele frequency data: ', round((nsnp*numpop*8+nsnp*112)/1e6), ' MB\n'))
    # 8, 112: estimated scaling factors for AF columns and annotation columns
  }

  rlen = file.info(fl)$size/(nsnp+1)
  conn = file(fl, 'rb')
  invisible(readBin(conn, 'raw', n = rlen))
  afmatrix = countmatrix = matrix(NA, nsnp, numpop)
  colnames(afmatrix) = colnames(countmatrix) = upops
  rownames(afmatrix) = rownames(countmatrix) = snpfile$SNP
  popind3 = c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`))
  popindmat = matrix(sapply(1:numpop, function(i) indfile$X3[popind2] == upops[i]), ncol=numpop)
  cnt = 1
  sumna = 0
  modnum = ifelse(shiny::isRunning(), 1e5, 1e3)
  while(cnt <= nsnp) {
    if(cnt+blocksize > nsnp) {
      blocksize = nsnp-cnt+1
      popind3 = sort(c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`)))
    }
    bitmat = matrix(as.integer(rawToBits(readBin(conn, 'raw', n = rlen*blocksize))), ncol=8, byrow = TRUE)
    gmat = matrix(c(t(bitmat[,c(8,6,4,2)]*2+bitmat[,c(7,5,3,1)]))[popind3], ncol=blocksize)
    gmat[gmat==3]=NA # assuming non-missing genotypes are 0, 1, 2 missing is 3
    if(cnt == 1) {
      ploidy = apply(gmat, 1, function(x) max(1, length(unique(na.omit(x)))-1))
      alert_info(paste0('Detected ', sum(ploidy == 2), ' diploid samples and ',
                        sum(ploidy == 1), ' pseudohaploid samples\n'))
    }
    popfreqs = popcounts = matrix(NA, blocksize, numpop)
    for(i in seq_len(numpop)) {
      gm = gmat[popindmat[,i],, drop=FALSE]
      popfreqs[,i] = colMeans(gm, na.rm=TRUE)/2
      popcounts[,i] = colSums((!is.na(gm))*ploidy[popindmat[,i]])
    }
    popfreqs[is.nan(popfreqs)] = NA
    sumna = sumna + sum(is.na(popfreqs))
    sq = cnt:(cnt+blocksize-1)
    afmatrix[sq,] = popfreqs
    countmatrix[sq,] = popcounts
    if(verbose && cnt %% modnum == 1) alert_info(paste0((cnt-1)/1e3, 'k SNPs read...\r'))
    cnt = cnt+blocksize
  }
  if(verbose) {
    alert_info('\n')
    alert_success(paste0(cnt-1, ' SNPs read in total\n'))
    alert_warning(paste0(sumna, ' allele frequencies are missing (on average ',
                         round(sumna/numpop), ' per population)\n'))
  }

  #outdat = treat_missing(afmatrix[keepsnps,], countmatrix[keepsnps,], snpfile[keepsnps,],
  #                       na.action = na.action, verbose = verbose)
  outlist = list(afs = afmatrix, counts = countmatrix, snpfile = snpfile)
  outlist
}


#' Read allele frequencies from ancestrymap files
#'
#' @export
#' @param pref prefix of ancestrymap files (files have to end in `.geno`, `.ind`, `.snp`).
#' @param inds vector of samples from which to compute allele frequencies.
#' @param pops vector of populations from which to compute allele frequencies. If `NULL` (default), populations will be extracted from the third column in the `.ind` file
#' @param verbose print progress updates
#' @return a list with two items: allele frequency data and allele counts.
#' @examples
#' \dontrun{
#' afdat = ancestrymap_to_aftable(prefix, pops = pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
ancestrymap_to_aftable = function(pref, inds = NULL, pops = NULL, blocksize = 1000, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # pops is vector of populations for which to calculate AFs
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population

  stopifnot(is.null(pops) || is.null(inds) || length(inds) == length(pops))

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = cols(), progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = cols(), progress = FALSE)
  if(is.null(inds)) {
    inds = indfile$X1
  } else {
    indfile$X1[!indfile$X1 %in% inds] = NA
  }
  if(is.null(pops)) {
    pops = indfile$X3
    indfile$X3[!indfile$X1 %in% inds] = NA
  } else {
    indfile$X3[!is.na(indfile$X1)] = pops
  }
  upops = unique(na.omit(pops))

  popind2 = which(!is.na(indfile$X3))
  indfile %<>% filter(!is.na(X3))
  fl = paste0(pref, '.geno')
  geno = apply(do.call(rbind, str_split(readLines(fl), '')), 2, as.numeric)
  nindall = ncol(geno)
  nsnp = nrow(geno)
  numpop = length(upops)
  geno = geno[,popind2]
  colnames(geno) = inds
  rownames(geno) = snpfile$SNP

  ploidy = apply(geno, 1, function(x) max(1, length(unique(na.omit(x)))-1))
  counts = t(rowsum((!is.na(t(geno)))*ploidy, pops))
  afs = t(rowsum(t(geno)/(3-ploidy), pops, na.rm=T))/counts

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnp, ' SNPs\n'))
    alert_info(paste0('Calculating allele frequencies from ', length(popind2), ' samples in ', numpop, ' populations\n'))
  }
  outlist = namedList(afs, counts, snpfile)
  outlist
}


discard_from_aftable = function(afdat, maxmiss = 1, minmaf = 0, maxmaf = 0.5, auto_only = TRUE,
                                transitions = TRUE, transversions = TRUE, keepsnps = NULL) {
  # afdat is list with 'snpfile', 'afs', 'counts'
  # returns same list with SNPs removed
  # keepsnps overrides maxmiss and auto_only
  # maxmiss = 1 is equivalent to na.action = 'none'
  # maxmiss = 0 is equivalent to na.action = 'remove'

  snpdat = afdat$snpfile
  if(maxmiss < 1) snpdat$miss = rowMeans(afdat$counts == 0)
  else snpdat$miss = 0
  if(minmaf > 0 | maxmaf < 0.5) snpdat %<>% mutate(af = rowMeans(afdat$afs, na.rm=TRUE)/2, maf = pmin(af, 1-af))
  else snpdat %<>% mutate(af = 0.2, maf = 0.2)

  remaining = discard_snps(snpdat, maxmiss = maxmiss, auto_only = auto_only, minmaf = minmaf, maxmaf = maxmaf,
                           transitions = transitions, transversions = transversions, keepsnps = keepsnps)
  #keeprows = match(remaining, snpdat[['SNP']])
  map(afdat, ~.[remaining,])
}


discard_from_geno = function(geno, maxmiss = 0.25, auto_only = TRUE, minmaf = 0, maxmaf = 0.5,
                             transitions = TRUE, transversions = TRUE, keepsnps = NULL) {
  # afdat is list with 'snpfile', 'afs', 'counts'
  # returns same list with SNPs removed
  # keepsnps overrides maxmiss and auto_only
  # maxmiss = 1 is equivalent to na.action = 'none'
  # maxmiss = 0 is equivalent to na.action = 'remove'

  bim = ifelse('bim' %in% names(geno), 'bim', 'snp')
  bed = ifelse('bed' %in% names(geno), 'bed', 'geno')
  fam = ifelse('fam' %in% names(geno), 'fam', 'ind')

  stopifnot(all(c('SNP', 'CHR', 'A1', 'A2') %in% colnames(geno[[bim]])))
  snpdat = geno[[bim]]

  if(minmaf > 0 | maxmaf < 0.5) {
    snpdat %<>% mutate(af = rowMeans(geno[[bed]], na.rm=TRUE)/2,
                                          maf = pmin(af, 1-af))
  } else {
    snpdat %<>% mutate(af = 0.25, maf = 0.25)
  }
  if(maxmiss < 1) snpdat %<>% mutate(miss = rowMeans(is.na(geno[[bed]])))
  else snpdat %<>% mutate(miss = 0)

  allsnps = snpdat[['SNP']]
  remaining = discard_snps(snpdat, maxmiss = maxmiss, auto_only = auto_only, minmaf = minmaf, maxmaf = maxmaf,
                           transitions = transitions, transversions = transversions, keepsnps = keepsnps)
  stopifnot(length(remaining) > 0)

  #keeprows = match(remaining, allsnps)
  geno[[bed]] = geno[[bed]][remaining,]
  geno[[bim]] = geno[[bim]][remaining,]
  geno
}

discard_snps = function(snpdat, maxmiss = 0.25, keepsnps = NULL, auto_only = TRUE,
                        minmaf = 0, maxmaf = 0.5, transitions = TRUE, transversions = TRUE) {
  # input is a data frame with columns 'SNP', 'CHR', 'A1', 'A2', 'miss', 'maf'
  # output is vector of remaining row indices

  stopifnot(all(c('SNP', 'CHR', 'A1', 'A2', 'miss', 'maf') %in% colnames(snpdat)))

  snpdat %<>% mutate(.snpindex = 1:n())

  if(!is.null(keepsnps)) {
    snpdat %<>% filter(SNP %in% keepsnps)
    return(snpdat$.snpindex)
  }

  snpdat %<>% mutate(mutation = NA)
  if(!transitions || !transversions) {
    snpdat %<>%
      mutate(a1 = toupper(A1), a2 = toupper(A2),
             aa1 = ifelse(a1 < a2, a1, a2),
             aa2 = ifelse(a1 < a2, a2, a1),
             a = paste0(aa1, aa2),
             mutation = ifelse(a %in% c('AG', 'CT'), 'transition',
                                        ifelse(a %in% c('AC', 'AT', 'CG', 'GT'),
                                               'transversion', NA))) %>%
      select(-a1, -a2, -aa1, -aa2, -a)
  }
  snpdat %>%
    filter(
    miss <= maxmiss,
    between(maf, minmaf, maxmaf),
    !auto_only | as.numeric(gsub('[a-zA-Z]+', '', CHR)) <= 22,
    transitions | mutation != 'transition',
    transversions | mutation != 'transversion'
  ) %$% .snpindex
}



#' Read genotype data from packedancestrymap files
#'
#' This is currently slower than `read_plink` because it is implemented only in `R`, not in `C++`.
#' @export
#' @param inds optional vector of samples to read in
#' @param blocksize number of SNPs read in each block.
#' @inheritParams packedancestrymap_to_aftable
#' @return a list with the genotype data matrix, the `.ind` file, and the `.snp` file
#' @examples
#' \dontrun{
#' samples = c('Ind1', 'Ind2', 'Ind3')
#' geno = read_packedancestrymap(prefix, samples)
#' }
read_packedancestrymap_old = function(pref, inds = NULL, blocksize = 1000, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # inds: optional vector of individual IDs
  # returns list with geno (genotype matrix), snp (snp metadata), ind (sample metadata).

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = cols(), progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = cols(), progress = FALSE)
  indfile$X3 = indfile$X1
  if(!is.null(inds)) {
    stopifnot(all(inds %in% indfile$X1))
    indfile$X3[!indfile$X3 %in% inds] = NA
    inds = intersect(inds, indfile$X1)
  } else {
    inds = indfile$X1
  }
  fl = paste0(pref, '.geno')
  conn = file(fl, 'rb')
  on.exit(close(conn))
  hd = strsplit(readBin(conn, 'character', n = 1), ' +')[[1]]
  close(conn)
  nindall = as.numeric(hd[2])
  nsnp = as.numeric(hd[3])
  nind = length(inds)

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnp, ' SNPs.\n'))
    alert_info(paste0('Reading data for ', nind, ' samples\n'))
    alert_info(paste0('Expected size of genotype data: ', round((nsnp*nind*8+nsnp*112)/1e6), ' MB\n'))
    # 8, 112: estimated scaling factors for AF columns and annotation columns
  }

  rlen = file.info(fl)$size/(nsnp+1)
  conn = file(fl, 'rb')
  invisible(readBin(conn, 'raw', n = rlen))
  afmatrix = matrix(NA, nsnp, nind)
  colnames(afmatrix) = inds
  rownames(afmatrix) = snpfile$SNP
  popind2 = which(!is.na(indfile$X3))
  popind3 = c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`))
  popindmat = matrix(sapply(1:nind, function(i) indfile$X3[popind2] == inds[i]), ncol=nind)
  indfile %<>% filter(!is.na(X3))
  stopifnot(nrow(indfile) == nind)
  cnt = 1
  sumna = 0
  count_nonmissing = function(x) sum(!is.na(x))
  modnum = ifelse(shiny::isRunning(), 1e5, 1e3)
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
    if(verbose && cnt %% modnum == 1) alert_info(paste0((cnt-1)/1e3, 'k SNPs read...\r'))
    cnt = cnt+blocksize
  }
  if(verbose) {
    alert_info('\n')
    alert_success(paste0(cnt-1, ' SNPs read in total\n'))
    alert_warning(paste0(sumna, ' genotypes are missing (on average ', round(sumna/nind), ' per sample)\n'))
  }
  #outdat = treat_missing(afmatrix[keepsnps,], NULL, snpfile[keepsnps,], na.action = na.action, verbose = verbose)
  nr = nrow(afmatrix)
  nc = ncol(afmatrix)
  stopifnot(nr == nrow(snpfile))
  stopifnot(nc == nrow(indfile))
  outlist = list(geno = afmatrix, ind = indfile, snp = snpfile)
  outlist
}

#' Read genotype data from packedancestrymap files
#'
#' This is currently slower than `read_plink` because it is implemented only in `R`, not in `C++`.
#' @export
#' @param pref prefix of the packedancestrymap files
#' @param inds optional vector of samples to read in
#' @param first index of first SNP to read
#' @param last index of last SNP to read
#' @param transpose transpose genotype matrix
#' @return a list with the genotype data matrix, the `.ind` file, and the `.snp` file
#' @examples
#' \dontrun{
#' samples = c('Ind1', 'Ind2', 'Ind3')
#' geno = read_packedancestrymap(prefix, samples)
#' }
read_packedancestrymap = function(pref, inds = NULL, first = 1, last = Inf,
                                  transpose = FALSE, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # inds: optional vector of individual IDs
  # returns list with geno (genotype matrix), snp (snp metadata), ind (sample metadata).

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = cols(), progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = cols(), skip = first-1,
                        n_max = last-first+1, progress = FALSE)

  indfile$.keep = indfile$X1
  if(!is.null(inds)) {
    stopifnot(all(inds %in% indfile$X1))
    indfile$.keep[!indfile$.keep %in% inds] = NA
    inds = intersect(inds, indfile$X1)
  } else {
    inds = indfile$X1
  }
  indvec = 1-is.na(indfile$.keep)
  indfile %<>% filter(!is.na(.keep)) %>% select(-.keep)

  fl = paste0(pref, '.geno')
  conn = file(fl, 'rb')
  hd = strsplit(readBin(conn, 'character', n = 1), ' +')[[1]]
  close(conn)
  nindall = as.numeric(hd[2])
  nsnpall = as.numeric(hd[3])
  nind = length(inds)
  nsnp = min(last, nsnpall) - first + 1

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnpall, ' SNPs.\n'))
    alert_info(paste0('Reading data for ', nind, ' samples and ', nsnp, ' SNPs\n'))
    alert_info(paste0('Expected size of genotype data: ', round((nsnp*nind*8+nsnp*112)/1e6), ' MB\n'))
    # 8, 112: estimated scaling factors for AF columns and annotation columns
  }

  #print(c(fl, nsnpall, nindall, indvec, first = first-1,
  #        last = min(last, nsnpall), transpose = transpose, verbose = verbose))
  geno = cpp_read_packedancestrymap(fl, nsnpall, nindall, indvec, first = first-1,
                                    last = min(last, nsnpall), transpose = transpose, verbose = verbose)

  if(!transpose) {
    colnames(geno) = inds
    rownames(geno) = snpfile$SNP
  } else {
    rownames(geno) = inds
    colnames(geno) = snpfile$SNP
  }

  if(verbose) {
    alert_success(paste0(length(snpfile$SNP), ' SNPs read in total\n'))
  }

  outlist = list(geno = geno, ind = indfile, snp = snpfile)
  outlist
}



#' Read genotype data from ancestrymap files
#'
#' @export
#' @param inds optional vector of samples to read in
#' @inheritParams packedancestrymap_to_aftable
#' @return a list with the genotype data matrix, the `.ind` file, and the `.snp` file
#' @examples
#' \dontrun{
#' samples = c('Ind1', 'Ind2', 'Ind3')
#' geno = read_packedancestrymap(prefix, samples)
#' }
read_ancestrymap = function(pref, inds = NULL, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # inds: optional vector of individual IDs
  # returns list with geno (genotype matrix), snp (snp metadata), ind (sample metadata).
  # currently doesn't handle missing data

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = cols(), progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = cols(), progress = FALSE)
  indfile$X3 = indfile$X1
  if(!is.null(inds)) {
    stopifnot(all(inds %in% indfile$X1))
    indfile$X3[!indfile$X3 %in% inds] = NA
    inds = intersect(inds, indfile$X1)
  } else {
    inds = indfile$X1
  }
  popind2 = which(!is.na(indfile$X3))
  indfile %<>% filter(!is.na(X3))
  fl = paste0(pref, '.geno')
  geno = apply(do.call(rbind, str_split(readLines(fl), '')), 2, as.numeric)
  nindall = nrow(geno)
  geno = geno[,popind2]
  colnames(geno) = inds
  rownames(geno) = snpfile$SNP

  nsnp = nrow(snpfile)
  nind = length(inds)

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnp, ' SNPs.\n'))
    alert_info(paste0('Reading data for ', nind, ' samples\n'))
  }
  outlist = list(geno = geno, ind = indfile, snp = snpfile)
  outlist
}



#' Read allele frequencies from `PLINK` files
#'
#' @export
#' @param pref prefix of `PLINK` files (files have to end in `.bed`, `.bim`, `.fam`).
#' @param inds vector of samples from which to compute allele frequencies.
#' @inheritParams packedancestrymap_to_aftable
#' @return a list with two items: allele frequency data and individual counts.
#' @examples
#' \dontrun{
#' afdat = plink_to_aftable(prefix, pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
plink_to_aftable = function(pref, inds = NULL, pops = NULL, verbose = FALSE) {
  # This is based on Gad Abraham's "plink2R" package
  # Modified to return per-group allele frequencies rather than raw genotypes.

  stopifnot(is.null(pops) || is.null(inds) || length(inds) == length(pops))
  if(verbose) alert_info('Reading allele frequencies from PLINK file...\n')

  bedfile = paste0(pref, '.bed')
  famfile = paste0(pref, '.fam')
  bimfile = paste0(pref, '.bim')
  nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
  bim = read_table2(bimfile, col_names = nam, progress = FALSE, col_types = cols())
  fam = read_table2(famfile, col_names = FALSE, progress = FALSE, col_types = cols())
  if(!is.null(inds) && !is.null(pops)) {
    stopifnot(nrow(fam) == length(inds)) # note: fix this, so that this condition is not necessary
    fam$X1 = pops
    fam$X2 = inds
  }
  if(is.null(inds)) inds = fam$X2
  if(is.null(pops)) pops = fam$X1

  indvec = pop_indices(fam, pops = pops, inds = inds)
  indvec2 = which(indvec > 0)
  keepinds = unique(fam[[is.null(pops)+1]][indvec > 0])
  afs = read_plink_afs_cpp(normalizePath(bedfile), indvec, indvec2, verbose = verbose)
  afmatrix = afs[[1]]
  countmatrix = afs[[2]]
  rownames(afmatrix) = rownames(countmatrix) = bim$SNP
  colnames(afmatrix) = colnames(countmatrix) = keepinds

  #outdat = treat_missing(afmatrix, countmatrix, bim[keepsnps,], na.action = na.action, verbose = verbose)
  outlist = list(afs = afmatrix, counts = countmatrix, snpfile = bim)
  outlist
}

pop_indices = function(famdat, inds = NULL, pops = NULL) {
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


#' Read genotype data from `PLINK` files
#'
#' This is based on Gad Abraham's \code{\href{https://github.com/gabraham/plink2R}{plink2R}} package,
#' but genotypes are `m` x `n`, not `n` x `m`.
#' See \code{\href{https://www.rdocumentation.org/packages/genio}{genio}} for a dedicated `R` package for
#' reading and writing `PLINK` files.
#' @export
#' @param inds optional vector of samples to read in
#' @param blocksize number of SNPs read in each block.
#' @inheritParams packedancestrymap_to_aftable
#' @return a list with the genotype data matrix, the `.ind` file, and the `.snp` file
#' @examples
#' \dontrun{
#' samples = c('Ind1', 'Ind2', 'Ind3')
#' geno = read_packedancestrymap(prefix, samples)
#' }
read_plink = function(pref, inds = NULL, auto_only = TRUE, verbose = FALSE) {
  # This is based on Gad Abraham's "plink2R" package, but genotypes are m x n, not n x m
  stopifnot(!is.null(auto_only))
  if(verbose) alert_info('Reading PLINK files...\n')

  bedfile = paste0(pref, '.bed')
  famfile = paste0(pref, '.fam')
  bimfile = paste0(pref, '.bim')
  nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
  bim = read_table2(bimfile, col_names = nam, col_types = cols(), progress = FALSE)
  fam = read_table2(famfile, col_names = FALSE, col_types = cols(), progress = FALSE)
  if(is.null(inds)) inds = fam[[2]]
  indvec = pop_indices(fam, pops = NULL, inds = inds)
  indvec2 = which(indvec > 0)
  #indvec2 = which(fam[[2]] %in% inds)
  keepinds = fam[[2]][indvec2]

  g = 2 * read_plink_afs_cpp(normalizePath(bedfile), indvec=indvec, indvec2, verbose = verbose)[[1]]
  rownames(g) = bim$SNP
  colnames(g) = keepinds

  #outdat = treat_missing(g, NULL, bim, na.action = na.action, verbose = verbose)
  outlist = list(bed = g, fam = fam[indvec2,], bim = bim)
  outlist
}



#' Write blocked f2 estimates to disk
#'
#' This function takes a 3d array of blocked f2 estimates, splits it by population pair,
#' and writes each pair to a separate `.rds` file.
#' @export
#' @param f2_block 3d array with f2 block jackknife estimates. The first two dimensions of
#' `f2_block` have to have population names. Each file will be stored under \code{{outdir}/{pop1}/{pop2}.rds}.
#' @param outdir directory into which to write the files.
#' @param overwrite should existing `.rds` files be overwritten?
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
#' This function reads f2 block jackknife estimates which were writtend to disk by \code{\link{write_f2}}
#' and returns a 3d array of f2 block jackknife estimates.
#' @export
#' @param f2_dir directory from which to read files
#' @param pops the populations for which f2 statistics should be read. Defaults to all populations,
#' which may require a lot of memory.
#' @return a 3d array of block jackknife estimates
#' @seealso \code{\link{write_f2}}
#' @examples
#' \dontrun{
#' read_f2(f2_dir, pops = c('pop1', 'pop2', 'pop3'))
#' }
read_f2 = function(f2_dir, pops = NULL) {
  # reads f2 block jackknife rds files and returns 3d array
  # pops is vector of populations which should be read. defaults to all populations.
  remove_na = is.null(pops)
  if(is.null(pops)) {
    pops = list.dirs(f2_dir, full.names = FALSE, recursive = FALSE)
  }
  block_lengths = readRDS(paste0(f2_dir, '/block_lengths.rds'))
  numblocks = length(block_lengths)
  numpops = length(pops)
  f2_blocks = array(NA, c(numpops, numpops, numblocks))
  dimnames(f2_blocks)[[1]] = dimnames(f2_blocks)[[2]] = pops
  dimnames(f2_blocks)[[3]] = paste0('l', block_lengths)
  for(pop1 in pops) {
    for(pop2 in pops) {
      if(pop1 <= pop2) {
        fl = paste0(f2_dir, '/', pop1, '/', pop2, '.rds')
        f2 = readRDS(fl)
        f2_blocks[pop1, pop2, ] = f2_blocks[pop2, pop1, ] = f2
        if(any(is.na(f2))) warning(paste0('missing values in ', pop1, ' - ', pop2, '!'))
      }
    }
  }
  if(remove_na) {
    keep = apply(f2_blocks, 1, function(x) sum(is.na(x)) == 0)
    f2_blocks = f2_blocks[keep, keep, ]
  }
  f2_blocks
}

write_indiv = function(data, ind, outdir, overwrite = FALSE) {
  fl = paste0(outdir, '/indivs/', ind, '.rds')
  if(!file.exists(fl) | overwrite) saveRDS(data, file=fl)
}

write_pairdat2 = function(data, ind1, ind2, outdir, overwrite = FALSE) {
  i1 = pmin(ind1, ind2)
  i2 = pmax(ind1, ind2)
  fl = paste0(outdir, '/pairs/', i1, '/', i2, '.rds')
  if(!file.exists(fl) | overwrite) saveRDS(data, file = fl)
}

write_pairdat = function(aa_arr, nn_arr, outdir, overwrite = FALSE) {

  if(!dir.exists(outdir)) dir.create(paste0(outdir, '/pairs'), recursive = TRUE, showWarnings = FALSE)
  stopifnot(!is.null(dimnames(aa_arr)[[1]]) && !is.null(dimnames(aa_arr)[[2]]))
  d1 = dim(aa_arr)[1]
  d2 = dim(aa_arr)[2]
  nam1 = dimnames(aa_arr)[[1]]
  nam2 = dimnames(aa_arr)[[2]]
  for(i in seq_len(d1)) {
    for(j in seq_len(d2)) {
      ind1 = min(nam1[i], nam2[j])
      ind2 = max(nam1[i], nam2[j])
      dir = paste0(outdir, '/pairs/', ind1, '/')
      fl = paste0(dir, ind2, '.rds')
      if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
      prods = cbind(aa = aa_arr[i, j, ], nn = nn_arr[i, j, ])
      if(!file.exists(fl) | overwrite) saveRDS(prods, file=fl)
    }
  }
}


#' Write allele frequency estimates to disk
#'
#' This function takes a 3d array of f2 block jackknife estimates, splits it by population pair,
#' and writes each pair to a separate \code{.rds} file.
#' @export
#' @param afmat matrix with allele frequency estimates. Has to have column names. Data for each population
#' will be stored in a variable called \code{allele_frequencies} under \code{{path}/{pop}.RData}
#' @param path directory into which to write the files.
#' @param countmat optional matrix with allele counts for each SNP and population.
#' Dimensions have to match \code{afmat}. Data will be stored in same file (\code{{path}/{pop}.rds})
#' in the variable \code{num_individuals}
#' @param overwrite should existing \code{.rds} files be overwritten?
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
    fl = paste0(path, '/', pop, '.rds')
    if(!file.exists(fl) | overwrite) save(allele_frequencies, num_alleles, population, file=fl)
  }
}


#' Split a matrix into blocks
#'
#' This function splits a large matrix into smaller blocks with `cols_per_chunk` columns per block,
#' and saves them as `.rds` files with prefix `prefix`
#' @export
#' @param mat the matrix to be split
#' @param cols_per_chunk the number of columns per block
#' @param prefix prefix of output files
#' @seealso \code{\link{packedancestrymap_to_aftable}}, \code{\link{write_split_f2_block}}
#' @examples
#' \dontrun{
#' afdat = packedancestrymap_to_aftable('path/to/packedancestrymap_prefix', allpopulations)
#' split_mat(afdat$afs, cols_per_chunk = 20, prefix = 'afdat_split_v42.1/afs')
#' }
split_mat = function(mat, cols_per_chunk, prefix, overwrite = FALSE, verbose = TRUE) {

  dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
  npops = ncol(mat)
  starts = seq(1, npops, cols_per_chunk)
  numparts = length(starts)
  ends = c(lead(starts)[-numparts]-1, npops)
  for(i in seq_len(numparts)) {
    if(verbose) cat(paste0('\rpart ', i, ' of ', numparts))
    spl = mat[, starts[i]:ends[i]]
    fl = paste0(prefix, i, '.rds')
    if(overwrite || !file.exists(fl)) saveRDS(spl, file = fl)
  }
  if(verbose) cat('\n')
}


#' Compute f2 blocks and write them to disk
#'
#' This is intended for computing f2-statistics for a large number of populations,
#' too many to do everything in working memory.
#' It assumes that the allele frequencies have already been computed and are stored in `.rds` files,
#' split into consecutive blocks for a set of populations. This function calls \code{\link{write_f2}},
#' which takes a (sub-)chunk of pairwise f2-statistics, and writes to disk one pair at a time.
#' @export
#' @param afmatprefix prefix of the allele frequency `.rds` files created by \code{\link{split_mat}}
#' @param countmatprefix prefix of the allele frequency `.rds` files created by \code{\link{split_mat}}
#' @param outdir directory where the f2 blocks will be stored
#' @param chunk1 index of the first chunk of populations
#' @param chunk2 index of the second chunk of populations
#' @param popcounts named vector with number of samples for each population.
#' @param block_lengths vector with lengths of each jackknife block. \code{sum(block_lengths)} has to
#' match the number of SNPs.
#' @param f2_denom scaling factor applied to f2-statistics. If set to 0.278, this will be approximately equal to Fst.
#' @param verbose print progress updates
#' @seealso \code{\link{split_mat}} for creating split allele frequency data,
#' \code{\link{write_f2}} for writing split f2 block jackknife estimates
#' @examples
#' \dontrun{
#' afdat = packedancestrymap_to_aftable('path/to/packedancestrymap_prefix', allpopulations)
#' split_mat(afdat$afs, cols_per_chunk = 20, prefix = 'afdat_split_v42.1/afs')
#' split_mat(afdat$counts, cols_per_chunk = 20, prefix = 'afdat_split_v42.1/counts')
#' numchunks = 185 # this should be the number of split allele frequency files
#' block_lengths = get_block_lengths(afdat$snpfile)
#' for(j in 1:numchunks) {
#'   for(j in i:numchunks) {
#'     write_split_f2_block('afmatall_split_v42.1/afs', 'afmatall_split_v42.1/counts', 'f2blocks_v42.1/',
#'                          chunk1 = i, chunk2 = j, block_lengths)
#'     }
#'   }
#' }
#' expand_grid(i = 1:numchunks, j = 1:numchunks) %>%
#'   filter(j >= i) %$%
#'   furrr::future_map2(i, j, ~write_split_f2_block('afmatall_split_v42.1/afs', 'afmatall_split_v42.1/counts', 'f2blocks_v42.1/',
#'                                                   chunk1 = .x, chunk2 = .y, block_lengths))
#' furrr::future_map(1:numchunks, ~{i=.; map(i:numchunks, ~{
#'   write_split_f2_block('afmatall_split_v42.1/afs', 'afmatall_split_v42.1/counts', 'f2blocks_v42.1/',
#'     chunk1 = i, chunk2 = ., block_lengths)
#'   })})
write_split_f2_block = function(afmatprefix, countmatprefix, outdir, chunk1, chunk2,
                                block_lengths, verbose = TRUE) {
  # reads two afmat blocks, computes f2 jackknife blocks, and writes output to outdir

  afmat1 = readRDS(paste0(afmatprefix, chunk1, '.rds'))
  afmat2 = readRDS(paste0(afmatprefix, chunk2, '.rds'))
  countmat1 = readRDS(paste0(countmatprefix, chunk1, '.rds'))
  countmat2 = readRDS(paste0(countmatprefix, chunk2, '.rds'))
  nam1 = colnames(afmat1)
  nam2 = colnames(afmat2)
  nsnp = nrow(afmat1)
  numblocks = length(block_lengths)
  filenames = expand_grid(nam1, nam2) %>%
    transmute(nam = paste0(outdir, '/', pmin(nam1, nam2), '/', pmax(nam1, nam2), '.rds')) %$% nam
  if(all(file.exists(filenames))) return()

  f2_subblock = mats_to_f2arr(afmat1, afmat2, countmat1, countmat2) %>%
   block_arr_mean(block_lengths) %>%
    replace_nan_with_na() %>%
   `dimnames<-`(list(nam1, nam2, paste0('l', block_lengths)))
  if(chunk1 == chunk2) for(i in 1:dim(f2_subblock)[1]) f2_subblock[i,i,] = 0
  write_f2(f2_subblock, outdir = outdir)
}

#' Compute and write block lengths
#'
#' @export
#' @examples
#' \dontrun{
#' pref = 'path/to/packedancestrymap_prefix'
#' genodir = 'split_v42.1'
#' outdir = 'indpairs_v42.1'
#' write_block_lengths(pref, outdir)
#' dat = read_packedancestrymap(pref, inds = inds)
#' split_mat(dat$geno, cols_per_chunk = 20, outdir = genodir)
#' write_split_inddat(genodir, outdir)
#' numchunks = 185 # this should be the number of split allele frequency files
#' for(j in 1:numchunks) {
#'   for(j in i:numchunks) {
#'     write_split_pairdat(genodir, outdir, chunk1 = i, chunk2 = j)
#'     }
#'   }
#' }
write_block_lengths = function(pref, outdir, dist = 0.05, distcol = 'cm', auto_only = TRUE) {

  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  plinkfile = paste0(pref, '.bim')
  pafile = paste0(pref, '.snp')
  plink = file.exists(plinkfile)
  if(plink) {
    nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
    snpfile = plinkfile
  } else {
    nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
    snpfile = pafile
  }
  snps = read_table2(snpfile, col_names = nam, col_types = cols(), progress = FALSE)
  if(auto_only) snps %<>% filter(between(CHR, 1, 22))
  block_lengths = get_block_lengths(snps, dist = dist, distcol = distcol)
  saveRDS(block_lengths, file = paste0(outdir, '/block_lengths.rds'))
}


#' @export
write_split_inddat = function(genodir, outdir, overwrite = FALSE, maxmem = 8000, verbose = TRUE) {
  # reads split genotype matrices, computes a and n, and writes output to outdir

  files = list.files(genodir, '\\.rds$', full.names = TRUE)
  nfiles = length(files)
  block_lengths = readRDS(paste0(outdir, '/block_lengths.rds'))
  dir.create(paste0(outdir, '/indivs'), showWarnings = FALSE)

  if(verbose) alert_info(paste0('Writing data for ', length(files), ' files\n'))
  for(i in seq_len(nfiles)) {
    if(verbose) alert_info(paste0('file ', i, ' out of ', nfiles, '\r'))
    xmat = readRDS(files[i])
    xmat_to_inddat(xmat, block_lengths, outdir = outdir, overwrite = overwrite,
                   maxmem = maxmem, verbose = FALSE)
  }
  if(verbose) alert_info(paste0('\n'))
}


#' @export
write_split_pairdat = function(genodir, outdir, chunk1, chunk2, overwrite = FALSE) {
  # reads two (possibly split) genotype matrices, computes aa and nn, and writes output to outdir

  block_lengths = readRDS(paste0(outdir, '/block_lengths.rds'))
  xmat1 = readRDS(paste0(genodir, '/part_', chunk1, '.rds'))
  xmat2 = readRDS(paste0(genodir, '/part_', chunk2, '.rds'))
  nam1 = colnames(xmat1)
  nam2 = colnames(xmat2)

  pairnames = expand_grid(nam1, nam2) %>%
    transmute(nam = paste0(outdir, '/pairs/', pmin(nam1, nam2), '/', pmax(nam1, nam2), '.rds')) %$% nam
  if(all(file.exists(pairnames))) return()

  arrs = xmats_to_pairarrs(xmat1, xmat2)

  aa_subblock = block_arr_mean(arrs$aa, block_lengths)
  nn_subblock = block_arr_mean(arrs$nn, block_lengths)
  write_pairdat(aa_subblock, nn_subblock, outdir = outdir, overwrite = overwrite)
}



#' Compute and store blocked f2 statistics
#'
#' Prepare data for various admixtools functions. Reads data from packedancestrymap or PLINK files,
#' and computes allele frequencies and f2 block jackknife statistics for selected populations.
#' Data is either written to disk in files for each population and population pair, or returned as list.
#' This function calls \code{\link{packedancestrymap_to_aftable}} or \code{\link{plink_to_aftable}}
#' and \code{\link{afs_to_f2_blocks}}.
#' @export
#' @param pref prefix of packedancestrymap or PLINK files. packedancestrymap has to end in `.geno`, `.snp`, `.ind`,
#' PLINK has to end in `.bed`, `.bim`, `.fam`
#' @param outdir the directory to which to write data
#' @param inds the individuals for which data should be extracted.
#' @param pops the populations for which data should be extracted. `pops` and `inds` cannot be specified
#' at the same time. If none are specified, all populations will be extracted.
#' @param dist genetic distance in Morgan. Default is 0.05 (50 cM).
#' @param maxmem split up allele frequency data into blocks, if memory requirements exceed `maxmem` MB.
#' @param maxmiss discard SNPs which are missing in a fraction of populations higher than `maxmiss`
#' @param minmaf discard SNPs with minor allele frequency less than `minmaf`
#' @param maxmaf discard SNPs which minor allele frequency greater than than `maxmaf`
#' @param transitions set to `FALSE` to exclude transition SNPs
#' @param transversions set to `FALSE` to exclude transversion SNPs
#' @param keepsnps SNP IDs of SNPs to keep. overrides other SNP filtering options
#' @param overwrite should existing files be overwritten?
#' @param format supply this if the prefix can refer to genotype data in different formats
#' and you want to choose which one to read
#' @param verbose print progress updates
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' extract_f2(pref, pops = c('popA', 'popB', 'popC'))
#' }
extract_f2 = function(pref, outdir, inds = NULL, pops = NULL, dist = 0.05, maxmem = 8000,
                      maxmiss = 0.25, minmaf = 0, maxmaf = 0.5, transitions = TRUE, transversions = TRUE,
                      keepsnps = NULL, overwrite = FALSE, format = NULL, verbose = TRUE) {

  outdir = normalizePath(outdir, mustWork = FALSE)
  if(length(list.files(outdir)) > 0 && verbose) alert_danger('output directory not empty!')
  afdat = anygeno_to_aftable(pref, inds = inds, pops = pops, format = format, verbose = verbose)
  afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf,
                                  transitions = transitions, transversions = transversions,
                                  keepsnps = keepsnps, auto_only = TRUE)
  block_lengths = get_block_lengths(afdat$snpfile, dist = dist)
  afs_to_f2_blocks(afdat$afs, afdat$counts, block_lengths,
                   outdir = outdir, overwrite = overwrite,
                   maxmem = maxmem, verbose = verbose)
  #allele_counts = colMeans(afdat$counts, na.rm = TRUE)
  saveRDS(block_lengths, file = paste0(outdir, '/block_lengths.rds'))
  #saveRDS(allele_counts, file = paste0(outdir, '/allele_counts.rds'))
  if(verbose) alert_info(paste0('Data written to ', outdir, '/\n'))
}


anygeno_to_aftable = function(pref, inds = NULL, pops = NULL, format = NULL, verbose = TRUE) {

  if(is.null(format)) {
    if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) format = 'plink'
    else if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
      if(is_binfile(paste0(pref, '.geno'))) format = 'packedancestrymap'
      else format = 'ancestrymap'
    }
    else stop('Genotype files not found!')
  }
  geno_to_aftable = paste0(format, '_to_aftable')
  stopifnot(exists(geno_to_aftable))
  afdat = get(geno_to_aftable)(pref, inds = inds, pops = pops, verbose = verbose)
  afdat
}


read_anygeno = function(pref, inds = NULL, format = format, verbose = TRUE) {

  if(is.null(format)) {
    if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) format = 'plink'
    else if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
      if(is_binfile(paste0(pref, '.geno'))) format = 'packedancestrymap'
      else format = 'ancestrymap'
    }
    else stop('Genotype files not found!')
  }
  if(tolower(format) == 'packedancestrymap') {
    read_geno = function(...) {g = read_packedancestrymap(...); list(bed = g$geno, fam = g$ind, bim = g$snp)}
  } else if(tolower(format) == 'ancestrymap') {
    read_geno = function(...) {g = read_ancestrymap(...); list(bed = g$geno, fam = g$ind, bim = g$snp)}
  } else if(tolower(format) == 'plink') {
    read_geno = read_plink
  } else stop('Genotype files not found!')
  if(verbose) alert_info(paste0('Reading genotypes in ', tolower(format), 'format...\n'))
  read_geno(pref, inds, verbose = verbose)

}

#' Extract and store data needed to compute blocked f2
#'
#' Prepare data for various admixtools functions. This function reads data from packedancestrymap or PLINK files,
#' and extracts data required to compute blocked f-statistics for any sets of samples. The data consists of
#' `.rds` files with total and alternative allele counts for each individual, and products of total
#' and alternative allele counts for each pair.
#' The function calls \code{\link{packedancestrymap_to_aftable}} or \code{\link{plink_to_aftable}}
#' and \code{\link{afs_to_f2_blocks}}.
#'
#' @export
#' @param inds the individuals for which to extract data
#' @param maxmiss discard SNPs which are missing in a fraction of individuals greater than `maxmiss`
#' @param minmaf discard SNPs with minor allele frequency less than `minmaf`
#' @param maxmaf discard SNPs which minor allele frequency greater than than `maxmaf`
#' @param transitions set to `FALSE` to exclude transition SNPs
#' @param transversions set to `FALSE` to exclude transversion SNPs
#' @param keepsnps SNP IDs of SNPs to keep. overrides other SNP filtering options
#' @inheritParams extract_f2
#' @examples
#' \dontrun{
#'
#' }
extract_counts = function(pref, outdir, inds = NULL, dist = 0.05,  maxmiss = 0.25, minmaf = 0, maxmaf = 0.5,
                          transitions = TRUE, transversions = TRUE, keepsnps = NULL,
                          maxmem = 8000, overwrite = FALSE, format = NULL, verbose = TRUE) {
  dir.create(outdir, showWarnings = FALSE)
  if(length(list.files(outdir)) > 0) stop('output directory not empty!')
  bfile = paste0(outdir, '/block_lengths.rds')
  if(file.exists(bfile)) {
    # todo: snp filters
    extract_more_counts(pref, outdir, inds = inds, maxmem = maxmem,
                        overwrite = overwrite, format = format, verbose = verbose)
    return()
  }

  g = read_anygeno(pref, inds, format = format, verbose = verbose)
  g %<>% discard_from_geno(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf,
                           transitions = transitions, transversions = transversions,
                           keepsnps = keepsnps, auto_only = TRUE)
  randsnps = sample(1:nrow(g$bed), floor(nrow(g$bed)/2))
  g$bed[randsnps, ] = 2 - g$bed[randsnps, ]

  if(verbose) alert_info(paste0(nrow(g$bed), ' SNPs remain.\nDetermining SNP blocks...\n'))
  block_lengths = get_block_lengths(g$bim, dist = dist)
  saveRDS(block_lengths, file = bfile)
  xmat_to_inddat(g$bed, block_lengths,
                 outdir = outdir, overwrite = overwrite,
                 maxmem = maxmem, verbose = verbose)
  xmat_to_pairdat(g$bed, block_lengths,
                  outdir = outdir, overwrite = overwrite,
                  maxmem = maxmem, verbose = verbose)
  if(verbose) alert_info(paste0('Data written to ', outdir, '/\n'))
}

extract_more_counts = function(pref, outdir, inds = NULL,
                               maxmem = 8000, overwrite = FALSE, format = NULL, verbose = TRUE) {

  bfile = paste0(outdir, '/block_lengths.rds')
  oldinds = list.files(paste0(outdir, '/indivs/')) %>% str_replace('\\.rds', '')
  if(is.null(format)) {
    if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) {
      format = 'plink'
    } else if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
      format = 'packedancestrymap'
    }
  }
  if(format == 'plink') {
    read_geno = read_plink
  } else if(format == 'packedancestrymap') {
    read_geno = function(...) {
      g = read_packedancestrymap(...)
      list(bed = g$geno, fam = g$ind, bim = g$snp)
    }
  } else stop('Genotype files not found!')

  if(verbose) alert_info(paste0('Adding pairs to existing directory...\n'))
  block_lengths = readRDS(bfile)

  gnew = read_anygeno(pref, inds, format = format, verbose = verbose)$bed
  if(!is.null(inds) && !isTRUE(all.equal(sort(inds), sort(oldinds)))) {
    stop('Robert needs to implement adding samples in case of different geno files')
    gold = read_anygeno(pref, oldinds, format = format, verbose = verbose)$bed
    gcomb = cbind(gold, gnew)
  } else {
    gcomb = gnew
  }
  # continue here: update xmat_to_pairdat for cases where gold and gnew differ
  xmat_to_inddat(gnew, block_lengths,
                  outdir = outdir, overwrite = overwrite,
                  maxmem = maxmem, verbose = verbose)
  xmat_to_pairdat(gcomb, block_lengths,
                  outdir = outdir, overwrite = overwrite,
                  maxmem = maxmem, verbose = verbose)
  if(verbose) alert_info(paste0('Data written to ', outdir, '/\n'))
}


#' Compute and store blocked allele frequency data
#'
#' Prepare data for various admixtools functions. Reads data from packedancestrymap or PLINK files,
#' and computes allele frequencies for selected populations and stores it as `.rds` files in outdir.
#' @export
#' @param pref prefix of packedancestrymap or PLINK files. packedancestrymap has to end in `.geno`, `.snp`, `.ind`,
#' PLINK has to end in `.bed`, `.bim`, `.fam`
#' @param outdir the directory to which to write data
#' @param inds the individuals for which data should be extracted.
#' @param pops the populations for which data should be extracted. `pops` and `inds` cannot be specified
#' at the same time. If none are specified, all populations will be extracted.
#' @param dist genetic distance in Morgan. Default is 0.05 (50 cM).
#' @param cols_per_chunk number of populations per chunk. higher value will lead to fewer,
#' but more memory intensive jobs when computing f2
#' @param maxmiss discard SNPs which are missing in a fraction of populations higher than `maxmiss`
#' @param minmaf discard SNPs with minor allele frequency less than `minmaf`
#' @param maxmaf discard SNPs which minor allele frequency greater than than `maxmaf`
#' @param transitions set to `FALSE` to exclude transition SNPs
#' @param transversions set to `FALSE` to exclude transversion SNPs
#' @param keepsnps SNP IDs of SNPs to keep. overrides other SNP filtering options
#' @param format supply this if the prefix can refer to genotype data in different formats
#' and you want to choose which one to read
#' @param verbose print progress updates
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' outdir = 'dir/for/afdata/'
#' extract_afs(pref, outdir)
#' }
extract_afs = function(pref, outdir, inds = NULL, pops = NULL, dist = 0.05, cols_per_chunk = 20,
                       maxmiss = 0.25, minmaf = 0, maxmaf = 0.5, transitions = TRUE, transversions = TRUE,
                       keepsnps = NULL, format = NULL, verbose = TRUE) {

  # read data and compute allele frequencies
  afdat = anygeno_to_aftable(pref, inds = inds, pops = pops, format = format, verbose = verbose)
  afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf,
                                  transitions = transitions, transversions = transversions,
                                  keepsnps = keepsnps, auto_only = TRUE)
  if(verbose) alert_info(paste0(nrow(afdat$afs), ' SNPs remain after filtering\n'))

  # split allele frequency data into chunks and write to disk
  split_mat(afdat$afs, cols_per_chunk = cols_per_chunk, prefix = paste0(outdir, '/afs'), verbose = verbose)
  split_mat(afdat$counts, cols_per_chunk = cols_per_chunk, prefix = paste0(outdir, '/counts'), verbose = verbose)
  # compute jackknife blocks
  block_lengths = get_block_lengths(afdat$snpfile, dist = dist, distcol = 'cm')
  saveRDS(block_lengths, file = paste0(outdir, '/block_lengths.rds'))
}



#' Copy a subset of f2-statistics to a new directory
#'
#' @export
#' @param from directory with f2-statistics
#' @param to target directory
#' @param pops the populations to copy
#' @param verbose print progress updates
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' outdir = 'dir/for/afdata/'
#' extract_afs(pref, outdir)
#' }
extract_f2_subset = function(from, to, pops) {

  if(!dir.exists(to)) dir.create(to)
  file.copy(paste0(from, '/block_lengths.rds'), paste0(to, '/block_lengths.rds'))
  for(p1 in pops) {
    dr = paste0(to, '/', p1)
    if(!dir.exists(dr)) dir.create(dr)
    for(p2 in pops) {
      if(p1 <= p2) {
        fn = paste0(p1, '/', p2, '.rds')
        file.copy(paste0(from, '/', fn), paste0(to, '/', fn))
      }
    }
  }
}



#' Read f2 jackknife blocks from genotype data
#'
#' Prepare data for various admixtools functions. Reads data from packedancestrymap or PLINK files,
#' and computes allele frequencies and f2 block jackknife statistics for selected populations.
#' Data is either written to disk in files for each population and population pair, or returned as list.
#' This function calls \code{\link{packedancestrymap_to_aftable}} / \code{\link{plink_to_aftable}}
#' and \code{\link{afs_to_f2_blocks}}.
#' @export
#' @inheritParams extract_f2
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' f2_blocks = f2_from_geno(pref, pops = c('Protoss', 'Zerg', 'Terran'))
#' }
f2_from_geno = function(pref, inds = NULL, pops = NULL,
                        maxmem = 8000, format = NULL, verbose = TRUE) {

  stopifnot(is.character(pref))
  afdat = anygeno_to_aftable(pref, inds = inds, pops = pops, format = format, verbose = verbose)
  afs = afdat$afs
  counts = afdat$counts
  block_lengths = get_block_lengths(afdat$snpfile)
  f2_blocks = afs_to_f2_blocks(afs, counts, block_lengths,
                               maxmem = maxmem, verbose = verbose)
  f2_blocks
}

f2_from_genomat = function(geno, snpfile, pops) {
  # geno is numeric matrix
  # pops is vector ncol(matrix) with population labels
  stopifnot(is.matrix(geno) && is.numeric(geno))
  stopifnot(length(pops) == ncol(geno))

  counts = t(rowsum((!is.na(t(geno)))+0, pops))
  afs = t(rowsum(t(geno), pops))/counts/2
  block_lengths = get_block_lengths(snpfile)
  f2_blocks = afs_to_f2_blocks(afs, counts, block_lengths)
  f2_blocks
}



# this should not be needed in practice; used for testing
f2_from_geno_indivs = function(pref, inds = NULL, pops = NULL,
                               format = NULL, maxmem = 8000, apply_corr = TRUE, verbose = TRUE) {

  if(is.null(format)) {
    if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) format = 'plink'
    else if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) format = 'packedancestrymap'
  }
  if(tolower(format) == 'packedancestrymap') {
    read_geno = function(...) {g = read_packedancestrymap(...); list(bed = g$geno, fam = g$ind, bim = g$snp)}
  } else if(tolower(format) == 'plink') {
    read_geno = read_plink
  } else stop('Genotype files not found!')

  g = read_geno(pref, inds, verbose = verbose)

  block_lengths = get_block_lengths(g$bim)
  indivs = xmat_to_inddat(g$bed, block_lengths, maxmem = maxmem, verbose = verbose)
  pairs = xmat_to_pairdat(g$bed, block_lengths, maxmem = maxmem, verbose = verbose)

  if(is.null(inds)) inds = g[[2]][[2]]
  if(is.null(pops)) pops = inds
  poplist = tibble(ind = inds, pop = pops)

  f2_blocks = indpairs_to_f2blocks(indivs, pairs, poplist, block_lengths, apply_corr = apply_corr)
  f2_blocks
}


#' Read blocked f2 statistics from disk
#'
#' @export
#' @param dir directory with precomputed f2 statistics, or precomputed individual pair data
#' @param inds the individuals for which data should be read. Defaults to all individuals,
#' which may require a lot of memory.
#' @param pops the populations for which data should be read. Defaults to all populations,
#' which may require a lot of memory.
#' @param return_array should a 3d array or a data frame be returned
#' @param apply_corr should f2 correction factor be subtracted. mostly used for testing, and only applicable for counts
#' @param verbose print progress updates
#' @return a 3d array of f2 statistics
#' @examples
#' \dontrun{
#' dir = 'my/f2/dir/'
#' f2_blocks = f2_from_precomp(dir, pops = c('pop1', 'pop2', 'pop3'))
#' }
f2_from_precomp = function(dir, inds = NULL, pops = NULL, pops2 = NULL, return_array = TRUE,
                           apply_corr = TRUE, verbose = TRUE) {

  if(!is.null(pops) && !is.null(inds) && length(pops) != length(inds)) stop('pops and inds are not the same length!')
  indpairs = dir.exists(paste0(dir, '/indivs'))
  if(!is.null(pops2)) return(f2_from_precomp_nonsquare(dir, pops, pops2, verbose = verbose))

  if(indpairs) {
    if(is.null(inds)) {
      #if(!is.null(pops)) stop('Individual pair data and pop IDs, but no individual IDs supplied!')
      inds = pops
      if(is.null(pops)) inds = pops = list.dirs(paste0(dir, '/pairs'), full.names = FALSE, recursive = FALSE)
    }
    if(is.null(pops)) pops = inds
    poplist = tibble(ind = inds, pop = pops)
    if(verbose) alert_info(paste0('Reading precomputed data for ', nrow(poplist), ' individuals grouped into ',
                                  length(unique(pops)), ' populations...\n'))
    f2_blocks = f2_from_precomp_indivs(dir, poplist = poplist, return_array = return_array,
                                       apply_corr = apply_corr, verbose = verbose)$f2_blocks
  } else {
    if(!is.null(inds)) stop('Individual IDs supplied, but no "indivs" directory found!')
    if(is.null(pops)) pops = list.dirs(dir, full.names = FALSE, recursive = FALSE)
    if(verbose) alert_info(paste0('Reading precomputed data for ', length(pops), ' populations...\n'))
    f2_blocks = read_f2(dir, pops)
  }

  f2_blocks
}


f2_from_precomp_indivs = function(dir, poplist = NULL, return_array = TRUE, apply_corr = TRUE, verbose = TRUE) {

  remove_na = is.null(poplist)
  if(is.null(poplist)) {
    ind = pop = list.dirs(paste0(dir, '/pairs'), full.names = FALSE, recursive = FALSE)
    poplist = tibble(ind, pop)
  }
  inds = poplist$ind
  block_lengths = readRDS(paste0(dir, '/block_lengths.rds'))
  indivs = read_indivs(dir, inds, block_lengths = block_lengths)
  pairs = read_pairs(dir, inds, block_lengths = block_lengths)

  if(verbose) alert_info('Data read. Computing f2...\n')
  f2_blocks = indpairs_to_f2blocks(indivs, pairs, poplist, block_lengths,
                                   apply_corr = apply_corr, return_array = return_array)
  if(remove_na) {
    keep = apply(f2_blocks, 1, function(x) sum(is.na(x)) == 0)
    f2_blocks %<>% `[`(keep, keep, )
  }
  namedList(f2_blocks, indivs, pairs, poplist)
}


f2_from_precomp_nonsquare = function(dir, pops1, pops2, verbose = TRUE) {

  block_lengths = readRDS(paste0(dir, '/block_lengths.rds'))

  # p1 = pmin(pops1, pops2)
  # p2 = pmax(pops1, pops2)
  # map2(p1, p2, ~{
  #   fl = paste0(dir, '/', .x, '/', .y, '.rds')
  #   tibble(pop1 = .x, pop2 = .y, f2 = readRDS(fl)) %>%
  #     mutate(block = 1:n(), block_lengths)
  #   }) %>% bind_rows %>%
  #   bind_rows(rename(., pop1 = pop2, pop2 = pop1)) %>%
  #   distinct

  arr = array(NA, c(length(pops1), length(pops2), length(block_lengths)),
              list(pops1, pops2, paste0('l', block_lengths)))
  for(pop1 in pops1) {
    for(pop2 in pops2) {
      p1 = min(pop1, pop2)
      p2 = max(pop1, pop2)
      dat = readRDS(paste0(dir, '/', p1, '/', p2, '.rds'))
      arr[pop1, pop2, ] = dat
    }
  }
  arr
}



read_indivs = function(dir, inds, block_lengths) {

  bl = seq_along(block_lengths)
  indivs = list()
  for(ind in inds) {
    data = readRDS(paste0(dir, '/indivs/', ind, '.rds'))
    if(any(is.na(data))) warning(paste0('missing values in ', ind, '!'))
    indivs %<>% append(list(tibble(ind, bl) %>% bind_cols(as_tibble(data))))
  }
  indivs %>% bind_rows
}


read_pairs = function(dir, inds, inds2 = inds, block_lengths) {

  bl = seq_along(block_lengths)
  same = isTRUE(all.equal(inds, inds2))
  pairs = list()
  for(ind1 in inds) {
    for(ind2 in inds2) {
      if(ind1 <= ind2 || !same) {
        fl = paste0(dir, '/pairs/', min(ind1, ind2), '/', max(ind1, ind2), '.rds')
        prods = readRDS(fl)
        if(any(is.na(prods))) warning(paste0('missing values in ', ind1, ' - ', ind2, '!'))
        pairs %<>% append(list(tibble(ind1, ind2, bl) %>% bind_cols(as_tibble(prods))))
      }
    }
  }
  pairs %>% bind_rows
}


# currently not used because f2_blocks is recomputed from scratch each time
add_indivs = function(indivs, dir, add) {
  # updates (returns) indivs data for selected individuals
  # not tested

  stopifnot(length(intersect(indivs$ind, add)) == 0)
  block_lengths = readRDS(paste0(dir, '/block_lengths.rds'))
  bl = seq_along(block_lengths)
  inds = union(indivs$ind, add)
  for(ind in inds) {
    data = readRDS(paste0(dir, '/indivs/', ind, '.rds'))
    if(any(is.na(data))) warning(paste0('missing values in ', ind, '!'))
    newdat = tibble(ind, bl) %>% bind_cols(as_tibble(data))
    indivs %<>% bind_rows(newdat)
  }
  indivs
}


# currently not used because f2_blocks is recomputed from scratch each time
add_pairs = function(pairs, dir, add) {
  # updates (returns) pairs data for selected individuals

  stopifnot(length(intersect(pairs$ind1, add)) == 0)
  block_lengths = readRDS(paste0(dir, '/block_lengths.rds'))
  bl = seq_along(block_lengths)
  inds = union(pairs$ind1, add)

  for(ind1 in inds) {
    for(ind2 in inds) {
      if(ind1 <= ind2 && (ind1 %in% add || ind2 %in% add)) {
        fl = paste0(dir, '/pairs/', ind1, '/', ind2, '.rds')
        prods = readRDS(fl)
        if(any(is.na(prods))) warning(paste0('missing values in ', ind1, ' - ', ind2, '!'))
        newdat = tibble(ind1, ind2, bl) %>% bind_cols(as_tibble(prods))
        pairs %<>% bind_rows(newdat)
      }
    }
  }
  pairs
}


xmat_to_inddat = function(xmat, block_lengths, f2_denom = 1, maxmem = 8000,
                          outdir = NULL, overwrite = FALSE, verbose = TRUE) {
  # input is genotype matrix
  # either writes data to outdir, or returns it
  # output is either a data frame with a and n for all indivs, or nothing

  nr = nrow(xmat)
  nc = ncol(xmat)
  stopifnot(nr == sum(block_lengths))

  if(verbose) alert_info(paste0('Writing summary data for ', nc, ' samples...\n'))

  xmat %<>% fix_ploidy
  ploidy = attr(xmat, 'ploidy')
  bids = rep(seq_along(block_lengths), block_lengths)

  indivs = as_tibble(xmat) %>%
    rownames_to_column(var='SNP') %>%
    pivot_longer(-SNP, names_to = 'ind', values_to = 'a') %>%
    mutate(n = 1-is.na(a), n = n*ploidy[ind], a = replace_na(a, 0)) %>%
    add_column(bl = rep(bids, each=nc), .before = 'ind') %>%
    group_by(bl, ind) %>% summarize(a = mean(a), n = mean(n)) %>% ungroup

  if(!is.null(outdir)) {
    inddir = paste0(outdir, '/indivs')
    if(!dir.exists(inddir)) dir.create(inddir, recursive = TRUE, showWarnings = FALSE)
    indivs %>% select(-bl) %>% split(.$ind) %>% map(~as.matrix(.[,-1])) %>%
      imap(write_indiv, outdir = outdir, overwrite = overwrite)
  } else {
    return(indivs)
  }
}


xmat_to_pairdat = function(xmat, block_lengths, f2_denom = 1, maxmem = 8000,
                           outdir = NULL, overwrite = FALSE, verbose = TRUE) {
  # input is genotype matrix
  # either writes data to outdir, or returns it
  # output is either a data frame with aa and nn for all pairs, or nothing

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
    reqmem = round(mem2/1e6)
    alert_info(paste0('Computing pairwise stats for all SNPs and sample pairs requires ',
                         reqmem, ' MB RAM without splitting\n'))
    if(numsplits2 > 1) alert_info(paste0('Splitting into ', numsplits2,
                                         ' blocks of ', width, ' samples and up to ', maxmem,
                                         ' MB (', choose(numsplits2+1,2), ' block pairs)\n'))
    else alert_info(paste0('Computing without splitting since ', reqmem, ' < ', maxmem, ' (maxmem)...\n'))
  }

  bids = rep(seq_along(block_lengths), block_lengths)

  cmb = combn(0:numsplits2, 2)+(1:0)
  aa_list = nn_list = replicate(numsplits2, list())

  for(i in 1:ncol(cmb)) {
    if(numsplits2 > 1 & verbose) alert_info(paste0('sample pair block ', i, ' out of ', ncol(cmb), '\r'))
    c1 = cmb[1,i]
    c2 = cmb[2,i]
    s1 = starts[c1]:ends[c1]
    s2 = starts[c2]:ends[c2]
    a1 = xmat[,s1, drop=F]
    a2 = xmat[,s2, drop=F]

    if(!overwrite) {
      nam = sort(unique(c(colnames(a1), colnames(a2))))
      files = expand_grid(n1=nam, n2=nam) %>% filter(n1 <= n2) %$% paste0(outdir, '/pairs/', n1, '/', n2, '.rds')
      if(all(file.exists(files))) next
    }
    arrs = xmats_to_pairarrs(a1, a2)
    arrs$pp = arrs$aa/arrs$nn
    arrs$pp[!is.finite(arrs$pp)] = NA
    aa_subblock = block_arr_mean(arrs$aa, block_lengths)
    nn_subblock = block_arr_mean(arrs$nn, block_lengths)
    if(!is.null(outdir)) {
      write_pairdat(aa_subblock, nn_subblock, outdir = outdir, overwrite = overwrite)
    } else {
      aa_list[[c1]][[c2]] = aa_subblock
      aa_list[[c2]][[c1]] = aperm(aa_list[[c1]][[c2]], c(2,1,3))
      nn_list[[c1]][[c2]] = nn_subblock
      nn_list[[c2]][[c1]] = aperm(nn_list[[c1]][[c2]], c(2,1,3))
    }
  }
  if(numsplits2 > 1 & verbose) alert_info('\n')

  if(is.null(outdir)) {
    assemble_arrays = function(l)
      do.call(abind, list(lapply(l, function(x)
        do.call(abind, list(x, along=2))), along=1))
    aa_arr_full = assemble_arrays(aa_list)
    nn_arr_full = assemble_arrays(nn_list)
    dimnames(aa_arr_full) = dimnames(nn_arr_full) = list(ind1 = colnames(xmat),
                                                         ind2 = colnames(xmat),
                                                         bl = paste0('l', block_lengths))
    pairs = aa_arr_full %>% as.array %>% cubelyr::as.tbl_cube(met_name = 'aa') %>% as_tibble %>%
      mutate(nn = nn_arr_full %>% as.array %>% cubelyr::as.tbl_cube(met_name = 'nn') %>% as_tibble %$% nn)

    return(pairs)
  }
}


#' Group precomputed data
#'
#' Computing f2 statistics for populations consisting of many individuals can be slow and require a lot of memory.
#' To speed this up, this function groups individuals into populations, computes allele counts and products of
#' pairwise allele counts with all individuals and groups, and writes the data to disk. f2 statistics for a
#' combination of grouped and ungrouped precomputed data can then be read using \code{\link{f2_from_precomp}},
#' replacing individual IDs of the grouped samples with the new group labels.
#' All groupings are listed in `{dir}/groups/{groupname}.rds`
#' @export
#' @param dir directory with precomputed individual pair data
#' @param inds vector of individuals to group
#' @param pops vector of group names, either length 1, or same length as `inds`
#' @param overwrite should existing groups be overwritten?
#' @param verbose print progress updates
#' @seealso \code{\link{delete_groups}}
#' @examples
#' \dontrun{
#' dir = 'my/f2/dir/'
#' inds = c('ind1', 'ind2', 'ind3', 'ind4', 'ind5')
#' pops = c('pop_A', 'pop_A', 'pop_A', 'pop_B', 'pop_B')
#' group_samples(dir, inds, pops)
#' }
group_samples = function(dir, inds, pops, overwrite = FALSE, verbose = TRUE) {
  # groups samples and writes a, n, aa, nn files for groups as if they were individuals
  # pops has to be length one, or same as inds
  # data is written to /pairs/, /groups/, and /indivs/
  # can be deleted with delete_groups

  stopifnot(length(pops) == 1 || length(pops) == length(inds))
  tibble(ind = inds, pop = pops) %>% split(.$pop) %>%
    imap(~group_samples_onepop(dir, .x$ind, .y, overwrite = overwrite, verbose = verbose)) %>%
    invisible
}


group_samples_onepop = function(dir, inds, pop, overwrite = FALSE, verbose = TRUE) {
  # groups samples and writes a, n, aa, nn files for groups as if they were individuals
  # pop has to be length one
  # data is written to /pairs/, /groups/, and /indivs/
  # can be deleted with delete_groups

  stopifnot(length(pop) == 1)
  groups = list.files(paste0(dir, '/groups/')) %>% str_replace('\\.rds$', '')
  if(overwrite && pop %in% groups) delete_groups(dir, pop, verbose = FALSE)
  allinds = list.dirs(paste0(dir, '/pairs'), full.names = F) %>% `[`(nchar(.) > 0)
  if(pop %in% allinds && !(overwrite && pop %in% groups)) stop(paste0('Group name ', pop, ' already exists!'))

  poplist = tibble(ind = inds, pop)
  block_lengths = readRDS(paste0(dir, '/block_lengths.rds'))
  numblocks = length(block_lengths)
  indivs = read_indivs(dir, inds, block_lengths = block_lengths)

  indsums = indivs %>% left_join(poplist, by='ind') %>%
    group_by(pop, bl) %>% summarize(a = sum(a), n = sum(n)) %>% ungroup

  indsums %>% select(-bl) %>% split(.$pop) %>% map(~as.matrix(.[,-1])) %>%
    imap(write_indiv, outdir = dir, overwrite = overwrite)

  rest = setdiff(allinds, inds)
  pairs = read_pairs(dir, inds, allinds, block_lengths = block_lengths)
  pairsums = pairs %>%
    left_join(poplist %>% transmute(ind1 = ind, pop1 = pop), by = 'ind1') %>%
    left_join(tibble(ind2 = c(rest, poplist$ind), pop2 = c(rest, poplist$pop)), by = 'ind2') %>%
    group_by(pop1, pop2, bl) %>%
    summarize(aa = sum(aa), nn = sum(nn)) %>% ungroup

  dir.create(paste0(dir, '/pairs/', pop), showWarnings = FALSE)
  pairsums %>%
    mutate(p1 = pmin(pop1, pop2), p2 = pmax(pop1, pop2), pop1 = p1, pop2 = p2) %>%
    select(-bl, -p1, -p2) %>%
    split(paste(.$pop1, .$pop2, sep = '_')) %>%
    map(~write_pairdat2(data = as.matrix(select(., aa, nn)), ind1 = .$pop1[1], ind2 = .$pop2[2],
        outdir = dir, overwrite = overwrite))

  dir.create(paste0(dir, '/groups'), showWarnings = FALSE)
  poplist %>%
    split(.$pop) %>%
    map(~saveRDS(.$ind, file = paste0(dir, '/groups/', .$pop[1], '.rds')))
  if(verbose) alert_info(paste0('Grouped ', length(inds),' samples into "', pop, '"\n'))
}


#' Delete groups
#'
#' This function deletes data for groups created by `group_samples`
#' @export
#' @param dir directory with precomputed individual pair data
#' @param groups the groups to delete. defaults to all groups
#' @param verbose print progress updates
#' @return invisibly returns sample IDs in deleted groups as character vector
#' @seealso \code{\link{group_samples}}
#' @examples
#' \dontrun{
#' dir = 'my/f2/dir/'
#' inds = c('ind1', 'ind2', 'ind3', 'ind4', 'ind5')
#' pops = c('pop_A', 'pop_A', 'pop_A', 'pop_B', 'pop_B')
#' group_samples(dir, inds, pops)
#' }
delete_groups = function(dir, groups = NULL, verbose = TRUE) {
  # delte groups created by group_samples
  # defaults to all groups
  # returns sample IDs in deleted groups

  if(is.null(groups)) groups = list.files(paste0(dir, '/groups/')) %>% str_replace('\\.rds$', '')
  stopifnot(length(groups) > 0)
  allinds = list.dirs(paste0(dir, '/pairs'), full.names = FALSE) %>% `[`(nchar(.) > 0)
  for(ind1 in allinds) {
    for(ind2 in groups) {
      fl = paste0(dir, '/pairs/', min(ind1, ind2), '/', max(ind1, ind2), '.rds')
      if(file.exists(fl)) file.remove(fl)
    }
  }
  ids = unique(unlist(map(groups, ~readRDS(paste0(dir, '/groups/', ., '.rds')))))
  file.remove(paste0(dir, '/pairs/', groups))
  file.remove(paste0(dir, '/indivs/', groups, '.rds'))
  file.remove(paste0(dir, '/groups/', groups, '.rds'))
  if(verbose) alert_info('Grouped data deleted\n')
  invisible(ids)
}

#' @export
is_group = function(dir, group) file.exists(paste0(dir, '/groups/', group, '.rds'))


is_binfile = function(filename) {
  # tries to infer if filename is binary packedancestrymap file or text ancestrymap file
  conn = file(filename, 'rb')
  on.exit(close(conn))
  dat = readBin(conn, 'int', 1e3, 1, signed=F)
  max(dat) > 128 || length(unique(dat)) > 5
}

is_packedancestrymap_prefix = function(input) {
  if(!is.character(input) || length(input) > 1) return(FALSE)
  filesexist = all(file.exists(paste0(input, c('.geno', '.ind', '.snp'))))
  filesexist && is_binfile(paste0(input, c('.geno')))
}


