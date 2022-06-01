
#' Read allele frequencies from packedancestrymap files
#'
#' @export
#' @param pref Prefix of packedancestrymap files (files have to end in `.geno`, `.ind`, `.snp`)
#' @param inds Individuals from which to compute allele frequencies
#' @param pops Populations from which to compute allele frequencies. If `NULL` (default), populations will be extracted from the third column in the `.ind` file. If population labels are provided, they should have the same length as `inds`, and will be matched to them by position
#' @param adjust_pseudohaploid Genotypes of pseudohaploid samples are usually coded as `0` or `2`, even though only one allele is observed. `adjust_pseudohaploid` ensures that the observed allele count increases only by `1` for each pseudohaploid sample. If `TRUE` (default), samples that don't have any genotypes coded as `1` among the first 1000 SNPs are automatically identified as pseudohaploid. This leads to slightly more accurate estimates of f-statistics. Setting this parameter to `FALSE` is equivalent to the ADMIXTOOLS `inbreed: NO` option. Setting `adjust_pseudohaploid` to an integer `n` will check the first `n` SNPs instead of the first 1000 SNPs.
#' @param verbose Print progress updates
#' @return A list with three data frames: allele frequency data, allele counts, and SNP metadata
#' @examples
#' \dontrun{
#' afdat = packedancestrymap_to_afs(prefix, pops = pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
packedancestrymap_to_afs = function(pref, inds = NULL, pops = NULL, adjust_pseudohaploid = TRUE,
                                    first = 1, last = NULL, verbose = TRUE) {


  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # pops is vector of populations for which to calculate AFs
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population

  if(verbose) alert_info('Reading allele frequencies from packedancestrymap files...\n')

  pref %<>% normalizePath(mustWork = FALSE)
  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = 'ccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = 'ccddcc', progress = FALSE)
  nindall = nrow(indfile)
  nsnpall = as.numeric(nrow(snpfile))
  first = max(1, first)
  last = if(is.null(last)) nsnpall else min(last, nsnpall)

  ip = match_samples(indfile$X1, indfile$X3, inds, pops)
  indvec = ip$indvec - 1
  upops = ip$upops

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnpall, ' SNPs\n'))
    alert_info(paste0('Calculating allele frequencies from ', sum(indvec != -1), ' samples in ', length(upops), ' populations\n'))
    alert_info(paste0('Expected size of allele frequency data: ', round((nsnpall*length(upops)*8+nsnpall*112)/1e6), ' MB\n'))
    # 8, 112: estimated scaling factors for AF columns and annotation columns
  }

  ntest = if(is.numeric(adjust_pseudohaploid)) adjust_pseudohaploid else 1000
  if(adjust_pseudohaploid) ploidy = cpp_packedancestrymap_ploidy(paste0(pref, '.geno'), nsnpall, nindall, indvec, ntest)
  else ploidy = rep(2, nindall)
  afdat = cpp_packedancestrymap_to_afs(paste0(pref, '.geno'), nsnpall, nindall, indvec, first = first-1,
                                       last = last, ploidy = ploidy,
                                       transpose = FALSE, verbose = verbose)

  if(verbose) {
    alert_success(paste0(nrow(afdat$afs), ' SNPs read in total\n'))
  }

  colnames(afdat$afs) = colnames(afdat$counts) = upops
  rownames(afdat$afs) = rownames(afdat$counts) = snpfile$SNP[first:last]

  afdat$snpfile = snpfile
  afdat
}



eigenstrat_to_afs_old = function(pref, inds = NULL, pops = NULL, adjust_pseudohaploid = TRUE, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # pops is vector of populations for which to calculate AFs
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population

  if(verbose) alert_info('Reading allele frequencies from EIGENSTRAT files...\n')

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = 'ccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = 'ccddcc', progress = FALSE)

  ip = match_samples(indfile$X1, indfile$X3, inds, pops)
  indvec = ip$indvec
  upops = ip$upops
  popind2 = which(indvec > 0)
  pops = upops[indvec]

  fl = paste0(pref, '.geno')
  geno = apply(do.call(rbind, str_split(read_lines(fl), '')), 2, as.numeric)
  nindall = ncol(geno)
  nsnp = nrow(geno)
  numpop = length(upops)
  geno = geno[,popind2]
  geno[geno == 9] = NA
  colnames(geno) = inds
  rownames(geno) = snpfile$SNP

  if(adjust_pseudohaploid) ploidy = apply(geno, 2, function(x) max(1, length(unique(na.omit(x)))-1))
  else ploidy = rep(2, ncol(geno))
  counts = t(rowsum((!is.na(t(geno)))*ploidy, pops))[,upops]
  afs = t(rowsum(t(geno)/(3-ploidy), pops, na.rm=T))[,upops]/counts


  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnp, ' SNPs\n'))
    alert_info(paste0('Calculating allele frequencies from ', length(popind2), ' samples in ', numpop, ' populations\n'))
  }

  outlist = namedList(afs, counts, snpfile)
  outlist
}


#' Read allele frequencies from *EIGENSTRAT* files
#'
#' @export
#' @param pref Prefix of *EIGENSTRAT* files (files have to end in `.geno`, `.ind`, `.snp`)
#' @param numparts Number of parts into which the genotype file is split. Lowering this number can speed things up,
#' but will take more memory.
#' @inheritParams packedancestrymap_to_afs
#' @return A list with three data frames: allele frequency data, allele counts, and SNP metadata
#' @examples
#' \dontrun{
#' afdat = eigenstrat_to_afs(prefix, pops = pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
eigenstrat_to_afs = function(pref, inds = NULL, pops = NULL, numparts = 100,
                             adjust_pseudohaploid = TRUE, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # pops is vector of populations for which to calculate AFs
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population

  if(verbose) alert_info('Reading allele frequencies from EIGENSTRAT files...\n')

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = 'ccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = 'ccddcc', progress = FALSE)

  ip = match_samples(indfile$X1, indfile$X3, inds, pops)
  indvec = ip$indvec
  upops = ip$upops
  popind2 = which(indvec > 0)
  pops = upops[indvec]
  nindall = nrow(indfile)
  nsnp = nrow(snpfile)
  numpop = length(upops)

  fl = normalizePath(paste0(pref, '.geno'))
  if(adjust_pseudohaploid) {
    ntest = if(is.numeric(adjust_pseudohaploid)) adjust_pseudohaploid else 1000
    geno = cpp_read_eigenstrat(fl, nsnp, nindall, (indvec != 0)+0, 0, min(nsnp, ntest), FALSE, FALSE)
    ploidy = apply(geno, 2, function(x) max(1, length(unique(na.omit(x)))-1))
  } else ploidy = rep(2, nindall)

  afs = counts = matrix(NA, nsnp, length(upops))
  start = ceiling(seq(1, nsnp+1, len = numparts+1))
  end = start[-1]-1

  for(i in 1:numparts) {

    if(verbose && numparts > 1) alert_info(paste0('Reading part ', i,' of ', numparts,'...\r'))
    geno = cpp_read_eigenstrat(fl, nsnp, nindall, (indvec != 0)+0, start[i]-1, end[i], FALSE, verbose && numparts == 1)
    colnames(geno) = inds
    counts[start[i]:end[i],] = t(rowsum((!is.na(t(geno)))*ploidy, pops))[,upops]
    afs[start[i]:end[i],] = t(rowsum(t(geno)/(3-ploidy), pops, na.rm=T))[,upops]/counts[start[i]:end[i],]
  }
  afs[!is.finite(afs)] = NA
  rownames(afs) = rownames(counts) = snpfile$SNP
  colnames(afs) = colnames(counts) = upops

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnp, ' SNPs\n'))
    alert_info(paste0('Calculating allele frequencies from ', length(popind2), ' samples in ', numpop, ' populations\n'))
  }

  outlist = namedList(afs, counts, snpfile)
  outlist
}



discard_from_aftable = function(afdat, maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE, outpop = NULL, auto_only = TRUE,
                                poly_only = FALSE, transitions = TRUE, transversions = TRUE, keepsnps = NULL) {
  # afdat is list with 'snpfile', 'afs', 'counts'
  # returns same list with SNPs removed
  # keepsnps overrides maxmiss and auto_only

  snpdat = afdat$snpfile
  if(maxmiss < 1) snpdat$miss = rowMeans(afdat$counts == 0)
  else snpdat %<>% mutate(miss = 0)
  if(minmaf > 0 | maxmaf < 0.5) snpdat %<>% mutate(af = weighted_row_means(afdat$afs, afdat$counts), maf = pmin(af, 1-af))
  else snpdat %<>% mutate(af = 0.2, maf = 0.2)
  if(minac2) {
    popvec = seq_len(ncol(afdat$counts))
    if(minac2 == 2) {
      # only consider non-singleton populations
      popvec = apply(afdat$counts, 2, max) > 1
    }
    snpdat %<>% mutate(minac = apply(afdat$counts[,popvec], 1, min))
  } else snpdat %<>% mutate(minac = 2)

  if(poly_only) snpdat %<>% mutate(poly = cpp_is_polymorphic(afdat$afs))
  else snpdat %<>% mutate(poly = TRUE)

  if(!is.null(outpop)) {
    if(!outpop %in% colnames(afdat$afs)) stop("'outpop' should be a population name!")
    snpdat %<>% mutate(outgroupaf = afdat$afs[,outpop])
  } else snpdat %<>% mutate(outgroupaf = 0.5)

  remaining = discard_snps(snpdat, maxmiss = maxmiss, auto_only = auto_only, poly_only = poly_only,
                           minmaf = minmaf, maxmaf = maxmaf, minac2 = minac2,
                           transitions = transitions, transversions = transversions, keepsnps = keepsnps)
  #keeprows = match(remaining, snpdat[['SNP']])
  if(length(remaining) == 0) stop("No SNPs remain! Select fewer populations (particularly with low coverage), or use the option 'allsnps'!")
  map(afdat, ~.[remaining,,drop = FALSE])
}


discard_from_geno = function(geno, maxmiss = 0, auto_only = TRUE, poly_only = FALSE,
                             minmaf = 0, maxmaf = 0.5,
                             transitions = TRUE, transversions = TRUE, keepsnps = NULL) {
  # geno is list with 'snpfile', 'afs', 'counts'
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

  if(poly_only) snpdat %<>% mutate(poly = cpp_is_polymorphic(geno[[bed]]))
  else snpdat %<>% mutate(poly = TRUE)

  snpdat %<>% mutate(outgroupaf = 0.5)

  allsnps = snpdat[['SNP']]
  remaining = discard_snps(snpdat, maxmiss = maxmiss, auto_only = auto_only, poly_only = poly_only,
                           minmaf = minmaf, maxmaf = maxmaf,
                           transitions = transitions, transversions = transversions, keepsnps = keepsnps)
  stopifnot(length(remaining) > 0)

  #keeprows = match(remaining, allsnps)
  geno[[bed]] = geno[[bed]][remaining,]
  geno[[bim]] = geno[[bim]][remaining,]
  geno
}

discard_snps = function(snpdat, maxmiss = 1, keepsnps = NULL, auto_only = TRUE, poly_only = FALSE,
                        minmaf = 0, maxmaf = 0.5, minac2 = FALSE, transitions = TRUE, transversions = TRUE) {
  # input is a data frame with columns 'SNP', 'CHR', 'A1', 'A2', 'miss', 'maf'
  # output is vector of remaining row indices

  stopifnot(all(c('SNP', 'CHR', 'A1', 'A2', 'miss', 'maf') %in% colnames(snpdat)))

  snpdat %<>% mutate(.snpindex = seq_len(n()))

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
  if(auto_only) {
    snpdat %<>% mutate(CHR = as.numeric(gsub('[a-zA-Z]+', '', CHR)))
    if(any(is.na(snpdat$CHR))) stop("Could not parse chromosome numbers! Set 'auto_only = FALSE' to ignore chromosome labels!")
  }
  snpdat %>%
    filter(
    miss <= maxmiss,
    between(maf, minmaf, maxmaf),
    outgroupaf > 0 & outgroupaf < 1,
    !minac2 | minac > 1,
    !auto_only | CHR <= 22,
    !poly_only | poly == 1,
    transitions | mutation != 'transition',
    transversions | mutation != 'transversion'
  ) %$% .snpindex
}



#' Read genotype data from packedancestrymap files
#'
#' @export
#' @param pref Prefix of the packedancestrymap files
#' @param inds Individuals for which data should be read. Defaults to all individuals
#' @param pops Populations for which data should be read. Cannot be provided together with 'inds'
#' @param first Index of first SNP to read
#' @param last Index of last SNP to read
#' @param transpose Transpose genotype matrix (default is `snps` x `individuals`)
#' @param verbose Print progress updates
#' @return A list with the genotype data matrix, the `.ind` file, and the `.snp` file
#' @examples
#' \dontrun{
#' samples = c('Ind1', 'Ind2', 'Ind3')
#' geno = read_packedancestrymap(prefix, samples)
#' }
read_packedancestrymap = function(pref, inds = NULL, pops = NULL, first = 1, last = Inf,
                                  transpose = FALSE, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # inds: optional vector of individual IDs
  # returns list with geno (genotype matrix), snp (snp metadata), ind (sample metadata).

  pref = normalizePath(pref, mustWork = FALSE)
  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = 'ccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = 'ccddcc', skip = first-1,
                        n_max = last-first+1, progress = FALSE)
  if(!is.null(pops)) {
    if(!is.null(inds)) stop("'inds' and 'pops' cannot both be provided!")
    inds = indfile %>% filter(X3 %in% pops) %>% pull(X1)
  }

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



#' Read genotype data from *EIGENSTRAT* files
#'
#' @export
#' @inheritParams read_packedancestrymap
#' @return A list with the genotype data matrix, the `.ind` file, and the `.snp` file
#' @examples
#' \dontrun{
#' samples = c('Ind1', 'Ind2', 'Ind3')
#' geno = read_packedancestrymap(prefix, samples)
#' }
read_eigenstrat = function(pref, inds = NULL, pops = NULL, first = 1, last = Inf, transpose = FALSE, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # inds: optional vector of individual IDs
  # returns list with geno (genotype matrix), snp (snp metadata), ind (sample metadata).

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = 'ccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = 'ccddcc', progress = FALSE)
  if(!is.null(pops)) {
    if(!is.null(inds)) stop("'inds' and 'pops' cannot both be provided!")
    inds = indfile %>% filter(X3 %in% pops) %>% pull(X1)
  }

  nsnpall = nrow(snpfile)
  nindall = nrow(indfile)
  indfile$X3 = indfile$X1
  if(!is.null(inds)) {
    stopifnot(all(inds %in% indfile$X1))
    indfile$X3[!indfile$X1 %in% inds] = NA
  } else {
    inds = indfile$X1
  }
  indvec = (!is.na(indfile$X3))+0
  indfile %<>% filter(!is.na(X3))
  nind = nrow(indfile)
  if(!is.finite(last)) last = nsnpall
  nsnp = last - first + 1
  snpfile %<>% slice(first:last)

  if(verbose) {
    alert_info(paste0('Reading data for ', nind, ' samples and ', nsnp, ' SNPs...\n'))
  }
  geno = cpp_read_eigenstrat(paste0(pref, '.geno'), nsnpall, nindall, indvec, first = first-1, last,
                             transpose = transpose, verbose = verbose)

  colnames(geno) = if(transpose) snpfile$SNP else inds
  rownames(geno) = if(transpose) inds else snpfile$SNP

  outlist = list(geno = geno, ind = indfile, snp = snpfile)
  outlist
}

# this exists just so there is an equivalent to cpp_read_plink and cpp_read_packedancestrymap
# cpp_read_eigenstrat = function(genofile, nsnp, nind, indvec, first, last, transpose = FALSE, verbose = TRUE) {
#
#   #geno = apply(do.call(rbind, str_split(readLines(genofile, last)[(first+1):last], ''))[,indvec,drop=FALSE], 2, as.numeric)
#
#   geno = genofile %>%
#     read_lines(first, n_max = last-first, progress = FALSE) %>%
#     str_split('') %>%
#     do.call(rbind, .) %>%
#     `[`(,which(indvec==1),drop=FALSE) %>%
#     apply(2, as.numeric) %>%
#     replace(. == 9, NA)
#   if(transpose) geno %<>% t
#   geno
# }

# this exists just so there is an equivalent to cpp_read_plink and cpp_read_packedancestrymap
cpp_eigenstrat_to_afs = function(genofile, nsnp, nind, indvec, first, last, ploidy, transpose, verbose) {
  geno = cpp_read_eigenstrat(genofile, nsnp = nsnp, nind = nind, indvec = (indvec > -1)+0 , first = first, last = last,
                             transpose = FALSE, verbose = verbose)
  pops = indvec[indvec > -1]
  counts = t(rowsum((!is.na(t(geno)))*ploidy, pops))
  afs = t(rowsum(t(geno)/(3-ploidy), pops, na.rm=T))/counts
  namedList(afs, counts)
}


eigenstrat_ploidy = function(genofile, nsnp, nind, indvec, ntest = 1000) {
  # assumes -1 in indvec is do-not-use, to be consistent with cpp_ploidy functions

  geno = apply(do.call(rbind, str_split(readLines(genofile, ntest), '')), 2, as.numeric)[,indvec != -1]
  geno[geno == 9] = NA
  apply(geno, 2, function(x) length(unique(na.omit(x)))-1)
}


#' Read allele frequencies from `PLINK` files
#'
#' @export
#' @param pref prefix of `PLINK` files (files have to end in `.bed`, `.bim`, `.fam`).
#' @param numblocks Number of blocks in which to read genotype file. Setting this to a number
#' greater than one is more memory efficient, but slower.
#' @param poly_only Only keep SNPs with mean allele frequency not equal to 0 or 1.
#' @inheritParams packedancestrymap_to_afs
#' @return A list with three items: Allele frequency matrix, allele count matrix, and SNP meta data.
#' @examples
#' \dontrun{
#' afdat = plink_to_afs(prefix, pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
plink_to_afs = function(pref, inds = NULL, pops = NULL, adjust_pseudohaploid = TRUE,
                        first = 1, last = NULL, numblocks = 1, poly_only = FALSE, verbose = TRUE) {
  # This is based on Gad Abraham's "plink2R" package
  # Modified to return per-group allele frequencies rather than raw genotypes.

  if(verbose) alert_info('Reading allele frequencies from PLINK files...\n')

  bedfile = paste0(pref, '.bed')
  famfile = paste0(pref, '.fam')
  bimfile = paste0(pref, '.bim')
  nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
  bim = read_table2(bimfile, col_names = nam, progress = FALSE, col_types = 'ccddcc')
  fam = read_table2(famfile, col_names = FALSE, progress = FALSE, col_types = 'cccccc')
  nsnpall = nrow(bim)
  nindall = nrow(fam)
  first = max(0, first)
  last = if(is.null(last)) nsnpall else min(last, nsnpall)

  ip = match_samples(fam$X2, fam$X1, inds, pops)
  indvec = ip$indvec
  upops = ip$upops

  indvec2 = which(indvec > 0)
  #keepinds = unique(fam[[is.null(pops)+1]][indvec > 0])

  if(verbose) {
    alert_info(paste0(basename(pref), '.bed has ', nindall, ' samples and ', nsnpall, ' SNPs\n'))
    alert_info(paste0('Calculating allele frequencies from ', sum(indvec != 0), ' samples in ', length(upops), ' populations\n'))
    alert_info(paste0('Expected size of allele frequency data: ', 2*round((nsnpall*length(upops)*8+nsnpall*112)/1e6), ' MB\n'))
    # 8, 112: estimated scaling factors for AF columns and annotation columns
  }

  afmatrix = countmatrix = matrix(NA, last - first + 1, length(upops))
  aflist = countlist = list()
  firsts = round(seq(first, last, length = numblocks))
  lasts = c(firsts[-1]-1, last)
  snp_indices = c()
  ntest = if(is.numeric(adjust_pseudohaploid)) adjust_pseudohaploid else 1000
  if(adjust_pseudohaploid) ploidy = cpp_plink_ploidy(normalizePath(bedfile), nsnpall, nindall, indvec, ntest)
  else ploidy = rep(2, nindall)

  for(i in seq_len(numblocks)) {
    if(numblocks > 1) alert_info(paste0('Reading block ', i,' of ', numblocks,'...\r'))
    dat = cpp_plink_to_afs(normalizePath(bedfile), nsnpall, nindall, indvec-1,
                           first = firsts[i]-1, last = lasts[i], ploidy, FALSE, verbose && numblocks == 1)
    if(poly_only) {
      keep = !rowMeans(dat$afs, na.rm=TRUE) %in% 0:1
      dat$afs %<>% `[`(keep,)
      dat$counts %<>% `[`(keep,)
      snp_indices %<>% c((firsts[i]:lasts[i])[keep])
    } else {
      snp_indices %<>% c((firsts[i]:lasts[i]))
    }
    aflist[[i]] = dat$afs
    countlist[[i]] = dat$counts
  }
  bim %<>% slice(snp_indices)

  afmatrix = do.call(rbind, aflist)
  countmatrix = do.call(rbind, countlist)
  rownames(afmatrix) = rownames(countmatrix) = bim$SNP
  colnames(afmatrix) = colnames(countmatrix) = upops

  if(verbose) {
    alert_success(paste0(nrow(afmatrix), ' SNPs read in total\n'))
  }

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


match_samples = function(haveinds, havepops, inds, pops) {
  # takes individuals in file, populations in file, requested individuals, requested populations
  # returns list with two items:
  # integer indvec with population number (0 means do not use, same length as haveinds/havepops)
  # unique population labels, with position corresponding to indvec and first occurence in pops

  if(!is.null(inds)) {
    if(any(duplicated(inds))) stop("Individual IDs are duplicated!")
    if(!all(inds %in% haveinds)) stop(paste0("Individuals missing in indfile:\n", paste(setdiff(inds, haveinds), collapse=', ')))
  } else if(!is.null(pops)) {
    if(!all(pops %in% havepops)) stop(paste0("Populations missing in indfile:\n", paste(setdiff(pops, havepops), collapse=', ')))
  }

  if(is.null(inds) && is.null(pops)) {
    pops = havepops
  } else if(!is.null(inds) && !is.null(pops)) {
    if(length(inds) != length(pops)) stop("'inds' and 'pops' should have the same length!")
    haveinds[!haveinds %in% inds] = NA
    havepops = NA
    havepops[!is.na(haveinds)] = pops[match(na.omit(haveinds), inds)]
  } else if(is.null(inds) && !is.null(pops)) {
    havepops[!havepops %in% pops] = NA
  } else {
    havepops[!haveinds %in% inds] = NA
    pops = havepops
  }

  upops = unique(na.omit(pops))

  indvec = as.numeric(factor(havepops, levels = upops))
  indvec[is.na(indvec)] = 0

  namedList(indvec, upops)
}



#' Read genotype data from `PLINK` files
#'
#' See \href{https://www.rdocumentation.org/packages/genio}{genio} for a dedicated `R` package for
#' reading and writing `PLINK` files. This function is based on a similar function in the `plink2R` package.
#' @export
#' @param inds Individuals for which data should be read. Defaults to all individuals
#' @param pops Populations for which data should be read. Cannot be provided together with 'inds'
#' @inheritParams packedancestrymap_to_afs
#' @return A list with the genotype data matrix, the `.ind` file, and the `.snp` file
#' @examples
#' \dontrun{
#' samples = c('Ind1', 'Ind2', 'Ind3')
#' geno = read_packedancestrymap(prefix, samples)
#' }
read_plink = function(pref, inds = NULL, pops = NULL, verbose = FALSE) {
  # This is based on Gad Abraham's "plink2R" package, but genotypes are m x n, not n x m
  if(verbose) alert_info('Reading PLINK files...\n')

  bedfile = paste0(pref, '.bed')
  famfile = paste0(pref, '.fam')
  bimfile = paste0(pref, '.bim')
  nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
  bim = read_table2(bimfile, col_names = nam, col_types = 'ccddcc', progress = FALSE)
  fam = read_table2(famfile, col_names = FALSE, col_types = 'cccccc', progress = FALSE)

  if(!is.null(pops)) {
    if(!is.null(inds)) stop("'inds' and 'pops' cannot both be provided!")
    inds = fam %>% filter(X1 %in% pops) %>% pull(X2)
  }
  if(is.null(inds)) inds = fam[[2]]
  indvec = pop_indices(fam, pops = NULL, inds = inds)
  indvec2 = which(indvec > 0)
  #indvec2 = which(fam[[2]] %in% inds)
  keepinds = fam[[2]][indvec2]

  # g = 2 * cpp_read_plink_afs(normalizePath(bedfile), indvec=indvec, indvec2,
  #                            adjust_pseudohaploid = TRUE, verbose = verbose)[[1]]
  g = cpp_read_plink(normalizePath(bedfile), nrow(bim), nrow(fam), (indvec > 0)+0, 0, nrow(bim), F, verbose)
  rownames(g) = bim$SNP
  colnames(g) = keepinds

  outlist = list(bed = g, fam = fam[indvec2,], bim = bim)
  outlist
}



#' Write blocked f2 estimates to disk
#'
#' This function takes 3d arrays of blocked f2, allele frequency products, and counts, splits them by population pair,
#' and writes each pair to a separate `.rds` file under \code{{outdir}/{pop1}/{pop2}.rds}.
#' @export
#' @param est_arr A 3d array with blocked f2, allele frequency products, or fst estimates for each population pair.
#' The first two dimensions of each array have to have population names.
#' @param count_arr A 3d array of the same dimension with counts
#' @param id Postfix showing the type of statistic ("f2", "ap", or "fst")
#' @param outdir Directory where data will be stored
#' @param overwrite Overwrite existing files in `outdir`
#' @seealso \code{\link{read_f2}}
#' @examples
#' \dontrun{
#' write_f2(f2_arr, count_arr, outdir = 'path/to/f2stats/')
#' }
write_f2 = function(est_arr, count_arr, outdir, id = 'f2', overwrite = FALSE) {

  if(!dir.exists(outdir)) dir.create(outdir)
  d1 = dim(est_arr)[1]
  d2 = dim(est_arr)[2]
  nam1 = dimnames(est_arr)[[1]]
  nam2 = dimnames(est_arr)[[2]]
  for(i in seq_len(d1)) {
    for(j in seq_len(d2)) {
      min1 = stringi::stri_cmp_le(nam1[i], nam2[j], locale = 'C')
      pop1 = if(min1) nam1[i] else nam2[j]
      pop2 = if(min1) nam2[j] else nam1[i]
      # pop1 = min(nam1[i], nam2[j])
      # pop2 = max(nam1[i], nam2[j])
      #if(pop1 <= pop2) {
        mat = cbind(as.vector(est_arr[i, j, ]), as.vector(count_arr[i, j, ]))
        colnames(mat) = c(id, 'counts')
        dir = paste0(outdir, '/', pop1, '/')
        fl = paste0(dir, pop2, '_', id,'.rds')
        if(!dir.exists(dir)) dir.create(dir)
        if(!file.exists(fl) || overwrite) saveRDS(mat, file = fl)
      #}
    }
  }
}


#' Read blocked f2 estimates from disk
#'
#' This function reads blocked f2 estimates (or allele frequency products) which were writtend to disk by \code{\link{write_f2}}
#' and returns them as a 3d array.
#' @export
#' @param f2_dir Directory from which to read files
#' @param pops Populations for which f2 statistics should be read. Defaults to all populations,
#' which may require a lot of memory.
#' @param pops2 Specify this if you only want to read a subset of all population pairs. The resulting array will differ on 1st and 2nd dimension and will not work with all functions.
#' @param afprod Return allele frequency products instead of f2 estimates
#' @param counts Return allele counts instead of f2 estimates
#' @param remove_na Remove blocks with missing values
#' @param verbose Print progress updates
#' @return A 3d array of block jackknife estimates
#' @seealso \code{\link{write_f2}}
#' @examples
#' \dontrun{
#' read_f2(f2_dir, pops = c('pop1', 'pop2', 'pop3'))
#' }
read_f2 = function(f2_dir, pops = NULL, pops2 = NULL, type = 'f2',
                   counts = FALSE, remove_na = TRUE, verbose = FALSE) {
  # assumes f2 is in first column, afprod in second column

  if(is.null(pops)) pops = list.dirs(f2_dir, full.names = FALSE, recursive = FALSE)
  if(is.null(pops2)) pops2 = pops
  stopifnot(!any(duplicated(pops)))
  stopifnot(!any(duplicated(pops2)))

  fl = paste0(f2_dir, '/block_lengths_',type,'.rds')
  if(!file.exists(fl)) stop('block_lengths file not found. Please run extract_f2() again, as the file format has recently changed.')
  block_lengths = readRDS(fl)

  f2_blocks = array(NA, c(length(pops), length(pops2), length(block_lengths)),
              list(pops, pops2, paste0('l', block_lengths)))

  popcomb = expand_grid(pops, pops2) %>%
    #mutate(p1 = pmin(pops, pops2), p2 = pmax(pops, pops2)) %>%
    mutate(min1 = stringi::stri_cmp_lt(pops, pops2, locale = 'C'),
           p1 = ifelse(min1, pops, pops2), p2 = ifelse(min1, pops2, pops)) %>%
    select(-min1) %>%
    add_count(p1, p2) %>%
    filter(!duplicated(paste(p1, p2)))

  col = if(counts) 2 else 1
  for(i in seq_len(nrow(popcomb))) {
    pop1 = popcomb$pops[i]
    pop2 = popcomb$pops2[i]
    if(verbose) alert_info(paste0('Reading ', type,
                                  ' data for pair ', i, ' out of ', nrow(popcomb),'...\r'))
    fl = paste0(f2_dir, '/', popcomb$p1[i], '/', popcomb$p2[i], '_', type, '.rds')
    if(!file.exists(fl)) stop(paste0('File ', fl, ' not found! You may have to recompute the f-statistics!'))
    dat = readRDS(fl)[,col]
    f2_blocks[pop1, pop2, ] = dat
    if(popcomb$n[i] == 2) f2_blocks[pop2, pop1, ] = dat
    if(type == 'fst' && pop1 == pop2) f2_blocks[pop1, pop1, ] = 0
    #if(any(is.na(dat))) warning(paste0('missing values in ', pop1, ' - ', pop2, '!'))
  }
  if(verbose) alert_info(paste0('\n'))
  if(remove_na || verbose) {
    keep = apply(f2_blocks, 3, function(x) sum(!is.finite(x)) == 0)
    if(!all(keep)) {
      if(remove_na) {
        if(sum(keep) == 0) stop("No blocks remain after discarding blocks with missing values! If you want to compute FST for pseudohaploid samples, set `adjust_pseudohaploid = FALSE`")
        warning(paste0('Discarding ', sum(!keep), ' block(s) due to missing values!\n',
                       'Discarded block(s): ', paste0(which(!keep), collapse = ', ')))
        f2_blocks = f2_blocks[,, keep, drop = FALSE]
      } else {
        warning(paste0('The following block(s) have missing values: ', paste0(which(!keep), collapse = ', ')))
      }
    }
  }
  if(counts) f2_blocks = f2_blocks * rep(block_lengths, each = prod(dim(f2_blocks)[1:2]))
  f2_blocks
}



write_indiv = function(data, ind, outdir, overwrite = FALSE) {
  fl = paste0(outdir, '/indivs/', ind, '.rds')
  if(!file.exists(fl) | overwrite) saveRDS(data, file=fl)
}

write_pairdat2 = function(data, ind1, ind2, outdir, overwrite = FALSE) {
  min1 = stringi::stri_cmp_lt(ind1, ind2, locale = 'C')
  i1 = ifelse(min1, ind1, ind2)
  i2 = ifelse(min1, ind2, ind1)
  # i1 = pmin(ind1, ind2)
  # i2 = pmax(ind1, ind2)
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




#' Split a matrix into blocks
#'
#' This function splits a large matrix into smaller blocks with `cols_per_chunk` columns per block,
#' and saves them as `.rds` files with prefix `prefix`
#' @export
#' @param mat The matrix to be split
#' @param cols_per_chunk Number of columns per block
#' @param prefix Prefix of output files
#' @param overwrite Overwrite existing files (default `TRUE`)
#' @param verbose Print progress updates
#' @seealso \code{\link{packedancestrymap_to_afs}}, \code{\link{afs_to_f2}}
#' @examples
#' \dontrun{
#' afdat = packedancestrymap_to_afs('path/to/packedancestrymap_prefix', allpopulations)
#' split_mat(afdat$afs, cols_per_chunk = 20, prefix = 'afdat_split_v42.1/afs')
#' }
split_mat = function(mat, cols_per_chunk, prefix, overwrite = TRUE, verbose = TRUE) {

  dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
  npops = ncol(mat)
  starts = seq(1, npops, cols_per_chunk)
  numparts = length(starts)
  ends = c(lead(starts)[-numparts]-1, npops)
  for(i in seq_len(numparts)) {
    if(verbose) cat(paste0('\rpart ', i, ' of ', numparts))
    spl = mat[, starts[i]:ends[i],drop=F]
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
#' which takes a (sub-)chunk of pairwise f2-statistics, and writes one pair at a time to disk.
#' @export
#' @param afdir Directory with allele frequency and counts `.rds` files created by \code{\link{split_mat}}
#' @param outdir Directory where data will be stored
#' @param chunk1 Index of the first chunk of populations
#' @param chunk2 Index of the second chunk of populations
#' @param blgsize SNP block size in Morgan. Default is 0.05 (50 cM). If `blgsize` is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
#' @param snpwt A vector of scaling factors applied to the f2-statistics for each SNP. The length has to match the number of SNPs.
#' @param overwrite Overwrite existing files (default `FALSE`)
#' @param verbose Print progress updates
#' @seealso \code{\link{extract_f2}} Does the same thing in one step for smaller data.
#' @examples
#' \dontrun{
#' afdir = 'tmp_af_dir/'
#' f2dir = 'f2_dir'
#' extract_afs('path/to/packedancestrymap_prefix', afdir)
#' numchunks = length(list.files(afdir, 'afs.+rds'))
#' # numchunks should be the number of split allele frequency files
#' for(i in 1:numchunks) {
#'   for(j in i:numchunks) {
#'     afs_to_f2(afdir, f2dir, chunk1 = i, chunk2 = j)
#'   }
#' }
#' }
#' # Alternatively, the following code will do the same, while submitting each chunk as a separate job.
#' # (if \code{\link[future]{plan}} has been set up appropriately)
#' \dontrun{
#' furrr::future_map(1:numchunks, ~{i=.; map(i:numchunks, ~{
#'   afs_to_f2(afdir, f2dir, chunk1 = i, chunk2 = .)
#'   })})
#'   }
afs_to_f2 = function(afdir, outdir, chunk1, chunk2, blgsize = 0.05, snpwt = NULL, overwrite = FALSE,
                     type = 'f2', poly_only = FALSE, snpdat = NULL, apply_corr = TRUE, verbose = TRUE) {
  # reads data from afdir, computes f2 jackknife blocks, and writes output to outdir

  if(is.null(snpdat)) {
    fl = paste0(afdir, '/snpdat.tsv.gz')
    nc = ncol(read_table2(fl, n_max = 0, col_types=cols()))
    snpdat = read_table2(fl, col_types = paste0('ccddcc', paste0(rep('?', nc-6), collapse='')), progress = FALSE)
  }

  am1 = readRDS(paste0(afdir, '/afs', chunk1, '.rds'))
  am2 = readRDS(paste0(afdir, '/afs', chunk2, '.rds'))
  cm1 = readRDS(paste0(afdir, '/counts', chunk1, '.rds'))
  cm2 = readRDS(paste0(afdir, '/counts', chunk2, '.rds'))
  nam1 = colnames(am1)
  nam2 = colnames(am2)
  nsnp = nrow(am1)

  if(poly_only) {
    snps = snpdat$poly
    am1 = am1[snps,,drop=F]
    am2 = am2[snps,,drop=F]
    cm1 = cm1[snps,,drop=F]
    cm2 = cm2[snps,,drop=F]
  } else snps = seq_len(nsnp)

  fl = paste0(outdir, '/block_lengths_',type,'.rds')
  if(!file.exists(fl)) {
    block_lengths = get_block_lengths(snpdat[snps,], blgsize = blgsize)
    saveRDS(block_lengths, file = fl)
  } else {
    block_lengths = readRDS(fl)
  }

  filenames = expand_grid(nam1, nam2) %>%
    mutate(min1 = stringi::stri_cmp_lt(nam1, nam2, locale = 'C'),
           p1 = ifelse(min1, nam1, nam2), p2 = ifelse(min1, nam2, nam1)) %>%
    select(-min1) %>%
    transmute(nam = paste0(outdir, '/', p1, '/', p2)) %$%
    nam %>% rep(each = 2) %>% paste0('_', type, '.rds')
  if(all(file.exists(filenames)) && !overwrite) return()

  countsap = aparr = fstarr = NULL

  fun = get(paste0('mats_to_', type, 'arr'))
  if(sum(block_lengths) != nrow(am1)) stop("block_lengths and am1 don't match!")

  arr = fun(am1, am2, cm1, cm2, block_lengths, snpwt, apply_corr = apply_corr)
  counts = mats_to_ctarr(am1, am2, cm1, cm2, block_lengths)
  if(chunk1 == chunk2) for(i in 1:dim(arr)[1]) arr[i, i, ] = 0
  write_f2(arr, counts, outdir = outdir, id = type, overwrite = overwrite)
}



write_split_inddat = function(genodir, outdir, overwrite = FALSE, maxmem = 8000, verbose = TRUE) {
  # reads split genotype matrices, computes a and n, and writes output to outdir

  files = list.files(genodir, '^afs.\\.rds$', full.names = TRUE)
  nfiles = length(files)
  block_lengths = readRDS(paste0(outdir, '/block_lengths.rds'))
  dir.create(paste0(outdir, '/indivs'), showWarnings = FALSE)

  if(verbose) alert_info(paste0('Writing individual data from ', length(files), ' chunks\n'))
  for(i in seq_len(nfiles)) {
    if(verbose) alert_info(paste0('Chunk ', i, ' out of ', nfiles, '\r'))
    xmat = 2*readRDS(files[i])
    xmat_to_inddat(xmat, block_lengths, outdir = outdir, overwrite = overwrite,
                   maxmem = maxmem, verbose = FALSE)
  }
  if(verbose) alert_info(paste0('\n'))
}

#' Compute count blocks and write them to disk
#'
#' This is intended for computing allele count data for a large number of individuals,
#' too many to do everything in working memory.
#' It assumes that the allele frequencies have already been computed and are stored in `.rds` files,
#' split into consecutive blocks for a set of individuals.
#' @export
#' @param genodir Directory with genotype files split in chunks created by \code{\link{extract_afs}} (with each individual its own population).
#' @param outdir Directory where allele count data will be stored
#' @param chunk1 Index of the first chunk of individuals
#' @param chunk2 Index of the second chunk of individuals
#' @param overwrite Overwrite existing files (default `FALSE`)
#' @param verbose Print progress updates
#' @seealso \code{\link{extract_counts}} Does the same thing in one step for smaller data.
afs_to_counts = function(genodir, outdir, chunk1, chunk2, overwrite = FALSE, verbose = TRUE) {
  # reads two (possibly split) genotype matrices, computes aa and nn, and writes output to outdir

  block_lengths = readRDS(paste0(outdir, '/block_lengths.rds'))
  xmat1 = 2*readRDS(paste0(genodir, '/afs', chunk1, '.rds'))
  xmat2 = 2*readRDS(paste0(genodir, '/afs', chunk2, '.rds'))
  nam1 = colnames(xmat1)
  nam2 = colnames(xmat2)

  pairnames = expand_grid(nam1, nam2) %>%
    mutate(min1 = stringi::stri_cmp_lt(nam1, nam2, locale = 'C'),
           p1 = ifelse(min1, nam1, nam2), p2 = ifelse(min1, nam2, nam1)) %>%
    select(-min1) %>%
    transmute(nam = paste0(outdir, '/pairs/', p1, '/', p2, '.rds')) %$% nam
  if(all(file.exists(pairnames))) return()

  arrs = xmats_to_pairarrs(xmat1, xmat2)
  arrs2 = xmats_to_pairarrs(2-xmat1, 2-xmat2)
  aa_subblock = (block_arr_mean(arrs$aa, block_lengths) + block_arr_mean(arrs2$aa, block_lengths))/2
  nn_subblock = (block_arr_mean(arrs$nn, block_lengths) + block_arr_mean(arrs2$nn, block_lengths))/2

  write_pairdat(aa_subblock, nn_subblock, outdir = outdir, overwrite = overwrite)
}



#' Compute and store blocked f2 statistics
#'
#' This function prepares data for various other *ADMIXTOOLS 2* functions. It reads data from genotype files,
#' computes allele frequencies and blocked f2-statistics for selected populations, and writes the results to `outdir`.
#' @export
#' @param pref Prefix of *PLINK/EIGENSTRAT/PACKEDANCESTRYMAP* files.
#' *EIGENSTRAT/PACKEDANCESTRYMAP* have to end in `.geno`, `.snp`, `.ind`, *PLINK* has to end in `.bed`, `.bim`, `.fam`
#' @param outdir Directory where data will be stored.
#' @param inds Individuals for which data should be extracted
#' @param pops Populations for which data should be extracted. If both `pops` and `inds` are provided, they should have the same length and will be matched by position. If only `pops` is provided, all individuals from the `.ind` or `.fam` file in those populations will be extracted. If only `inds` is provided, each indivdual will be assigned to its own population of the same name. If neither `pops` nor `inds` is provided, all individuals and populations in the `.ind` or `.fam` file will be extracted.
#' @param blgsize SNP block size in Morgan. Default is 0.05 (50 cM). If `blgsize` is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
#' @param maxmem Maximum amount of memory to be used. If the required amount of memory exceeds `maxmem`, allele frequency data will be split into blocks, and the computation will be performed separately on each block pair.
#' @param maxmiss Discard SNPs which are missing in a fraction of populations higher than `maxmiss`
#' @param minmaf Discard SNPs with minor allele frequency less than `minmaf`
#' @param maxmaf Discard SNPs with minor allele frequency greater than than `maxmaf`
#' @param minac2 Discard SNPs with allele count lower than 2 in any population (default `FALSE`). This option should be set to `TRUE` when computing f3-statistics where one population consists mostly of pseudohaploid samples. Otherwise heterozygosity estimates and thus f3-estimates can be biased. `minac2 == 2` will discard SNPs with allele count lower than 2 in any non-singleton population (this option is experimental and is based on the hypothesis that using SNPs with allele count lower than 2 only leads to biases in non-singleton populations). While the `minac2` option discards SNPs with allele count lower than 2 in any population, the \code{\link{qp3pop}} function will only discard SNPs with allele count lower than 2 in the first (target) population (when the first argument is the prefix of a genotype file).
#' @param pops2 If specified, only a pairs between `pops` and `pops2` will be computed
#' @param outpop Keep only SNPs which are heterozygous in this population
#' @param outpop_scale Scale f2-statistics by the inverse `outpop` heteroygosity (`1/(p*(1-p))`). Providing `outpop` and setting `outpop_scale` to `TRUE` will give the same results as the original *qpGraph* when the `outpop` parameter has been set, but it has the disadvantage of treating one population different from the others. This may limit the use of these f2-statistics for other models.
#' @param transitions Set this to `FALSE` to exclude transition SNPs
#' @param transversions Set this to `FALSE` to exclude transversion SNPs
#' @param auto_only Keep only SNPs on chromosomes 1 to 22
#' @param keepsnps SNP IDs of SNPs to keep. Overrides other SNP filtering options
#' @param overwrite Overwrite existing files in `outdir`
#' @param format Supply this if the prefix can refer to genotype data in different formats
#' and you want to choose which one to read. Should be `plink` to read `.bed`, `.bim`, `.fam` files, or `eigenstrat`, or `packedancestrymap` to read `.geno`, `.snp`, `.ind` files.
#' @param adjust_pseudohaploid Genotypes of pseudohaploid samples are usually coded as `0` or `2`, even though only one allele is observed. `adjust_pseudohaploid` ensures that the observed allele count increases only by `1` for each pseudohaploid sample. If `TRUE` (default), samples that don't have any genotypes coded as `1` among the first 1000 SNPs are automatically identified as pseudohaploid. This leads to slightly more accurate estimates of f-statistics. Setting this parameter to `FALSE` treats all samples as diploid and is equivalent to the *ADMIXTOOLS* `inbreed: NO` option. Setting `adjust_pseudohaploid` to an integer `n` will check the first `n` SNPs instead of the first 1000 SNPs.
#' @param cols_per_chunk Number of allele frequency chunks to store on disk. Setting this to a positive integer makes the function slower, but requires less memory. The default value for `cols_per_chunk` in \code{\link{extract_afs}} is 10. Lower numbers will lower the memory requirement but increase the time it takes.
#' @param fst Write files with pairwise FST for every population pair. Setting this to FALSE can make `extract_f2` faster and will require less memory.
#' @param afprod  Write files with allele frequency products for every population pair. Setting this to FALSE can make `extract_f2` faster and will require less memory.
#' @param poly_only Specify whether SNPs with identical allele frequencies in every population should be discarded (`poly_only = TRUE`), or whether they should be used (`poly_only = FALSE`). By default (`poly_only = c("f2")`), these SNPs will be used to compute FST and allele frequency products, but not to compute f2 (this is the default option in the original ADMIXTOOLS).
#' @param apply_corr Apply small-sample-size correction when computing f2-statistics (default `TRUE`)
#' @param qpfstats Compute smoothed f2-statistics (default `FALSE`). In the presence
#' of large amounts of missing data, this option can be used to retain information
#' from all SNPs while introducing less bias than setting `maxmiss` to values greater
#' than 0. When setting `qpfstats = TRUE`, most other options to `extract_f2` will
#' be ignored. See \code{\link{qpfstats}} for more information. Arguments to
#' \code{\link{qpfstats}} can be passed via `...`
#' @param n_cores Parallelize computation across `n_cores` cores via the `doParallel` package.
#' @param verbose Print progress updates
#' @param ... Pass arguments to \code{\link{qpfstats}}
#' @return SNP metadata (invisibly)
#' @seealso \code{\link{f2_from_precomp}} for reading the stored f2-statistics back
#' into R, \code{\link{f2_from_geno}} to skip writting f2-statistics to disk and return them directly
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' f2dir = 'my/f2dir/'
#' extract_f2(pref, f2dir, pops = c('popA', 'popB', 'popC'))
#' }
extract_f2 = function(pref, outdir, inds = NULL, pops = NULL, blgsize = 0.05, maxmem = 8000,
                      maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE, pops2 = NULL, outpop = NULL, outpop_scale = TRUE,
                      transitions = TRUE, transversions = TRUE,
                      auto_only = TRUE, keepsnps = NULL, overwrite = FALSE, format = NULL,
                      adjust_pseudohaploid = TRUE, cols_per_chunk = NULL, fst = TRUE, afprod = TRUE,
                      poly_only = c('f2'), apply_corr = TRUE, qpfstats = FALSE, n_cores = 1, verbose = TRUE, ...) {

  if(!is.null(cols_per_chunk)) {
    stopifnot(is.null(pops2))
    snpfile = extract_f2_large(pref, outdir, inds = inds, pops = pops, blgsize = blgsize,
                               cols_per_chunk = cols_per_chunk, maxmiss = maxmiss,
                               minmaf = minmaf, maxmaf = maxmaf, minac2 = minac2, outpop = outpop, outpop_scale = outpop_scale,
                               transitions = transitions, transversions = transversions,
                               keepsnps = keepsnps, snpblocks = snpblocks, overwrite = overwrite,
                               format = format, adjust_pseudohaploid = adjust_pseudohaploid,
                               afprod = afprod, fst = fst, poly_only = poly_only, verbose = verbose)
    return(invisible(snpfile))
  }
  if(!is.null(outdir)) {
    outdir = normalizePath(outdir, mustWork = FALSE)
    if(length(list.files(outdir)) > 0 && !overwrite) stop('Output directory not empty! Set overwrite to TRUE if you want to overwrite files!')
  }
  if(is.null(inds) && is.null(pops) && verbose && max(file.info(paste0(pref, '.geno'))$size, file.info(paste0(pref, '.bed'))$size, na.rm = T)/1e9 > 1) alert_danger('No poplations or individuals provided. Extracting f2-stats for all population pairs. If that takes too long, you can either specify the "pops" or "inds" parameter, or follow the example in "afs_to_f2".\n')

  if(qpfstats) {
    #if(!is.null(inds)) stop('"qpfstats" option is incompatible with "inds" option!')
    warning("Most `extract_f2` options will be ignored when using the `qpfstats` option!")
    if(is.null(pops)) {
      fi = format_info(pref)
      pops = unique(read_table2(paste0(pref, fi$indend), col_names=F, col_types = fi$indtype)[[which(fi$indnam == 'pop')]])
    }
    f2blocks = qpfstats(pref, pops, ...)
    if(is.null(outdir)) {
      return(list(f2_blocks = f2blocks))
    }
    else {
      block_lengths = parse_number(dimnames(f2blocks)[[3]])
      counts = array(1, dim(f2blocks))
      write_f2(f2blocks, counts, outdir, overwrite = overwrite)
      write_f2(f2blocks, counts, outdir, overwrite = overwrite, id = 'ap')
      saveRDS(block_lengths, paste0(outdir, '/block_lengths_f2.rds'))
      saveRDS(block_lengths, paste0(outdir, '/block_lengths_ap.rds'))
      return()
    }
  }

  if(is.null(inds)) pops = union(pops, pops2)
  afdat = anygeno_to_aftable(pref, inds = inds, pops = pops, format = format,
                             adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)
  afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf, minac2 = minac2, outpop = outpop,
                                  transitions = transitions, transversions = transversions,
                                  keepsnps = keepsnps, auto_only = auto_only, poly_only = FALSE)
  afdat$snpfile %<>% mutate(poly = as.logical(cpp_is_polymorphic(afdat$afs)))
  if(sum(afdat$snpfile$poly) == 0) stop('There are no informative SNPs!')

  if(verbose) alert_warning(paste0(nrow(afdat$afs), ' SNPs remain after filtering. ',
                                   sum(afdat$snpfile$poly),' are polymorphic.\n'))

  if(isTRUE(poly_only)) poly_only = c('f2', 'ap', 'fst')
  arrs = afs_to_f2_blocks(afdat, outdir = outdir, overwrite = overwrite, maxmem = maxmem, poly_only = poly_only,
                          pops1 = pops, pops2 = pops2, outpop = if(outpop_scale) outpop else NULL,
                          blgsize = blgsize, afprod = afprod, fst = fst, apply_corr = apply_corr,
                          n_cores = n_cores, verbose = verbose)

  if(is.null(outdir)) return(arrs)

  if(verbose) alert_info(paste0('Data written to ', outdir, '/\n'))
  invisible(afdat$snpfile)
}

#' Compute blocked f2 statistics
#'
#' This function prepares data for various other *ADMIXTOOLS 2* functions. It reads data from genotype files,
#' computes allele frequencies and blocked f2-statistics for selected populations, and returns them as a 3d array.
#' @export
#' @inheritParams extract_f2
#' @param afprod Return negative average allele frequency products instead of f2-statistics.
#' Setting `afprod = TRUE` will result in more precise f4-statistics when the original data
#' had large amounts of missingness, and should be used in that case for \code{\link{qpdstat}}
#' and \code{\link{qpadm}}. It can also be used for outgroup f3-statistics with a fixed outgroup
#' (for example for \code{\link{qpgraph}}); values will be shifted by a constant amount compared
#' to regular f3-statistics. This shift affects the fit of a graph only by small amounts, possibly
#' less than bias in regular f3-statistics introduced by large amounts of missing data.
#' @return A 3d array of f2-statistics (or scaled allele frequency products if `afprod = TRUE`)
#' @seealso \code{\link{f2_from_precomp}} for reading previously stored f2-statistics into R, \code{\link{extract_f2}} for storing f2-statistics on disk
f2_from_geno = function(pref, inds = NULL, pops = NULL, blgsize = 0.05, maxmem = 8000,
                        maxmiss = 0, minmaf = 0, maxmaf = 0.5, pops2 = NULL, outpop = NULL, outpop_scale = TRUE,
                        transitions = TRUE, transversions = TRUE,
                        auto_only = TRUE, keepsnps = NULL, afprod = FALSE, fst = FALSE, poly_only = c("f2"),
                        format = NULL,
                        adjust_pseudohaploid = TRUE, remove_na = TRUE, apply_corr = TRUE, qpfstats = FALSE, verbose = TRUE, ...) {

  arrs = extract_f2(pref, outdir = NULL, inds = inds, pops = pops, blgsize = blgsize, maxmem = maxmem,
             maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf, pops2 = pops2, outpop = outpop,
             outpop_scale = outpop_scale, transitions = transitions, transversions = transversions,
             auto_only = auto_only, keepsnps = keepsnps, format = format,
             adjust_pseudohaploid = adjust_pseudohaploid, fst = fst, afprod = afprod,
             poly_only = poly_only, apply_corr = apply_corr, qpfstats = qpfstats, verbose = verbose, ...)
  if(afprod) {
    if(verbose) alert_info(paste0('Returning allele frequency product blocks\n'))
    #blocks = scale_ap_blocks(arrs$ap_blocks, from = min(arrs$f2_blocks, na.rm=T), to = max(arrs$f2_blocks, na.rm=T))
    blocks = scale_ap_blocks(arrs$ap_blocks, from = 0)
  } else if(fst) {
    if(verbose) alert_info(paste0('Returning fst blocks\n'))
    blocks = arrs$fst_blocks
  } else {
    if(verbose) alert_info(paste0('Returning f2 blocks\n'))
    blocks = arrs$f2_blocks
  }
  if(remove_na) {
    keep = apply(blocks, 3, function(x) sum(!is.finite(x))==0)
    blocks = blocks[,,keep, drop = FALSE]
    if(sum(keep) == 0) stop("No blocks remain after discarding blocks with missing values! If you want to compute FST for pseudohaploid samples, set `adjust_pseudohaploid = FALSE`")
    if(sum(!keep) > 0) warning(paste0('Discarding ', sum(!keep), ' block(s) due to missing values!\n',
                                      'Discarded block(s): ', paste0(which(!keep), collapse = ', ')))
  }
  blocks
}


scale_ap_blocks = function(ap_blocks, from = NULL, to = NULL) {
  # does the following things to ap_blocks:
  # 1. multiply by -2
  # 2. shift and scale so that all values are between from and to
  # 3. set diag to 0

  ap_blocks = -2*ap_blocks
  if(!is.null(from)) {
    if(!is.null(to)) {
      ap_blocks = (ap_blocks - min(ap_blocks, na.rm=T)) * (to - from)/diff(range(ap_blocks, na.rm = TRUE)) + from
    } else {
      ap_blocks = (ap_blocks - min(ap_blocks, na.rm=T)) + from
    }
  }
  if(isTRUE(all.equal(dimnames(ap_blocks)[[1]], dimnames(ap_blocks)[[2]]))) {
    d1 = dim(ap_blocks)[1]
    d3 = dim(ap_blocks)[3]
    dg = rep(seq_len(d1) + d1*(seq_len(d1)-1),d3) + rep(d1*d1*(seq_len(d3)-1), each = d1)
    ap_blocks[dg] = 0
  }
  ap_blocks
}



#' Compute and store blocked f2 statistics
#'
#' `extract_f2_large` does the same as \code{\link{extract_f2}}, but it requires less memory and is slower. `outdir` has to be set in `extract_f2_large`.
#' @inheritParams extract_f2
#' @param cols_per_chunk Number of populations per chunk. Lowering this number will lower the memory requirements when running
#' @details `extract_f2_large` requires less memory because it writes allele frequency data to disk, and doesn't store the allele frequency matrix for all populations and SNPs in memory. If you still run out of memory, reduce `cols_per_chunk`. This function is a wrapper around \code{\link{extract_afs}} and \code{\link{afs_to_f2}}, and is slower than \code{\link{extract_f2}}. It may be faster to call \code{\link{extract_afs}} and \code{\link{afs_to_f2}} directly, parallelizing over the different calls to \code{\link{afs_to_f2}}.
#' @return SNP metadata (invisibly)
#' @seealso \code{\link{extract_f2}}
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' f2dir = 'my/f2dir/'
#' extract_f2_large(pref, f2dir, pops = c('popA', 'popB', 'popC'))
#' }
extract_f2_large = function(pref, outdir, inds = NULL, pops = NULL, blgsize = 0.05, cols_per_chunk = 10,
                            maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE, outpop = NULL, outpop_scale = TRUE,
                            transitions = TRUE, transversions = TRUE,
                            keepsnps = NULL, snpblocks = NULL, overwrite = FALSE, format = NULL,
                            adjust_pseudohaploid = TRUE, afprod = TRUE, fst = TRUE, poly_only = c('f2'),
                            apply_corr = TRUE, verbose = TRUE) {

  if(verbose) alert_info(paste0('Extracting allele frequencies...\n'))
  snpdat = extract_afs(pref, outdir, inds = inds, pops = pops, cols_per_chunk = cols_per_chunk, numparts = 100,
              maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf, minac2 = minac2, outpop = outpop,
              transitions = transitions, transversions = transversions,
              keepsnps = keepsnps, format = format, poly_only = FALSE,
              adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)

  numchunks = length(list.files(outdir, 'afs.+rds'))

  if(!is.null(outpop) && outpop_scale) {
    p = snpdat$outpopaf
    snpwt = 1/(p*(1-p))
  } else snpwt = NULL
  if(isTRUE(poly_only)) poly_only = c('f2', 'ap', 'fst')

  if(verbose) alert_warning(paste0('Computing ', choose(numchunks+1, 2), ' chunk pairs. If this takes too long,
  consider running "extract_afs" and then paralellizing over "afs_to_f2".\n'))
  for(i in 1:numchunks) {
    for(j in i:numchunks) {
      if(verbose) alert_info(paste0('Writing pair ', i, ' - ', j, '...\r'))
      afs_to_f2(outdir, outdir, chunk1 = i, chunk2 = j, blgsize = blgsize, snpdat = snpdat,
                snpwt = snpwt, overwrite = overwrite, type = 'f2', poly_only = 'f2' %in% poly_only,
                apply_corr = apply_corr)
      if(afprod) afs_to_f2(outdir, outdir, chunk1 = i, chunk2 = j, blgsize = blgsize, snpdat = snpdat,
                           snpwt = snpwt, overwrite = overwrite, type = 'ap', poly_only = 'ap' %in% poly_only)
      if(fst) afs_to_f2(outdir, outdir, chunk1 = i, chunk2 = j, blgsize = blgsize, snpdat = snpdat,
                        snpwt = snpwt, overwrite = overwrite, type = 'fst', poly_only = 'fst' %in% poly_only)
    }
  }
  if(verbose) alert_info('\n')
  if(verbose) alert_info(paste0('Deleting allele frequency files...\n'))
  unlink(paste0(outdir, c('/afs*.rds', '/counts*.rds')))
}


anygeno_to_aftable = function(pref, inds = NULL, pops = NULL, format = NULL, adjust_pseudohaploid = TRUE, verbose = TRUE) {

  if(is.null(format)) {
    if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) {
      format = 'plink'
      geno_to_aftable = plink_to_afs
    } else if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
      if(is_packed(paste0(pref, '.geno'))) {
        format = 'packedancestrymap'
        geno_to_aftable = packedancestrymap_to_afs
      } else {
        format = 'eigenstrat'
        geno_to_aftable = eigenstrat_to_afs
      }
    } else stop('Genotype files not found!')
  }
  geno_to_aftable = switch(format,
                           'plink' = plink_to_afs,
                           'packedancestrymap' = packedancestrymap_to_afs,
                           'eigenstrat' = eigenstrat_to_afs)
  if(is.null(geno_to_aftable)) stop('Invalid format!')

  afdat = geno_to_aftable(pref, inds = inds, pops = pops, adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)
  afdat
}


read_anygeno = function(pref, inds = NULL, format = NULL, verbose = TRUE) {

  if(is.null(format)) {
    if(all(file.exists(paste0(pref, c('.bed', '.bim', '.fam'))))) format = 'plink'
    else if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
      if(is_packed(paste0(pref, '.geno'))) format = 'packedancestrymap'
      else format = 'eigenstrat'
    }
    else stop('Genotype files not found!')
  }
  if(tolower(format) == 'packedancestrymap') {
    read_geno = function(...) {g = read_packedancestrymap(...); list(bed = g$geno, fam = g$ind, bim = g$snp)}
  } else if(tolower(format) == 'eigenstrat') {
    read_geno = function(...) {g = read_eigenstrat(...); list(bed = g$geno, fam = g$ind, bim = g$snp)}
  } else if(tolower(format) == 'plink') {
    read_geno = read_plink
  } else stop('Genotype files not found!')
  if(verbose) alert_info(paste0('Reading genotypes in ', tolower(format), ' format...\n'))
  read_geno(pref, inds, verbose = verbose)
}

#' Extract and store data needed to compute blocked f2
#'
#' Prepare data for various *ADMIXTOOLS 2* functions. This function reads data from genotype files,
#' and extracts data required to compute blocked f-statistics for any sets of samples. The data consists of
#' `.rds` files with total and alternative allele counts for each individual, and products of total
#' and alternative allele counts for each pair.
#' The function calls \code{\link{packedancestrymap_to_afs}} or \code{\link{plink_to_afs}}
#' and \code{\link{afs_to_f2_blocks}}.
#'
#' @export
#' @param inds Individuals for which data should be read. Defaults to all individuals
#' @param maxmiss Discard SNPs which are missing in a fraction of individuals greater than `maxmiss`
#' @param cols_per_chunk Number of genotype chunks to store on disk. Setting this to a positive integer makes the function slower, but requires less memory. The default value for `cols_per_chunk` in \code{\link{extract_afs}} is 10. Lower numbers will lower the memory requirement but increase the time it takes.
#' @inheritParams extract_f2
extract_counts = function(pref, outdir, inds = NULL, blgsize = 0.05,  maxmiss = 0, minmaf = 0, maxmaf = 0.5,
                          transitions = TRUE, transversions = TRUE, auto_only = TRUE, keepsnps = NULL,
                          maxmem = 8000, overwrite = FALSE, format = NULL, cols_per_chunk = NULL, verbose = TRUE) {

  if(!is.null(cols_per_chunk)) {
    snpfile = extract_counts_large(pref, outdir, inds = inds, blgsize = blgsize,
                               cols_per_chunk = cols_per_chunk, maxmiss = maxmiss,
                               minmaf = minmaf, maxmaf = maxmaf,
                               transitions = transitions, transversions = transversions, auto_only = auto_only,
                               keepsnps = keepsnps, overwrite = overwrite,
                               format = format, verbose = verbose)
    return(invisible(snpfile))
  }

  dir.create(outdir, showWarnings = FALSE)
  if(length(list.files(outdir)) > 0 && !overwrite) stop('Output directory not empty! Set overwrite to TRUE if you want to overwrite files!')
  bfile = paste0(outdir, '/block_lengths.rds')
  # if(file.exists(bfile)) {
  #   todo: snp filters
  #   extract_more_counts(pref, outdir, inds = inds, maxmem = maxmem,
  #                      overwrite = overwrite, format = format, verbose = verbose)
  #   return()
  # }

  g = read_anygeno(pref, inds, format = format, verbose = verbose)
  g %<>% discard_from_geno(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf,
                           transitions = transitions, transversions = transversions,
                           keepsnps = keepsnps, auto_only = TRUE, poly_only = FALSE)

  if(verbose) alert_info(paste0(nrow(g$bed), ' SNPs remain.\nDetermining SNP blocks...\n'))
  block_lengths = get_block_lengths(g$bim, blgsize = blgsize)
  saveRDS(block_lengths, file = bfile)
  xmat_to_inddat(g$bed, block_lengths,
                 outdir = outdir, overwrite = overwrite,
                 maxmem = maxmem, verbose = verbose)
  xmat_to_pairdat(g$bed, block_lengths,
                  outdir = outdir, overwrite = overwrite,
                  maxmem = maxmem, verbose = verbose)
  if(verbose) alert_info(paste0('Data written to ', outdir, '/\n'))
}

extract_counts_large = function(pref, outdir, inds = NULL, blgsize = 0.05, cols_per_chunk = 10,
                                maxmiss = 0, minmaf = 0, maxmaf = 0.5,
                                transitions = TRUE, transversions = TRUE, auto_only = TRUE, keepsnps = NULL,
                                overwrite = FALSE, format = NULL, verbose = TRUE) {

  # same as extract_counts, but split across chunks of individuals
  pref %<>% normalizePath(mustWork = FALSE)
  l = format_info(pref, format)
  indfile = read_table2(paste0(pref, l$indend), col_names = l$indnam, col_types = l$indtype, progress = FALSE)
  if(!is.null(inds)) indfile %<>% filter(ind %in% inds)

  snpdat = extract_afs(pref, outdir, inds = indfile$ind, pops = indfile$ind, cols_per_chunk = cols_per_chunk, numparts = 100,
                       maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf,
                       transitions = transitions, transversions = transversions, auto_only = auto_only,
                       keepsnps = keepsnps, format = format, verbose = verbose)
  numchunks = length(list.files(outdir, 'afs.+rds'))

  block_lengths = get_block_lengths(snpdat, blgsize = blgsize)
  saveRDS(block_lengths, file = paste0(outdir, '/block_lengths.rds'))

  write_split_inddat(outdir, outdir, overwrite = overwrite, verbose = verbose)

  if(verbose) alert_warning(paste0('Computing ', choose(numchunks, 2), ' chunk pairs. If this takes too long,
  consider running "extract_afs" (with each individual its own population) and then paralellizing over "afs_to_counts".\n'))
  for(i in 1:numchunks) {
    for(j in i:numchunks) {
      if(verbose) alert_info(paste0('Writing pair ', i, ' - ', j, '...\r'))
      afs_to_counts(outdir, outdir, chunk1 = i, chunk2 = j, overwrite = overwrite)
    }
  }
  if(verbose) alert_info('\n')
  if(verbose) alert_info(paste0('Deleting allele frequency files...\n'))
  unlink(paste0(outdir, c('/afs*.rds', '/counts*.rds')))
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
#' Prepare data for various *ADMIXTOOLS 2* functions. Reads data from packedancestrymap or PLINK files,
#' and computes allele frequencies for selected populations and stores it as `.rds` files in outdir.
#' @export
#' @inheritParams extract_f2
#' @param cols_per_chunk Number of populations per chunk. Lowering this number will lower the memory requirements when running \code{\link{afs_to_f2}}, but more chunk pairs will have to be computed.
#' @return SNP metadata (invisibly)
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' outdir = 'dir/for/afdata/'
#' extract_afs(pref, outdir)
#' }
extract_afs_simple = function(pref, outdir, inds = NULL, pops = NULL, blgsize = 0.05, cols_per_chunk = 10,
                       maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE, outpop = NULL, transitions = TRUE, transversions = TRUE,
                       keepsnps = NULL, format = NULL, poly_only = FALSE, adjust_pseudohaploid = TRUE, verbose = TRUE) {

  # read data and compute allele frequencies
  afdat = anygeno_to_aftable(pref, inds = inds, pops = pops, format = format,
                             adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)
  afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf, minac2 = minac2, outpop = outpop,
                                  transitions = transitions, transversions = transversions,
                                  keepsnps = keepsnps, auto_only = TRUE, poly_only = poly_only)

  afdat$snpfile %<>% mutate(poly = as.logical(cpp_is_polymorphic(afdat$afs)))
  if(verbose) alert_warning(paste0(nrow(afdat$afs), ' SNPs remain after filtering. ',
                                   sum(afdat$snpfile$poly),' are polymorphic.\n'))

  # split allele frequency data into chunks and write to disk
  split_mat(afdat$afs, cols_per_chunk = cols_per_chunk, prefix = paste0(outdir, '/afs'), verbose = verbose)
  split_mat(afdat$counts, cols_per_chunk = cols_per_chunk, prefix = paste0(outdir, '/counts'), verbose = verbose)
  # compute jackknife blocks
  block_lengths = get_block_lengths(afdat$snpfile %>% filter(poly), blgsize = blgsize)
  block_lengths_a = get_block_lengths(afdat$snpfile, blgsize = blgsize)
  saveRDS(block_lengths, file = paste0(outdir, '/block_lengths.rds'))
  saveRDS(block_lengths_a, file = paste0(outdir, '/block_lengths_a.rds'))
  write_tsv(afdat$snpfile, paste0(outdir, '/snpdat.tsv.gz'))
  invisible(afdat$snpfile)
}


#' Compute and store blocked allele frequency data
#'
#' Prepare data for various *ADMIXTOOLS 2* functions. Reads data from packedancestrymap or PLINK files,
#' and computes allele frequencies for selected populations and stores it as `.rds` files in outdir.
#' @export
#' @inheritParams extract_f2
#' @param cols_per_chunk Number of populations per chunk. Lowering this number will lower the memory requirements when running \code{\link{afs_to_f2}}, but more chunk pairs will have to be computed.
#' @param numparts Number of parts in which genotype data will be read for computing allele frequencies
#' @return SNP metadata (invisibly)
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' outdir = 'dir/for/afdata/'
#' extract_afs(pref, outdir)
#' }
extract_afs = function(pref, outdir, inds = NULL, pops = NULL, cols_per_chunk = 10, numparts = 100,
                       maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE, outpop = NULL, auto_only = TRUE,
                       transitions = TRUE, transversions = TRUE, keepsnps = NULL, format = NULL,
                       poly_only = FALSE, adjust_pseudohaploid = TRUE, verbose = TRUE) {

  pref %<>% normalizePath(mustWork = FALSE)
  l = format_info(pref, format)

  if(verbose) alert_info(paste0('Reading metadata...\n'))
  indfile = read_table2(paste0(pref, l$indend), col_names = l$indnam, col_types = l$indtype, progress = FALSE)
  snpfile = read_table2(paste0(pref, l$snpend), col_names = l$snpnam, col_types = 'ccddcc', progress = FALSE)
  nindall = nrow(indfile)
  nsnpall = nrow(snpfile)

  ip = match_samples(indfile$ind, indfile$pop, inds, pops)
  indvec = ip$indvec - 1

  if(auto_only) snpfile %<>% mutate(CHR = as.numeric(gsub('[a-zA-Z]+', '', CHR))) %>% filter(CHR <= 22)

  starts = seq(0, nrow(snpfile), length.out = numparts+1) %>% round %>% head(-1)
  ends = c(lead(starts)[-numparts], nrow(snpfile))

  snpparts = list()
  ntest = if(is.numeric(adjust_pseudohaploid)) adjust_pseudohaploid else 1000
  if(adjust_pseudohaploid) ploidy = l$cpp_geno_ploidy(paste0(pref, l$genoend), nsnpall, nindall, indvec, ntest)
  else ploidy = rep(2, nindall)

  for(i in 1:numparts) {
    if(verbose) alert_info(paste0('Reading part ', i, ' out of ', numparts, '...\r'))
    # read data and compute allele frequencies
    afdat = l$cpp_geno_to_afs(paste0(pref, l$genoend), nsnpall, nindall, indvec, first = starts[i],
                                last = ends[i], ploidy = ploidy, transpose = FALSE, verbose = FALSE)
    afdat$snpfile = snpfile %>% slice((starts[i]+1):(ends[i]))

    colnames(afdat$afs) = colnames(afdat$counts) = ip$upops
    rownames(afdat$afs) = rownames(afdat$counts) = afdat$snpfile$SNP

    afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf, minac2 = minac2, outpop = outpop,
                                    transitions = transitions, transversions = transversions,
                                    keepsnps = keepsnps, auto_only = TRUE, poly_only = poly_only)

    afdat$snpfile %<>% mutate(poly = as.logical(cpp_is_polymorphic(afdat$afs)))
    if(!is.null(outpop)) afdat$snpfile %<>% mutate(outpopaf = afdat$afs[,outpop])
    snpparts[[i]] = afdat$snpfile

    # split allele frequency data into chunks and write to disk
    partdir = paste0(outdir, '/part',i,'/')
    dir.create(partdir, recursive = TRUE, showWarnings = FALSE)
    split_mat(afdat$afs, cols_per_chunk = cols_per_chunk, prefix = paste0(partdir,'/afs'), verbose = FALSE)
    split_mat(afdat$counts, cols_per_chunk = cols_per_chunk, prefix = paste0(partdir, '/counts'), verbose = FALSE)
  }
  if(verbose) alert_info('\n')
  snpparts %<>% bind_rows %>% mutate(CHR = as.character(CHR))
  if(verbose) alert_warning(paste0(nrow(snpparts), ' SNPs remain after filtering. ',
                                   sum(snpparts$poly),' are polymorphic.\n'))

  numchunks = length(list.files(partdir, 'afs.+rds'))
  for(j in seq_len(numchunks)) {
    if(verbose) alert_info(paste0('Merging chunk ', j, ' out of ', numchunks, '...\r'))
    am = cm = list()
    for(i in seq_len(numparts)) {
      partdir = paste0(outdir, '/part', i, '/')
      am[[i]] = readRDS(paste0(partdir, '/afs', j, '.rds'))
      cm[[i]] = readRDS(paste0(partdir, '/counts', j, '.rds'))
    }
    saveRDS(do.call(rbind, am), file = paste0(outdir,'/afs', j, '.rds'))
    saveRDS(do.call(rbind, cm), file = paste0(outdir,'/counts', j, '.rds'))
  }
  unlink(paste0(outdir, '/part', seq_len(numparts), '/'), recursive = TRUE)
  if(verbose) alert_info('\n')

  write_tsv(snpparts, paste0(outdir, '/snpdat.tsv.gz'))
  invisible(snpparts)
}


format_info = function(pref, format = NULL) {
  # returns named list with genotype format specific information

  l = list()
  if(is_ancestrymap_prefix(pref) && is_packed(paste0(pref, '.geno')) || isTRUE(format == 'packedancestrymap')) {
    l$format = 'packedancestrymap'
    l$snpnam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
    l$indnam = c('ind', 'sex', 'pop')
    l$indtype = 'ccc'
    l$snpend = '.snp'
    l$indend = '.ind'
    l$genoend = '.geno'
    l$cpp_geno_ploidy = cpp_packedancestrymap_ploidy
    l$cpp_geno_to_afs = cpp_packedancestrymap_to_afs
    l$cpp_read_geno = cpp_read_packedancestrymap
  } else if(is_ancestrymap_prefix(pref) && !is_packed(paste0(pref, '.geno')) || isTRUE(format == 'eigenstrat')) {
    l$format = 'eigenstrat'
    l$snpnam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
    l$indnam = c('ind', 'sex', 'pop')
    l$indtype = 'ccc'
    l$snpend = '.snp'
    l$indend = '.ind'
    l$genoend = '.geno'
    l$cpp_geno_ploidy = eigenstrat_ploidy
    l$cpp_geno_to_afs = cpp_eigenstrat_to_afs
    l$cpp_read_geno = cpp_read_eigenstrat
  } else if(is_plink_prefix(pref) || isTRUE(format == 'plink')) {
    l$format = 'plink'
    l$snpnam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
    l$indnam = c('pop', 'ind', 'p1', 'p2', 'sex', 'pheno')
    l$indtype = 'cccccc'
    l$snpend = '.bim'
    l$indend = '.fam'
    l$genoend = '.bed'
    l$cpp_geno_ploidy = cpp_plink_ploidy
    l$cpp_geno_to_afs = cpp_plink_to_afs
    l$cpp_read_geno = cpp_read_plink
  } else stop('Genotype files not found!')
  l
}


#' Copy f2-statistics
#'
#' Copy a subset of f2-statistics to a new directory
#' @export
#' @param from Directory with f2-statistics
#' @param to Target directory
#' @param pops Populations to copy
#' @param verbose Print progress updates
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' outdir = 'dir/for/afdata/'
#' extract_f2_subset(pref, outdir)
#' }
extract_f2_subset = function(from, to, pops, verbose = TRUE) {

  if(!dir.exists(to)) dir.create(to)
  file.copy(paste0(from, '/block_lengths.rds'), paste0(to, '/block_lengths.rds'))
  file.copy(paste0(from, '/block_lengths_a.rds'), paste0(to, '/block_lengths_a.rds'))
  for(p1 in pops) {
    dir = paste0(to, '/', p1)
    if(!dir.exists(dir)) dir.create(dir)
    for(p2 in pops) {
      if(stringi::stri_cmp_le(p1, p2, locale = 'C')) {
        fn1 = paste0(p1, '/', p2, '_f2.rds')
        fn2 = paste0(p1, '/', p2, '_ap.rds')
        file.copy(paste0(from, '/', fn1), paste0(to, '/', fn1))
        file.copy(paste0(from, '/', fn2), paste0(to, '/', fn2))
      }
    }
  }
}



# this should not be needed in practice; used for testing
f2_from_geno_indivs = function(pref, inds = NULL, pops = NULL, format = NULL, maxmem = 8000,
                               apply_corr = TRUE, verbose = TRUE) {

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
#' @param dir Directory with precomputed f2 statistics, or precomputed individual pair data
#' @param inds Individuals for which data should be read. Defaults to all individuals,
#' which may require a lot of memory.
#' @param pops Populations for which data should be read. Defaults to all populations,
#' which may require a lot of memory.
#' @param pops2 Optional second vector of populations. Useful if a f4 statistics of a few against many populations should be computed. `pops2` should not be specified in other cases, as most functions depend on f2-statistics for all population pairs in `pops`.
#' @param afprod Return negative average allele frequency products instead of f2 estimates. This will result in more precise f4-statistics when the original data had large amounts of missingness, and should be used in that case for \code{\link{qpdstat}} and \code{\link{qpadm}}. It can also be used for outgroup f3-statistics with a fixed outgroup (for example for \code{\link{qpgraph}}); values will be shifted by a constant amount compared to regular f3-statistics. This shift affects the fit of a graph only by small amounts, possibly less than bias in regular f3-statistics introduced by large amounts of missing data.
#' This option is currently ineffective when reading data extracted with \code{\link{extract_counts}}.
#' @param return_array Return a 3d array (default). If false, a data frame will be returned.
#' @param apply_corr Subtract the f2 correction factor. Setting this to `FALSE` can occasionally be useful
#' @param remove_na Remove blocks with missing values
#' @param verbose Print progress updates
#' @return A 3d array of f2 statistics
#' @examples
#' \dontrun{
#' dir = 'my/f2/dir/'
#' f2_blocks = f2_from_precomp(dir, pops = c('pop1', 'pop2', 'pop3'))
#' }
f2_from_precomp = function(dir, inds = NULL, pops = NULL, pops2 = NULL, afprod = FALSE, fst = FALSE,
                           return_array = TRUE, apply_corr = TRUE, remove_na = TRUE, verbose = TRUE) {

  if(!is.null(pops) && !is.null(inds) && length(pops) != length(inds)) stop("'pops' and 'inds' are not the same length!")
  if(!dir.exists(dir)) stop(paste0("Directory ", normalizePath(dir, mustWork = FALSE), " doesn't exist!"))
  indpairs = dir.exists(paste0(dir, '/indivs'))

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
    f2_blocks = f2_from_precomp_indivs(dir, poplist = poplist, afprod = afprod, return_array = return_array,
                                       apply_corr = apply_corr, verbose = verbose)$f2_blocks
  } else {
    if(!is.null(inds)) stop('Individual IDs supplied, but no "indivs" directory found!')
    if(is.null(pops)) pops = list.dirs(dir, full.names = FALSE, recursive = FALSE)
    if(is.null(pops2)) pops2 = pops
    if(verbose) alert_info(paste0('Reading precomputed data for ', length(union(pops, pops2)), ' populations...\n'))
    type = if(fst) 'fst' else if(afprod) 'ap' else 'f2'
    f2_blocks = read_f2(dir, pops, pops2, type = type, remove_na = remove_na, verbose = verbose)
  }
  if(afprod) f2_blocks = scale_ap_blocks(f2_blocks, from = 0)
  else if(any(apply(f2_blocks,1:2,mean,na.rm=T)<0) && !fst) warning(paste('Some f2-statistic estimates are',
         'negative across blocks. This is probably caused by too many missing or rare SNPs in',
         'populations with low sample size, which makes the f2 bias correction unreliable.',
         "Consider running 'f2_from_precomp' with 'afprod = TRUE'."))
  f2_blocks
}


f2_from_precomp_indivs = function(dir, poplist = NULL, afprod = FALSE, return_array = TRUE, apply_corr = TRUE, verbose = TRUE) {

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
  f2_blocks = indpairs_to_f2blocks(indivs, pairs, poplist, block_lengths, afprod = afprod,
                                   apply_corr = apply_corr, return_array = return_array)
  if(remove_na) {
    keep = apply(f2_blocks, 3, function(x) sum(is.na(x)) == 0)
    f2_blocks = f2_blocks[,,keep]
  }
  namedList(f2_blocks, indivs, pairs, poplist)
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
      if(stringi::stri_cmp_le(ind1, ind2, locale = 'C') || !same) {
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
      if(stringi::stri_cmp_le(ind1, ind2, locale = 'C') && (ind1 %in% add || ind2 %in% add)) {
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


xmat_to_inddat = function(xmat, block_lengths, maxmem = 8000,
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
    group_by(bl, ind) %>%
    #summarize(a = mean(a), n = mean(n)) %>%
    summarize(a = (mean(a)+mean(n-a))/2, n = mean(n)) %>%
    ungroup

  if(!is.null(outdir)) {
    inddir = paste0(outdir, '/indivs')
    if(!dir.exists(inddir)) dir.create(inddir, recursive = TRUE, showWarnings = FALSE)
    indivs %>% select(-bl) %>% split(.$ind) %>% map(~as.matrix(.[,-1])) %>%
      imap(write_indiv, outdir = outdir, overwrite = overwrite)
  } else {
    return(indivs)
  }
}


xmat_to_pairdat = function(xmat, block_lengths, maxmem = 8000,
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
                                         ' chunks of ', width, ' samples and up to ', maxmem,
                                         ' MB (', choose(numsplits2+1,2), ' chunk pairs)\n'))
    else alert_info(paste0('Computing without splitting (', reqmem, ' < maxmem)...\n'))
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
      files = expand_grid(n1=nam, n2=nam) %>%
        filter(stringi::stri_cmp_le(n1, n2, locale = 'C')) %$% paste0(outdir, '/pairs/', n1, '/', n2, '.rds')
      if(all(file.exists(files))) next
    }
    arrs = xmats_to_pairarrs(a1, a2)
    arrs2 = xmats_to_pairarrs(2-a1, 2-a2)
    aa_subblock = (block_arr_mean(arrs$aa, block_lengths) + block_arr_mean(arrs2$aa, block_lengths))/2
    nn_subblock = (block_arr_mean(arrs$nn, block_lengths) + block_arr_mean(arrs2$nn, block_lengths))/2
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
#' @param dir Directory with precomputed individual pair data
#' @param inds Individuals to group
#' @param pops Group names, either length 1, or same length as `inds`
#' @param overwrite Overwrite existing files in `outdir`
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
    mutate(min1 = stringi::stri_cmp_lt(pop1, pop2, locale = 'C'),
           p1 = ifelse(min1, pop1, pop2), p2 = ifelse(min1, pop2, pop1)) %>%
    select(-min1) %>% mutate(pop1 = p1, pop2 = p2) %>%
    #mutate(p1 = pmin(pop1, pop2), p2 = pmax(pop1, pop2), pop1 = p1, pop2 = p2) %>%
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
#' @param dir Directory with precomputed individual pair data
#' @param groups Groups to delete. Defaults to all groups
#' @param verbose Print progress updates
#' @return Invisibly returns sample IDs in deleted groups as character vector
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

is_group = function(dir, group) file.exists(paste0(dir, '/groups/', group, '.rds'))


is_packed = function(filename) {
  # tries to infer if filename is binary packedancestrymap file or text eigenstrat file
  conn = file(filename, 'rb')
  on.exit(close(conn))
  dat = readBin(conn, 'int', 1e3, 1, signed=F)
  max(dat) > 128 || length(unique(dat)) > 5
}

is_ancestrymap_prefix = function(input) {
  if(!is.character(input) || length(input) > 1) return(FALSE)
  all(file.exists(paste0(input, c('.geno', '.ind', '.snp'))))
}

is_plink_prefix = function(input) {
  if(!is.character(input) || length(input) > 1) return(FALSE)
  filesexist = all(file.exists(paste0(input, c('.bed', '.bim', '.fam'))))
  filesexist
}

is_geno_prefix = function(input) {
  is_ancestrymap_prefix(input) || is_plink_prefix(input)
}

is_precomp_dir = function(input) {
  is.character(input) && dir.exists(input) && file.exists(paste0(input, '/block_lengths_f2.rds'))
}

#' Convert *EIGENSTRAT* or *PACKEDANCESTRYMAP* to *PLINK*
#'
#' This function converts *EIGENSTRAT/PACKEDANCESTRYMAP* format files to *PLINK* files, using \code{\link[genio]{write_plink}}.
#' When `inds` or `pops` is provided, only a subset of samples will be extracted.
#' This function can have a high memory footprint, because data for all SNPs will be read before writing the *PLINK* files.
#' @export
#' @param inpref Prefix of the input files
#' @param outpref Prefix of the *PLINK* output files
#' @param inds Individuals which should be extracted
#' @param pops Populations which should be extracted. Can not be provided together with `inds`
#' @param verbose Print progress updates
#' @aliases eigenstrat_to_plink
#' @section Alias:
#' `eigenstrat_to_plink`
packedancestrymap_to_plink = function(inpref, outpref, inds = NULL, pops = NULL, verbose = TRUE) {
  # extracts samples from geno file and writes new, smaller PLINK file using genio

  stopifnot(is.null(inds) || is.null(pops))
  stopifnot(is_ancestrymap_prefix(inpref))
  if(system.file(package = 'genio') == '') stop('Please install the "genio" package!')
  inpref = normalizePath(inpref, mustWork = F)
  if(!is.null(pops)) {
    inds = read_table2(paste0(inpref, '.ind'), col_names = F, col_types = 'ccc') %>%
      filter(X3 %in% pops) %$% X1
  }

  read_geno = ifelse(is_packed(paste0(inpref, '.geno')), read_packedancestrymap, read_eigenstrat)
  dat = read_geno(inpref, inds, verbose = verbose)

  genio::write_plink(normalizePath(outpref, mustWork = F),
                     dat$geno[,dat$ind$X1,drop=FALSE],
                     bim = dat$snp %>% set_colnames(c('id', 'chr', 'posg', 'pos', 'ref', 'alt')),
                     fam = dat$ind %>% transmute(fam = X3, id = X1, pat = 0, mat = 0, sex = X2, pheno = -9),
                     verbose = verbose)
}
#' @export
eigenstrat_to_plink = packedancestrymap_to_plink

#' Extract samples from PLINK files
#'
#' This function reads *PLINK* files, extracts a subset of samples, and writes new *PLINK* files using \code{\link[genio]{write_plink}}. It's probably slower than running the equivalent command in *PLINK* directly, but it can be useful to do this from within R.
#' When `inds` or `pops` is provided, only a subset of samples will be extracted.
#' @export
#' @param inpref Prefix of the *PLINK* input files
#' @param outpref Prefix of the *PLINK* output files
#' @param inds Individuals which should be extracted
#' @param pops Populations which should be extracted. Can not be provided together with `inds`
#' @param overwrite Set this to `TRUE` if `inpref == outpref` and you really want to overwrite the input files.
#' @param verbose Print progress updates
extract_samples = function(inpref, outpref, inds = NULL, pops = NULL, overwrite = FALSE, verbose = TRUE) {

  stopifnot(is.null(inds) || is.null(pops))
  if(system.file(package = 'genio') == '') stop('Please install the "genio" package!')
  if(inpref == outpref && !overwrite)
    stop("If you really want to overwrite the input files, set 'overwrite = TRUE'")
  if(!is.null(pops)) {
    inds = read_table2(paste0(inpref, '.fam'), col_names = F, col_types = 'cccccc') %>%
      filter(X1 %in% pops) %$% X2
  }
  dat = read_plink(inpref, inds, verbose = verbose)

  genio::write_plink(normalizePath(outpref, mustWork = F),
                     dat$bed[,dat$fam$X2],
                     bim = dat$bim %>% set_colnames(c('chr', 'id', 'posg', 'pos', 'ref', 'alt')),
                     fam = dat$fam %>% transmute(fam = X1, id = X2, pat = 0, mat = 0, sex = X2, pheno = -9),
                     verbose = verbose)
}


#' f4 from genotype data
#'
#' Compute per-block f4-statistics directly from genotype data
#' @export
#' @param pref Prefix of genotype files
#' @param popcombs A data frame with one population combination per row, and columns `pop1`, `pop2`, `pop3`, `pop4`. If there is an additional integer column named `model` and `allsnps = FALSE`, only SNPs present in every population in any given model will be used to compute f4-statistics for that model.
#' @param left Populations on the left side of f4 (`pop1` and `pop2`). Can be provided together with `right` in place of `popcombs`.
#' @param right Populations on the right side of f4 (`pop3` and `pop4`). Can be provided together with `left` in place of `popcombs`.
#' @param auto_only Use only chromosomes 1 to 22.
#' @param blgsize SNP block size in Morgan. Default is 0.05 (50 cM). If `blgsize` is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
#' @param block_lengths An optional vector with block lengths. If `NULL`, block lengths will be computed.
#' @param f4mode If `TRUE`: f4 is computed from allele frequencies `a`, `b`, `c`, and `d` as `(a-b)*(c-d)`. if `FALSE`, D-statistics are computed instead, defined as `(a-b)*(c-d) / ((a + b - 2*a*b) * (c + d - 2*c*d))`, which is the same as `(P(ABBA) - P(BABA)) / (P(ABBA) + P(BABA))`.
#' @param allsnps Use all SNPs with allele frequency estimates in every population of any given population quadruple. If `FALSE` (the default) only SNPs which are present in all populations in `popcombs` (or any given model in it) will be used. Setting `allsnps = TRUE` in the presence of large amounts of missing data might lead to false positive results.
#' @param verbose Print progress updates
#' @return A data frame with per-block f4-statistics for each population quadruple.
f4blockdat_from_geno = function(pref, popcombs = NULL, left = NULL, right = NULL, auto_only = TRUE,
                                blgsize = 0.05,
                                block_lengths = NULL, f4mode = TRUE, allsnps = FALSE,
                                poly_only = FALSE, verbose = TRUE, snpwt = NULL, ...) {

  stopifnot(!is.null(popcombs) || (!is.null(left) && !is.null(right)))
  stopifnot(is.null(popcombs) || is.null(left) && is.null(right))

  if(allsnps == 'qpfs') {
    return(f4blockdat_from_geno_qpfs(pref, popcombs, auto_only = auto_only,
                                     blgsize = blgsize,
                                     block_lengths = block_lengths, f4mode = f4mode, verbose = verbose))
  }

  if(is.null(popcombs)) popcombs = tibble(pop1 = left[1], pop2 = left[-1]) %>%
    expand_grid(tibble(pop3 = right[1], pop4 = right[-1]))
  if(is.matrix(popcombs)) {
    if(ncol(popcombs) != 4) stop("'popcombs' is a matrix but doens't have four columns!")
    popcombs %<>% as_tibble(.name_repair = ~paste0('pop', 1:4))
  }

  if('model' %in% names(popcombs)) {
    hasmodels = TRUE
  } else {
    hasmodels = FALSE
    popcombs %<>% mutate(model = 1)
  }
  if(allsnps) {
    modelvec = 0
    pc = popcombs %>% select(pop1:pop4) %>% distinct
  } else {
    modelvec = popcombs$model
    pc = popcombs
  }

  if(verbose) alert_info('Reading metadata...\n')
  pref = normalizePath(pref, mustWork = FALSE)
  l = format_info(pref)

  indfile = read_table2(paste0(pref, l$indend), col_names = l$indnam, col_types = l$indtype, progress = FALSE)
  snpfile = read_table2(paste0(pref, l$snpend), col_names = l$snpnam, col_types = 'ccddcc', progress = FALSE)
  cpp_read_geno = l$cpp_read_geno
  fl = paste0(pref, l$genoend)

  nsnpall = nrow(snpfile)
  nindall = nrow(indfile)

  snpfile$keep = TRUE
  if(auto_only) snpfile %<>% mutate(keep = as.numeric(gsub('[a-zA-Z]+', '', CHR)) <= 22)
  if('keepsnps' %in% names(list(...))) snpfile %<>% mutate(keep = keep & SNP %in% list(...)$keepsnps)
  nsnpaut = sum(snpfile$keep)
  pops = unique(c(popcombs$pop1, popcombs$pop2, popcombs$pop3, popcombs$pop4))

  if(!all(pops %in% indfile$pop))
    stop(paste0('Populations missing from indfile: ', paste0(setdiff(pops, indfile$pop), collapse = ', ')))
  if(!is.null(block_lengths) && sum(block_lengths) != nsnpaut)
    stop(paste0('block_lengths should sum to ', nsnpaut,' (the number of autosomal SNPs)'))
  if(any(duplicated(indfile$ind)))
    stop('Duplicate individual IDs are not allowed!')

  allinds = indfile$ind
  allpops = indfile$pop
  indfile %<>% filter(pop %in% pops)
  indvec = (allinds %in% indfile$ind)+0
  popvec = match(indfile$pop, pops)
  p1 = match(pc$pop1, pops)
  p2 = match(pc$pop2, pops)
  p3 = match(pc$pop3, pops)
  p4 = match(pc$pop4, pops)

  if(verbose) alert_info(paste0('Computing block lengths for ', sum(snpfile$keep),' SNPs...\n'))
  if(is.null(block_lengths)) block_lengths = get_block_lengths(snpfile %>% filter(keep), blgsize = blgsize)
  numblocks = length(block_lengths)

  snpfile$block = NA; snpfile$block[snpfile$keep] = rep(1:length(block_lengths), block_lengths)
  snpfile %<>% fill(block, .direction = 'updown')
  snpind = split(snpfile$keep, snpfile$block)
  bl = rle(snpfile$block)$lengths
  start = lag(cumsum(bl), default = 0)
  end = cumsum(bl)

  popind = popcombs %>%
    group_by(model) %>%
    summarize(pp = list(unique(c(pop1, pop2, pop3, pop4)))) %>%
    mutate(popind = map(pp, ~match(., pops))) %$% popind
  usesnps = matrix(0)

  numer = denom = cnt = matrix(NA, numblocks, nrow(pc))
  for(i in 1:numblocks) {
    if(verbose) alert_info(paste0('Computing ', nrow(pc),' f4-statistics for block ',
                                  i, ' out of ', numblocks, '...\r'))
    # replace following two lines with cpp_geno_to_afs?
    gmat = cpp_read_geno(fl, nsnpall, nindall, indvec, start[i], end[i], T, F)[,snpind[[i]]]
    at = gmat_to_aftable(gmat, popvec)
    if(!allsnps) {
      usesnps = popind %>% map(~(colSums(!is.finite(at[.,,drop=FALSE])) == 0)+0) %>% do.call(rbind, .)
      if(poly_only) {
        #fn = function(mat) apply(mat, 2, function(x) length(unique(na.omit(x))) > 1)+0
        fn = function(mat) apply(mat, 2, function(x) length(unique(na.omit(x))) > 1 | !max(na.omit(x)) %in% c(0,1))+0
        usesnps = (usesnps & (popind %>% map(~fn(at[.,,drop=FALSE])) %>% do.call(rbind, .)))+0
      }
    }
    num = cpp_aftable_to_dstatnum(at, p1, p2, p3, p4, modelvec, usesnps, allsnps, poly_only)
    if(!is.null(snpwt)) {
      num$num = t(t(num$num) * snpwt[(start[i]+1):end[i]])
    }
    numer[i,] = unname(rowMeans(num$num, na.rm = TRUE))
    cnt[i,] = c(num$cnt)
    if(!f4mode) {
      den = cpp_aftable_to_dstatden(at, p1, p2, p3, p4, modelvec, usesnps, allsnps, poly_only)
      denom[i,] = unname(rowMeans(den, na.rm = TRUE))
    }
  }
  if(verbose) cat('\n')

  out = pc %>%
    expand_grid(block = 1:numblocks) %>%
    mutate(est = c(numer), n = c(cnt))
  #if(!f4mode) out %<>% mutate(est = est/c(denom))
  if(!f4mode) out %<>% mutate(den = c(denom))
  if(allsnps) out = popcombs %>% left_join(out, by = paste0('pop', 1:4))
  if(!hasmodels) out$model = NULL

  out %>%
    mutate(length = block_lengths[block], est = if_else(is.finite(est) & n > 0, est, NA_real_))
}

#' f3 from genotype data
#'
#' Compute per-block f3-statistics directly from genotype data
#' @export
#' @param pref Prefix of genotype files
#' @param popcombs A data frame with one population combination per row, and columns `pop1`, `pop2`, `pop3`, `pop4`. If there is an additional integer column named `model` and `allsnps = FALSE`, only SNPs present in every population in any given model will be used to compute f4-statistics for that model.
#' @param auto_only Use only chromosomes 1 to 22.
#' @param blgsize SNP block size in Morgan. Default is 0.05 (50 cM). If `blgsize` is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.
#' @param block_lengths An optional vector with block lengths. If `NULL`, block lengths will be computed.
#' @param allsnps Use all SNPs with allele frequency estimates in every population of any given population quadruple. If `FALSE` (the default) only SNPs which are present in all populations in `popcombs` (or any given model in it) will be used. Setting `allsnps = TRUE` in the presence of large amounts of missing data might lead to false positive results.
#' @param adjust_pseudohaploid Genotypes of pseudohaploid samples are usually coded as `0` or `2`, even though only one allele is observed. `adjust_pseudohaploid` ensures that the observed allele count increases only by `1` for each pseudohaploid sample. If `TRUE` (default), samples that don't have any genotypes coded as `1` among the first 1000 SNPs are automatically identified as pseudohaploid. This leads to slightly more accurate estimates of f-statistics. Setting this parameter to `FALSE` is equivalent to the ADMIXTOOLS `inbreed: NO` option. Setting `adjust_pseudohaploid` to an integer `n` will check the first `n` SNPs instead of the first 1000 SNPs.
#' @param verbose Print progress updates
#' @return A data frame with per-block f4-statistics for each population quadruple.
f3blockdat_from_geno = function(pref, popcombs, auto_only = TRUE,
                                blgsize = 0.05,
                                block_lengths = NULL, allsnps = FALSE, adjust_pseudohaploid = TRUE,
                                poly_only = FALSE, apply_corr = TRUE, outgroupmode = FALSE, verbose = TRUE) {

  if(is.matrix(popcombs)) {
    if(ncol(popcombs) != 3) stop("'popcombs' is a matrix but doens't have three columns!")
    popcombs %<>% as_tibble(.name_repair = ~paste0('pop', 1:3))
  }

  if('model' %in% names(popcombs)) {
    hasmodels = TRUE
  } else {
    hasmodels = FALSE
    popcombs %<>% mutate(model = 1)
  }
  if(allsnps) {
    modelvec = 0
    pc = popcombs %>% select(pop1:pop3) %>% distinct
  } else {
    modelvec = popcombs$model
    pc = popcombs
  }

  if(verbose) alert_info('Reading metadata...\n')
  pref = normalizePath(pref, mustWork = FALSE)
  l = format_info(pref)

  indfile = read_table2(paste0(pref, l$indend), col_names = l$indnam, col_types = l$indtype, progress = FALSE)
  snpfile = read_table2(paste0(pref, l$snpend), col_names = l$snpnam, col_types = 'ccddcc', progress = FALSE)
  cpp_read_geno = l$cpp_read_geno
  fl = paste0(pref, l$genoend)

  nsnpall = nrow(snpfile)
  nindall = nrow(indfile)
  snpfile$keep = TRUE
  if(auto_only) snpfile %<>% mutate(keep = as.numeric(gsub('[a-zA-Z]+', '', CHR)) <= 22)
  #if('keepsnps' %in% names(list(...))) snpfile %<>% mutate(keep = keep & SNP %in% list(...)$keepsnps)
  #if(!is.null(keepsnps)) snpfile %<>% mutate(keep = keep & SNP %in% keepsnps)
  nsnpaut = sum(snpfile$keep)
  pops = unique(c(popcombs$pop1, popcombs$pop2, popcombs$pop3))

  if(!all(pops %in% indfile$pop))
    stop(paste0('Populations missing from indfile: ', paste0(setdiff(pops, indfile$pop), collapse = ', ')))
  if(!is.null(block_lengths) && sum(block_lengths) != nsnpaut)
    stop(paste0('block_lengths should sum to ', nsnpaut,' (the number of autosomal SNPs)'))
  if(any(duplicated(indfile$ind)))
    stop('Duplicate individual IDs are not allowed!')

  allinds = indfile$ind
  allpops = indfile$pop
  indfile %<>% filter(pop %in% pops)
  indvec = (allinds %in% indfile$ind)+0
  popvec = match(indfile$pop, pops)
  p1 = match(pc$pop1, pops)
  p2 = match(pc$pop2, pops)
  p3 = match(pc$pop3, pops)
  ntest = if(is.numeric(adjust_pseudohaploid)) adjust_pseudohaploid else 1000
  if(adjust_pseudohaploid) ploidy = l$cpp_geno_ploidy(fl, nsnpall, nindall, indvec, ntest)

  if(verbose) alert_info(paste0('Computing block lengths for ', sum(snpfile$keep),' SNPs...\n'))
  if(is.null(block_lengths)) block_lengths = get_block_lengths(snpfile %>% filter(keep), blgsize = blgsize)
  numblocks = length(block_lengths)

  snpfile$block = NA; snpfile$block[snpfile$keep] = rep(1:length(block_lengths), block_lengths)
  snpfile %<>% fill(block, .direction = 'updown')
  snpind = split(snpfile$keep, snpfile$block)
  bl = rle(snpfile$block)$lengths
  start = lag(cumsum(bl), default = 0)
  end = cumsum(bl)

  popind = popcombs %>%
    group_by(model) %>%
    summarize(pp = list(unique(c(pop1, pop2, pop3)))) %>%
    mutate(popind = map(pp, ~match(., pops))) %$% popind
  usesnps = matrix(0)

  numer = denom = cnt = matrix(NA, numblocks, nrow(pc))
  for(i in 1:numblocks) {
    if(verbose) alert_info(paste0('Computing ', nrow(pc),' f3-statistics for block ',
                                  i, ' out of ', numblocks, '...\r'))
    # replace following two lines with cpp_geno_to_afs?
    gmat = cpp_read_geno(fl, nsnpall, nindall, indvec, start[i], end[i], T, F)[,snpind[[i]],drop=F]
    ref = rowsum(gmat, popvec, na.rm = TRUE)
    at = ref / rowsum((!is.na(gmat))+0, popvec) / 2
    if(!allsnps) {
      # get SNP index matrix which will be used for all f-stats
      # otherwise, delegate decision to cpp_aftable_to_dstatnum and choose different snps for each f-stat
      usesnps = popind %>% map(~(colSums(!is.finite(at[.,,drop=FALSE])) == 0)+0) %>% do.call(rbind, .)
      if(poly_only) {
        fun = function(mat) {
          #apply(mat, 2, function(x) length(na.omit(unique(x))) > 1)+0
          apply(mat, 2, function(x) length(na.omit(unique(x))) > 1 | !suppressWarnings(max(na.omit(x))) %in% c(0,1))+0
        }
        usesnps2 = popind %>% map(~fun(at)) %>% do.call(rbind, .)
        usesnps = usesnps & usesnps2
      }
    }
    num = cpp_aftable_to_dstatnum(at, p1, p2, p1, p3, modelvec, usesnps, allsnps, poly_only)
    if(apply_corr || !outgroupmode) {
      gmatinv = 2-gmat
      gmatplo = gmat
      gmatploinv = gmatinv
      if(adjust_pseudohaploid) {
        gmatplo[ploidy[indvec == 1] == 1, ] = gmatplo[ploidy[indvec == 1] == 1, ]/2
        gmatploinv[ploidy[indvec == 1] == 1, ] = gmatploinv[ploidy[indvec == 1] == 1, ]/2
      }
      ref = rowsum(gmatplo, popvec, na.rm = TRUE)
      alt = rowsum(gmatploinv, popvec, na.rm = TRUE)
      tot = ref+alt
      h = (ref*alt)/(tot*(tot-1))
      h1 = h[p1,,drop=F]
      if(apply_corr) {
        corr1 = h1/tot[p1,,drop=F]
        corr1[p1 == p2 | p1 == p3,] = 0
        num$num = num$num - corr1
      }
      if(apply_corr == 2) {
        corr2 = h[p2,,drop=F]/tot[p2,,drop=F]
        corr2[p2 != p3,] = 0
        num$num = num$num - corr2
      }
      if(!outgroupmode) {
        h1[c(!is.finite(num$num))] = NA
        denom[i,] = unname(rowSums(2*h1, na.rm = TRUE))
      }
    }
    num$num[!is.finite(num$num)] = NA
    if(outgroupmode) denom[i,] = rowSums(is.finite(num$num))
    numer[i,] = unname(rowSums(num$num, na.rm = TRUE))
    cnt[i,] = c(num$cnt)
  }
  if(verbose) cat('\n')
  out = pc %>%
    expand_grid(block = 1:numblocks) %>%
    mutate(numer = c(numer/cnt), denom = c(denom/cnt), est = numer/denom, n = c(cnt))
  if(allsnps) out = popcombs %>% left_join(out, by = paste0('pop', 1:3))
  if(!hasmodels) out$model = NULL

  out %>%
    mutate(length = block_lengths[block], est = if_else(is.finite(est) & n > 0, est, NA_real_))
}



construct_fstat_matrix = function(popcomb) {

  pops = unique(unlist(popcomb))
  npop = length(pops)
  n2 = choose(npop, 2)

  indmat = matrix(NA, npop, npop)
  indmat[lower.tri(indmat)] = 1:choose(npop,2)
  indmat = t(indmat)
  indmat[lower.tri(indmat)] = 1:choose(npop,2)

  out = matrix(0, nrow(popcomb), n2)
  ind4 = apply(as.matrix(popcomb), 2, function(x) match(x, pops))
  for(i in 1:nrow(popcomb)) {
    x = ind4[i,]
    out[i, indmat[x[1],x[4]]] = 1
    out[i, indmat[x[2],x[3]]] = 1
    out[i, indmat[x[1],x[3]]] = -1
    out[i, indmat[x[2],x[4]]] = -1
    if(popcomb[i,1] == popcomb[i,3] && popcomb[i,2] == popcomb[i,4]) out[i,] = out[i,]*2
  }
  out/2
}


#' Get smoothed f2-statistics
#'
#' This function returns an array of (pseudo-) f2-statistics which are computed by
#' taking into account other f2-, f3-, and f4-statistics. The advantage of doing that
#' is that f3- and f4-statistics computed from these smoothed f2-statistics can be
#' more accurate for populations with large amounts of missing data. The function
#' uses SNPs which are missing in some populations in a manner which tends to introduce
#' less bias than setting `maxmiss` to values greater than 0.
#' @export
#' @param pref Prefix of genotype files
#' @param pops Populations for which to compute f2-statistics
#' @param include_f2 Should f2-statistics be used to get smoothed f2-statistics?
#' If `include_f2` is a positive integer, it specifies how many randomly chosen f2-statistics should be used.
#' @param include_f3 Should f3-statistics be used to get smoothed f2-statistics?
#' If `include_f3` is a positive integer, it specifies how many randomly chosen f3-statistics should be used.
#' @param include_f4 Should f4-statistics be used to get smoothed f2-statistics?
#' If `include_f4` is a positive integer, it specifies how many randomly chosen f4-statistics should be used.
#' @return A 3d-array of smoothed f2-statistics
#' @examples
#' \dontrun{
#' f2_blocks = qpfstats(geno_prefix, mypops)
#' }
qpfstats = function(pref, pops, include_f2 = TRUE, include_f3 = TRUE, include_f4 = TRUE, verbose = TRUE) {

  popcomb = tibble()
  sp = sort(pops)
  npop = length(pops)
  if(include_f2) {
    pcf2 = expand_grid(pop1 = sp, pop2 = sp) %>% filter(pop1 < pop2) %>%
      mutate(pop3 = pop1, pop4 = pop2)
    if(is.numeric(include_f2)) pcf2 %<>% slice_sample(n = include_f2)
    popcomb %<>% bind_rows(pcf2)
  }
  if(include_f3) {
    pcf3 = fstat_get_popcombs(NULL, pops, fnum = 3) %>%
      transmute(pop1, pop2, pop4 = pop3, pop3 = pop1)
    if(is.numeric(include_f3)) pcf3 %<>% slice_sample(n = include_f3)
    popcomb %<>% bind_rows(pcf3)
  }
  if(include_f4) {
    pcf4 = fstat_get_popcombs(NULL, pops, fnum = 4)
    if(is.numeric(include_f4)) pcf4 %<>% slice_sample(n = include_f4)
    popcomb %<>% bind_rows(pcf4)
  }
  f4blockdat = f4blockdat_from_geno(pref, popcomb, allsnps = TRUE)
  f4pass1 = f4blockdat %>% f4blockdat_to_f4out(FALSE)

  if(verbose) alert_info(paste0('Constructing matrix...\n'))
  x = construct_fstat_matrix(popcomb)
  ymat = f4blockdat %>%
    select(pop1:pop4, block, est) %>%
    pivot_wider(c(1:4), block, values_from = est) %>% select(-1:-4) %>% as.matrix
  y = f4pass1$est
  #ymat = ymat/f4pass1$se
  #x = x/f4pass1$se
  if(verbose) alert_info(paste0('Running regression...\n'))
  lh = solve((t(x) %*% x) + diag(ncol(x))*0.00001) %*% t(x)
  b = lh %*% ymat
  bglob = lh %*% y

  nblocks = ncol(b)
  bl = f4blockdat %>% slice(1:nblocks) %>% pull(length)
  f2blocks = array(0, c(npop, npop, nblocks), list(sp, sp, paste0('l', bl)))
  m = matrix(1:npop^2, npop, npop)
  m2 = matrix(1:npop^2, npop, npop, byrow = T)
  add = rep((0:(nblocks-1))*npop^2, each = choose(npop,2))
  ind1 = rep(m2[lower.tri(m2)], nblocks) + add
  ind2 = rep(m[lower.tri(m)], nblocks) + add
  f2blocks[ind1] = f2blocks[ind2] = c(b)

  f2mat = f2mat2 = matrix(0, npop, npop)
  f2mat[m2[lower.tri(m2)]] = f2mat[m[lower.tri(m)]] = bglob
  f2mat2[m2[lower.tri(m2)]] = f2mat2[m[lower.tri(m)]] = f2(f2blocks)$est
  f2blocks2 = f2blocks - c(f2mat2) + c(f2mat)

  f2blocks2
}



#' Convert graph to dot format
#' @export
#' @param graph Graph as igraph object or edge list (columns labelled 'from', 'to', 'weight')
#' @param outfile Output file name
#' @examples
#' \dontrun{
#' results = qpgraph(example_f2_blocks, example_graph)
#' write_dot(results$edges)
#' }
write_dot = function(graph, outfile = stdout(), size1 = 7.5, size2 = 10,
                     title = '', dot2pdf = FALSE) {
  # writes qpgraph output to a dot format file

  if('igraph' %in% class(graph)) {
    edges = graph %>% as_edgelist %>% as_tibble(.name_repair = ~c('from', 'to')) %>%
      add_count(to) %>% mutate(type = ifelse(n == 1, 'edge', 'admix')) %>% select(-n)
  } else edges = graph
  if(!'weight' %in% names(edges)) edges %<>% mutate(weight = 0)

  leaves = setdiff(edges$to, edges$from)
  root = setdiff(edges$from, edges$to)
  internal = setdiff(c(edges$from, edges$to), c(leaves))
  edges = mutate(edges, lab = ifelse(type == 'edge',
                                     paste0(' [ label = "', round(weight * 1000), '" ];'),
                                     paste0(' [ style=dotted, label = "', round(weight * 100), '%" ];')),
                 from = str_replace_all(from, '[\\.-]', ''),
                 to = str_replace_all(to, '[\\.-]', ''))
  nodes = paste0(internal, ' [shape = point];', collapse = '\n')

  out = paste0('digraph G {\nlabel = "',title,'";\nlabelloc=t;\nlabeljust=l;\n')
  out = paste0(out, 'size = "',size1,',',size2,'";\n')
  out = paste0(out, nodes, '\n')
  out = paste0(out, paste(edges$from, ' -> ', edges$to, edges$lab, collapse = '\n'))
  out = paste0(out, '\n}')

  writeLines(out, outfile)
  if(dot2pdf) {
    outpdf = str_replace(outfile, '\\..*$', '.pdf')
    system(paste0('dot -Tpdf < ', outfile,' > ', outpdf))
  }
}


geno_to_treemix = function(pref, outfile = paste0(pref, '.txt.gz'), pops = NULL, verbose = TRUE) {

  afs = anygeno_to_aftable(pref, verbose = verbose)
  if(!is.null(pops)) {
    afs$afs = afs$afs[,pops]
    afs$counts = afs$counts[,pops]
  }
  cnt = rowSums(afs$counts > 0)
  nomiss = cnt == ncol(afs$counts)
  m1 = replace_na(afs$counts * afs$afs, 0)[nomiss,]
  m2 = replace_na(afs$counts * (1-afs$afs), 0)[nomiss,]

  if(verbose) alert_info('Writing data in treemix format...\n')
  paste0(m1, ',', m2) %>% matrix(nrow(m1)) %>% as_tibble(.name_repair = ~colnames(m1)) %>%
    write_delim(outfile, delim=' ')
}


graph_to_treemix = function(edges, outpref = NULL) {

  # doesn't work for complex graphs. have to define tree from graph first; test how much choice of tree affects fit

  # simplify admixture nodes
  root = setdiff(edges$from, edges$to)
  connectors = setdiff(names(which(table(c(edges$from, edges$to)) == 2)), root)
  remove = edges %>% filter(to %in% connectors) %>% select(to, from) %>% deframe
  edges %<>% mutate(from = ifelse(from %in% connectors, remove[from], from)) %>% filter(!to %in% connectors)

  # simplify admixed leaves
  leaves = setdiff(edges$to, edges$from)
  for(i in 1:10) {
    admnodes = intersect(names(which(table(edges$to) == 2)), names(which(table(edges$from) == 1)))
    remove = edges %>% filter(to %in% leaves, from %in% admnodes) %>% select(from, to) %>% deframe
    edges %<>% mutate(to = ifelse(to %in% names(remove), remove[to], to)) %>% filter(!from %in% names(remove))
  }

  nodes = union(edges$from, edges$to)
  edges %<>% add_count(to) %>% group_by(from) %>%
    mutate(weight = ifelse(type == 'admix', weight, 1), mxw = weight == max(weight, na.rm=T)) %>% ungroup %>%
    group_by(to) %>% mutate(mig = weight != max(weight) & n > 1 & !mxw) %>% ungroup %>% mutate(weight = weight/sum(weight)) %>% ungroup

  mig = edges %>% filter(mig) %>% pull(from)

  edges = edges %>%
    transmute(from = match(from, nodes), to = match(to, nodes),
              l = ifelse(mig, 0, 1),
              w = ifelse(n > 1, weight, 1),
              m = ifelse(mig, 'MIG', 'NOT_MIG'))

  #tree = edges %>% filter(m == 'NOT_MIG') %>% mutate(from = paste0('n', from), to = paste0('n', to)) %>%
  #  edges_to_igraph()
  #parents = map(V(tree), ~names(neighbors(tree, ., mode = 'in')))
  nodes = tibble(i = (1:length(nodes)), nodes,
         r = ifelse(nodes == root, 'ROOT', 'NOT_ROOT'),
         m = ifelse(nodes %in% mig, 'MIG', 'NOT_MIG'),
         t = ifelse(nodes %in% leaves, 'TIP', 'NOT_TIP')) %>%
    mutate(nodes = ifelse(nodes %in% leaves, nodes, NA))

  if(!is.null(outpref)) {
    write_delim(edges, paste0(outpref, '.edges.gz'), col_names = FALSE, delim = ' ')
    write_delim(nodes, paste0(outpref, '.vertices.gz'), col_names = FALSE, delim = ' ')
  } else {
    return(namedList(edges, nodes))
  }
}


edges_to_treemix = function(edges, outpref) {


  graph = edges %>% select(1:2) %>% as.matrix %>% graph_from_edgelist()
  graph %<>% simplify_graph
  spl = split_graph(edges)
  del = paste(spl$deleted[,1], ' ', spl$deleted[,2])
  kep = paste(spl$kept[,1], ' ', spl$kept[,2])

  leaves = setdiff(edges$to, edges$from)
  root = setdiff(edges$from, edges$to)
  edges %<>% mutate(ee = paste(from, ' ', to)) %>% mutate(mig = ee %in% del)

  # simplify admixture nodes
  root = setdiff(edges$from, edges$to)
  connectors = setdiff(names(which(table(c(edges$from, edges$to)) == 2)), root)
  remove = edges %>% filter(to %in% connectors) %>% select(to, from) %>% deframe
  edges %<>% mutate(from = ifelse(from %in% connectors, remove[from], from)) %>% filter(!to %in% connectors)

  # simplify admixed leaves
  for(i in 1:10) {
    admnodes = intersect(names(which(table(edges$to) == 2)), names(which(table(edges$from) == 1)))
    remove = edges %>% filter(to %in% leaves, from %in% admnodes) %>% select(from, to) %>% deframe
    edges %<>% mutate(to = ifelse(to %in% names(remove), remove[to], to)) %>% filter(!from %in% names(remove))
  }

  # H to X
  hnodes1 = names(which(table(edges$from) == 1))
  hnodes2 = edges %>% filter(from %in% hnodes1) %>% pull(to)
  names(hnodes2) = hnodes1
  edges %<>% filter(!from %in% hnodes1) %>%
    mutate(to = ifelse(to %in% hnodes1, hnodes2[to], to))

  mignodes = edges %>% filter(mig) %>% pull(from)

  # shift mig targets up to next non-mig node
  while(edges %>% filter(to %in% mignodes, mig) %>% nrow > 0) {
    nd = edges %>% filter(to %in% mignodes, mig) %>% pull(to)
    prev = edges %>% filter(!mig, to %in% nd) %>% select(to, from) %>% deframe
    edges %<>% mutate(to = ifelse(to %in% nd & mig, prev[nd], to))
  }

  e = edges

  nodes = union(e$from, e$to)

  e %>%
    transmute(from = match(from, nodes), to = match(to, nodes),
              l = ifelse(mig, 0, 1),
              w = ifelse(type == 'admix', weight, 1),
              m = ifelse(mig, 'MIG', 'NOT_MIG')) %>%
    write_delim(paste0(outpref, '.edges.gz'), col_names = FALSE, delim = ' ')

  tibble(i = (1:length(nodes)), nodes,
         r = ifelse(nodes == root, 'ROOT', 'NOT_ROOT'),
         m = ifelse(nodes %in% mignodes, 'MIG', 'NOT_MIG'),
         t = ifelse(nodes %in% leaves, 'TIP', 'NOT_TIP')) %>%
    mutate(nodes = ifelse(nodes %in% leaves, nodes, NA)) %>%
    write_delim(paste0(outpref, '.vertices.gz'), col_names = FALSE, delim = ' ')
}

#' Read graph in dot format
#' @export
#' @param dotfile Name of a file with a dot formatted admixture graph
#' @examples
#' \dontrun{
#' graph = parse_dot("/my/graph.dot")
#' }
parse_dot = function(dotfile) {
  dotfile %>%
    read_lines %>%
    str_subset('->') %>%
    str_split('->') %>%
    map(~str_squish(.) %>% str_replace_all(' .+', '') %>% str_squish()) %>%
    do.call(rbind, .) %>%
    igraph::graph_from_edgelist()
}



