
#' Read allele frequencies from packedancestrymap files
#'
#' @export
#' @param pref Prefix of packedancestrymap files (files have to end in `.geno`, `.ind`, `.snp`)
#' @param inds Individuals from which to compute allele frequencies
#' @param pops Populations from which to compute allele frequencies. If `NULL` (default), populations will be extracted from the third column in the `.ind` file. If population labels are provided, they should have the same length as `inds`, and will be matched to them by position
#' @param adjust_pseudohaploid Genotypes of pseudohaploid samples are usually coded as `0` or `2`, even though only one allele is observed. `adjust_pseudohaploid` ensures that the observed allele count increases only by `1` for each pseudohaploid sample. If `TRUE` (default), samples that don't have any genotypes coded as `1` among the first 1000 SNPs are automatically identified as pseudohaploid. This leads to slightly more accurate estimates of f-statistics. Setting this parameter to `FALSE` is equivalent to the ADMIXTOOLS `inbreed: NO` option.
#' @param verbose Print progress updates
#' @return A list with three data frames: allele frequency data, allele counts, and SNP metadata
#' @examples
#' \dontrun{
#' afdat = packedancestrymap_to_aftable(prefix, pops = pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
packedancestrymap_to_aftable = function(pref, inds = NULL, pops = NULL, adjust_pseudohaploid = TRUE, verbose = TRUE) {


  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # pops is vector of populations for which to calculate AFs
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population

  if(verbose) alert_info('Reading allele frequencies from packedancestrymap files...\n')

  pref %<>% normalizePath(mustWork = FALSE)
  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = 'cccccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = 'ccnncc', progress = FALSE)
  nindall = nrow(indfile)
  nsnpall = as.numeric(nrow(snpfile))

  ip = match_samples(indfile$X1, indfile$X3, inds, pops)
  indvec = ip$indvec - 1
  upops = ip$upops

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnpall, ' SNPs\n'))
    alert_info(paste0('Calculating allele frequencies from ', sum(indvec != -1), ' samples in ', length(upops), ' populations\n'))
    alert_info(paste0('Expected size of allele frequency data: ', round((nsnpall*length(upops)*8+nsnpall*112)/1e6), ' MB\n'))
    # 8, 112: estimated scaling factors for AF columns and annotation columns
  }

  if(adjust_pseudohaploid) ploidy = cpp_packedancestrymap_ploidy(paste0(pref, genoend), nsnpall, nindall, indvec)
  else ploidy = rep(2, nindall)
  afdat = cpp_packedancestrymap_to_aftable(paste0(pref, '.geno'), nsnpall, nindall, indvec, first = 0,
                                           last = nsnpall, ploidy = ploidy,
                                           transpose = FALSE, verbose = verbose)

  if(verbose) {
    alert_success(paste0(nrow(afdat$afs), ' SNPs read in total\n'))
  }

  colnames(afdat$afs) = colnames(afdat$counts) = upops
  rownames(afdat$afs) = rownames(afdat$counts) = snpfile$SNP

  afdat$snpfile = snpfile
  afdat
}



#' Read allele frequencies from ancestrymap files
#'
#' @export
#' @param pref Prefix of ancestrymap files (files have to end in `.geno`, `.ind`, `.snp`)
#' @inheritParams packedancestrymap_to_aftable
#' @return A list with three data frames: allele frequency data, allele counts, and SNP metadata
#' @examples
#' \dontrun{
#' afdat = ancestrymap_to_aftable(prefix, pops = pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
ancestrymap_to_aftable = function(pref, inds = NULL, pops = NULL, adjust_pseudohaploid = TRUE, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # pops is vector of populations for which to calculate AFs
  # defaults to third column in ind file
  # inds: instead of specifying a list of populations for which to calculate AFs, you can specify a list of individuals
  # returns data.frame; first 6 columns: snpfile; remaining columns: AF for each population

  if(verbose) alert_info('Reading allele frequencies from ancestrymap files...\n')

  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = 'cccccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = 'ccnncc', progress = FALSE)

  ip = match_samples(indfile$X1, indfile$X3, inds, pops)
  indvec = ip$indvec
  upops = ip$upops
  popind2 = which(indvec > 0)
  pops = upops[indvec]

  fl = paste0(pref, '.geno')
  geno = apply(do.call(rbind, str_split(readLines(fl), '')), 2, as.numeric)
  nindall = ncol(geno)
  nsnp = nrow(geno)
  numpop = length(upops)
  geno = geno[,popind2]
  colnames(geno) = inds
  rownames(geno) = snpfile$SNP

  ploidy = apply(geno, 1, function(x) max(1, length(unique(na.omit(x)))-1))
  if(adjust_pseudohaploid) counts = t(rowsum((!is.na(t(geno)))*ploidy, pops))
  else counts = t(rowsum((!is.na(t(geno)))*2, pops))
  afs = t(rowsum(t(geno)/(3-ploidy), pops, na.rm=T))/counts

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnp, ' SNPs\n'))
    alert_info(paste0('Calculating allele frequencies from ', length(popind2), ' samples in ', numpop, ' populations\n'))
  }

  outlist = namedList(afs, counts, snpfile)
  outlist
}


discard_from_aftable = function(afdat, maxmiss = 1, minmaf = 0, maxmaf = 0.5, auto_only = TRUE,
                                poly_only = FALSE, transitions = TRUE, transversions = TRUE, keepsnps = NULL) {
  # afdat is list with 'snpfile', 'afs', 'counts'
  # returns same list with SNPs removed
  # keepsnps overrides maxmiss and auto_only
  # maxmiss = 1 is equivalent to na.action = 'none'
  # maxmiss = 0 is equivalent to na.action = 'remove'

  snpdat = afdat$snpfile
  if(maxmiss < 1) snpdat$miss = rowMeans(afdat$counts == 0)
  else snpdat %<>% mutate(miss = 0)
  if(minmaf > 0 | maxmaf < 0.5) snpdat %<>% mutate(af = rowMeans(afdat$afs, na.rm=TRUE)/2, maf = pmin(af, 1-af))
  else snpdat %<>% mutate(af = 0.2, maf = 0.2)

  if(poly_only) snpdat %<>% mutate(poly = cpp_is_polymorphic(afdat$afs))
  else snpdat %<>% mutate(poly = TRUE)


  remaining = discard_snps(snpdat, maxmiss = maxmiss, auto_only = auto_only, poly_only = poly_only,
                           minmaf = minmaf, maxmaf = maxmaf,
                           transitions = transitions, transversions = transversions, keepsnps = keepsnps)
  #keeprows = match(remaining, snpdat[['SNP']])
  map(afdat, ~.[remaining,,drop = FALSE])
}


discard_from_geno = function(geno, maxmiss = 1, auto_only = TRUE, poly_only = TRUE,
                             minmaf = 0, maxmaf = 0.5,
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

  if(poly_only) snpdat %<>% mutate(poly = cpp_is_polymorphic(afdat$afs))
  else snpdat %<>% mutate(poly = TRUE)

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
                        minmaf = 0, maxmaf = 0.5, transitions = TRUE, transversions = TRUE) {
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
  snpdat %>%
    filter(
    miss <= maxmiss,
    between(maf, minmaf, maxmaf),
    !auto_only | as.numeric(gsub('[a-zA-Z]+', '', CHR)) <= 22,
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
read_packedancestrymap = function(pref, inds = NULL, first = 1, last = Inf,
                                  transpose = FALSE, verbose = TRUE) {
  # pref is the prefix for packedancestrymap files (ending in .geno, .snp, .ind)
  # inds: optional vector of individual IDs
  # returns list with geno (genotype matrix), snp (snp metadata), ind (sample metadata).

  pref = normalizePath(pref, mustWork = FALSE)
  nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = 'cccccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = 'ccnncc', skip = first-1,
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
#' @param inds Individuals for which data should be read. Defaults to all individuals
#' @inheritParams packedancestrymap_to_aftable
#' @return A list with the genotype data matrix, the `.ind` file, and the `.snp` file
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
  indfile = read_table2(paste0(pref, '.ind'), col_names = FALSE, col_types = 'cccccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, '.snp'), col_names = nam, col_types = 'ccnncc', progress = FALSE)
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
#' @inheritParams packedancestrymap_to_aftable
#' @return a list with two items: allele frequency data and individual counts.
#' @examples
#' \dontrun{
#' afdat = plink_to_aftable(prefix, pops)
#' afs = afdat$afs
#' counts = afdat$counts
#' }
plink_to_aftable = function(pref, inds = NULL, pops = NULL, adjust_pseudohaploid = TRUE, verbose = TRUE) {
  # This is based on Gad Abraham's "plink2R" package
  # Modified to return per-group allele frequencies rather than raw genotypes.

  if(verbose) alert_info('Reading allele frequencies from PLINK files...\n')

  bedfile = paste0(pref, '.bed')
  famfile = paste0(pref, '.fam')
  bimfile = paste0(pref, '.bim')
  nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
  bim = read_table2(bimfile, col_names = nam, progress = FALSE, col_types = 'ccnncc')
  fam = read_table2(famfile, col_names = FALSE, progress = FALSE, col_types = 'cccccc')
  nsnpall = nrow(bim)
  nindall = nrow(fam)

  ip = match_samples(fam$X2, fam$X1, inds, pops)
  indvec = ip$indvec
  upops = ip$upops

  # if(!is.null(inds) && !is.null(pops)) {
  #   stopifnot(nrow(fam) == length(inds)) # note: fix this, so that this condition is not necessary
  #   fam$X1 = pops
  #   fam$X2 = inds
  # }
  # if(is.null(inds)) inds = fam$X2
  # if(is.null(pops)) pops = fam$X1
  # indvec = pop_indices(fam, pops = pops, inds = inds)

  indvec2 = which(indvec > 0)
  #keepinds = unique(fam[[is.null(pops)+1]][indvec > 0])

  if(verbose) {
    alert_info(paste0(basename(pref), '.geno has ', nindall, ' samples and ', nsnpall, ' SNPs\n'))
    alert_info(paste0('Calculating allele frequencies from ', sum(indvec != 0), ' samples in ', length(upops), ' populations\n'))
    alert_info(paste0('Expected size of allele frequency data: ', round((nsnpall*length(upops)*8+nsnpall*112)/1e6), ' MB\n'))
    # 8, 112: estimated scaling factors for AF columns and annotation columns
  }

  afs = cpp_read_plink_afs(normalizePath(bedfile), indvec, indvec2, adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)

  afmatrix = afs[[1]]
  countmatrix = afs[[2]]
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
#' This is based on Gad Abraham's \href{https://github.com/gabraham/plink2R}{plink2R} package,
#' but genotypes are `m` x `n`, not `n` x `m`.
#' See \href{https://www.rdocumentation.org/packages/genio}{genio} for a dedicated `R` package for
#' reading and writing `PLINK` files.
#' @export
#' @param inds Individuals for which data should be read. Defaults to all individuals
#' @inheritParams packedancestrymap_to_aftable
#' @return A list with the genotype data matrix, the `.ind` file, and the `.snp` file
#' @examples
#' \dontrun{
#' samples = c('Ind1', 'Ind2', 'Ind3')
#' geno = read_packedancestrymap(prefix, samples)
#' }
read_plink = function(pref, inds = NULL, verbose = FALSE) {
  # This is based on Gad Abraham's "plink2R" package, but genotypes are m x n, not n x m
  if(verbose) alert_info('Reading PLINK files...\n')

  bedfile = paste0(pref, '.bed')
  famfile = paste0(pref, '.fam')
  bimfile = paste0(pref, '.bim')
  nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
  bim = read_table2(bimfile, col_names = nam, col_types = 'ccnncc', progress = FALSE)
  fam = read_table2(famfile, col_names = FALSE, col_types = 'cccccc', progress = FALSE)
  if(is.null(inds)) inds = fam[[2]]
  indvec = pop_indices(fam, pops = NULL, inds = inds)
  indvec2 = which(indvec > 0)
  #indvec2 = which(fam[[2]] %in% inds)
  keepinds = fam[[2]][indvec2]

  g = 2 * cpp_read_plink_afs(normalizePath(bedfile), indvec=indvec, indvec2,
                             adjust_pseudohaploid = TRUE, verbose = verbose)[[1]]
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
#' @param f2_arrs 3d arrays with blocked f2, allele frequency products, and counts for each population pair.
#' The first two dimensions of each array have to have population names.
#' @param outdir Directory where data will be stored
#' @param overwrite Overwrite existing files in `outdir`
#' @seealso \code{\link{read_f2}}
#' @examples
#' \dontrun{
#' write_f2(f2_arrs, outdir = 'path/to/f2stats/')
#' }
write_f2 = function(f2_arrs, outdir, overwrite = FALSE) {

  if(!dir.exists(outdir)) dir.create(outdir)
  d1 = dim(f2_arrs[[1]])[1]
  d2 = dim(f2_arrs[[1]])[2]
  nam1 = dimnames(f2_arrs[[1]])[[1]]
  nam2 = dimnames(f2_arrs[[1]])[[2]]
  for(i in seq_len(d1)) {
    for(j in seq_len(d2)) {
      pop1 = min(nam1[i], nam2[j])
      pop2 = max(nam1[i], nam2[j])
      if(pop1 <= pop2) {
        mat1 = cbind(f2 = as.vector(f2_arrs$f2[i, j, ]), counts = as.vector(f2_arrs$counts[i, j, ]))
        mat2 = cbind(afprod = as.vector(f2_arrs$afprod[i, j, ]), countsap = as.vector(f2_arrs$countsap[i, j, ]))
        dir = paste0(outdir, '/', pop1, '/')
        fl1 = paste0(dir, pop2, '_f2.rds')
        fl2 = paste0(dir, pop2, '_ap.rds')
        if(!dir.exists(dir)) dir.create(dir)
        if(!file.exists(fl1) || overwrite) saveRDS(mat1, file = fl1)
        if(!file.exists(fl2) || overwrite) saveRDS(mat2, file = fl2)
      }
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
#' @param afprod Return allele frequency products instead of f2 estimates
#' @param remove_na Remove blocks with missing values
#' @param verbose Print progress updates
#' @return A 3d array of block jackknife estimates
#' @seealso \code{\link{write_f2}}
#' @examples
#' \dontrun{
#' read_f2(f2_dir, pops = c('pop1', 'pop2', 'pop3'))
#' }
read_f2 = function(f2_dir, pops = NULL, afprod = FALSE, remove_na = TRUE, verbose = FALSE) {
  # assumes f2 is in first column, afprod in second column

  if(is.null(pops)) {
    pops = list.dirs(f2_dir, full.names = FALSE, recursive = FALSE)
  }
  if(afprod) block_lengths = readRDS(paste0(f2_dir, '/block_lengths_a.rds'))
  else block_lengths = readRDS(paste0(f2_dir, '/block_lengths.rds'))
  numblocks = length(block_lengths)
  numpops = length(pops)
  f2_blocks = array(NA, c(numpops, numpops, numblocks))
  dimnames(f2_blocks)[[1]] = dimnames(f2_blocks)[[2]] = pops
  dimnames(f2_blocks)[[3]] = paste0('l', block_lengths)
  for(i in seq_along(pops)) {
    pop1 = pops[i]
    if(verbose) alert_info(paste0('reading ', ifelse(afprod, 'afprod', 'f2'),
                                  ' data for pop ', i, ' out of ', length(pops),'...\r'))
    for(pop2 in pops) {
      if(pop1 <= pop2) {
        pref = paste0(f2_dir, '/', pop1, '/', pop2)
        if(afprod) dat = readRDS(paste0(pref, '_ap.rds'))[,1]
        else dat = readRDS(paste0(pref, '_f2.rds'))[,1]
        f2_blocks[pop1, pop2, ] = f2_blocks[pop2, pop1, ] = dat
        if(any(is.na(dat))) warning(paste0('missing values in ', pop1, ' - ', pop2, '!'))
      }
    }
  }
  if(verbose) alert_info(paste0('\n'))
  if(remove_na) {
    keep = apply(f2_blocks, 3, function(x) sum(is.na(x)) == 0)
    if(!all(keep)) warning(paste0('Discarding ', sum(!keep), ' blocks due to missing values!'))
    f2_blocks = f2_blocks[,, keep]
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
#' @seealso \code{\link{packedancestrymap_to_aftable}}, \code{\link{write_split_f2_block}}
#' @examples
#' \dontrun{
#' afdat = packedancestrymap_to_aftable('path/to/packedancestrymap_prefix', allpopulations)
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
#' @param blgsize SNP block size in Morgan. Default is 0.05 (50 cM).
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
#'     write_split_f2_block(afdir, f2dir, chunk1 = i, chunk2 = j)
#'   }
#' }
#' }
#' # Alternatively, the following code will do the same, while submitting each chunk as a separate job.
#' # (if \code{\link[future]{plan}} has been set up appropriately)
#' \dontrun{
#' furrr::future_map(1:numchunks, ~{i=.; map(i:numchunks, ~{
#'   write_split_f2_block(afdir, f2dir, chunk1 = i, chunk2 = .)
#'   })})
#'   }
write_split_f2_block = function(afdir, outdir, chunk1, chunk2, blgsize = 0.05, overwrite = FALSE, verbose = TRUE) {
  # reads data from afdir, computes f2 jackknife blocks, and writes output to outdir

  snpdat = read_table2(paste0(afdir, '/snpdat.tsv.gz'), col_types = 'ccnncc?', progress = FALSE)
  poly = snpdat$poly
  am1 = readRDS(paste0(afdir, '/afs', chunk1, '.rds'))
  am2 = readRDS(paste0(afdir, '/afs', chunk2, '.rds'))
  cm1 = readRDS(paste0(afdir, '/counts', chunk1, '.rds'))
  cm2 = readRDS(paste0(afdir, '/counts', chunk2, '.rds'))
  nam1 = colnames(am1)
  nam2 = colnames(am2)
  nsnp = nrow(am1)

  fl = paste0(outdir, '/block_lengths.rds')
  fla = paste0(outdir, '/block_lengths_a.rds')
  if(!file.exists(fl)) {
    block_lengths = get_block_lengths(snpdat[poly,], blgsize = blgsize)
    saveRDS(block_lengths, file = fl)
  } else {
    block_lengths = readRDS(fl)
  }
  if(!file.exists(fla)) {
    block_lengths_a = get_block_lengths(snpdat, blgsize = blgsize)
    saveRDS(block_lengths_a, file = fla)
  } else {
    block_lengths_a = readRDS(fla)
  }

  filenames = expand_grid(nam1, nam2) %>%
    transmute(nam = paste0(outdir, '/', pmin(nam1, nam2), '/', pmax(nam1, nam2))) %$%
    nam %>% rep(each = 2) %>% paste0(c('_f2.rds', '_ap.rds'))
  if(all(file.exists(filenames)) && !overwrite) return()

  f2 = mats_to_f2arr(am1[poly,,drop=F], am2[poly,,drop=F], cm1[poly,,drop=F], cm2[poly,,drop=F], block_lengths)
  counts = mats_to_ctarr(am1[poly,,drop=F], am2[poly,,drop=F], cm1[poly,,drop=F], cm2[poly,,drop=F], block_lengths)
  afprod = mats_to_aparr(am1, am2, cm1, cm2, block_lengths_a)
  countsap = mats_to_ctarr(am1, am2, cm1, cm2, block_lengths_a)
  if(chunk1 == chunk2) for(i in 1:dim(f2)[1]) f2[i, i, ] = 0
  write_f2(namedList(f2, counts, afprod, countsap), outdir = outdir)
}



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
#' @param pref Prefix of packedancestrymap or PLINK files. packedancestrymap has to end in `.geno`, `.snp`, `.ind`,
#' PLINK has to end in `.bed`, `.bim`, `.fam`
#' @param outdir Directory where data will be stored
#' @param inds Individuals for which data should be extracted
#' @param pops Populations for which data should be extracted. If both `pops` and `inds` are provided, they should have the same length and will be matched by position. If only `pops` is provided, all individuals from the `.ind` or `.fam` file in those populations will be extracted. If only `inds` is provided, each indivdual will be assigned to its own population of the same name. If neither `pops` nor `inds` is provided, all individuals and populations in the `.ind` or `.fam` file will be extracted.
#' @param blgsize SNP block size in Morgan. Default is 0.05 (50 cM).
#' @param maxmem Maximum amount of memory to be used. If the required amount of memory exceeds `maxmem`, allele frequency data will be split into blocks, and the computation will be performed separately on each block pair.
#' @param maxmiss Discard SNPs which are missing in a fraction of populations higher than `maxmiss`
#' @param minmaf Discard SNPs with minor allele frequency less than `minmaf`
#' @param maxmaf Discard SNPs with minor allele frequency greater than than `maxmaf`
#' @param transitions Set this to `FALSE` to exclude transition SNPs
#' @param transversions Set this to `FALSE` to exclude transversion SNPs
#' @param keepsnps SNP IDs of SNPs to keep. Overrides other SNP filtering options
#' @param snpblocks Optional data frame with pre-assigned SNP blocks. Should have columns 'SNP' and 'block'. Can also have a column 'weight'. If present, the average 'weight' in each SNP block will be used for jackknife weights, instead of the number of SNPs in that block.
#' @param overwrite Overwrite existing files in `outdir`
#' @param format Supply this if the prefix can refer to genotype data in different formats
#' and you want to choose which one to read. Should be either `plink`, `ancestrymap`, or `packedancestrymap`
#' @param poly_only Exclude sites with identical allele frequencies in all populations
#' @param adjust_pseudohaploid Genotypes of pseudohaploid samples are usually coded as `0` or `2`, even though only one allele is observed. `adjust_pseudohaploid` ensures that the observed allele count increases only by `1` for each pseudohaploid sample. If `TRUE` (default), samples that don't have any genotypes coded as `1` among the first 1000 SNPs are automatically identified as pseudohaploid. This leads to slightly more accurate estimates of f-statistics. Setting this parameter to `FALSE` is equivalent to the ADMIXTOOLS `inbreed: NO` option.
#' @param verbose Print progress updates
#' @return SNP metadata (invisibly)
#' @seealso \code{\link{extract_f2_large}}
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' f2dir = 'my/f2dir/'
#' extract_f2(pref, f2dir, pops = c('popA', 'popB', 'popC'))
#' }
extract_f2 = function(pref, outdir, inds = NULL, pops = NULL, blgsize = 0.05, maxmem = 8000,
                      maxmiss = 1, minmaf = 0, maxmaf = 0.5, transitions = TRUE, transversions = TRUE,
                      keepsnps = NULL, snpblocks = NULL, overwrite = FALSE, format = NULL, poly_only = TRUE,
                      adjust_pseudohaploid = TRUE, verbose = TRUE) {

  outdir = normalizePath(outdir, mustWork = FALSE)
  if(length(list.files(outdir)) > 0 && !overwrite) stop('Output directory not empty! Set overwrite to TRUE if you want to overwrite files!')
  if(is.null(inds) && is.null(pops) && verbose && max(file.info(paste0(pref, '.geno'))$size, file.info(paste0(pref, '.bed'))$size, na.rm = T)/1e9 > 1) alert_danger('No poplations or individuals provided. Extracting f2-stats for all population pairs. If that takes too long, you can either specify the "pops" or "inds" parameter, or follow the example in "write_split_f2_block".')


  afdat = anygeno_to_aftable(pref, inds = inds, pops = pops, format = format,
                             adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)
  afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf,
                                  transitions = transitions, transversions = transversions,
                                  keepsnps = keepsnps, auto_only = TRUE, poly_only = FALSE)
  afdat$snpfile %<>% mutate(poly = as.logical(cpp_is_polymorphic(afdat$afs)))

  if(verbose) alert_warning(paste0(nrow(afdat$afs), ' SNPs remain after filtering. ',
                                   sum(afdat$snpfile$poly),' are polymorphic.\n'))

  afs_to_f2_blocks(afdat, outdir = outdir, overwrite = overwrite,
                   maxmem = maxmem, poly_only = poly_only, blgsize = blgsize, verbose = verbose)

  if(verbose) alert_info(paste0('Data written to ', outdir, '/\n'))
  invisible(afdat$snpfile)
}

#' Compute and store blocked f2 statistics
#'
#' `extract_f2_large` does the same as \code{\link{extract_f2}}, but it requires less memory.
#' @export
#' @inheritParams extract_f2
#' @param cols_per_chunk Number of populations per chunk. Lowering this number will lower the memory requirements when running
#' @details `extract_f2_large` requires less memory because it writes allele frequency data to disk, and doesn't store the allele frequency matrix for all populations and SNPs in memory. If you still run out of memory, reduce `cols_per_chunk`. This function is a wrapper around \code{\link{extract_afs}} and \code{\link{write_split_f2_block}}, and is slower than \code{\link{extract_f2}}. It may be faster to call \code{\link{extract_afs}} and \code{\link{write_split_f2_block}} directly, parallelizing over the different calls to \code{\link{write_split_f2_block}}.
#' @return SNP metadata (invisibly)
#' @seealso \code{\link{extract_f2}}
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' f2dir = 'my/f2dir/'
#' extract_f2_large(pref, f2dir, pops = c('popA', 'popB', 'popC'))
#' }
extract_f2_large = function(pref, outdir, inds = NULL, pops = NULL, blgsize = 0.05, cols_per_chunk = 10,
                            maxmiss = 1, minmaf = 0, maxmaf = 0.5, transitions = TRUE, transversions = TRUE,
                            keepsnps = NULL, snpblocks = NULL, overwrite = FALSE, format = NULL, poly_only = TRUE,
                            adjust_pseudohaploid = TRUE, verbose = TRUE) {

  if(verbose) alert_info(paste0('Extracting allele frequencies...\n'))
  extract_afs(pref, outdir, inds = inds, pops = pops, cols_per_chunk = cols_per_chunk, numparts = 100,
              maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf, transitions = transitions, transversions = transversions,
              keepsnps = keepsnps, format = format, poly_only = FALSE,
              adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)
  numchunks = length(list.files(outdir, 'afs.+rds'))

  if(verbose) alert_warning(paste0('Computing ', choose(numchunks, 2), ' chunk pairs. If this takes too long,
  consider running "extract_afs" and then paralellizing over "write_split_f2_block".\n'))
  for(i in 1:numchunks) {
    for(j in i:numchunks) {
      if(verbose) alert_info(paste0('Writing pair ', i, ' - ', j, '...\r'))
      write_split_f2_block(outdir, outdir, chunk1 = i, chunk2 = j, blgsize = blgsize, overwrite = overwrite)
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
      geno_to_aftable = plink_to_aftable
    } else if(all(file.exists(paste0(pref, c('.geno', '.snp', '.ind'))))) {
      if(is_binfile(paste0(pref, '.geno'))) {
        format = 'packedancestrymap'
        geno_to_aftable = packedancestrymap_to_aftable
      } else {
        format = 'ancestrymap'
        geno_to_aftable = ancestrymap_to_aftable
      }
    } else stop('Genotype files not found!')
  }
  geno_to_aftable = switch(format,
                           'plink' = plink_to_aftable,
                           'packedancestrymap' = packedancestrymap_to_aftable,
                           'ancestrymap' = ancestrymap_to_aftable)
  if(is.null(geno_to_aftable)) stop('Invalid format!')

  afdat = geno_to_aftable(pref, inds = inds, pops = pops, adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)
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
#' @param inds Individuals for which data should be read. Defaults to all individuals
#' @param maxmiss Discard SNPs which are missing in a fraction of individuals greater than `maxmiss`
#' @inheritParams extract_f2
extract_counts = function(pref, outdir, inds = NULL, blgsize = 0.05,  maxmiss = 1, minmaf = 0, maxmaf = 0.5,
                          transitions = TRUE, transversions = TRUE, keepsnps = NULL,
                          maxmem = 8000, overwrite = FALSE, format = NULL, verbose = TRUE) {
  dir.create(outdir, showWarnings = FALSE)
  if(length(list.files(outdir)) > 0 && !overwrite) stop('Output directory not empty! Set overwrite to TRUE if you want to overwrite files!')
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
                           keepsnps = keepsnps, auto_only = TRUE, poly_only = TRUE)
  randsnps = sample(1:nrow(g$bed), floor(nrow(g$bed)/2))
  g$bed[randsnps, ] = 2 - g$bed[randsnps, ]

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
#' @inheritParams extract_f2
#' @param cols_per_chunk Number of populations per chunk. Lowering this number will lower the memory requirements when running \link{\code{write_split_f2_block}}, but more chunk pairs will have to be computed.
#' @return SNP metadata (invisibly)
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' outdir = 'dir/for/afdata/'
#' extract_afs(pref, outdir)
#' }
extract_afs_old = function(pref, outdir, inds = NULL, pops = NULL, blgsize = 0.05, cols_per_chunk = 10,
                       maxmiss = 1, minmaf = 0, maxmaf = 0.5, transitions = TRUE, transversions = TRUE,
                       keepsnps = NULL, format = NULL, poly_only = FALSE, adjust_pseudohaploid = TRUE,
                       verbose = TRUE) {

  # read data and compute allele frequencies
  afdat = anygeno_to_aftable(pref, inds = inds, pops = pops, format = format,
                             adjust_pseudohaploid = adjust_pseudohaploid, verbose = verbose)
  afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf,
                                  transitions = transitions, transversions = transversions,
                                  keepsnps = keepsnps, auto_only = TRUE, poly_only = poly_only)

  afdat$snpfile %<>% mutate(poly = as.logical(cpp_is_polymorphic(afdat$afs)))
  if(verbose) alert_warning(paste0(nrow(afdat$afs), ' SNPs remain after filtering. ',
                                   sum(afdat$snpfile$poly),' are polymorphic.\n'))

  # split allele frequency data into chunks and write to disk
  split_mat(afdat$afs, cols_per_chunk = cols_per_chunk, prefix = paste0(outdir, '/afs'), verbose = verbose)
  split_mat(afdat$counts, cols_per_chunk = cols_per_chunk, prefix = paste0(outdir, '/counts'), verbose = verbose)
  # compute jackknife blocks
  block_lengths = get_block_lengths(afdat$snpfile %>% filter(poly), blgsize = blgsize, distcol = 'cm')
  block_lengths_a = get_block_lengths(afdat$snpfile, blgsize = blgsize, distcol = 'cm')
  saveRDS(block_lengths, file = paste0(outdir, '/block_lengths.rds'))
  saveRDS(block_lengths_a, file = paste0(outdir, '/block_lengths_a.rds'))
  write_tsv(afdat$snpfile, paste0(outdir, '/snpdat.tsv.gz'))
  invisible(afdat$snpfile)
}


#' Compute and store blocked allele frequency data
#'
#' Prepare data for various admixtools functions. Reads data from packedancestrymap or PLINK files,
#' and computes allele frequencies for selected populations and stores it as `.rds` files in outdir.
#' @export
#' @inheritParams extract_f2
#' @param cols_per_chunk Number of populations per chunk. Lowering this number will lower the memory requirements when running \link{\code{write_split_f2_block}}, but more chunk pairs will have to be computed.
#' @param numparts Number of parts in which genotype data will be read for computing allele frequencies
#' @return SNP metadata (invisibly)
#' @examples
#' \dontrun{
#' pref = 'my/genofiles/prefix'
#' outdir = 'dir/for/afdata/'
#' extract_afs(pref, outdir)
#' }
extract_afs = function(pref, outdir, inds = NULL, pops = NULL, cols_per_chunk = 10, numparts = 100,
                       maxmiss = 1, minmaf = 0, maxmaf = 0.5, transitions = TRUE, transversions = TRUE,
                       keepsnps = NULL, format = NULL, poly_only = FALSE, adjust_pseudohaploid = TRUE,
                       verbose = TRUE) {

  pref %<>% normalizePath(mustWork = FALSE)
  if(is_packedancestrymap_prefix(pref) || isTRUE(format == 'packedancestrymap')) {
    format = 'packedancestrymap'
    nam = c('SNP', 'CHR', 'cm', 'POS', 'A1', 'A2')
    indnam = c('ind', 'sex', 'pop')
    snpend = '.snp'
    indend = '.ind'
    genoend = '.geno'
    cpp_geno_ploidy = cpp_packedancestrymap_ploidy
    cpp_geno_to_aftable = cpp_packedancestrymap_to_aftable
  } else if(is_plink_prefix(pref) || isTRUE(format == 'plink')) {
    format = 'plink'
    nam = c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2')
    indnam = c('pop', 'ind', 'p1', 'p2', 'sex', 'pheno')
    snpend = '.bim'
    indend = '.fam'
    genoend = '.bed'
    cpp_geno_ploidy = cpp_plink_ploidy
    cpp_geno_to_aftable = cpp_plink_to_aftable
  } else stop('Genotype files not found!')

  if(verbose) alert_info(paste0('Reading metadata...\n'))
  indfile = read_table2(paste0(pref, indend), col_names = indnam, col_types = 'cccccc', progress = FALSE)
  snpfile = read_table2(paste0(pref, snpend), col_names = nam, col_types = 'ccnncc', progress = FALSE)
  nindall = nrow(indfile)
  nsnpall = nrow(snpfile)

  ip = match_samples(indfile$ind, indfile$pop, inds, pops)
  indvec = ip$indvec - 1

  snpfile %<>% mutate(CHR = as.numeric(gsub('[a-zA-Z]+', '', CHR))) %>% filter(CHR <= 22)

  starts = seq(0, nrow(snpfile), length.out = numparts+1) %>% round %>% head(-1)
  ends = c(lead(starts)[-numparts], nrow(snpfile))

  snpparts = list()
  if(adjust_pseudohaploid) ploidy = cpp_geno_ploidy(paste0(pref, genoend), nsnpall, nindall, indvec)
  else ploidy = rep(2, nindall)

  for(i in 1:numparts) {
    if(verbose) alert_info(paste0('Reading part ', i, ' out of ', numparts, '...\r'))
    # read data and compute allele frequencies
    afdat = cpp_geno_to_aftable(paste0(pref, genoend), nsnpall, nindall, indvec, first = starts[i],
                                last = ends[i], ploidy = ploidy, transpose = FALSE, verbose = FALSE)
    afdat$snpfile = snpfile %>% slice((starts[i]+1):(ends[i]))

    afdat %<>% discard_from_aftable(maxmiss = maxmiss, minmaf = minmaf, maxmaf = maxmaf,
                                    transitions = transitions, transversions = transversions,
                                    keepsnps = keepsnps, auto_only = TRUE, poly_only = poly_only)

    afdat$snpfile %<>% mutate(poly = as.logical(cpp_is_polymorphic(afdat$afs)))
    snpparts[[i]] = afdat$snpfile
    colnames(afdat$afs) = colnames(afdat$counts) = ip$upops
    rownames(afdat$afs) = rownames(afdat$counts) = afdat$snpfile$SNP

    # split allele frequency data into chunks and write to disk
    partdir = paste0(outdir, '/part',i,'/')
    dir.create(partdir, recursive = TRUE, showWarnings = FALSE)
    split_mat(afdat$afs, cols_per_chunk = cols_per_chunk, prefix = paste0(partdir,'/afs'), verbose = FALSE)
    split_mat(afdat$counts, cols_per_chunk = cols_per_chunk, prefix = paste0(partdir, '/counts'), verbose = FALSE)
  }
  if(verbose) alert_info('\n')
  snpparts %<>% bind_rows
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
      if(p1 <= p2) {
        fn1 = paste0(p1, '/', p2, '_f2.rds')
        fn2 = paste0(p1, '/', p2, '_ap.rds')
        file.copy(paste0(from, '/', fn1), paste0(to, '/', fn1))
        file.copy(paste0(from, '/', fn2), paste0(to, '/', fn2))
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
#' f2_blocks = f2_from_geno(pref, pops = c('pop1', 'pop2', 'pop3'))
#' }
f2_from_geno = function(pref, inds = NULL, pops = NULL,
                        maxmem = 8000, blgsize = 0.05, format = NULL, verbose = TRUE) {

  stopifnot(is.character(pref))
  afdat = anygeno_to_aftable(pref, inds = inds, pops = pops, format = format, verbose = verbose)
  f2_blocks = afs_to_f2_blocks(afdat, maxmem = maxmem, blgsize = blgsize, verbose = verbose)
  f2_blocks
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
#' @param afprod Return negative average allele frequency products instead of f2 estimates. This will result in more precise f4-statistics when the original data had large amounts of missingness, and should be used in that case for qpdstat and qpadm. It can also be used for outgroup f3-statistics with a fixed outgroup (for example for `qpgraph`); values will be shifted by a constant amount compared to regular f3-statistics. This shift affects the fit of a graph only by small amounts, possibly less than bias in regular f3-statistics introduced by large amounts of missing data.
#' This option is currently ineffective when reading data extracted with `extract_counts`.
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
f2_from_precomp = function(dir, inds = NULL, pops = NULL, pops2 = NULL, afprod = FALSE, return_array = TRUE,
                           apply_corr = TRUE, remove_na = TRUE, verbose = TRUE) {

  if(!is.null(pops) && !is.null(inds) && length(pops) != length(inds)) stop("'pops' and 'inds' are not the same length!")
  indpairs = dir.exists(paste0(dir, '/indivs'))
  if(!is.null(pops2)) return(f2_from_precomp_nonsquare(dir, pops, pops2, afprod = afprod,
                                                       remove_na = remove_na, verbose = verbose))

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
    f2_blocks = read_f2(dir, pops, afprod = FALSE, remove_na = remove_na, verbose = verbose)
    if(afprod) {
      f2_blocks_ap = -2*read_f2(dir, pops, afprod = TRUE, remove_na = remove_na, verbose = verbose)
      f2_blocks = (f2_blocks_ap - min(f2_blocks_ap, na.rm=T)) * diff(range(f2_blocks, na.rm = TRUE))/diff(range(f2_blocks_ap, na.rm = TRUE)) + min(f2_blocks, na.rm=T)
    }
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
    keep = apply(f2_blocks, 3, function(x) sum(is.na(x)) == 0)
    f2_blocks = f2_blocks[,,keep]
  }
  namedList(f2_blocks, indivs, pairs, poplist)
}


f2_from_precomp_nonsquare = function(dir, pops1, pops2, afprod = FALSE, remove_na = TRUE, verbose = TRUE) {

  if(afprod) block_lengths = readRDS(paste0(dir, '/block_lengths_a.rds'))
  else block_lengths = readRDS(paste0(dir, '/block_lengths.rds'))
  arr = array(NA, c(length(pops1), length(pops2), length(block_lengths)),
              list(pops1, pops2, paste0('l', block_lengths)))

  pops = expand_grid(pops1, pops2) %>%
    mutate(p1 = pmin(pops1, pops2), p2 = pmax(pops1, pops2))
  # putting the if statment inside the for loop makes it slower
  # check if there is an efficient tidyverse solution for assigning values to external objects (arr columns)
  if(afprod) {
    for(i in 1:nrow(pops)) {
      mat = readRDS(paste0(dir, '/', pops$p1[i], '/', pops$p2[i], '_ap.rds'))
      arr[pops$pops1[i], pops$pops2[i], ] = -2*mat[,1]
    }
  } else {
    for(i in 1:nrow(pops)) {
      mat = readRDS(paste0(dir, '/', pops$p1[i], '/', pops$p2[i], '_f2.rds'))
      arr[pops$pops1[i], pops$pops2[i], ] = mat[,1]
    }
  }
  if(remove_na) {
    keep = apply(arr, 3, function(x) sum(is.na(x)) == 0)
    arr = arr[,,keep, drop = FALSE]
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

is_plink_prefix = function(input) {
  if(!is.character(input) || length(input) > 1) return(FALSE)
  filesexist = all(file.exists(paste0(input, c('.bed', '.bim', '.fam'))))
  filesexist
}

is_geno_prefix = function(input) {
  is_packedancestrymap_prefix(input) || is_plink_prefix(input)
}

#' Extract samples from packedancestrymap files
#'
#' This function extracts samples from packedancestrymap files and saves them as PLINK files, using the package 'genio'
#' @export
#' @param inpref Prefix of the packedancestrymap input files
#' @param outpref Prefix of the PLINK output files
#' @param inds Individuals which should be extracted
#' @param pops Populations which should be extracted. Can not be provided together with 'inds'
extract_samples = function(inpref, outpref, inds = NULL, pops = NULL) {
  # extracts samples from geno file and writes new, smaller PLINK file using genio

  stopifnot(is.null(inds) || is.null(pops))
  if(!is.null(pops)) {
    inds = read_table2(paste0(inpref, '.ind'), col_names = F, col_types = 'cccccc') %>%
      filter(X3 %in% pops) %$% X1
  }
  dat = read_packedancestrymap(inpref, inds)
  genodat = dat[[1]]
  inddat = dat[[2]]
  snpdat = dat[[3]]

  nam = c('id', 'chr', 'posg', 'pos', 'ref', 'alt')
  bim = snpdat %>% set_colnames(nam)
  fam = inddat %>%
    transmute(fam = X3, id = X1, pat = 0, mat = 0, sex = X2, pheno = -9)
  genio::write_plink(normalizePath(outpref, mustWork = F), genodat, bim = bim, fam = fam)
}

