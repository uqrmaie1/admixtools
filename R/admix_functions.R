

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# general functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# number of all binary trees with a given number of leaf nodes
numtrees = function(n) factorial(2*n-2)/(2^(n-1)*factorial(n-1))

# number of all binary tree topologies
numtreestop = function(n) factorial(2*n)/factorial(n+1)/factorial(n)

# number of possible DAGs
numdags = function(n) {
  if(n <= 1) return(1)
  sum(sapply(1:n, function(k) (-1)^(k+1) * choose(n, k) * 2^(k*(n-k)) * numdags(n-k)))
}

namedList = function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}

mean_impute = function(mat, by=2) {
  out = apply(mat, by, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
  if(by == 1) out = t(out)
  out
}

multistart = function (parmat, fn, gr = NULL, lower = -Inf, upper = Inf, method = NULL, 
          hessian = FALSE, control = list(), verbose=TRUE, ...) 
{
  # same as from optimx library, but can be quiet
  nset <- nrow(parmat)
  npar <- ncol(parmat)
  if (nset < 1) 
    stop("multistart: no starting parameters!")
  ans.ret <- matrix(NA, nrow = nset, ncol = npar + 4)
  ans.ret <- data.frame(ans.ret)
  pstring <- colnames(parmat)
  if (is.null(pstring)) {
    pstring <- NULL
    for (j in 1:npar) {
      pstring[[j]] <- paste("p", j, sep = "")
    }
  }
  cnames <- c(pstring, "value", "fevals", "gevals", "convergence")
  colnames(ans.ret) <- cnames
  row.names(ans.ret) <- 1:nset
  for (imeth in 1:nset) {
    start <- parmat[imeth, ]
    ans <- optimr(par = start, fn = fn, gr = gr, lower = lower, 
                  upper = upper, method = method, hessian = hessian, 
                  control = control, ...)
    addvec <- c(ans$par, ans$value, ans$counts[1], ans$counts[2], 
                ans$convergence)
    ans.ret[imeth, ] <- addvec
    if(verbose) cat(paste0('\r', imeth, ' out of ', nset))
  }
  if(verbose) cat('\n')
  ans.ret
}


packedancestrymap_to_aftable = function(pref, pops=NULL, inds=NULL, blocksize=10000, na.action='none') {
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
  
  cat(blue(paste0(basename(pref), '.geno has ', bold(nind), ' samples and ', bold(nsnp), ' SNPs.\n')))
  cat(green(paste0('Calculating allele frequencies from ', bold(length(popind)), ' samples in ', bold(numpop), ' populations...\n')))
  cat(blue(paste0('Expected size of allele frequency data: ', bold(round((nsnp*numpop*8+nsnp*112)/1e6)), ' MB\n')))
  # 8, 112: estimated scaling factors for AF columns and annotation columns
  
  rlen = file.info(fl)$size/(nsnp+1)
  conn = file(fl, 'rb')
  invisible(readBin(conn, 'raw', n = rlen))
  afmatrix = matrix(NA, nsnp, numpop)
  colnames(afmatrix) = pops
  popind3 = c(outer(popind2, ((1:blocksize)-1)*rlen*4, `+`))
  popindmat = matrix(sapply(1:numpop, function(i) indfile$X3[popind2] == pops[i]), ncol=numpop)
  cnt=1
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
    afmatrix[cnt:(cnt+blocksize-1),] = popfreqs
    cat(paste0('\r', bold((cnt-1)/1e3), 'k SNPs read...'))
    cnt = cnt+blocksize
  }
  cat(paste0('\r', bold(cnt-1, ' SNPs'), ' read\n'))
  sumna = sum(is.na(afmatrix))
  cat(blue(paste0(bold(sumna), ' allele frequencies are missing (', bold(round(sumna/numpop)), ' per population)\n')))

  if(na.action == 'impute') {
    afmatrix = mean_impute(afmatrix, by=1)
    isnan = is.nan(afmatrix)
    afmatrix[isnan] = NA
    allna = apply(afmatrix, 1, function(x) all(is.na(x)))
    afmatrix = afmatrix[!allna,]
    snpfile = snpfile[!allna,]
    cat(green(paste0(bold(sumna-sum(isnan)), ' allele frequencies were imputed, ', bold(sum(allna)), ' SNPs were removed\n')))
  }
  if(na.action == 'remove') {
    anyna = apply(afmatrix, 1, function(x) any(is.na(x)))
    afmatrix = afmatrix[!anyna,]
    snpfile = snpfile[!anyna,]
    cat(green(paste0(bold(sum(anyna)), ' SNPs were removed\n')))
  }
  snpfile %>% bind_cols(as_tibble(afmatrix))
}


plink_to_aftable = function(pref, pops, outdir='./freqs/', plink='/Users/robert/Downloads/plink_mac/plink') {
  # pref is plink prefix
  # pops is vector of populations, with length equal to nrow(fam); should not have whitespace
  
  system(paste0('mkdir ', outdir))
  bim = read_table2(paste0(pref, '.bim'), col_names=F)
  bim %>% transmute(X2, X5) %>% write_tsv(paste0(outdir, 'refalleles.txt'), col_names=F)
  fam = read_table2(paste0(pref, '.fam'), col_names=F) %>% bind_cols(p = pops)
  freqs = counts = bim %>% set_colnames(c('CHR', 'SNP', 'cm', 'POS', 'A1', 'A2'))
  for(pop in unique(na.omit(pops))) {
    print(pop)
    fam %>% filter(p == pop) %>% write_tsv(paste0(outdir, 'samples.txt'), col_names=F)
    cmd = paste0(plink, " --bfile ", pref, " --keep ", outdir, "samples.txt --freq --keep-allele-order --out ", outdir, "freqs --reference-allele ", outdir, "refalleles.txt --allow-no-sex")
    system(cmd, ignore.stdout=T)
    dat = read_table2(paste0(outdir, 'freqs.frq'), col_types = cols())
    freqs %<>% bind_cols(!!pop:=dat$MAF)
    counts %<>% bind_cols(!!pop:=dat$NCHROBS)
  }
  list(freqs, counts)
}


cor.test.jack = function(x, y, blocks=1000, blockids=NA, covariance=F) {
  # assumes x and y are ordered (probably genomic position)
  # either blocks or blockids has to be specified; blocks is ignored if blockids is specified
  keep = is.finite(x) & is.finite(y)
  x = x[keep]
  y = y[keep]
  if(!is.na(blockids[1])) {
    blockids = as.numeric(as.factor(blockids[keep]))
    blocks = max(blockids)
  } else {
    blockids = floor(seq(1, blocks+1-1e-3, len=length(x)))
  }
  mat = cbind(x, y, blockids)
  #part_ests = sapply(1:blocks, function(i) cor(mat[mat[,3] != i, 1], mat[mat[,3] != i, 2]))
  lengths = rle(mat[,3])$lengths
  parts_sm1 = sapply(split(mat[,1], mat[,3]), sum)
  parts_sm2 = sapply(split(mat[,2], mat[,3]), sum)
  parts_mn1 = parts_sm1/lengths
  parts_mn2 = parts_sm2/lengths
  parts_dp = sapply(split((mat[,1]-parts_mn1[mat[,3]])*(mat[,2]-parts_mn2[mat[,3]]), mat[,3]), sum)
  parts_dp1 = sapply(split((mat[,1]-parts_mn1[mat[,3]])^2, mat[,3]), sum)
  parts_dp2 = sapply(split((mat[,2]-parts_mn2[mat[,3]])^2, mat[,3]), sum)
  tot_dp = sum((mat[,1]-mean(mat[,1]))*(mat[,2]-mean(mat[,2])))
  tot_dp1 = sum((mat[,1]-mean(mat[,1]))^2)
  tot_dp2 = sum((mat[,2]-mean(mat[,2]))^2)
  
  lo_dp = rep(tot_dp, blocks) - parts_dp
  lo_dp1 = rep(tot_dp1, blocks) - parts_dp1
  lo_dp2 = rep(tot_dp2, blocks) - parts_dp2
  
  if(covariance) lo_corr = lo_dp/(nrow(mat)-lengths)
  else lo_corr = lo_dp/sqrt(lo_dp1*lo_dp2)
  jack_est = mean(lo_corr)
  jack_se = sqrt((blocks-1)/blocks * sum((jack_est - lo_corr)^2))
  norm_est = cor(x, y)
  norm_se = cor.test.plus(x, y)$se
  list(estimate=jack_est, se=jack_se, p.value=2*pnorm(-abs(jack_est/jack_se)), estimate_normal=norm_est, se_normal=norm_se, p.value_normal=2*pnorm(-abs(norm_est/norm_se)))
}

get_cov_blocks = function(mat, blockids) {
  # assumes x and y are ordered (probably genomic position)
  # returns k * ncol * ncol array
  keep = apply(mat, 1, function(x) all(is.finite(x)))
  mat = mat[keep,]
  
  blockids = as.numeric(as.factor(blockids[keep]))
  numblocks = length(unique(blockids))
  
  #mat = cbind(x, y, blockids)
  splmat = matrix(blockids, nrow(mat), ncol(mat))
  makemat = function(x) matrix(x, ncol=ncol(mat))
  sumfun = function(x) colSums(makemat(x))
  lengths = rle(blockids)$lengths
  parts_mn = t(sapply(split(mat, splmat), sumfun))/lengths
  
  matcentered = mat - parts_mn[rep(1:length(lengths), lengths),]
  parts_dp = lapply(split(matcentered, splmat), function(x) crossprod(makemat(x)))
  tot_dp = crossprod(scale(mat, scale=F))
  lo_cov = lapply(1:length(parts_dp), function(i) (tot_dp - parts_dp[[i]])/(nrow(mat) - lengths[i]))
  array(unlist(lo_cov), c(ncol(mat), ncol(mat), numblocks))
  
}

get_cov_blocks2 = function(m1, m2, blockids) {
  # assumes x and y are ordered (probably genomic position)
  # returns k * ncol * ncol array
  keep = apply(cbind(m1, m2), 1, function(x) all(is.finite(x)))
  m1 = m1[keep,]
  m2 = m2[keep,]
  
  blockids = as.numeric(as.factor(blockids[keep]))
  numblocks = length(unique(blockids))
  
  splmat1 = matrix(blockids, nrow(m1), ncol(m1))
  splmat2 = matrix(blockids, nrow(m2), ncol(m2))
  makemat = function(x, ncol) matrix(x, ncol=ncol)
  sumfun = function(x, ncol) colSums(makemat(x, ncol))
  lengths = rle(blockids)$lengths
  parts_mn1 = t(sapply(split(m1, splmat1), sumfun, ncol(m1)))/lengths
  parts_mn2 = t(sapply(split(m2, splmat2), sumfun, ncol(m2)))/lengths
  
  matcentered1 = m1 - parts_mn1[rep(1:length(lengths), lengths),]
  matcentered2 = m2 - parts_mn2[rep(1:length(lengths), lengths),]
  parts_dp1 = lapply(split(matcentered1, splmat1), function(x) crossprod(makemat(x, ncol(m1))))
  parts_dp2 = lapply(split(matcentered2, splmat2), function(x) crossprod(makemat(x, ncol(m2))))
  spl1 = lapply(split(matcentered1, splmat1), makemat, ncol(m1))
  spl2 = lapply(split(matcentered2, splmat2), makemat, ncol(m2))
  parts_dp = lapply(1:length(spl1), function(i) t(spl1[[i]]) %*% spl2[[i]])
  tot_dp = t(scale(m1, scale=F)) %*% scale(m2, scale=F)
  lo_cov = lapply(1:length(parts_dp), function(i) (tot_dp - parts_dp[[i]])/(nrow(m1) - lengths[i]))
  array(unlist(lo_cov), c(ncol(m1), ncol(m2), numblocks))
  
}

get_dp_blocks2 = function(m1, m2, blockids) {
  # same as get_cov_blocks2, but without centering
  # assumes x and y are ordered (probably genomic position)
  # returns k * ncol(m1) * ncol(m2) array
  keep = apply(cbind(m1, m2), 1, function(x) all(is.finite(x)))
  m1 = m1[keep,]
  m2 = m2[keep,]
  
  blockids = as.numeric(as.factor(blockids[keep]))
  numblocks = length(unique(blockids))
  
  splmat1 = matrix(blockids, nrow(m1), ncol(m1))
  splmat2 = matrix(blockids, nrow(m2), ncol(m2))
  makemat = function(x, ncol) matrix(x, ncol=ncol)
  sumfun = function(x, ncol) colSums(makemat(x, ncol))
  lengths = rle(blockids)$lengths
  parts_mn1 = t(sapply(split(m1, splmat1), sumfun, ncol(m1)))/lengths
  parts_mn2 = t(sapply(split(m2, splmat2), sumfun, ncol(m2)))/lengths
  
  parts_dp1 = lapply(split(m1, splmat1), function(x) crossprod(makemat(x, ncol(m1))))
  parts_dp2 = lapply(split(m2, splmat2), function(x) crossprod(makemat(x, ncol(m2))))
  spl1 = lapply(split(m1, splmat1), makemat, ncol(m1))
  spl2 = lapply(split(m2, splmat2), makemat, ncol(m2))
  parts_dp = lapply(1:length(spl1), function(i) t(spl1[[i]]) %*% spl2[[i]])
  tot_dp = t(m1) %*% m2
  lo_cov = lapply(1:length(parts_dp), function(i) (tot_dp - parts_dp[[i]])/(nrow(m1) - lengths[i]))
  array(unlist(lo_cov), c(ncol(m1), ncol(m2), numblocks))
  
}


get_dp_blocks_colsq = function(mat, blockids) {
  # get block jackknife dot products of euclidean distance of each column (=f2 if input is afdiff mat)
  # assumes x and y are ordered (probably genomic position)
  # returns k * ncol(m) matrix
  keep = apply(mat, 1, function(x) all(is.finite(x)))
  mat = mat[keep,]
  
  blockids = as.numeric(as.factor(blockids[keep]))
  numblocks = length(unique(blockids))
  lengths = rle(blockids)$lengths
  
  splmat = matrix(blockids, nrow(mat), ncol(mat))
  makemat = function(x) matrix(x, ncol=ncol(mat))
  
  parts_dp = sapply(split(mat, splmat), function(x) colSums(makemat(x)^2))
  tot_dp = colSums(mat^2)
  lo_dp = rep(tot_dp, numblocks) - parts_dp
  t(t(lo_dp)/(nrow(mat) - lengths))
}

# setblocks = function(dat, dist=0.05) {
#   # dat must have SNP, CHR, POS
#   # must be arranged by CHR, POS and maybe grouping factor
#   # still not exactly equal. probably something at chromosome boundaries
#   dat %<>% mutate(dff = pmax(0, c(0, diff(POS))), cumsum = 0, newblock=FALSE)
#   # find way to vectorize this
#   for(i in 1:nrow(dat)) {
#     if(i %% 10000 == 0) print(i)
#     if(i < nrow(dat) && dat$newblock[i]) next
#     if(i == 1 || i < nrow(dat) && dat$CHR[i] != dat$CHR[i-1]) dat$cumsum[i] = 0
#     else dat$cumsum[i] = dat$cumsum[i-1] + dat$POS[i] - dat$POS[i-1]
#     if(i != nrow(dat) && (dat$cumsum[i] > dist ||  dat$CHR[i+1] != dat$CHR[i])) {
#       dat$newblock[i+1] = TRUE
#       if(dat$CHR[i+1] != dat$CHR[i]) {
#         dat$cumsum[i+1] = dat$POS[i+1]
#       } else {
#         dat$cumsum[i+1] = dat$POS[i+1] - dat$POS[i]
#       }
#     }
#   }
#   dat %>% mutate(block = cumsum(newblock)+1)
# }

setblocks4 = function(dat, dist=0.05, distcol='cm') {
  # requires 'CHR' and either 'POS' or 'cm'
  # dat needs to be ordered first by 'CHR', then by 'POS' or 'cm'
  # starts new block at SNP after the first SNP which is not within dist of the last block start
  
  sb = function(cumpos, CHR, POS) {o = cumpos; cumpos[3]=POS; cumpos[1] = o[2]; cumpos[2] = pmax(0, POS-o[3]+o[2]); if(o[2] %% dist < o[1] %% dist && o[2] > dist) cumpos[2] = POS-o[3]; cumpos}
  
  newdat = do.call(rbind, (accumulate2(.x=dat$CHR, .y=dat[[distcol]], .f=sb, .init=c(0, 0, dat[[distcol]][1])))) %>% as_tibble()
  
  dat %<>% bind_cols(newdat %>% slice(-1))
  dat %>% mutate(newblock = V2 > lead(V2, default=0) & CHR == lead(CHR, default=0) | CHR > lag(CHR, default=0), block=cumsum(newblock))
  
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# qpDstat
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


qpDstat = function(pops, bin='/Users/robert/Downloads/AdmixTools-master/bin/qpDstat', pref='/Users/robert/Downloads/v37.2_HO_subs3', outdir='.', blgsize='100', env='') {
  # pops is matrix with 4 columns
  if(is.null(dim(pops))) pops = t(pops)
  
  parfile = paste0('genotypename: ', pref, '.geno\n',
                   'snpname: ', pref, '.snp\n',
                   'indivname: ', pref, '.ind\n',
                   'popfilename: ',outdir,'/popfile\n',
                   'f4mode: YES\n',
                   'printsd:  YES')
  
  write(parfile, paste0(outdir, '/parfile'))
  print(pops)
  write(pops, paste0(outdir, '/popfile'), ncolumns = ncol(pops))
  
  cmd = paste0(env,' ', bin, ' -p ', outdir, '/parfile > ', outdir, '/out')
  
  system(cmd)
  system('cat ./out')
  
}

dstat = function(afs, blocks, pops) {
  est = f4(afs[,pops[1]], afs[,pops[2]], afs[,pops[3]], afs[,pops[4]])
  se = cor.test.jack(afs[,pops[1]] - afs[,pops[2]], afs[,pops[3]] - afs[,pops[4]], blockids = blocks, covariance=T)$se
  z = est/se
  p = ztop(z)
  c(est, se, z, p)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# F4-ratio test
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




qpf4ratio = function(A, O, X, B, C, bin='/Users/robert/Downloads/AdmixTools-master/bin/qpF4ratio', pref='/Users/robert/Downloads/AdmixTools-master/convertf/example', outdir='.', env='', blgsize='100') {
  # p1 has to be in genotype file
  # p2,p3,p4 can be vector of population labels
  
  parfile = paste0('genotypename: ', pref, '.geno\n',
                   'snpname: ', pref, '.snp\n',
                   'indivname: ', pref, '.ind\n',
                   'popfilename: ',outdir,'/f4list\n',
                   'blgsize: ', blgsize)
  
  write(parfile, paste0(outdir, '/parfile'))
  write(paste(X, C, O, A, B, C, O, A), paste0(outdir, '/f4list'))
  # 'Ami Eskimo_Naukan : Mbuti Even :: French Eskimo_Naukan : Mbuti Even'
  
  cmd = paste0(env,' ', bin, ' -p ', outdir, '/parfile > ', outdir, '/out')
  
  system(cmd)
  vals = as.numeric(unlist(quietly(read_table2)('/Users/robert/Dropbox/postdoc2/projects/admixprograms/out', skip=21, n_max = 1, col_names = F))[11:13])
  names(vals) = c('alpha', 'se', 'z')
  vals
}



f4 = function(p1, p2, p3, p4) {
  cov(p1-p2, p3-p4)
}

f4ratio = function(afs, A, O, X, B, C) {
  f41 = f4(afs[,X], afs[,C], afs[,O], afs[,A])
  f42 = f4(afs[,B], afs[,C], afs[,O], afs[,A])
  f41/f42
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# qpAdm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

qpAdm = function(leftpops, rightpops, bin='/Users/robert/Downloads/AdmixTools-master/bin/qpAdm', pref='/Users/robert/Downloads/v37.2_HO_subs3', outdir='.', env='') {
  
  parfile = paste0('genotypename: ', pref, '.geno\n',
                   'snpname: ', pref, '.snp\n',
                   'indivname: ', pref, '.ind\n',
                   'popleft: ',outdir,'/leftlist\n',
                   'popright: ',outdir,'/rightlist\n',
                   'details: YES\n',
                   'hashcheck: NO')
  
  write(parfile, paste0(outdir, '/parfile'))
  write(leftpops, paste0(outdir, '/leftlist'))
  write(rightpops, paste0(outdir, '/rightlist'))
  
  cmd = paste0(env,' ', bin, ' -p ', outdir, '/parfile > ', outdir, '/out')
  
  system(cmd)
  #system('cat ./out')
  
}

getf4x = function(afs, leftpops, rightpops) {
  # needs to be fixed; compare to jackknife function
  pops = union(leftpops, rightpops)
  tot = length(pops)
  f2mat = matrix(NA, tot, tot)
  for(i in 1:tot) {
    for(j in (i+1):tot) {
      if(j <= tot & j > i) {
        f2mat[i,j] = f2mat[j,i] = mean((afs[,i] - afs[,j])^2)
      }
    }
  }
  colnames(f2mat) = rownames(f2mat) = pops
  #print(f2mat)
  
  u0 = match(leftpops[1], pops)
  v0 = match(rightpops[1], pops)
  f4mat = matrix(NA, length(leftpops)-1, length(rightpops)-1)
  for(i in 2:length(leftpops)) {
    lpos = match(leftpops[i], pops)
    for(j in 2:length(rightpops)) {
      rpos = match(rightpops[j], pops)
      # Eq. 24 b; f4(p1, p2; p3, p4) = (f2(p1, p4) + f2(p2, p3) - f2(p1, p3) - f2(p2, p4))
      f4mat[i-1, j-1] = (f2mat[u0, rpos] + f2mat[lpos, v0] - f2mat[u0, v0] - f2mat[lpos, rpos])/2
    }
  }
  rownames(f4mat) = leftpops[-1]
  colnames(f4mat) = rightpops[-1]
  f4mat
}


get_weights_covariance = function(f4x, qinv, likfun) {
  
  rnk = dim(f4x)[1]-1
  numblocks = dim(f4x)[3]
  X = apply(f4x, 1:2, mean)
  
  wmat = matrix(NA, numblocks, dim(f4x)[1])
  for(i in 1:numblocks) {
    f4svd = svd(f4x[,,i])
    A = f4svd$u[, 1:rnk, drop=F]
    # B = t(f4svd$v[, 1:rnk, drop=F])
    # d = f4svd$d[1:rnk, drop=F]
    # B = diag(d) %*% B
    # nr = nrow(A)
    # nc = ncol(B)
    # xx = optim(par = c(A, B), likfun, args=list(d, X, qinv, nr, nc, rnk), method='L-BFGS-B', control=list(maxit=1e4))
    # Aopt = matrix(xx$par[1:(nr*rnk)], nr, rnk)
    # Bopt = matrix(xx$par[-(1:(nr*rnk))], rnk, nc)
    # w = MASS::Null(Aopt)[,1]
    w = MASS::Null(A)[,1]
    wmat[i,] = w/sum(w)
  }
  jackmeans = colMeans(wmat)
  mnc = jackmeans - t(wmat)
  (numblocks-1)/numblocks * (mnc %*% t(mnc))
  #sqrt(diag((numblocks-1)/numblocks * (mnc %*% t(mnc))))
}

getQ = function(f4x) {
  # assuming Q is (nl*nr) * (nl*nr) matrix of f4 SEs
  numblocks = length(unique(blockids))
  f4xflat = matrix(f4x, (dim(f4x)[1])*(dim(f4x)[2]), numblocks)
  bdiff = t(rowMeans(f4xflat) - f4xflat)
  (numblocks-1)/numblocks * crossprod(bdiff)
}

likfun = function(mats, args) {
  
  d = args[[1]]
  X = args[[2]]
  qinv = args[[3]]
  nr = args[[4]]
  nc = args[[5]]
  rnk = args[[6]]
  
  A = matrix(mats[1:(nr*rnk)], nr, rnk)
  B = matrix(mats[-(1:(nr*rnk))], rnk, nc)
  A = t(t(A)/colSums(A^2))
  B = t(t(B)/rowSums(B^2))
  #E = X - (A %*% diag(d, nrow=length(d)) %*% B)
  E = X - (A %*% B)
  lik = -sum(qinv * (E %x% t(E)))/2
  if(!is.finite(lik)) {
    print('Inf!')
    return(.Machine$double.xmin)
  }
  lik
}

scalefun = function(mat) {
  # continue here: look at 'printstrmat' in 'f4rank.c'
  y = sqrt(colSums(mat^2))
  y = y/sqrt(nrow(mat))
  1/y  
}

adm = function(leftpops, rightpops, afs=NULL, blockids=NA, f2blockarr=NULL, opt=TRUE) {
  stopifnot(!is.null(afs) | !is.null(f2blockarr))
  
  if(is.na(blockids[1])) blockids = afs %>% setblocks4(dist=0.05) %$% block
  numblocks = length(unique(blockids))
  
  if(!is.null(f2blockarr)) {
    f4x = get_f4_from_f2blockarr(f2blockarr, leftpops, rightpops)
  } else {
    lmat = as.matrix(afs[,leftpops])
    rmat = as.matrix(afs[,rightpops])
    ldiff = lmat[,-1] - lmat[,1]
    rdiff = rmat[,-1] - rmat[,1]
    f4x = get_dp_blocks2(ldiff, rdiff, blockids)
  }
  
  X = apply(f4x, 1:2, mean)
  Q = getQ(f4x)
  qinv = pinv(Q)
  
  rnk = dim(f4x)[1]-1
  f4svd = svd(X)
  A = f4svd$u[, 1:rnk, drop=FALSE]
  B = t(f4svd$v[, 1:rnk, drop=FALSE])
  d = f4svd$d[1:rnk, drop=FALSE]
  B = diag(d) %*% B
  
  if(opt) {
    nr = nrow(A)
    nc = ncol(B)
    xx = optim(par = c(A, B), likfun, args=list(d, X, qinv, nr, nc, rnk), method='L-BFGS-B', control=list(maxit=1e4, fnscale=1))
    A = matrix(xx$par[1:(nr*rnk)], nr, rnk)
    B = matrix(xx$par[-(1:(nr*rnk))], rnk, nc)
  }
  
  w = MASS::Null(A)[,1]
  w = w/sum(w)
  
  covmat = get_weights_covariance(f4x, qinv, likfun)
  list(weights = w, se = sqrt(diag(covmat)))
}

get_f2blockarr = function(afmat, blockids) {
  nc = ncol(afmat)
  stopifnot(nc <= 50)
  diffmat = afmat[,rep(1:nc,nc)] - afmat[,rep(1:nc,each=nc)]
  dp_blocks = get_dp_blocks_colsq(diffmat, blockids)
  dp_blocks = array(dp_blocks, dim = c(nc, nc, length(unique(blockids))))
  dimnames(dp_blocks) = list(colnames(afmat), colnames(afmat), NULL)
  dp_blocks
}

get_f4_from_f2blockarr = function(f2blockarr, leftpops, rightpops) {
  leftind = match(leftpops, dimnames(f2blockarr)[[1]])
  rightind = match(rightpops, dimnames(f2blockarr)[[2]])
  nl = length(leftpops)
  nr = length(rightpops)
  # Eq. 24 b; f4(p1, p2; p3, p4) = (f2(p1, p4) + f2(p2, p3) - f2(p1, p3) - f2(p2, p4))
  f4arr =
    f2blockarr[rep(leftind[ 1], nl-1), rightind[-1],,drop=F] +
    f2blockarr[leftind[-1], rep(rightind[ 1], nr-1),,drop=F] -
    f2blockarr[rep(leftind[ 1], nl-1), rep(rightind[ 1], nr-1),,drop=F] -
    f2blockarr[leftind[-1], rightind[-1],,drop=F]
  f4arr/2
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# qpGraph
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

qpGraph2 = function(edges, outpop='NULL', bin='/Users/robert/Downloads/AdmixTools-master/bin/qpGraph', pref='/Users/robert/Downloads/v37.2_HO_subs3', outdir='.', printonly=FALSE, lambdascale=-1, lsqmode='NO', diag=0.0001, hires='NO', forcezmode='NO', allsnps='NO', bigiter=100, env='') {
  # wrapper around AdmixTools qpGraph
  # makes parfile and graphfile
  
  parfile = paste0('indivname:       ', pref, '.ind\n',
                   'snpname:         ', pref, '.snp\n',   
                   'genotypename:    ', pref, '.geno\n',
                   'outpop:         ', outpop, '\n',
                   'blgsize: 0.05\n',
                   'details: YES\n',
                   'fstdetails: YES\n',
                   'diag: ', diag, '\n',
                   'lsqmode: ', lsqmode, '\n',
                   'hires: ', hires, '\n',
                   'forcezmode: ', forcezmode, '\n',
                   'allsnps: ', allsnps, '\n',
                   'lambdascale: ', lambdascale, '\n',
                   'bigiter: ', bigiter, '\n')
  
  edg = as_tibble(edges) %>% set_colnames(c('V1', 'V2')) %>% group_by(V2) %>% mutate(type = ifelse(n()==1, 'edge', 'admix')) %>% ungroup
  e1 = edg %>% filter(type == 'edge') %$% V1
  e2 = edg %>% filter(type == 'edge') %$% V2
  a1 = edg %>% filter(type == 'admix') %$% V1
  a2 = edg %>% filter(type == 'admix') %$% V2
  leaves = setdiff(edg$V2, edg$V1)
  admix = tibble()
  for(m in unique(a2)) {
    admix %<>% bind_rows(tibble(v1='admix', v2=m, v3=edg %>% filter(V2 == m) %$% V1[1], v4=edg %>% filter(V2 == m) %$% V1[2]))
  }
  
  pops = union(edg[[1]], edg[[2]])
  simfile = tibble(v1 = c('root'), v2 = c('R'), v3='', v4='') %>%
    bind_rows(tibble(v1 = 'label', v2=leaves, v3=leaves, v4='')) %>%
    bind_rows(tibble(v1='edge', v2=paste0('e', 1:length(e1)), v3=e1, v4=e2)) %>%
    bind_rows(admix)

  pf = paste0(outdir, '/parfile')
  gf = paste0(outdir, '/graphfile')
  of = paste0(outdir, '/out')
  
  write(parfile, pf)
  simfile %>% write_tsv(gf, col_names=F)

  qpGraph(bin=bin, parfile=pf, graphfile=gf, outfile=of, printonly=printonly, env=env)
}

qpGraph = function(bin='./qpGraph', parfile='./parfile', graphfile='./graphfile', outfile='./out', printonly=FALSE, env='') {
  # wrapper around AdmixTools qpGraph
  # input is locations of parfile and graphfile
  # output is parsed output
  
  cmd = paste0(env,' ', bin, ' -p ', parfile, ' -g ', graphfile, ' > ', outfile)
  if(printonly) {
    print(cmd)
  } else{
    system(cmd)
    return(parse_qpGraph_output(outfile))
  }
}

parse_qpGraph_output = function(outfile) {
  # reads qpGraph output file
  # returns list of three objects:
  # 'edgeweights': data.frame of branch lengths and admixture weights
  # 'score': best fit score
  # 'f2fit': data.frame of estimated and fitted f2 values
  
  dat = read_table(outfile, col_names=F, col_types = cols(), guess_max = 1e6)
  
  edgeweights = dat %>% filter(grepl('^ledge|^redge|^admix', X1)) %>% separate(X1, c('type', 'name', 'from', 'to', 'w', 'w2'), sep=' +', convert = T, extra='drop', fill='right') %>%
    bind_rows(filter(., type=='admix') %>% mutate(type='aedge', to=name, name='', w2=NA)) %>% 
    bind_rows(filter(., type=='admix') %>% mutate(type='aedge', from=to, to=name, name='', w=w2, w2=NA)) %>%
    filter(!type == 'admix') %>% mutate(type = ifelse(type=='aedge', 'admix', 'edge')) %>% select(-w2, -name)
  
  score = dat %>% filter(grepl('^final score', X1)) %>% separate(X1, c('a', 'b', 'score'), sep=' +', convert = T, extra='drop', fill='right') %$% score
  
  f2fit = dat %>% filter(grepl(' f2: ', X1)) %>% separate(X1, c('pop1', 'pop2', 'fst','f2fit','f2est','diff','se','z'), sep=' +', convert = T) %>% select(-fst)
  
  namedList(edgeweights, score, f2fit)
}

parse_qpGraph_graphfile = function(graphfile) {
  # reads graph in qpGraph format
  # returns edge matrix (adjacency list)
  
  read_lines(graphfile) %>% tibble %>% filter(grepl('^edge|admix', .)) %>% separate('.', c('type', 'name', 'from', 'to'), sep = '\\s+') %>%
    bind_rows(filter(., type=='admix') %>% mutate(type='edge', to=name, name='')) %>%
    bind_rows(filter(., type=='admix') %>% mutate(type='edge', from=to, to=name, name='')) %>%
    filter(type != 'admix') %>% select(from, to) %>% as.matrix
}

parse_qpGraph_parfile = function(parfile) {
  # reads qpGraph parfile
  # returns named list of parameters
  # all genotype files have to have same prefix
  
  read_table2(parfile, comment = '#', col_names = c('par', 'value'), col_types = cols()) %>% mutate(par=gsub(':$', '', gsub('genotypename', 'pref', par))) %>% mutate(value=gsub('\\.geno$', '', gsub('S1', filter(., par=='S1')$value[1], gsub('DIR', filter(., par=='DIR')$value[1], value))), value=ifelse(value=='YES', TRUE, ifelse(value=='NO', FALSE, value))) %>% filter(!par %in% c('DIR', 'S1', 'indivname', 'snpname')) %>% t %>% as_tibble() %>% set_colnames(slice(., 1)) %>% slice(-1) %>% type_convert(col_type=cols()) %>% as.list
  
}



permute_leaves = function(edges, fix_outgroup=FALSE) {
  # shuffles the leaf nodes of a tree
  # assumes outgroup is at edges[1,2]
  
  edges = as.matrix(edges)
  leafindex = which(!edges[,2] %in% edges[,1])
  if(fix_outgroup) leafindex = setdiff(leafindex, 1)
  edges[leafindex,2] = sample(edges[leafindex,2])
  edges
}

permute_topology = function(edges, fix_outgroup=FALSE) {
  # creates a random new tree with the same leaf nodes and number of admixture events
  # assumes outgroup is at edges[1,2]
  
  edges = as.matrix(edges)
  outpop = edges[1,2]
  leafindex = which(!edges[,2] %in% edges[,1])
  leaves = setdiff(edges[,2], edges[,1])
  numadmix = sum(table(edges[,2]) == 2)
  if(fix_outgroup) leaves = setdiff(leaves, outpop)
  newick = random_newick(leaves)
  if(fix_outgroup) newick = paste0('(', outpop, ',', newick, ')')
  edges = newick_to_edges(newick)
  edges = insert_admix(edges, numadmix)
  edges
}



fst_qpgraph = function(p1, p2, c1, c2) {
  # c1, c2: counts of individuals
  # is the same as fst_reich
  num = (p1-p2)^2 - p1*(1-p1)/(2*c1-1) - p2*(1-p2)/(2*c2-1)
  denom = p1 + p2 - 2*p1*p2
  num/denom
}



bj_arr_lo_sum = function(arr, blockids) {
  # returns leave-one-block-out array sums
  # group over last dimension
  blockids = as.numeric(as.factor(blockids))
  numblocks = length(unique(blockids))
  tot = apply(arr, 1:2, sum)
  # might be faster with split
  parts = sapply(1:numblocks, function(i) apply(arr[,,blockids==i], 1:2, sum), simplify = 'array')
  rray(tot) - rray(parts)
}

bj_arr_lo_mean = function(arr, blockids) {
  # returns leave-one-block-out array means
  # group over last dimension
  
  #blockids = as.numeric(as.factor(blockids))
  lengths = rle(blockids)$lengths
  numblocks = length(unique(blockids))
  tot = apply(arr, 1:2, sum)
  parts = array(NA, c(dim(arr)[1], dim(arr)[2], numblocks))
  for(i in 1:numblocks) {
    if(i %% 10 == 0 || i == numblocks) cat(paste0('\r',  bold(i), ' out of ', bold(numblocks), ' SNP blocks processed...'))
    parts[,,i] = apply(arr[,,blockids==i], 1:2, sum)
  }
  cat('\n')
  sums = rray(tot) - rray(parts)
  out = sums / rray(length(blockids)-lengths, c(1,1,numblocks))
  if(!is.null(dimnames(arr)[[1]])) dimnames(out)[[1]] = dimnames(arr)[[1]]
  if(!is.null(dimnames(arr)[[2]])) dimnames(out)[[2]] = dimnames(arr)[[2]]
  out
}

bj_mat_stats = function(bj_lo_mat, blockids) {
  # input is matrix (one block per column)
  # output is list with vector of jackknife means and matrix of pairwise jackknife covariances
  # uses mean jackknife estimate instead of overall mean; probably makes very little difference
  # lengths normalization used to have -1, removed this because it resulted in negative values; where did this come from?
  lengths = rle(blockids)$lengths
  jest = rowMeans(bj_lo_mat)
  mnc = t(jest - bj_lo_mat) * sqrt((sum(lengths)/lengths-1)/length(lengths))
  jvar = (t(mnc) %*% mnc) 
  namedList(jest, jvar)
}

bj_arr_stats = function(bj_lo_arr, blockids) {
  # input is 3d array
  # output is list with jackknife means and jackknife variances
  # uses mean jackknife estimate instead of overall mean; probably makes very little difference
  lengths = rle(blockids)$lengths
  jest = apply(bj_lo_arr, 1:2, mean)
  xtau = rray(sum(lengths)/lengths-1, c(1,1,dim(bj_lo_arr)[3])) * (jest - bj_lo_arr)^2
  jvar = apply(xtau, 1:2, mean)
  namedList(jest, jvar)
}

arr3d_to_pairmat = function(arr, diag=TRUE) {
  # input is arr of dimension n x n x m
  # output is choose(n+1, 2) x m
  indx1 = which(lower.tri(arr[,,1], diag = diag))
  d1 = dim(arr)[1]
  d3 = dim(arr)[3]
  npair = choose(d1+diag, 2)
  indx = rep(indx1, d3) + rep((0:(d3-1))*d1^2, each=npair)
  matrix(as.array(arr)[indx], npair, d3)
}

get_yscal = function(fstest_blocks, fstnum_blocks, fstdenom_blocks, blockids) {
  # returns fst scaling factor 'yscal'
  fstest_blocks_2d = arr3d_to_pairmat(fstest_blocks)
  fst_blocks_2d = arr3d_to_pairmat(fstnum_blocks/(fstdenom_blocks+1e-10))
  f2_jest =  bj_mat_stats(fstest_blocks_2d, blockids)[[1]]
  fst_jest = apply(fstnum_blocks, 1:2, sum)/(apply(fstdenom_blocks, 1:2, sum)+1e-10)
  fst_jest = c(fst_jest[lower.tri(fst_jest, diag=T)])
  fst_jvar =  bj_mat_stats(fst_blocks_2d, blockids)[[2]]
  f2_jz = f2_jest/(sqrt(diag(fst_jvar)) + 1e-10)
  fst_jz = fst_jest/(sqrt(diag(fst_jvar)) + 1e-10)
  sum(c(fst_jz)*c(f2_jz)) / sum(c(f2_jz)*c(f2_jz))
}


plot_graph = function(edges, layout='tree') {
  # edges: 2 column data.frame or matrix of edges specifying an admixture graph (should be DAG)
  require(ggdag)
  edges = as.tibble(edges)
  names(edges)[1:2] = c('V1', 'V2')
  x = dag(paste(edges[[1]], edges[[2]], sep='->', collapse=' '))
  admixnodes = unique(edges[[2]][edges[[2]] %in% names(which(table(edges[[2]]) > 1))])

  dat = tidy_dagitty(x, layout=layout) %>% dplyr::mutate(xmean = (x+xend)/2, ymean=(y+yend)/2) %>% group_by(to) %>% dplyr::mutate(cnt = n()) %>% ungroup %>% dplyr::mutate(admix = (cnt > 1)*2+1) %>% ungroup %>% mutate(type = ifelse(name == 'R', 'root', ifelse(!name %in% edges[[1]], 'leaf', ifelse(name %in% admixnodes, 'admix', 'normal'))))

  if(!'label' %in% names(edges)) edges %<>% mutate(label='')
  dat %<>% dplyr::left_join(x=.$dat, y=edges, by=c('name'='V1','to'='V2'))
  plt = dat %>% ggplot(aes(x = x, y = y, xend = xend, yend = yend))  + geom_dag_edges(aes(edge_linetype=admix, label=label), edge_colour='grey') + geom_dag_label_repel(aes(fill=type, label=name), force=0) + theme_dag(legend.position='none') + scale_fill_manual(values=c('admix'=gg_color_hue(3)[1], 'leaf'=gg_color_hue(3)[2], 'normal'=gg_color_hue(3)[3], 'root'='#FFFFFF'))
  plt
}

plot_comparison = function(refout, myout) {
  # plots a comparison of qpGraph output and qpg output
  
  f2comp = refout$f2fit %>% select(-f2fit, -diff, -z) %>% mutate(prog = 'official') %>% bind_rows(myout$f2out %>% mutate(pop1=substr(pop1, 1, 3), pop2=substr(pop2, 1, 3), prog = 'my')) %>% gather(type, v, f2est, se) %>% spread(prog, v) %>% rename(from=pop1, to=pop2)
  
  myout[[1]] %>% left_join(refout$edgeweights %>% select(-type), by=c('from', 'to')) %>% dplyr::rename(official=w, my=weight) %>% bind_rows(f2comp) %>% ggplot(aes(official, my)) + geom_point() + facet_wrap(~ type, scales='free') + geom_abline() + xlab(paste0('official stats (score: ', refout$score,')')) + ylab(paste0('my stats (score: ', round(myout$score, 2),')')) + theme(panel.background = element_blank(), axis.line = element_line()) + geom_smooth(method='lm', se=F, formula=y~x-1)
  
}

edges_to_adjmat = function(edges) {
  # given a data.frame of directed edges, return nvertex*nvertex adjacency matrix
  # not tested / not needed?
  nam1 = unique(edges$V1)
  nam2 = unique(edges$V2)
  ind1 = match(edges$V1, nam1)
  ind2 = match(edges$V2, nam2)
  n1 = length(nam1)
  n2 = length(nam2)
  adjmat = matrix(0, n1, n2)
  rownames(adjmat) = nam1
  colnames(adjmat) = nam2
  adjmat[cbind(ind1, ind2)] = 1
  adjmat
}


graph_to_pwts = function(grph) {
  # input igraph object
  # output: (numpops-1)*numedges matrix 'pwts' which indicates all paths from pop to outpop; admixture edges are mapped onto parent edges and weighted
  # assumes first pop is outpop, and 'R' is root
  
  leaves = V(grph)$name[degree(grph, v = V(grph), mode = c('out')) == 0]
  pwts = matrix(0, length(E(grph)), length(leaves))
  colnames(pwts) = leaves
  rownames(pwts) = attr(E(grph), 'vnames')
  
  admixnodes = which(degree(grph, mode='in')==2)
  admixedges = unlist(incident_edges(grph, admixnodes, mode='in'))

  allpaths = all_simple_paths(grph, 'R', leaves, mode='out')
  pathcounts = table(names(sapply(allpaths, function(x) tail(x,1))))
  for(i in seq_len(length(allpaths))) {
    target = names(tail(allpaths[[i]],1))
    ln = length(allpaths[[i]])
    pth2 = allpaths[[i]][c(1, 1+rep(seq_len(ln-2), each=2), ln)]
    rowind = as.vector(E(grph)[get.edge.ids(grph, pth2)])
    pwts[rowind,target] = pwts[rowind,target] + 1/pathcounts[target]
  }
  
  if(!is.null(admixedges)) pwts = pwts[-admixedges,]
  pwts[,-1] - pwts[,1]
}

graph_to_weightind = function(grph) {
  # input igraph object
  # output: matrix with 8 columns: admixture node, admixture node index, edge id, edge id after removing actual admix edges (for pwts), left/right (0/1), all leave nodes connected to admix node, leaf indices (for pwts)
  # nrow = num leave nodes x (num edges left + num edges right), summed across all admixture nodes
  # doesn't include actual admixture edges (only normal edges between split nodes and admixture nodes)
  # assumes first pop is outpop, and 'R' is root
  
  leaves = which(degree(grph, v = V(grph), mode = c('out')) == 0)
  admixnodes = which(degree(grph, mode='in')==2)
  admixedges = unlist(incident_edges(grph, admixnodes, mode='in'))
  normedges = setdiff(1:length(E(grph)), admixedges)
  
  nam = c('admixnode', 'admixnodeindex', 'edgeid', 'edgeid2', 'lr', 'leaf', 'leafindex')
  out = matrix(0, 0, length(nam))
  colnames(out) = nam
  
  srt = topo_sort(grph, mode='out')
  admixnodes_sorted = srt[match(admixnodes, srt)]
  
  for(i in seq_len(length(admixnodes))) {
    thisleaves = intersect(leaves, subcomponent(grph, admixnodes[i], mode='out'))
    paths = all_simple_paths(grph, admixnodes[i], 'R', mode='in')
    nghb = neighbors(grph, admixnodes[i], mode='in')
    for(j in 1:length(paths)) {
      p = paths[[j]]
      ln = length(p)
      pth = p[c(1, 1+rep(seq_len(ln-2), each=2), ln)]
      e = setdiff(get.edge.ids(grph, rev(pth)), admixedges)
      lr = (nghb[1] != p[2]) + 0
      newind = cbind(admixnodes[i], i, e, match(e, normedges), lr, rep(thisleaves, each=length(e)), rep(match(thisleaves, leaves[-1]), each=length(e)))
      out = rbind(out, newind)
    }
  }
  out[!duplicated(out),]
}

fill_pwts = function(pwts, wim, weights) {
  # puts weights onto pwts, using index matrix
  # todo: make this faster by adding two columns to weightindmat with aggregation indices; first a combination of admixnodeindex, edgeid2, leafindex, for sums, then edgeid2, leafindex for products; second one has to index into the first one; then use tapply to aggregate;
  # problem: second index is shorter, can't put it in weightindmat
  
  if(length(weights)==0) return(pwts)
  wm = wim[,5] + c(1,-1)[wim[,5]+1]*weights[wim[,2]]
  agg = aggregate(wm, list(wim[,2], wim[,4], wim[,7]), sum)
  agg2 = aggregate(agg[,4], list(agg[,2], agg[,3]), prod)
  pwts[as.matrix(agg2[,1:2])] = agg2[,3]
  pwts
}

fill_pwts = function(pwts, wim, weights) {
  # puts weights onto pwts, using index matrix
  # both tibble and data.table implementation are about equally slow as aggregate
  # 5x faster than fill_pwts; todo don't compute index vectors each time, instead pass them as arguments
  
  if(length(weights)==0) return(pwts)
  #agg2 = wim %>% group_by(admixnodeindex, edgeid2, leafindex) %>% summarize(sm = sum(lr+lr*c(1,-1)[lr+1]*weights[admixnodeindex])) %>% group_by(edgeid2, leafindex) %>% summarize(prd = prod(sm))
  # wim2 = data.table::data.table(cbind(wim, wim[,5] + c(1,-1)[wim[,5]+1]*weights[wim[,2]]))
  # wim2 = wim2[,lapply(.SD, sum, .SDcols=V8), by=list(admixnodeindex, edgeid2, leafindex), ]
  # agg2 = wim2[,lapply(.SD, prod), by=list(edgeid2, leafindex), ]
  
  wm = wim[,5] + c(1,-1)[wim[,5]+1]*weights[wim[,2]]
  x = vapply(split(wm, paste(wim[,2], wim[,4], wim[,7])), sum, numeric(1))
  x2 = vapply(split(x, gsub('^.+? ', '', names(x))), prod, numeric(1))
  m = matrix(as.numeric(do.call(rbind, strsplit(names(x2), ' '))), ncol=2)
  
  pwts[m] = x2
  pwts
}



optimweightsfun = function(weights, args) {
  # likelihood function used in optimizing admixture weights
  # weights is vector of admixture weights to be optmized; only values for first incoming edge; 2nd is 1 - first
  
  pwts = args[[1]]
  weightindmat = args[[2]] # indices into pwts with weight positions
  ppinv = args[[3]]
  f3_jest = args[[4]]

  pwts = fill_pwts(pwts, weightindmat, weights)
  cmb = combn(0:ncol(pwts),2)
  ppwts_2d = t(pwts[,cmb[1,]+1]*pwts[,cmb[2,]])

  q2 = opt_edge_lengths(ppwts_2d, ppinv, f3_jest) # function call makes this ~ 3% slower
  w2 = (ppwts_2d %*% q2) - f3_jest
  lik = t(w2) %*% ppinv %*% w2
  #cat(paste0(lik[1,1], ' XXX ', paste(weights, collapse=' '), '\n'))
  lik[1,1]
}



opt_edge_lengths = function(ppwts_2d, ppinv, f3_jest) {
  # finds optimal edge lengths
  # pwts: nedge x (npop-1) design matrix with paths to outpop
  # ppinv: inverse of npair x npair matrix of varianc-covariance matrix of jackknife f3 stats
  # f3_jest: estimated f3 stats
  
  pppp = t(ppwts_2d) %*% ppinv
  cc = pppp %*% ppwts_2d
  diag(cc) = diag(cc) + mean(diag(cc))*0.0001
  cc = (cc+t(cc))/2
  q1 = -(pppp %*% f3_jest)[,1]
  q2 = pracma::quadprog(cc, q1, lb=0)$xmin
  q2
}

get_score = function(ppwts_2d, ppinv, f3_jest, q2) {
  
  w2 = (ppwts_2d %*% q2) - f3_jest
  lik = t(w2) %*% ppinv %*% w2
  lik[1,1]
}

random_newick = function(n, start='', end='') {
  # recursive function which returns topology of a random, binary tree in newick format
  # redirects to 'random_newick_named' when called with character vector of lables
  if(is.character(n)[1]) return(random_newick_named(n, start, end))
  if(n == 1) return('')
  n1 = sample(n-1,1)
  n2 = n-n1
  return(paste0(start, '(',random_newick(n1, ''),',',random_newick(n2, ''),')', end))
}

random_newick_named = function(names, start='', end='') {
  # recursive function which returns a labelled random, binary tree in newick format with named leaves in order of input
  n = length(names)
  if(n == 1) return(names)
  n1 = sample(n-1,1)
  return(paste0(start, '(',random_newick_named(names[1:n1]),',',random_newick_named(names[(n1+1):n]),')', end))
}

newick_to_edges = function(newick, node='R', edgemat=matrix(NA,0,2)) {
  # turns binary tree in newick format into matrix of edges (adjacency list)
  
  newick = gsub('^\\(', '', gsub('\\)$', '', gsub(';$', '', newick)))
  opencount = 0
  for(i in 1:nchar(newick)) {
    char = substr(newick, i,i)
    if(char == '(') opencount = opencount+1
    if(char == ')') opencount = opencount-1
    if(char == ',' && opencount == 0) break
  }
  left = substr(newick, 1, i-1)
  right = substr(newick, i+1, nchar(newick))
  stopifnot(str_count(newick, ',') >= 1)
  if(str_count(left, ',') == 0) {
    nodel = left
    edgesleft = matrix(NA, 0, 2)
  } else {
    nodel = paste0(node, 'l')
    edgesleft = newick_to_edges(left, nodel)
  }
  if(str_count(right, ',') == 0) {
    noder = right
    edgesright = matrix(NA, 0, 2)
  } else {
    noder = paste0(node, 'r')
    edgesright = newick_to_edges(right, noder)
  }
  rbind(c(node, nodel), edgesleft, c(node, noder), edgesright, edgemat)
}


insert_admix = function(edges, n=1) {
  # inserts admixture node and edges at random location
  # currently doesn't always working properly when n > 1
  # assumes earlier edges are not downstream of later edges
  # assumes first edge is R -> outpop
  
  if(n < 1) return(edges)
  edges = as.matrix(edges)
  admixedges = which(edges[,2] %in% names(which(table(edges[,2]) > 1)))
  cnt = length(admixedges)/2+1
  eg = sort(sample(setdiff(2:nrow(edges), admixedges), 2))
  o1 = edges[eg[1], ]
  o2 = edges[eg[2], ]
  newv11 = paste0(o1[1], 'x', cnt)
  newv12 = paste0(o1[1], 'xx', cnt)
  newv21 = paste0(o2[1], 'x', cnt)
  newv22 = paste0('admix', cnt)
  enew = rbind(edges[1:(eg[1]-1),], c(o1[1], newv11), c(newv11, newv12), c(newv11, o1[2]))
  if(eg[2]-eg[1] > 1) enew = rbind(enew, edges[(eg[1]+1):(eg[2]-1),])
  enew = rbind(enew, c(o2[1], newv21), c(newv21, newv22), c(newv22, o2[2]), c(newv12, newv22))
  if(nrow(edges) - eg[2] > 1) enew = rbind(enew, edges[(eg[2]+1):nrow(edges),])
  if(n == 1) return(enew)
  return(insert_admix(enew, n-1))
}


afs_to_jackknife_blocks = function(afs, popcounts, outpop=NULL, fstscale=FALSE, blgsize=0.05, maxmem=1e3) {
  
  # input:
  # data.frame with allele frequencies for each population
  # vector with pop counts
  # maxmem in MB
  # column 'cm' is required to compute jackknife blocks
  # if outpop is not null, SNPs are restricted to those which are heterozygous in the outgroup
  # outpop is moved to first position, order of other pops remain unchanged
  # output: list with two elements:
  # 3d array 'fstest_blocks' (npop x npop x nblocks) with all pairwise leave-one-out f2 stats, leaving each block out at a time
  # blockids: vector of length nsnp with jackknife block numbers
  # yscal: scaling factor applied to f2 stats
  
  infocols = 6
  snpweights = 1
  if(!is.null(outpop)) {
    stopifnot(outpop %in% names(afs))
    eps = 1e-10
    afs %<>% filter(between(!!sym(outpop), eps, 1-eps)) %>% select(1:infocols, outpop, everything())
    cat(blue(paste0(bold(nrow(afs)), ' SNPs remain which are polymorphic in outgroup ', bold(outpop), '\n')))
    stopifnot(nrow(afs) > 0)
    snpweights = 1/(afs[[outpop]] * (1-afs[[outpop]]))
    snpweights = rray(snpweights, c(1,1,length(snpweights)))
  }
  blockids = afs %>% setblocks4(dist=blgsize) %$% block
  snpinfo = afs %>% select(1:infocols)
  afs = afs %>% select(-1:-infocols) %>% as.matrix
  popcounts = as.vector(popcounts[colnames(afs)])
  
  mem1 = lobstr::obj_size(afs)
  mem2 = mem1*ncol(afs)
  cat(green(paste0('allele frequency matrix for ', bold(nrow(afs)), ' SNPs and ', bold(ncol(afs)), ' populations is ', bold(round(mem1/1e6)), ' MB\n')))
  cat(blue(paste0('matrix of pairwise f2 for all SNPs and population pairs will require ', bold(round(mem2/1e6)), ' MB\n')))
  
  if(mem2/1e6 > maxmem) {
    numsplits = ceiling(mem2/1e6/maxmem)
    cat(blue(paste0('splitting into ', bold(numsplits), ' blocks of up to ', maxmem, ' MB (',choose(numsplits,2),' block pairs)\n')))
    fstest_blocks = get_split_fstest_blocks(afs, blockids, snpweights, popcounts, numsplits=numsplits)
  } else {
    afrr1 = rray(t(afs), dim=c(ncol(afs), 1, nrow(afs)))
    afrr2 = rray_transpose(afrr1, c(2,1,3))
    pqarr = afrr1*(1-afrr1)/(2*popcounts-1) + afrr2*(1-afrr2)/t(2*popcounts-1)
    numer = (afrr1 - afrr2)^2 - pqarr
    dimnames(numer)[[1]] = dimnames(numer)[[2]] = colnames(afs)
    fstest_blocks = bj_arr_lo_mean(numer*snpweights, blockids)
    for(i in 1:dim(fstest_blocks)[1]) fstest_blocks[i,i,] = 0
  }
  
  yscal = 1
  if(fstscale) {
    denom = afrr1 + afrr2 - 2*afrr1*afrr2
    fstnum_blocks = bj_arr_lo_mean(numer, blockids)
    fstdenom_blocks = bj_arr_lo_mean(denom, blockids)
    for(i in 1:dim(fstnum_blocks)[1]) {
      fstnum_blocks[i,i,] = 0
      fstdenom_blocks[i,i,] = 0
    }
    # fix this; Nick's code calculates scaling factor based only on the populations used; does it matter if using all pops instead?
    yscal = get_yscal(fstest_blocks, fstnum_blocks, fstdenom_blocks, blockids)
    fstest_blocks = fstest_blocks*yscal
    cat(('f2 stats scaled by factor of ' %+% bold(round(yscal, 2))))
  }
  namedList(fstest_blocks, blockids, yscal)
}

get_split_fstest_blocks = function(afmat, blockids, snpweights, popcounts, numsplits) {
  # splits afmat into blocks by column, computes lo jackknife blocks on each pair of blocks, and combines into 3d array
  
  width = ceiling(ncol(afmat)/numsplits)
  starts = seq(1, ncol(afmat), width)
  ends = c(na.omit(lead(starts)-1), ncol(afmat))
  cmb = combn(0:numsplits, 2)
  cmb[1,] = cmb[1,]+1
  arrlist = replicate(numsplits, list())
  nsnp = nrow(afmat)
  
  for(i in 1:ncol(cmb)) {
    cat(red(paste0('pop combination ', i, ' out of ', ncol(cmb), '\n')))
    c1 = cmb[1,i]
    c2 = cmb[2,i]
    from1 = starts[c1]
    from2 = starts[c2]
    to1 = ends[c1]
    to2 = ends[c2]
    b1 = afmat[,from1:to1]
    b2 = afmat[,from2:to2]
    afrr1 = rray(t(b1), dim=c(ncol(b1), 1, nsnp))
    afrr2 = rray(t(b2), dim=c(1, ncol(b2), nsnp))
    pqarr = afrr1*(1-afrr1)/(2*popcounts[from1:to1]-1) + afrr2*(1-afrr2)/t(2*popcounts[from2:to2]-1)
    numer = (afrr1 - afrr2)^2 - pqarr
    dimnames(numer)[[1]] = colnames(b1)
    dimnames(numer)[[2]] = colnames(b2)
    arrlist[[c1]][[c2]] = bj_arr_lo_mean(numer*snpweights, blockids)
    arrlist[[c2]][[c1]] = aperm(arrlist[[c1]][[c2]], c(2,1,3))
  }
  fstest_blocks = do.call(abind, list(lapply(arrlist, function(x) do.call(abind, list(x, along=2))), along=1))
  for(i in 1:dim(fstest_blocks)[1]) fstest_blocks[i,i,] = 0
  fstest_blocks
}

qpGraphR_slow = function(parfile, graphfile, na.action = 'remove', blocksize=1000) {
  # runs qpGraphR for given parfile and graphfile (file locations)
  # running this is not recommended, as it reads in the genotype data and computes block-jackknife f-statistics each time
  # it will be much faster to instead run 'packedancestrymap_to_aftable' and 'afs_to_jackknife_blocks' first, and then 'qpGraphR'
  
  pars = parse_qpGraph_parfile(parfile)
  implemented = names(as.list(args(qpGraphR)))
  cat(red(paste0('The following parameters are not implemented:\n', paste0(setdiff(names(pars), c(implemented, 'pref', 'outpop')), collapse=', '), '\n')))
  indfile = read_table2(paste0(pars$pref, '.ind'), col_names = FALSE, col_types = cols())
  
  edges = parse_qpGraph_graphfile(graphfile)
  pops = setdiff(edges[,2], edges[,1])
  popcounts = indfile %>% filter(X3 %in% pops) %$% table(X3)

  afs = packedancestrymap_to_aftable(pars$pref, pops=pops, blocksize=1000, na.action = na.action)
  bj = afs_to_jackknife_blocks(afs=afs, popcounts=popcounts, outpop=pars$outpop, blgsize=pars$blgsize)
  pars$fstest_blocks = bj[[1]]
  pars$blockids = bj[[2]]
  pars$graph = graphfile
  
  myout = do.call(qpg, pars[names(pars) %in% implemented])

}


qpGraphR = function(fstest_blocks, blockids, graph, lsqmode=FALSE, fnscale=1e-6, fudge=1e-5, numstart=NULL, seed=123, verbose=TRUE) {
  # modelled after AdmixTools qpGraph
  
  # input:
  # fstest_blocks: block-jackknife leave-one-block-out estimates of f2 statistics (weighed by inverse of outgroup heterozygosity)
  # blockids: vector of length nsnps which was used in creating fstest_blocks
  # graph: either 2 column data.frame or matrix of edges specifying an admixture graph (should be DAG), or the location of a qpGraph format graph file
  # numstart: number of different random combinations of starting values; defaults to 10*nadmix
  # first edge has to be root -> outpop; need to implement topology check

  # output list:
  # edges: input graph with added branch lengths / admixture weights
  # score: likelihood of this tree
  # f3_jest: estimated f3 statistics

  # make graph
  if(length(graph) == 1) edges = parse_qpGraph_graphfile(graph)
  else edges = as.matrix(graph)
  grph = graph_from_edgelist(edges)
  nedges = length(E(grph))
  
  admixnodes = which(degree(grph, mode='in')==2)
  nadmix = length(admixnodes)
  admixedgesfull = sapply(seq_len(nadmix), function(i) incident_edges(grph, admixnodes, mode='in')[[i]][1:2])
  normedges = setdiff(1:nedges, admixedgesfull)
  weightindmat = graph_to_weightind(grph)
  pwts = graph_to_pwts(grph)

  pops = V(grph)$name[degree(grph, v = V(grph), mode = c('out')) == 0]
  stopifnot(all(pops %in% dimnames(fstest_blocks)[[1]]))
  popind = match(pops, dimnames(fstest_blocks)[[1]])
  
  npop = length(pops)
  npair = choose(npop, 2)

  # process f-stats
  f2 = bj_arr_stats(fstest_blocks[popind,popind,], blockids)
  f2out = tibble(pop1=combn(pops, 2)[1,], pop2=combn(pops, 2)[2,], f2est = f2[[1]][lower.tri(f2[[1]])], se = sqrt(f2[[2]][lower.tri(f2[[2]])]))
  # todo: write function to compute all pairwise fitted f2 stats from edge weights (if that's an interesting output)

  f3_blocks = (fstest_blocks[,popind[1],] + fstest_blocks[popind[1],,] - fstest_blocks)/2
  f3_blocks_2d = arr3d_to_pairmat(f3_blocks[popind[-1],popind[-1],])
  sts = bj_mat_stats(f3_blocks_2d, blockids)
  f3_jest = sts[[1]]
  f3_jvar = sts[[2]]
  f3out = tibble(pop1=pops[combn(0:(npop-1), 2)[1,]+1], pop2=pops[combn(0:(npop-1), 2)[2,]], f3est = f3_jest, se = sqrt(diag(f3_jvar)))
  diag(f3_jvar) = diag(f3_jvar) + sum(diag(f3_jvar))*fudge
  # in qpGraph fudge is 1e-5; sometimes quadprog doesn't converge unless this is larger; has large effect on magnitude of likelihood score
  
  if(lsqmode) ppinv = diag(1/diag(f3_jvar))
  else ppinv = solve(f3_jvar)
  
  if(nadmix > 0) {
    #parmat = as.matrix(expand.grid(replicate(nadmix, seq(0,1,len=numgrid), simplify = FALSE)))
    if(is.null(numstart)) numstart = 10*nadmix
    set.seed(seed)
    parmat = matrix(runif(numstart*nadmix), numstart)
    if(verbose) cat(blue(paste0('testing ', nrow(parmat), ' combinations of admixture weight starting values\n')))
    opt = multistart(parmat, optimweightsfun, args=list(pwts, weightindmat, ppinv, f3_jest), method='L-BFGS-B', lower=0, upper=1, control=list(maxit=1e4, fnscale=fnscale), verbose=verbose)
    best = opt %>% top_n(1, -value)
    wts = as.matrix(best[,1:nadmix])[1,]
    score0 = best$value
  } else {
    wts = numeric()
  }
  
  # startval = runif(nadmix)
  # startval = rep(0.9, nadmix)
  # opt = optim(par = startval, optimweightsfun, args=list(pwts, weightindmat, ppinv, f3_jest), method='L-BFGS-B', lower=0, upper=1, control=list(maxit=1e4, fnscale=fnscale))
  # wts = opt$par
  # score = opt$value

  pwts = fill_pwts(pwts, weightindmat, wts)
  cmb = combn(0:ncol(pwts),2)
  ppwts_2d = t(pwts[,cmb[1,]+1]*pwts[,cmb[2,]])
  q2 = opt_edge_lengths(ppwts_2d, ppinv, f3_jest)
  score = get_score(ppwts_2d, ppinv, f3_jest, q2)
  
  weight = rep(NA, nedges)
  if(nadmix > 0) {
    weight[admixedgesfull[1,]] = wts
    weight[admixedgesfull[2,]] = 1-wts
  }
  weight[normedges] = q2
  edgesout = as_tibble(edges) %>% set_colnames(c('from', 'to')) %>% mutate(type = ifelse(1:n() %in% normedges, 'edge', 'admix'), weight=weight) 
  
  namedList(edgesout, score, f2out, f3out, opt)
}






