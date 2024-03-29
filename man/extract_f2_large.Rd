% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{extract_f2_large}
\alias{extract_f2_large}
\title{Compute and store blocked f2 statistics}
\usage{
extract_f2_large(
  pref,
  outdir,
  inds = NULL,
  pops = NULL,
  blgsize = 0.05,
  cols_per_chunk = 10,
  maxmiss = 0,
  minmaf = 0,
  maxmaf = 0.5,
  minac2 = FALSE,
  outpop = NULL,
  outpop_scale = TRUE,
  transitions = TRUE,
  transversions = TRUE,
  keepsnps = NULL,
  snpblocks = NULL,
  overwrite = FALSE,
  format = NULL,
  adjust_pseudohaploid = TRUE,
  afprod = TRUE,
  fst = TRUE,
  poly_only = c("f2"),
  apply_corr = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{pref}{Prefix of \emph{PLINK/EIGENSTRAT/PACKEDANCESTRYMAP} files.
\emph{EIGENSTRAT/PACKEDANCESTRYMAP} have to end in \code{.geno}, \code{.snp}, \code{.ind}, \emph{PLINK} has to end in \code{.bed}, \code{.bim}, \code{.fam}}

\item{outdir}{Directory where data will be stored.}

\item{inds}{Individuals for which data should be extracted}

\item{pops}{Populations for which data should be extracted. If both \code{pops} and \code{inds} are provided, they should have the same length and will be matched by position. If only \code{pops} is provided, all individuals from the \code{.ind} or \code{.fam} file in those populations will be extracted. If only \code{inds} is provided, each indivdual will be assigned to its own population of the same name. If neither \code{pops} nor \code{inds} is provided, all individuals and populations in the \code{.ind} or \code{.fam} file will be extracted.}

\item{blgsize}{SNP block size in Morgan. Default is 0.05 (5 cM). If \code{blgsize} is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.}

\item{cols_per_chunk}{Number of populations per chunk. Lowering this number will lower the memory requirements when running}

\item{maxmiss}{Discard SNPs which are missing in a fraction of populations higher than \code{maxmiss}}

\item{minmaf}{Discard SNPs with minor allele frequency less than \code{minmaf}}

\item{maxmaf}{Discard SNPs with minor allele frequency greater than than \code{maxmaf}}

\item{minac2}{Discard SNPs with allele count lower than 2 in any population (default \code{FALSE}). This option should be set to \code{TRUE} when computing f3-statistics where one population consists mostly of pseudohaploid samples. Otherwise heterozygosity estimates and thus f3-estimates can be biased. \code{minac2 == 2} will discard SNPs with allele count lower than 2 in any non-singleton population (this option is experimental and is based on the hypothesis that using SNPs with allele count lower than 2 only leads to biases in non-singleton populations). While the \code{minac2} option discards SNPs with allele count lower than 2 in any population, the \code{\link{qp3pop}} function will only discard SNPs with allele count lower than 2 in the first (target) population (when the first argument is the prefix of a genotype file).}

\item{outpop}{Keep only SNPs which are heterozygous in this population}

\item{outpop_scale}{Scale f2-statistics by the inverse \code{outpop} heteroygosity (\code{1/(p*(1-p))}). Providing \code{outpop} and setting \code{outpop_scale} to \code{TRUE} will give the same results as the original \emph{qpGraph} when the \code{outpop} parameter has been set, but it has the disadvantage of treating one population different from the others. This may limit the use of these f2-statistics for other models.}

\item{transitions}{Set this to \code{FALSE} to exclude transition SNPs}

\item{transversions}{Set this to \code{FALSE} to exclude transversion SNPs}

\item{keepsnps}{SNP IDs of SNPs to keep. Overrides other SNP filtering options}

\item{overwrite}{Overwrite existing files in \code{outdir}}

\item{format}{Supply this if the prefix can refer to genotype data in different formats
and you want to choose which one to read. Should be \code{plink} to read \code{.bed}, \code{.bim}, \code{.fam} files, or \code{eigenstrat}, or \code{packedancestrymap} to read \code{.geno}, \code{.snp}, \code{.ind} files.}

\item{adjust_pseudohaploid}{Genotypes of pseudohaploid samples are usually coded as \code{0} or \code{2}, even though only one allele is observed. \code{adjust_pseudohaploid} ensures that the observed allele count increases only by \code{1} for each pseudohaploid sample. If \code{TRUE} (default), samples that don't have any genotypes coded as \code{1} among the first 1000 SNPs are automatically identified as pseudohaploid. This leads to slightly more accurate estimates of f-statistics. Setting this parameter to \code{FALSE} treats all samples as diploid and is equivalent to the \emph{ADMIXTOOLS} \code{inbreed: NO} option. Setting \code{adjust_pseudohaploid} to an integer \code{n} will check the first \code{n} SNPs instead of the first 1000 SNPs.}

\item{afprod}{Write files with allele frequency products for every population pair. Setting this to FALSE can make \code{extract_f2} faster and will require less memory.}

\item{fst}{Write files with pairwise FST for every population pair. Setting this to FALSE can make \code{extract_f2} faster and will require less memory.}

\item{poly_only}{Specify whether SNPs with identical allele frequencies in every population should be discarded (\code{poly_only = TRUE}), or whether they should be used (\code{poly_only = FALSE}). By default (\code{poly_only = c("f2")}), these SNPs will be used to compute FST and allele frequency products, but not to compute f2 (this is the default option in the original ADMIXTOOLS).}

\item{apply_corr}{Apply small-sample-size correction when computing f2-statistics (default \code{TRUE})}

\item{verbose}{Print progress updates}
}
\value{
SNP metadata (invisibly)
}
\description{
\code{extract_f2_large} does the same as \code{\link{extract_f2}}, but it requires less memory and is slower. \code{outdir} has to be set in \code{extract_f2_large}.
}
\details{
\code{extract_f2_large} requires less memory because it writes allele frequency data to disk, and doesn't store the allele frequency matrix for all populations and SNPs in memory. If you still run out of memory, reduce \code{cols_per_chunk}. This function is a wrapper around \code{\link{extract_afs}} and \code{\link{afs_to_f2}}, and is slower than \code{\link{extract_f2}}. It may be faster to call \code{\link{extract_afs}} and \code{\link{afs_to_f2}} directly, parallelizing over the different calls to \code{\link{afs_to_f2}}.
}
\examples{
\dontrun{
pref = 'my/genofiles/prefix'
f2dir = 'my/f2dir/'
extract_f2_large(pref, f2dir, pops = c('popA', 'popB', 'popC'))
}
}
\seealso{
\code{\link{extract_f2}}
}
