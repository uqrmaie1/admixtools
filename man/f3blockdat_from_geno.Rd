% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{f3blockdat_from_geno}
\alias{f3blockdat_from_geno}
\title{f3 from genotype data}
\usage{
f3blockdat_from_geno(
  pref,
  popcombs,
  auto_only = TRUE,
  blgsize = 0.05,
  block_lengths = NULL,
  allsnps = FALSE,
  poly_only = FALSE,
  apply_corr = TRUE,
  outgroupmode = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{pref}{Prefix of genotype files}

\item{popcombs}{A data frame with one population combination per row, and columns \code{pop1}, \code{pop2}, \code{pop3}, \code{pop4}. If there is an additional integer column named \code{model} and \code{allsnps = FALSE}, only SNPs present in every population in any given model will be used to compute f4-statistics for that model.}

\item{auto_only}{Use only chromosomes 1 to 22.}

\item{blgsize}{SNP block size in Morgan. Default is 0.05 (50 cM). If \code{blgsize} is 100 or greater, if will be interpreted as base pair distance rather than centimorgan distance.}

\item{block_lengths}{An optional vector with block lengths. If \code{NULL}, block lengths will be computed.}

\item{allsnps}{Use all SNPs with allele frequency estimates in every population of any given population quadruple. If \code{FALSE} (the default) only SNPs which are present in all populations in \code{popcombs} (or any given model in it) will be used. Setting \code{allsnps = TRUE} in the presence of large amounts of missing data might lead to false positive results.}

\item{verbose}{Print progress updates}
}
\value{
A data frame with per-block f4-statistics for each population quadruple.
}
\description{
Compute per-block f3-statistics directly from genotype data
}
