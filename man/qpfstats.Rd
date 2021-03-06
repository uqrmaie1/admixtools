% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{qpfstats}
\alias{qpfstats}
\title{Get smoothed f2-statistics}
\usage{
qpfstats(
  pref,
  pops,
  include_f2 = TRUE,
  include_f3 = TRUE,
  include_f4 = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{pref}{Prefix of genotype files}

\item{pops}{Populations for which to compute f2-statistics}

\item{include_f2}{Should f2-statistics be used to get smoothed f2-statistics?
If \code{include_f2} is a positive integer, it specifies how many randomly chosen f2-statistics should be used.}

\item{include_f3}{Should f3-statistics be used to get smoothed f2-statistics?
If \code{include_f3} is a positive integer, it specifies how many randomly chosen f3-statistics should be used.}

\item{include_f4}{Should f4-statistics be used to get smoothed f2-statistics?
If \code{include_f4} is a positive integer, it specifies how many randomly chosen f4-statistics should be used.}
}
\value{
A 3d-array of smoothed f2-statistics
}
\description{
This function returns an array of (pseudo-) f2-statistics which are computed by
taking into account other f2-, f3-, and f4-statistics. The advantage of doing that
is that f3- and f4-statistics computed from these smoothed f2-statistics can be
more accurate for populations with large amounts of missing data. The function
uses SNPs which are missing in some populations in a manner which tends to introduce
less bias than setting \code{maxmiss} to values greater than 0.
}
\examples{
\dontrun{
f2_blocks = qpfstats(geno_prefix, mypops)
}
}
