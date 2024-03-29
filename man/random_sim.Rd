% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{random_sim}
\alias{random_sim}
\title{Generate a random graph and simulate it in msprime v1.x}
\usage{
random_sim(
  nleaf,
  nadmix,
  outpref = "random_sim",
  max_depth = NULL,
  ind_per_pop = 1,
  mutation_rate = 1.25e-08,
  admix_weights = 0.5,
  neff = 1000,
  time = 1000,
  fix_leaf = FALSE,
  outpop = NULL,
  nchr = 1,
  recomb_rate = 2e-08,
  seq_length = 1000,
  ghost_lineages = TRUE,
  run = FALSE
)
}
\arguments{
\item{nleaf}{The number of leaf nodes}

\item{nadmix}{The number of admixture events}

\item{outpref}{A prefix of output files}

\item{max_depth}{A constraint specifying the maximum time depth of the admixture graph (in generations)}

\item{ind_per_pop}{The number of individuals to simulate for each population. If a scalar value, it will be constant across all populations. Alternatively, it can be a named vector with a different value for each population.}

\item{mutation_rate}{Mutation rate per site per generation. The default is \code{1.25e-8} per base pair per generation.}

\item{admix_weights}{Admixture weights. If a float value (0 < value < 1), admixture weights for each admixture event will be (value, 1-value).
Alternatively, it can be a range, i.e., \code{c(0.1, 0.4)} specifying lower and upper limits of a uniform distribution from which the admixture weight value will be drawn. By default, all admixture edges have a weight of 0.5.}

\item{neff}{Effective population size (in diploid individuals). If a scalar value, it will be constant across all populations. Alternatively, it can be a range, i.e., \code{c(500, 1000)} specifying lower and upper limits of an uniform distribution from which values will be drawn}

\item{time}{Time between nodes. Either a scalar value (1000 by default) with the dates generated by \code{\link{pseudo_dates}}, or a range, i.e., \code{c(500, 1000)} specifying lower and upper limits of a uniform distribution from which values will be drawn (see \code{\link{random_dates}})}

\item{fix_leaf}{A boolean specifying if the dates of the leaf nodes will be fixed at time 0. If \code{TRUE}, all samples will be drawn at the end of the simulation (i.e., from “today”).}

\item{outpop}{A name of the (optional) outgroup population.}

\item{nchr}{The number of chromosomes to simulate}

\item{recomb_rate}{A float value specifying recombination rate along the chromosomes. The default is \code{2e-8} per base pair per generation.}

\item{seq_length}{The sequence length of the chromosomes. If it is a scalar value, the sequence length will be constant for all chromosomes.
Alternatively, it can be a vector with a length equal to number of chromosomes (i.e., \code{c(100,50)} to simulate 2 chromosomes with the lengths of 100 and 50 base pairs).}

\item{ghost_lineages}{A boolean value specifying whether ghost lineages will be allowed.
If \code{TRUE}, admixture happens at the time points defined by the y-axis generated while plotting the graph by \code{\link{plot_graph}}
If \code{FALSE} (default), admixture occurs at the time of the previous split event}

\item{run}{If \code{FALSE}, the function will terminate after writing the msprime script. If \code{TRUE}, it will try and execute the function with the default python installation.
If you want to use some other python installation, you can set \verb{run = /my/python}.}
}
\value{
A list with the path of simulation script, a data frame of graph edges, dates and population sizes:
\itemize{
\item{out - }{A file name and path of the simulation script}
\item{edges - }{An edge dataframe with admixture weights}
\item{dates - }{A named vector with a date for each node}
\item{neffs - }{A named vector with an effective population size for each node}
}
}
\description{
This is basically a wrapper function around the \code{\link{msprime_genome}} that allows user to create a random graph and simulate it in msprime v1.x.
}
\examples{
# Create simulation script that simulates 2 chromosomes that are 50base long
# where maximum depth of the tree is 5000 generations, and plot the output graph
\dontrun{
out = random_sim(nleaf=4, nadmix=0, max_depth=5000, nchr=2, seq_length=50)
plot_graph(out$edges, dates = out$dates, neff = out$neff, hide_weights = TRUE)
}
}
