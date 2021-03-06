% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toposearch.R
\name{random_admixturegraph}
\alias{random_admixturegraph}
\title{Generate a random admixture graph}
\usage{
random_admixturegraph(
  leaves,
  numadmix = 0,
  simple = TRUE,
  outpop = NULL,
  nonzero_f4 = NULL,
  admix_constraints = NULL,
  event_order = NULL,
  ntry = 100
)
}
\arguments{
\item{leaves}{Names of the leaf nodes, or a number specifying how many leaf nodes there should be}

\item{numadmix}{Number of admixture events}

\item{simple}{Should edges leading to admixture nodes consist of separate admix edges and normal edges}

\item{outpop}{Outgroup population}
}
\description{
This function randomly generates an admixture graph for a given set of leaf nodes
}
\examples{
rand_graph = random_admixturegraph(10, numadmix = 5)
plot_graph(rand_graph)
}
