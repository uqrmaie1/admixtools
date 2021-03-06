% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toposearch.R
\name{graph_plusone}
\alias{graph_plusone}
\title{Find all graphs which result from adding one admixture edge}
\usage{
graph_plusone(graph, ntry = Inf)
}
\arguments{
\item{graph}{Admixture graph in \code{igraph} format}

\item{ntry}{Specify this to return only a subset of all possible graphs with one more edge}
}
\value{
A data frame with columns \code{from}, \code{to}, and \code{graph}
}
\description{
Find all graphs which result from adding one admixture edge
}
\examples{
\dontrun{
newgraphs = graph_plusone(example_igraph)
# now evaluate the new graphs
newgraphs \%>\%
  rowwise \%>\%
  mutate(res = list(qpgraph(example_f2_blocks, graph))) \%>\%
  unnest_wider(res)
}
}
