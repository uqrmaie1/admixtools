% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toposearch.R
\name{find_normedges}
\alias{find_normedges}
\title{Find drift edges}
\usage{
find_normedges(graph, exclude_first = FALSE)
}
\arguments{
\item{graph}{An admixture graph}

\item{exclude_first}{Do not return edge from root to outgroup}
}
\value{
A data frame with columns \code{from} and \code{to} with drift edges
}
\description{
Find drift edges
}
\seealso{
\code{\link{find_newedges}} \code{\link{find_admixedges}}
}
