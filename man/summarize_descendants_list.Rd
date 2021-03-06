% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toposearch.R
\name{summarize_descendants_list}
\alias{summarize_descendants_list}
\title{List leaf nodes for all internal nodes in a list of graphs}
\usage{
summarize_descendants_list(graphlist, rename = FALSE)
}
\arguments{
\item{graphlist}{A list of admixture graphs}

\item{rename}{If \code{FALSE} (default), the output will be a data frame indicating how often each node in each graph is observed in all other graphs. If \code{TRUE}, the output will be a list, where the inner nodes will be renamed to have these percentages as part of their name. \code{\link{plot_graph}} will print the percentages of graphs renamed in this way.}
}
\value{
A data frame with columns \code{graph}, \code{from}, \code{n}, \code{frac}
}
\description{
List leaf nodes for all internal nodes in a list of graphs
}
\examples{
\dontrun{
summarize_descendants_list(graphlist)
}
}
