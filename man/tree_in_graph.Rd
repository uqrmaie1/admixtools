% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toposearch.R
\name{tree_in_graph}
\alias{tree_in_graph}
\title{Test if a tree is part of a graph}
\usage{
tree_in_graph(tree, graph)
}
\arguments{
\item{tree}{An admixture graph without admixture event}

\item{graph}{An admixture graph}
}
\value{
\code{TRUE} if all admixture constraints are satisfied, else \code{FALSE}
\dontrun{
tree = graph_splittrees(example_igraph) %>% pull(graph) %>% pluck(1)
tree_in_graph(tree, example_igraph)
}
}
\description{
This function tests whether a tree is part of a graph. This is useful for testing whether a Y-chromosome tree is consistent with an autosomal admixture graph. Leaf node names matter, but internal node names are ignored.
}
