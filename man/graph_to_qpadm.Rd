% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toposearch.R
\name{graph_to_qpadm}
\alias{graph_to_qpadm}
\title{Get all qpadm models for a graph}
\usage{
graph_to_qpadm(
  graph,
  target,
  left = NULL,
  right = NULL,
  models = NULL,
  weights = TRUE,
  f4dat = NULL,
  allpops = TRUE,
  more_right = TRUE,
  return_f4 = FALSE,
  eps = 1e-10
)
}
\arguments{
\item{graph}{An admixture graph}

\item{target}{Name of the target population.}

\item{left}{Left populations (provide this optionally if you want to test only a single qpadm model)}

\item{right}{Right populations (provide this optionally if you want to test only a single qpadm model)}

\item{models}{A two column nested data frame with models to be evaluated, one model per row. The first column, \code{l}, should contain the left populations, the second column, \code{r}, should contain the right populations. The target population is provided separately in the \code{target} argument.}

\item{weights}{Set this to \code{FALSE} to return only information on the ranks, not the weights, of each qpadm model. The ranks should depend only on the graph topology, while the weights and weight-validity (all weights for left populations between 0 and 1) can depend on the branch lengths of the graph. By default f4-statistics are based on equal branch lengths and admixture weights of 0.5. This can be overridden by providing \code{f4dat}.}

\item{f4dat}{A data frame of f4-statistics which can be provided to override the default branch lengths.}

\item{allpops}{Evaluate only models which use all populations in the admixture graph. See \code{\link{qpadm_models}}}

\item{return_f4}{Include f4 statistic matrices in the results (default \code{FALSE})}

\item{eps}{Epsilon value close to zero which is used for determining which f4 matrix elements should be considered non-zero, and which weights are strictly between 0 and 1.}
}
\value{
A data frame with one qpadm model per row and columns \code{valid_rank} and \code{valid_weights} indicating whether a model should be valid under the graph.
}
\description{
This function tests which qpadm models should be valid for an admixture graph and a target population. By default, all models returned by \code{\link{qpadm_models}} are tested. For large graphs this will be too slow, and you may want test only some models by providing the \code{models} argument, or only a single model by providing the \code{left} and \code{right} arguments.
}
\details{
Two validity criteria are tested for each qpadm model: Rank validity and weight validity. Rank validity means that the rank of the f4 statistic matrix for only left and right populations is the same as the rank of the f4 statistic matrix that includes the target population among the left populations. Weight validity means that the estimated admixture weights for all left populations are between 0 and 1.
}
\note{
An earlier version of this function tried to use the graph topology for identifying valid qpadm models, but this didn't work reliably. Christian Huber had the idea of using the ranks of expected f4 statistic matrices instead.
}
\examples{
\dontrun{
graph2 = example_igraph \%>\% simplify_graph() \%>\%
  delete_admix('N3N0', 'N3N4') \%>\% delete_admix('N3N1', 'N3N8')
graph_to_qpadm(graph2, 'Mbuti.DG') \%>\% filter(valid_rank, valid_weights)
}
}
\seealso{
\code{\link{qpadm_models}}, \code{\link{graph_f2_function}}
}
