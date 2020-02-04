#' @import dplyr
#' @import readr
#' @import purrr
#' @import tidyr
#' @import stringr
#' @import ggplot2
#' @import utils
#' @importFrom magrittr set_colnames "%$%" "%<>%"
#' @importFrom rray rray rray_transpose
#' @importFrom lobstr obj_size
#' @importFrom abind abind
#' @importFrom crayon blue red green bold italic
#' @importFrom tibble as_tibble deframe enframe add_column
#' @importFrom stats na.omit setNames runif
#' @importFrom grDevices hcl
#' @importFrom rlang .data
#' @importFrom plotly plot_ly add_trace add_markers layout
#' @importFrom Rcpp cppFunction
#' @importFrom igraph V E neighbors subcomponent get.edge.ids degree incident_edges all_simple_paths
#' graph_from_edgelist as_edgelist shortest_paths
#' @useDynLib admixtools

TRUE

#' Lightweight R implementations of some AdmixTools programs.
#'
#' The admixtools package provides implementations of some AdmixTools programs
#' which make use of pre-computed f2 block jackknife statistics.
#' Aside from the core functions \code{\link{qp3pop}}, \code{\link{qpdstat}}, \code{\link{qpadm}},
#' and \code{\link{qpgraph}}, there are also wrapper functions and parsing functions around the original
#' AdmixTools software, as well as functions to read genotype files and precompute the data necessary
#' to quickly compute f statistics.
#'
#' @author Robert Maier \email{<rmaier@@broadinstitute.org>}
#' @author Nick Patterson
#'
#' @references Patterson, N. et al. (2012) \emph{Ancient admixture in human history.} Genetics
#'
#' @seealso
#' \url{https://uqrmaie1.github.io/admixtools/index.html}
#'
#' @docType package
#' @name admixtools
NULL

