#' @import dplyr
#' @import readr
#' @import purrr
#' @import tidyr
#' @import stringr
#' @import ggplot2
#' @import utils
#' @importFrom shiny isRunning
#' @importFrom magrittr set_colnames "%$%" "%<>%"
#' @importFrom lobstr obj_size
#' @importFrom abind abind
#' @importFrom crayon blue red green bold italic
#' @importFrom tibble as_tibble deframe enframe add_column rownames_to_column
#' @importFrom stats na.omit setNames runif
#' @importFrom grDevices hcl
#' @importFrom rlang .data
#' @importFrom plotly plot_ly add_trace add_markers layout
#' @importFrom Rcpp cppFunction
#' @importFrom igraph V E neighbors subcomponent get.edge.ids degree incident_edges all_simple_paths
#' graph_from_edgelist as_edgelist shortest_paths adjacent_vertices is.dag as_ids add_vertices
#' add_edges delete_edges difference set_vertex_attr
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


#a = devtools::check()
#a$notes[[2]] %>% str_split('’\n') %>% `[[`(1) %>% str_extract_all('‘.+') %>% unlist %>% str_sub(2) %>% unique %>% paste0(., collapse = '",\n"') %>% paste0('globalVariables(c("', ., '"))') %>% write_lines('r/admixtools.R', append = T)


