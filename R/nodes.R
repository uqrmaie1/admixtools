# admixtools typed condition classes ----------------------------------
# This file establishes the typed error / warning classes used by the
# sampled-internal-nodes feature. Error classes are not "registered"
# anywhere; they are conventions on the `class =` argument to
# rlang::abort() and rlang::warn(). Documented here as the single
# source of truth.
#
# Errors (rlang::abort):
#   admixtools_invalid_graph
#       strict-validation failure: orphan rows, duplicate name,
#       malformed nodes tibble, time-invariant break.
#   admixtools_pop_missing_in_f2_cache (reserved; not yet thrown)
#       qpgraph entry: sampled population absent from f2 cache.
#
# Warnings (rlang::warn):
#   admixtools_duplicate_node_name
#       non-strict prune of duplicate name (last wins).
#   admixtools_ambiguous_nodes_metadata
#       igraph input with both vertex_attr and R attr("nodes").
#   admixtools_orphan_nodes_pruned
#       non-strict prune of nodes-tibble rows without edges.
#   admixtools_dropping_node_attrs
#       as_edge_tibble discards a non-empty nodes tibble.
#   admixtools_samples_arg_overridden
#       graph_to_lgo samples= conflicts with nodes attr (nodes attr wins).
#   admixtools_nodes_twoN_overridden
#       graph_to_lgo twoN= overrides a captured nodes-attr twoN (arg wins).
# ---------------------------------------------------------------------

#' Return the nodes attribute of an edge tibble.
#'
#' Returns the `nodes` tibble carried as an attribute on an edge tibble.
#' Returns an empty zero-row tibble (with the canonical seven columns)
#' when no attribute is attached. Errors on igraph or other
#' non-data.frame inputs.
#'
#' @param graph An edge tibble.
#' @return A tibble with columns `name`, `samples`, `twoN_param`,
#'   `twoN`, `time_param`, `time`, `admix_event_time`.
#' @export
graph_nodes <- function(graph) {
  if (!is.data.frame(graph))
    rlang::abort(
      "`graph` must be an edge tibble (data.frame).",
      class = "admixtools_invalid_graph")
  nt <- attr(graph, "nodes")
  if (is.null(nt))
    return(tibble::tibble(
      name             = character(0),
      samples          = integer(0),
      twoN_param       = character(0),
      twoN             = numeric(0),
      time_param       = character(0),
      time             = numeric(0),
      admix_event_time = numeric(0)
    ))
  nt
}

#' Set or update node attributes.
#'
#' @param graph An edge tibble.
#' @param name Character vector of node names. Must be unique, must
#'   exist in the canonical node set.
#' @param samples Integer scalar or vector of length(name);
#'   NULL leaves unchanged.
#' @param twoN_param Character scalar or vector; NULL leaves unchanged.
#' @param twoN Numeric scalar or vector; NULL leaves unchanged.
#' @param time_param Character scalar or vector; NULL leaves unchanged.
#' @param time Numeric scalar or vector; NULL leaves unchanged.
#' @param admix_event_time Numeric scalar or vector; NULL leaves unchanged.
#'
#' @details
#' Scalar arguments recycle to length(name). Vector arguments must
#' match length(name) exactly. Unknown names error with
#' `admixtools_invalid_graph`. Duplicate names error with
#' `admixtools_invalid_graph`. NA in a scalar or vector sets the cell
#' to NA.
#'
#' @return The graph with updated nodes attribute. Refreshes edge time
#'   views if `time` or `admix_event_time` was set.
#' @export
set_node_attrs <- function(graph, name,
                           samples = NULL, twoN_param = NULL,
                           twoN = NULL, time_param = NULL, time = NULL,
                           admix_event_time = NULL) {
  # 1. Validate name
  if (!is.character(name) || any(is.na(name)) || length(name) == 0)
    rlang::abort(
      "`name` must be a non-empty character vector with no NAs.",
      class = "admixtools_invalid_graph")
  if (anyDuplicated(name))
    rlang::abort(
      sprintf("`name` contains duplicates: %s",
              paste(name[duplicated(name)], collapse = ", ")),
      class = "admixtools_invalid_graph")
  canonical <- unique(c(graph$from, graph$to))
  unknown <- setdiff(name, canonical)
  if (length(unknown) > 0)
    rlang::abort(
      sprintf("Unknown node names not in graph: %s",
              paste(unknown, collapse = ", ")),
      class = "admixtools_invalid_graph")

  # 2. Validate length of each non-NULL update
  validate_len <- function(arg, argname) {
    if (is.null(arg)) return(NULL)
    if (length(arg) == 1) return(rep(arg, length(name)))
    if (length(arg) != length(name))
      rlang::abort(
        sprintf("`%s` has length %d; must be 1 or length(name)=%d.",
                argname, length(arg), length(name)),
        class = "admixtools_invalid_graph")
    arg
  }
  samples          <- validate_len(samples,          "samples")
  twoN_param       <- validate_len(twoN_param,       "twoN_param")
  twoN             <- validate_len(twoN,             "twoN")
  time_param       <- validate_len(time_param,       "time_param")
  time             <- validate_len(time,             "time")
  admix_event_time <- validate_len(admix_event_time, "admix_event_time")

  # 3. Coerce types
  if (!is.null(samples) && !is.integer(samples))
    samples <- as.integer(samples)

  # 4. Merge into nodes tibble
  nt <- graph_nodes(graph)
  existing <- match(name, nt$name)
  new_rows <- is.na(existing)
  if (any(new_rows)) {
    # Seed new rows with typed NAs so column types do not get clobbered
    add <- tibble::tibble(
      name             = name[new_rows],
      samples          = NA_integer_,
      twoN_param       = NA_character_,
      twoN             = NA_real_,
      time_param       = NA_character_,
      time             = NA_real_,
      admix_event_time = NA_real_
    )
    nt <- dplyr::bind_rows(nt, add)
    existing <- match(name, nt$name)
  }
  if (!is.null(samples))          nt$samples[existing]          <- samples
  if (!is.null(twoN_param))       nt$twoN_param[existing]       <- twoN_param
  if (!is.null(twoN))             nt$twoN[existing]             <- twoN
  if (!is.null(time_param))       nt$time_param[existing]       <- time_param
  if (!is.null(time))             nt$time[existing]             <- time
  if (!is.null(admix_event_time)) nt$admix_event_time[existing] <- admix_event_time

  attr(graph, "nodes") <- nt

  # 5. Write denormalized edge time views for only the nodes named in this
  # call.  Targeted (not full-rebuild) so that:
  #   (a) nodes tracked for non-time attributes do not clobber pre-existing
  #       branch-length edge$time values;
  #   (b) an explicit time=NA_real_ correctly clears the edge column entry.
  if (!is.null(time)) {
    time_vals <- setNames(time, name)
    affected  <- !is.na(graph$to) & graph$to %in% name
    if (any(affected)) {
      if (!"time" %in% names(graph)) graph$time <- NA_real_
      graph$time[affected] <- unname(time_vals[graph$to[affected]])
    }
  }
  if (!is.null(admix_event_time)) {
    aet_vals <- setNames(admix_event_time, name)
    affected  <- !is.na(graph$to) & graph$to %in% name
    if (any(affected)) {
      if (!"admix_event_time" %in% names(graph)) graph$admix_event_time <- NA_real_
      graph$admix_event_time[affected] <- unname(aet_vals[graph$to[affected]])
    }
  }
  graph
}

#' Return the union of topological leaves and internal-sampled nodes.
#'
#' Used by qpgraph and other consumers that need the full set of
#' populations with genotype data, regardless of whether they sit at
#' the topology's leaves. Orphan rows (nodes-tibble rows whose name has
#' no incident edge) are silently ignored.
#'
#' @param graph An edge tibble or igraph object.
#' @return Character vector of node names (sampled populations).
#' @keywords internal
get_sampled_nodes <- function(graph) {
  if (inherits(graph, "igraph")) {
    ig <- graph
  } else {
    ig <- edges_to_igraph(graph)
  }
  leaves <- get_leafnames(ig)
  if (!is.data.frame(graph)) return(leaves)
  nt <- graph_nodes(graph)
  canonical <- unique(c(graph$from, graph$to))
  in_set <- nt$name %in% canonical
  internal_sampled <- nt$name[in_set & !is.na(nt$samples)]
  union(leaves, internal_sampled)
}

#' Refresh the denormalized edge time views from the nodes tibble.
#'
#' Rebuilds edges$time (from nodes$time) and edges$admix_event_time
#' (from nodes$admix_event_time) for every edge whose destination node
#' appears in the nodes tibble. Edges whose destination is not recorded
#' in the nodes tibble are left unchanged. Columns are only created
#' when at least one node has a non-NA value.
#'
#' This is a full rebuild over every tracked node. `set_node_attrs` does
#' not call it; that function performs its own targeted edge-time write so
#' it can also clear an edge entry when a node time is set to NA. Call
#' `refresh_edge_times` after mutating the nodes-tibble columns directly,
#' before reading the edge views.
#'
#' @param graph An edge tibble.
#' @return The graph with edges$time and edges$admix_event_time refreshed.
#' @keywords internal
refresh_edge_times <- function(graph) {
  nt <- graph_nodes(graph)
  if (nrow(nt) == 0) return(graph)
  if ("time" %in% names(nt) && any(!is.na(nt$time))) {
    time_map  <- setNames(nt$time, nt$name)
    new_times <- time_map[graph$to]          # NA for untracked or NA-time nodes
    idx       <- !is.na(new_times)
    if (any(idx)) {
      if (!"time" %in% names(graph)) graph$time <- NA_real_
      graph$time[idx] <- unname(new_times[idx])
    }
  }
  if ("admix_event_time" %in% names(nt) && any(!is.na(nt$admix_event_time))) {
    aet_map  <- setNames(nt$admix_event_time, nt$name)
    new_aet  <- aet_map[graph$to]
    idx      <- !is.na(new_aet)
    if (any(idx)) {
      if (!"admix_event_time" %in% names(graph)) graph$admix_event_time <- NA_real_
      graph$admix_event_time[idx] <- unname(new_aet[idx])
    }
  }
  graph
}

#' Drop nodes-tibble rows whose name has no incident edge.
#'
#' Called once per emitted candidate by find_graphs and topology-
#' mutating helpers. Per-call cost O(N) where N is nrow(nodes_tbl).
#'
#' @param graph An edge tibble.
#' @return The graph with stale nodes rows removed.
#' @keywords internal
prune_nodes_attr <- function(graph) {
  nt <- attr(graph, "nodes")
  if (is.null(nt) || nrow(nt) == 0) return(graph)
  canonical <- unique(c(graph$from, graph$to))
  keep <- nt$name %in% canonical
  if (all(keep)) return(graph)
  attr(graph, "nodes") <- nt[keep, ]
  graph
}

#' Coerce to a plain edge tibble (drop the nodes attribute).
#'
#' Use this to opt out of the nodes-aware API and revert to the legacy
#' edge-only representation. Warns if the dropped nodes tibble carried
#' any non-NA data.
#'
#' @param graph An edge tibble (potentially with a nodes attribute).
#' @return The edge tibble with `attr(., "nodes")` cleared.
#' @export
as_edge_tibble <- function(graph) {
  nt <- attr(graph, "nodes")
  if (!is.null(nt) && nrow(nt) > 0) {
    data_cols <- setdiff(names(nt), "name")
    non_na_data <- length(data_cols) > 0 &&
      any(vapply(nt[data_cols], function(col) any(!is.na(col)), logical(1)))
    if (non_na_data)
      rlang::warn(
        "Dropping a non-empty nodes attribute.",
        class = "admixtools_dropping_node_attrs")
  }
  attr(graph, "nodes") <- NULL
  graph
}

#' Get fitted times as a named numeric vector.
#'
#' Backward-compatibility shim for callers that used the legacy
#' `attr(result, "node_times")` attribute set by `read_legofit_output`.
#' Resolves to `nodes$time` from the new nodes tibble.
#'
#' @param graph An edge tibble.
#' @return Named numeric vector keyed by node name.
#' @export
node_times <- function(graph) {
  nt <- graph_nodes(graph)
  times <- if ("time" %in% names(nt)) nt$time else rep(NA_real_, nrow(nt))
  setNames(times, nt$name)
}
