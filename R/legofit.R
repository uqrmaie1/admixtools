# dplyr verbs reference column names as bare symbols; declare them here so
# R CMD check does not warn about "no visible binding for global variable".
utils::globalVariables(c("from", "to", "type", "weight", "n",
                         "family", "parameter", "point_estimate", "lo", "hi",
                         "est", "low", "high", "name"))

#' Default drift-to-time conversion (identity)
#'
#' The default `drift_to_time` function for [graph_to_lgo()]. Returns
#' `drift_vec` unchanged (treats drift values directly as LEGOFIT time
#' starting values).
#'
#' Supply a custom function here when you want to convert drift values to
#' absolute generation times before writing the `.lgo`. The function must
#' accept two arguments: `drift_vec` (numeric) and `segment_vec` (character),
#' and return a numeric vector of the same length as `drift_vec`.
#'
#' @param drift_vec numeric vector of per-edge drift values.
#' @param segment_vec character vector of segment names (the `to` endpoint of
#'   each edge); ignored by this default but required by the contract so
#'   user-supplied alternatives can index per-segment.
#' @return Numeric vector, same length as `drift_vec`.
#' @export
default_drift_to_time <- function(drift_vec, segment_vec = NULL) {
  drift_vec
}

#' @keywords internal
coerce_to_edge_tibble <- function(x) {
  if (inherits(x, "igraph")) {
    edges <- igraph::as_edgelist(x) %>%
      tibble::as_tibble(.name_repair = ~ c("from", "to"))
    # Edge schema: `weight` (numeric), `type` (character), `time` (numeric).
    # `samples` and `twoN` are segment-level concepts (downstream they're
    # resolved against node names, not edges); not copied off the igraph here.
    edge_attr_types <- list(weight = "numeric", type = "character",
                            time = "numeric")
    for (attr in names(edge_attr_types)) {
      if (attr %in% igraph::edge_attr_names(x)) {
        val <- igraph::edge_attr(x, attr)
        if (length(val) != nrow(edges)) {
          rlang::abort(
            sprintf("igraph edge attribute `%s` has length %d but graph has %d edges.",
                    attr, length(val), nrow(edges)),
            class = "legofit_invalid_input"
          )
        }
        expected <- edge_attr_types[[attr]]
        type_ok <- switch(expected,
          numeric   = is.numeric(val),
          character = is.character(val)
        )
        if (!type_ok) {
          rlang::abort(
            sprintf("igraph edge attribute `%s` must be %s; got %s.",
                    attr, expected, typeof(val)),
            class = "legofit_invalid_input"
          )
        }
        edges[[attr]] <- val
      }
    }
    # Extract vertex attributes into the nodes tibble.
    # Vertex names are required; absent vertex names are an error.
    if (is.null(igraph::V(x)$name))
      rlang::abort(
        "igraph vertices have no `name` attribute; cannot build nodes tibble.",
        class = "admixtools_invalid_graph")
    va <- igraph::vertex_attr(x)
    schema_cols <- c("samples", "twoN_param", "twoN",
                     "time_param", "time", "admix_event_time")
    nt <- tibble::tibble(
      name             = va$name,
      samples          = NA_integer_,
      twoN_param       = NA_character_,
      twoN             = NA_real_,
      time_param       = NA_character_,
      time             = NA_real_,
      admix_event_time = NA_real_
    )
    for (col in schema_cols) {
      # Skip vertex-level time/admix_event_time when the same column was already
      # extracted from edge attributes above — they carry different semantics
      # (branch lengths vs absolute times) and would cause false-positive errors
      # from validate_edge_tibble(strict=TRUE).
      if ((col == "time" || col == "admix_event_time") && col %in% names(edges) &&
          any(!is.na(edges[[col]]))) next
      if (!is.null(va[[col]])) nt[[col]] <- va[[col]]
    }
    # Enforce the canonical schema types regardless of how igraph stored each
    # vertex attribute (e.g. an integer-valued `time` attr must not flip the
    # column from double to integer). `samples` is integer, the numeric columns
    # are double, and the *_param columns are character.
    nt$samples          <- as.integer(nt$samples)
    nt$twoN             <- as.double(nt$twoN)
    nt$time             <- as.double(nt$time)
    nt$admix_event_time <- as.double(nt$admix_event_time)
    nt$twoN_param       <- as.character(nt$twoN_param)
    nt$time_param       <- as.character(nt$time_param)
    # If R attr("nodes") is ALSO present, warn and prefer vertex_attr
    r_attr_nodes <- attr(x, "nodes")
    if (!is.null(r_attr_nodes))
      rlang::warn(
        "igraph input carries both vertex attrs and an R `nodes` attribute; using vertex attrs.",
        class = "admixtools_ambiguous_nodes_metadata")
    # Attach to the edges tibble we built earlier in this branch
    attr(edges, "nodes") <- nt
  } else if (is.data.frame(x)) {
    edges <- tibble::as_tibble(x)
  } else {
    rlang::abort(
      paste0("`graph` must be an igraph or edge tibble; got: ", class(x)[1]),
      class = "legofit_invalid_input"
    )
  }

  # Derive 'type' from in-degree if not present. Use "normal" / "admix"
  # to match the convention emitted by random_sim() and consumed by
  # msprime_sim() elsewhere in admixtools.
  if (!"type" %in% names(edges)) {
    edges <- edges %>%
      dplyr::add_count(to) %>%
      dplyr::mutate(type = ifelse(n == 1, "normal", "admix")) %>%
      dplyr::select(-n)
  }

  # Default weight when absent
  if (!"weight" %in% names(edges)) edges$weight <- 0

  edges
}

# Tolerance for the denormalized time-consistency check below. Intentionally
# tighter than topology_walk_bottom_up's 1e-6: this compares an exact
# denormalized copy of nodes$<col>, not an accumulated branch-length sum, so
# any drift beyond floating-point rounding is a genuine inconsistency.
TIME_CONSISTENCY_TOL <- 1e-10

# Verify that the denormalized edge view of a per-node time column agrees with
# the nodes tibble. `col` is "time" or "admix_event_time". A node carrying a
# non-NA value in nt[[col]] requires every incoming edge (rows of x whose
# `to` is that node) to carry the same value in x[[col]]; both a wrong value
# and a missing (NA) edge value are inconsistencies. The freshly-coerced
# igraph case (vertex time set, edge time all-NA) is the motivating example.
# Aborts with `admixtools_invalid_graph` on any violation, else returns
# invisibly.
#' @keywords internal
check_denormalized_time <- function(x, nt, col, tol = TIME_CONSISTENCY_TOL) {
  if (!col %in% names(nt) || all(is.na(nt[[col]]))) return(invisible(NULL))
  if (!col %in% names(x))
    rlang::abort(
      sprintf("nodes$%s is set but edges$%s column is absent; call refresh_edge_times().",
              col, col),
      class = "admixtools_invalid_graph")
  expected <- setNames(nt[[col]], nt$name)[x$to]
  edge_val <- x[[col]]
  # Construct `mismatch` so it is pure logical (never NA) for any NA, NaN, or
  # +/-Inf input: `!is.na(expected)` gates out untracked nodes (FALSE & NA =
  # FALSE in R), and the both-finite split avoids Inf - Inf = NaN poisoning
  # `any()`. Finite values must agree within `tol`; non-finite values must be
  # identical (so Inf == Inf is consistent, Inf vs 5 is not).
  both_finite <- is.finite(expected) & is.finite(edge_val)
  matches <- (both_finite & abs(edge_val - expected) <= tol) |
             (!both_finite & !is.na(expected) & !is.na(edge_val) &
              expected == edge_val)
  mismatch <- !is.na(expected) & !matches
  if (any(mismatch))
    rlang::abort(
      sprintf("edges$%s disagrees with nodes$%s; call refresh_edge_times().",
              col, col),
      class = "admixtools_invalid_graph")
  invisible(NULL)
}

#' @keywords internal
validate_edge_tibble <- function(x, strict = FALSE) {
  required <- c("from", "to", "type", "weight")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    rlang::abort(
      c("Missing required columns in edge tibble.",
        "x" = paste("Missing:", paste(missing_cols, collapse = ", ")),
        "i" = "Required columns: from, to, type, weight."),
      class = "legofit_invalid_input"
    )
  }
  type_ok <- c(
    from   = is.character(x$from),
    to     = is.character(x$to),
    type   = is.character(x$type),
    weight = is.numeric(x$weight)
  )
  if (!all(type_ok)) {
    rlang::abort(
      c("Edge tibble has invalid column types.",
        "x" = paste("Wrong type:", paste(names(type_ok)[!type_ok], collapse = ", ")),
        "i" = "from, to, type must be character; weight must be numeric."),
      class = "legofit_invalid_input"
    )
  }
  # Accept both "normal" (random_sim / msprime_sim convention) and "edge"
  # (qpgraph convention). Admixtools is internally inconsistent on this;
  # we treat any non-"admix" type as a derive-edge downstream.
  if (!all(x$type %in% c("normal", "edge", "admix"))) {
    rlang::abort(
      c("Invalid `type` values.",
        "x" = paste("Got:", paste(unique(x$type), collapse = ", ")),
        "i" = "Allowed: 'normal' or 'edge' for non-admix edges, 'admix' for admix edges."),
      class = "legofit_invalid_input"
    )
  }
  # admix invariant: each 'to' with type='admix' has exactly 2 incoming edges
  admix_groups <- table(x$to[x$type == "admix"])
  bad <- names(admix_groups)[admix_groups != 2]
  if (length(bad) > 0) {
    rlang::abort(
      c("admix-typed edges must come in pairs sharing `to`.",
        "x" = paste("Offending node(s):", paste(bad, collapse = ", ")),
        "i" = "Each admixture event needs exactly 2 incoming admix edges."),
      class = "legofit_invalid_input"
    )
  }
  # strict-mode checks against the nodes attribute
  nt <- attr(x, "nodes")
  if (!is.null(nt) && nrow(nt) > 0) {
    canonical <- unique(c(x$from, x$to))
    orphans <- setdiff(nt$name, canonical)
    if (length(orphans) > 0) {
      msg <- sprintf(
        "Nodes-tibble rows reference nodes not in the graph: %s",
        paste(orphans, collapse = ", "))
      if (strict) {
        rlang::abort(msg, class = "admixtools_invalid_graph")
      } else {
        rlang::warn(msg, class = "admixtools_orphan_nodes_pruned")
        attr(x, "nodes") <- nt[!nt$name %in% orphans, ]
        nt <- attr(x, "nodes")
      }
    }
    if (anyDuplicated(nt$name)) {
      msg <- sprintf(
        "Duplicate node names in nodes tibble: %s",
        paste(nt$name[duplicated(nt$name)], collapse = ", "))
      if (strict) {
        rlang::abort(msg, class = "admixtools_invalid_graph")
      } else {
        rlang::warn(msg, class = "admixtools_duplicate_node_name")
        attr(x, "nodes") <- nt[!duplicated(nt$name, fromLast = TRUE), ]
      }
    }
    # Strict mode: verify the denormalized edge views agree with the nodes
    # tibble. Production callers (graph_to_lgo and the two graph-alignment
    # paths) use strict = FALSE, where edges$time is a branch length rather
    # than a denormalized absolute node time and must not be checked here.
    if (strict) {
      check_denormalized_time(x, nt, "time")
      check_denormalized_time(x, nt, "admix_event_time")
    }
  }
  invisible(x)
}

#' @keywords internal
generate_param_names <- function(edges) {
  all_nodes <- unique(c(edges$from, edges$to))
  admix_dests <- unique(edges$to[edges$type == "admix"])
  # Collision guard: a node named admix_<X> where <X> is also an admix
  # destination would cause T_admix_<X> to be ambiguous between the per-event
  # admix time and the per-node time for node admix_<X>.
  # Mirrors the ghost-segment collision check in assemble_lgo (commit 50a88d3).
  ambiguous_nodes <- grep("^admix_(.+)$", all_nodes, value = TRUE, perl = TRUE)
  if (length(ambiguous_nodes) > 0) {
    suffix <- sub("^admix_(.+)$", "\\1", ambiguous_nodes, perl = TRUE)
    bad    <- suffix[suffix %in% admix_dests]
    if (length(bad) > 0) {
      rlang::abort(
        c("Node names collide with admix-event parameter naming.",
          "x" = paste("Offending nodes:",
                      paste(ambiguous_nodes[suffix %in% admix_dests],
                            collapse = ", ")),
          "i" = "T_admix_<X> would be ambiguous between event time and node time."),
        class = "legofit_invalid_input"
      )
    }
  }
  # Each admix event also gets its own "admix time" parameter — the time
  # at which the event happens. LEGOFIT requires both parent segments of
  # an admix to reference a SINGLE shared time parameter (not merely the
  # same value); per-admix-event params satisfy this constraint.
  list(
    time       = setNames(paste0("T_", all_nodes), all_nodes),
    twoN       = setNames(paste0("twoN_", all_nodes), all_nodes),
    mixFrac    = setNames(paste0("m_", admix_dests), admix_dests),
    admix_time = setNames(paste0("T_admix_", admix_dests), admix_dests)
  )
}

# Generate ghost-segment names for each admix parent. The LEGOFIT
# convention is to introduce a segment per admix-parent lineage at the
# admix event time. We name them <dest>_<original> so the round-trip
# parser can recognize them by structure.
#' @keywords internal
admix_ghost_name <- function(dest, parent) {
  paste0(dest, "_", parent)
}

#' @keywords internal
resolve_twoN_decls <- function(edges, twoN) {
  if (!is.data.frame(edges) || !all(c("from", "to") %in% names(edges))) {
    rlang::abort(
      "resolve_twoN_decls() requires an edge tibble with `from` and `to`.",
      class = "legofit_invalid_input"
    )
  }
  all_nodes <- unique(c(edges$from, edges$to))

  if (is.null(twoN)) {
    return(list(
      declarations = "twoN fixed one=1",
      per_segment  = setNames(rep("one", length(all_nodes)), all_nodes)
    ))
  }
  if (length(twoN) == 1 && is.numeric(twoN)) {
    return(list(
      declarations = sprintf("twoN fixed shared=%g", twoN),
      per_segment  = setNames(rep("shared", length(all_nodes)), all_nodes)
    ))
  }
  if (is.numeric(twoN) && !is.null(names(twoN))) {
    missing_nodes <- setdiff(all_nodes, names(twoN))
    if (length(missing_nodes) > 0) {
      rlang::abort(
        c("twoN named vector missing entries for some segments.",
          "x" = paste("Missing:", paste(missing_nodes, collapse = ", "))),
        class = "legofit_invalid_input"
      )
    }
    extra <- setdiff(names(twoN), all_nodes)
    if (length(extra) > 0) {
      rlang::abort(
        c("twoN named vector has entries that do not match any node.",
          "x" = paste("Unknown:", paste(extra, collapse = ", ")),
          "i" = paste("Known:", paste(all_nodes, collapse = ", "))),
        class = "legofit_invalid_input"
      )
    }
    return(list(
      declarations = paste(sprintf("twoN free twoN_%s=%g", names(twoN), twoN),
                           collapse = "\n"),
      per_segment  = setNames(paste0("twoN_", names(twoN)), names(twoN))
    ))
  }
  rlang::abort("twoN must be NULL, a scalar numeric, or a named numeric vector.",
               class = "legofit_invalid_input")
}

#' @keywords internal
strip_outgroup <- function(edges, outpop) {
  if (!outpop %in% edges$to) {
    rlang::abort(
      c("`outpop` not found among the graph's leaf nodes.",
        "x" = paste("Got:", outpop),
        "i" = paste("Available leaves:",
                    paste(setdiff(edges$to, edges$from), collapse = ", "))),
      class = "legofit_invalid_input"
    )
  }
  outpop_parent <- edges$from[edges$to == outpop]
  if (length(outpop_parent) != 1) {
    rlang::abort("`outpop` must have exactly one parent (root split).",
                 class = "legofit_invalid_input")
  }
  root <- setdiff(edges$from, edges$to)
  if (length(root) != 1 || outpop_parent != root) {
    rlang::abort(
      c("`outpop` is not adjacent to the topology root.",
        "i" = "Only direct outgroup leaves can be stripped via `outpop`."),
      class = "legofit_invalid_input"
    )
  }
  n_root_children <- sum(edges$from == root)
  if (n_root_children != 2) {
    rlang::abort(
      c("Cannot strip `outpop`: root must be strictly bifurcating.",
        "x" = sprintf("Root `%s` has %d children, not 2.",
                      root, n_root_children),
        "i" = "Re-root the graph before calling graph_to_lgo, or omit `outpop`."),
      class = "legofit_invalid_input"
    )
  }
  edges %>% dplyr::filter(to != outpop, from != root)
}

#' @keywords internal
resolve_scalar_or_named <- function(arg, nodes, default = NULL) {
  if (!is.character(nodes) || length(nodes) < 1) {
    rlang::abort(
      "resolve_scalar_or_named() requires a non-empty character vector of node names.",
      class = "legofit_invalid_input"
    )
  }
  if (is.null(arg)) {
    if (is.null(default)) {
      rlang::abort(
        "`arg` is NULL and no default supplied.",
        class = "legofit_invalid_input"
      )
    }
    return(setNames(rep(default, length(nodes)), nodes))
  }
  if (is.null(names(arg))) {
    if (length(arg) != 1) {
      rlang::abort(
        "Unnamed numeric must be length 1 (used as scalar default).",
        class = "legofit_invalid_input"
      )
    }
    return(setNames(rep(arg, length(nodes)), nodes))
  }
  # Named vector: reject unknown names (usually a typo'd leaf name).
  unknown <- setdiff(names(arg), nodes)
  if (length(unknown) > 0) {
    rlang::abort(
      c("Named numeric has entries that do not match any node.",
        "x" = paste("Unknown:", paste(unknown, collapse = ", ")),
        "i" = paste("Known:", paste(nodes, collapse = ", "))),
      class = "legofit_invalid_input"
    )
  }
  # Validate coverage
  missing_nodes <- setdiff(nodes, names(arg))
  if (length(missing_nodes) > 0) {
    if (is.null(default)) {
      rlang::abort(
        c("Named numeric missing entries (no default to fall back on).",
          "x" = paste("Missing:", paste(missing_nodes, collapse = ", "))),
        class = "legofit_invalid_input"
      )
    }
    out <- setNames(rep(default, length(nodes)), nodes)
  } else {
    out <- setNames(rep(if (!is.null(default)) default else NA_real_,
                        length(nodes)), nodes)
  }
  overlap <- intersect(nodes, names(arg))
  out[overlap] <- arg[overlap]
  out
}

#' @keywords internal
topology_walk_bottom_up <- function(edges, branch_lengths, leaf_times,
                                    tol = 1e-6) {
  if (length(branch_lengths) != nrow(edges)) {
    rlang::abort(
      sprintf("branch_lengths length (%d) must equal nrow(edges) (%d).",
              length(branch_lengths), nrow(edges)),
      class = "legofit_invalid_input"
    )
  }
  all_nodes <- unique(c(edges$from, edges$to))
  node_time <- setNames(rep(NA_real_, length(all_nodes)), all_nodes)

  leaves <- setdiff(edges$to, edges$from)
  if (!setequal(names(leaf_times), leaves)) {
    rlang::abort(
      "`leaf_times` must name exactly the leaf nodes.",
      class = "legofit_invalid_input"
    )
  }
  node_time[leaves] <- leaf_times[leaves]

  progress <- TRUE
  while (progress) {
    progress <- FALSE
    for (i in seq_len(nrow(edges))) {
      u <- edges$from[[i]]
      v <- edges$to[[i]]
      if (is.na(node_time[u]) && !is.na(node_time[v])) {
        node_time[u] <- node_time[v] + branch_lengths[[i]]
        progress <- TRUE
      } else if (!is.na(node_time[u]) && !is.na(node_time[v])) {
        proposed <- node_time[v] + branch_lengths[[i]]
        if (abs(node_time[u] - proposed) > tol) {
          rlang::abort(
            c(sprintf("Inconsistent branch lengths at node `%s`.", u),
              "x" = sprintf("Computed t=%g via one path and t=%g via edge %s->%s.",
                            node_time[u], proposed, u, v),
              "i" = "Input drifts imply inconsistent absolute times. This is common: random and real qpgraph drifts are rarely additively consistent across paths.",
              "i" = "Use `time_handling = \"free\"` (depth-seeded free times; the LEGOFIT optimizer fits them, but the input absolute times are not preserved), or `time_handling = \"init\"` with an explicit `time` column to supply consistent times yourself."),
            class = "legofit_invalid_input"
          )
        }
      }
    }
  }

  if (any(is.na(node_time))) {
    unreachable <- names(node_time)[is.na(node_time)]
    rlang::abort(
      c("Some nodes unreachable from leaves during topology walk.",
        "x" = paste("Unreachable:", paste(unreachable, collapse = ", ")),
        "i" = "Input graph may be disconnected, may contain cycles, or the leaf set may be incomplete."),
      class = "legofit_invalid_input"
    )
  }
  node_time
}

#' Longest-path depth from any leaf to each node.
#'
#' Used to seed starting values for `time_handling = "free"`. LEGOFIT
#' requires `free` declarations to carry a starting value; depth gives
#' a monotonic deterministic series (parent > child) without requiring
#' consistent drift values across the topology.
#' @keywords internal
compute_node_depths <- function(edges) {
  all_nodes <- unique(c(edges$from, edges$to))
  leaves <- setdiff(edges$to, edges$from)
  depth <- setNames(rep(NA_integer_, length(all_nodes)), all_nodes)
  depth[leaves] <- 0L
  # Cap iterations at V: in any acyclic graph the longest path is at
  # most V-1, so depth converges within V passes. A cycle would let
  # the longest-path algorithm grow without bound, so cap and error.
  max_iter <- length(all_nodes)
  progress <- TRUE
  iter <- 0L
  while (progress && iter < max_iter) {
    progress <- FALSE
    iter <- iter + 1L
    for (i in seq_len(nrow(edges))) {
      u <- edges$from[[i]]
      v <- edges$to[[i]]
      if (!is.na(depth[v])) {
        candidate <- depth[v] + 1L
        if (is.na(depth[u]) || candidate > depth[u]) {
          depth[u] <- candidate
          progress <- TRUE
        }
      }
    }
  }
  if (progress && iter >= max_iter) {
    rlang::abort(
      c("Cycle detected during depth computation.",
        "x" = sprintf("Longest path did not converge in %d iterations.", max_iter),
        "i" = "Input graph must be a DAG; cycles are not supported."),
      class = "legofit_invalid_input"
    )
  }
  if (any(is.na(depth))) {
    unreachable <- names(depth)[is.na(depth)]
    rlang::abort(
      c("Some nodes are unreachable from any leaf during depth computation.",
        "x" = paste("Unreachable:", paste(unreachable, collapse = ", ")),
        "i" = "Input graph may contain a disconnected component."),
      class = "legofit_invalid_input"
    )
  }
  depth
}

#' @keywords internal
compute_times <- function(edges, time_handling, drift_to_time,
                          dates_terminal, fix_times = FALSE) {
  required_cols <- c("from", "to", "type", "weight")
  if (!is.data.frame(edges) || !all(required_cols %in% names(edges))) {
    rlang::abort(
      "compute_times() requires an edge tibble with from/to/type/weight columns.",
      class = "legofit_invalid_input"
    )
  }
  if (!time_handling %in% c("fix_admix", "init", "free")) {
    rlang::abort(
      sprintf("Invalid `time_handling`: %s (allowed: fix_admix, init, free).",
              time_handling),
      class = "legofit_invalid_input"
    )
  }
  if (!is.function(drift_to_time)) {
    rlang::abort(
      "`drift_to_time` must be a function.",
      class = "legofit_invalid_input"
    )
  }
  all_nodes <- unique(c(edges$from, edges$to))
  leaves    <- setdiff(edges$to, edges$from)

  if (length(leaves) == 0) {
    rlang::abort(
      c("Graph has no leaves.",
        "x" = "Every node appears as both `from` and `to`.",
        "i" = "Likely a cyclic input or malformed edges tibble; a finite DAG always has at least one leaf."),
      class = "legofit_invalid_input"
    )
  }

  admix_dests <- unique(edges$to[edges$type == "admix"])
  terminal_t  <- resolve_scalar_or_named(dates_terminal, leaves, default = 0)

  if (time_handling == "free") {
    # LEGOFIT requires `free` declarations to carry a starting value;
    # emit topological-depth-based values so parent > child holds.
    # The optimizer will adjust from there.
    non_leaf <- setdiff(all_nodes, leaves)
    depths   <- compute_node_depths(edges)
    init_t   <- setNames(as.numeric(depths[non_leaf]), non_leaf)
    free_vec <- setNames(!(all_nodes %in% leaves), all_nodes)
    if (isTRUE(fix_times)) free_vec[] <- FALSE   # fix_times: pin internal times
    return(list(
      value = c(terminal_t, init_t)[all_nodes],
      free  = free_vec,
      omit  = character(0)
    ))
  }

  # "init" or "fix_admix": derive branch lengths, then topology-walk
  bl <- if ("time" %in% names(edges)) {
    ifelse(!is.na(edges$time),
           edges$time,
           drift_to_time(edges$weight, edges$to))
  } else {
    drift_to_time(edges$weight, edges$to)
  }
  if (time_handling == "fix_admix") {
    bl <- ifelse(edges$type == "admix", 0, bl)
  }

  node_time <- topology_walk_bottom_up(edges, bl, leaf_times = terminal_t)

  # Path B (ghost segments) means admix dests get their own t= field
  # like any other node — they aren't omitted from time declarations.
  # The LEGOFIT param-sharing constraint is handled by the ghost
  # segments anchoring at T_admix_<dest> instead.
  free <- setNames(!(all_nodes %in% leaves), all_nodes)
  # fix_times: emit internal node times as `time fixed` at their computed
  # absolute values, so a model with free `twoN` is no longer degenerate up to a
  # global scale and absolute twoN becomes recoverable.
  if (isTRUE(fix_times)) free[] <- FALSE

  list(value = node_time[all_nodes], free = free, omit = character(0))
}

# Single source of truth for (high, low) parent-pair ordering in admix events.
# Sort: desc(weight) primary, ascending `from` as tiebreaker.
admix_parent_pairs <- function(edges) {
  edges %>%
    dplyr::filter(type == "admix") %>%
    dplyr::group_by(to) %>%
    dplyr::arrange(desc(weight), from, .by_group = TRUE) %>%
    dplyr::summarise(
      high          = dplyr::first(from),
      low           = dplyr::last(from),
      mixfrac_value = dplyr::last(weight),
      .groups       = "drop"
    )
}

format_time_declarations <- function(time_param_names, times,
                                     admix_time_param_names = NULL,
                                     admix_pairs = NULL,
                                     fix_times = FALSE) {
  # Per-node time declarations (every node gets one — admix dests are no
  # longer omitted; ghost segments handle the param-sharing constraint).
  keep <- names(time_param_names)
  fixed   <- keep[!times$free[keep]]
  freed   <- keep[ times$free[keep]]
  ordered <- c(fixed, freed)
  free_str <- ifelse(times$free[ordered], "free", "fixed")
  node_lines <- sprintf("time %s %s=%g",
                        free_str,
                        time_param_names[ordered],
                        times$value[ordered])
  # Captured source params (round-trip) may be shared across nodes; emit one
  # declaration per distinct param name. No-op for generated `T_<node>` names,
  # which are unique per node.
  node_lines <- node_lines[!duplicated(time_param_names[ordered])]

  # Per-admix-event time declarations (T_admix_<M>). Free by default (the admix
  # event time is what LEGOFIT optimizes), but fixed under `fix_times`: leaving it
  # free while twoN is free reintroduces the scale degeneracy and wrecks absolute
  # twoN recovery. The admix-event time is structurally non-identifiable either
  # way (its segment carries a single lineage, so no coalescence), so pinning it
  # is safe.
  admix_lines <- character(0)
  if (!is.null(admix_pairs) && nrow(admix_pairs) > 0) {
    admix_lines <- sprintf("time %s %s=%g",
                           if (isTRUE(fix_times)) "fixed" else "free",
                           admix_time_param_names[admix_pairs$to],
                           times$value[admix_pairs$high])
  }

  out <- c(node_lines, admix_lines)
  if (length(out) == 0) return(character(0))
  paste(out, collapse = "\n")
}

format_mixfrac_declarations <- function(mixfrac_param_names, edges) {
  pairs <- admix_parent_pairs(edges)
  if (nrow(pairs) == 0) return(character(0))
  paste(sprintf("mixFrac free %s=%g",
                mixfrac_param_names[pairs$to],
                pairs$mixfrac_value),
        collapse = "\n")
}

format_segment_declarations <- function(nodes, time_param_names,
                                        twoN_param_per_node, samples_vec) {
  # Every segment carries an explicit `t=<param>`. (Ghost segments for
  # admix events are emitted separately by format_ghost_segments.)
  # The t= and twoN= tokens come from time_param_names / twoN_param_per_node.
  # graph_to_lgo overrides those maps (and the matching declarations) with the
  # source-file tokens when round tripping a graph that carries a nodes tibble,
  # so segments and the Parameters block always reference the same declared
  # names.
  lines <- vapply(nodes, function(n) {
    parts <- paste0("segment ", n, " t=", time_param_names[n],
                    " twoN=", twoN_param_per_node[n])
    if (!is.na(samples_vec[n])) parts <- paste0(parts, " samples=", samples_vec[n])
    parts
  }, character(1))
  paste(lines, collapse = "\n")
}

# Emit the "ghost" segments LEGOFIT requires for each admix event. Both
# parents of an admix must reference the SAME time parameter; we anchor
# them at T_admix_<dest>. The ghost is named <dest>_<parent>.
format_ghost_segments <- function(admix_pairs, admix_time_param_names,
                                  twoN_param_per_node) {
  if (nrow(admix_pairs) == 0) return(character(0))
  lines <- character(0)
  for (i in seq_len(nrow(admix_pairs))) {
    m    <- admix_pairs$to[[i]]
    high <- admix_pairs$high[[i]]
    low  <- admix_pairs$low[[i]]
    anchor <- admix_time_param_names[[m]]
    lines <- c(lines,
      sprintf("segment %s t=%s twoN=%s",
              admix_ghost_name(m, high), anchor, twoN_param_per_node[[high]]),
      sprintf("segment %s t=%s twoN=%s",
              admix_ghost_name(m, low),  anchor, twoN_param_per_node[[low]])
    )
  }
  paste(lines, collapse = "\n")
}

# Each ghost derives from its corresponding real parent population.
format_ghost_derives <- function(admix_pairs) {
  if (nrow(admix_pairs) == 0) return(character(0))
  lines <- character(0)
  for (i in seq_len(nrow(admix_pairs))) {
    m    <- admix_pairs$to[[i]]
    high <- admix_pairs$high[[i]]
    low  <- admix_pairs$low[[i]]
    lines <- c(lines,
      sprintf("derive %s from %s", admix_ghost_name(m, high), high),
      sprintf("derive %s from %s", admix_ghost_name(m, low),  low)
    )
  }
  paste(lines, collapse = "\n")
}

format_derive_statements <- function(edges_non_admix) {
  paste(sprintf("derive %s from %s",
                edges_non_admix$to,
                edges_non_admix$from),
        collapse = "\n")
}

format_mix_statements <- function(edges_admix, mixfrac_param_names) {
  pairs <- admix_parent_pairs(edges_admix)
  if (nrow(pairs) == 0) return(character(0))
  # Mix references the ghost segments (named <dest>_<parent>), not the
  # real parents, because LEGOFIT requires both mix parents to share a
  # time parameter and the ghosts anchor at T_admix_<dest>.
  paste(sprintf("mix %s from %s + %s * %s",
                pairs$to,
                admix_ghost_name(pairs$to, pairs$high),
                mixfrac_param_names[pairs$to],
                admix_ghost_name(pairs$to, pairs$low)),
        collapse = "\n")
}

#' @keywords internal
assemble_lgo <- function(edges, params, times, twoN_decls, samples_vec,
                         fix_times = FALSE) {
  all_nodes   <- unique(c(edges$from, edges$to))
  admix_pairs <- admix_parent_pairs(edges)

  # Guard against user-supplied node names colliding with the synthetic
  # ghost-segment names we will emit (<dest>_<parent>). A collision would
  # create duplicate `segment` declarations and break LEGOFIT parse.
  if (nrow(admix_pairs) > 0) {
    proposed_ghosts <- c(
      admix_ghost_name(admix_pairs$to, admix_pairs$high),
      admix_ghost_name(admix_pairs$to, admix_pairs$low)
    )
    collisions <- intersect(proposed_ghosts, all_nodes)
    if (length(collisions) > 0) {
      rlang::abort(
        c("Node names collide with synthetic ghost-segment names.",
          "x" = paste("Collisions:", paste(collisions, collapse = ", ")),
          "i" = "Ghost names follow <dest>_<parent>; rename the offending nodes."),
        class = "legofit_invalid_input"
      )
    }
  }

  time_decls   <- format_time_declarations(params$time, times,
                                           params$admix_time, admix_pairs,
                                           fix_times = fix_times)
  mfrac_decls  <- format_mixfrac_declarations(params$mixFrac, edges)
  seg_decls    <- format_segment_declarations(all_nodes, params$time,
                                              twoN_decls$per_segment,
                                              samples_vec)
  ghost_segs   <- format_ghost_segments(admix_pairs, params$admix_time,
                                        twoN_decls$per_segment)
  derive_stmts <- format_derive_statements(
    edges %>% dplyr::filter(type != "admix")
  )
  ghost_derives <- format_ghost_derives(admix_pairs)
  mix_stmts     <- format_mix_statements(
    edges %>% dplyr::filter(type == "admix"),
    params$mixFrac
  )

  # Combine real and ghost segment / derive blocks, dropping empty pieces
  # so we don't emit blank inner lines.
  seg_block    <- paste(c(seg_decls, ghost_segs)[nzchar(c(seg_decls, ghost_segs))],
                        collapse = "\n")
  derive_block <- paste(c(derive_stmts, ghost_derives)[nzchar(c(derive_stmts, ghost_derives))],
                        collapse = "\n")

  paste(
    "## .lgo generated by admixtools::graph_to_lgo()",
    "",
    "### Parameters",
    time_decls,
    twoN_decls$declarations,
    mfrac_decls,
    "",
    "### Segments",
    seg_block,
    "",
    "### Relationships",
    derive_block,
    mix_stmts,
    sep = "\n"
  )
}

# Unified regex for time/twoN/mixFrac/param declarations.
# Optional bounds group (non-capturing) captures lo (group 3) and hi (group 4).
# Name rule per parse.c:248-254: letter followed by letters, digits,
# underscores, colons, or periods.
# Must be used with perl = TRUE.
PARAM_RE <- paste0(
  "^\\s*",
  "(time|twoN|mixFrac|param)",            # G1: type keyword
  "\\s+",
  "(fixed|free|constrained)",             # G2: subtype
  "\\s+",
  "(?:\\[\\s*([^,\\]]+?)\\s*,",          # G3: lo (optional; non-capturing outer)
  "\\s*([^\\]]+?)\\s*\\]\\s+)?",          # G4: hi (optional; non-capturing outer)
  "([A-Za-z][A-Za-z0-9_.:]*)",            # G5: param name
  "\\s*=\\s*",
  "(.*?)",                                # G6: value / constrained expression
  "\\s*$"
)

# Pre-process lines: join a line ending with one of `+ - * /` (per
# parse.c:643-668 get_one_line continuation logic) with the next line.
# Comment stripping must have already been applied — this function sees
# post-strip lines. Trailing whitespace is trimmed defensively (matching
# LEGOFIT's get_one_line whitespace-strip before the suffix check).
#' @keywords internal
join_continuations <- function(lines) {
  lines  <- trimws(lines, which = "right")
  joined <- character(0)
  buffer <- ""
  for (ln in lines) {
    buffer <- if (nzchar(buffer)) paste(buffer, ln) else ln
    buffer <- trimws(buffer, which = "right")
    last_char <- substring(buffer, nchar(buffer))
    if (last_char %in% c("+", "-", "*", "/")) next
    joined <- c(joined, buffer)
    buffer <- ""
  }
  if (nzchar(buffer)) {
    rlang::abort(
      sprintf("Unexpected end of input after continuation line: %s", buffer),
      class = "legofit_lgo_unsupported"
    )
  }
  joined
}

# Safely evaluate an arithmetic expression from a `time constrained` RHS.
# Only allows numeric literals, identifier lookups against `env`, and the
# operators + - * / ^ (). Anything else aborts with
# `legofit_lgo_unsupported`. Avoids the security risk of plain
# `eval(parse(text=...))` on third-party .lgo files.
#' @keywords internal
safe_eval_arith <- function(expr, env) {
  if (is.numeric(expr)) return(as.numeric(expr))
  if (is.name(expr)) {
    name <- as.character(expr)
    if (!exists(name, envir = env, inherits = FALSE)) {
      rlang::abort(
        sprintf("Unknown identifier `%s` in constrained expression.", name),
        class = "legofit_lgo_unsupported"
      )
    }
    return(get(name, envir = env, inherits = FALSE))
  }
  if (is.call(expr)) {
    op <- as.character(expr[[1]])
    if (!op %in% c("+", "-", "*", "/", "^", "(")) {
      rlang::abort(
        sprintf("Unsafe operator `%s` in constrained expression.", op),
        class = "legofit_lgo_unsupported"
      )
    }
    args <- lapply(as.list(expr)[-1], safe_eval_arith, env = env)
    return(do.call(op, args))
  }
  rlang::abort(
    "Unsupported expression syntax in constrained declaration.",
    class = "legofit_lgo_unsupported"
  )
}

#' Parse a LEGOFIT .lgo file
#'
#' Parses the LEGOFIT 1.87 `.lgo` grammar produced by [graph_to_lgo()] and
#' most third-party `.lgo` files. Supported grammar: line continuation
#' (`+`, `-`, `*`, `/` trailing operators per parse.c:643-668), bounded
#' `free` declarations (`[lo, hi]` before name per parse.c:211-235),
#' `constrained` arithmetic expressions (safe evaluator), and both `"edges"`
#' and `"igraph"` return types.
#'
#' Ghost segments written by [graph_to_lgo()] (named `<dest>_<parent>`) are
#' recognised by the narrow 3-condition rule and collapsed back into
#' admixtools admix edges. Third-party segments that happen to appear as mix
#' parents but do NOT match the `<dest>_<parent>` naming convention (e.g.,
#' Rogers 2020 rha20.lgo's `d2`, `s2`) are preserved as real edges.
#'
#' Exactly one of `path` or `text` must be supplied.
#'
#' @param path Path to a `.lgo` file (UTF-8 encoded).
#' @param text A character scalar (full file text) or character vector
#'   (one element per line).
#' @param as Output format: `"edges"` (default) returns a tibble with columns
#'   `from`, `to`, `type`, `weight`; `"igraph"` returns an igraph with `type`
#'   and `weight` edge attributes.  Both carry a `param_bounds` attribute
#'   (list keyed by param name, value `c(lo, hi)`) when bounded `free`
#'   declarations are present.
#' @return A tibble (or igraph) representing the graph topology.
#' @examples
#' lgo <- "
#' time fixed T_x = 0
#' time fixed T_y = 0
#' time free  T_R = 2
#' twoN fixed one = 1
#' segment R t=T_R twoN=one
#' segment x t=T_x twoN=one samples=1
#' segment y t=T_y twoN=one samples=1
#' derive x from R
#' derive y from R
#' "
#' read_lgo(text = lgo)
#' @export
read_lgo <- function(path = NULL, text = NULL, as = c("edges", "igraph")) {
  as <- match.arg(as)

  if (is.null(path) == is.null(text)) {
    rlang::abort(
      "Exactly one of `path` or `text` must be supplied.",
      class = "legofit_invalid_input"
    )
  }

  lines <- if (!is.null(path)) {
    readLines(path, encoding = "UTF-8")
  } else if (is.character(text) && length(text) == 1) {
    strsplit(text, "\n", fixed = TRUE)[[1]]
  } else if (is.character(text)) {
    text
  } else {
    rlang::abort(
      "`text` must be a character scalar or vector.",
      class = "legofit_invalid_input"
    )
  }

  # Strip comments (#...) and trailing whitespace; skip blank lines.
  # Comment stripping is per-physical-line, BEFORE continuation joining
  # (matching parse.c:590 get_one_line behavior).
  lines <- gsub("#.*$", "", lines)
  lines <- trimws(lines, which = "right")
  lines <- lines[nzchar(lines)]

  # Join lines ending with a trailing binary operator (parse.c:643-668).
  lines <- join_continuations(lines)

  # Accumulators
  params_rows   <- list()
  derive_rows   <- list()
  mix_rows      <- list()
  param_bounds  <- list()   # keyed by param name; value = c(lo, hi)
  segments_meta <- list()   # name -> list(samples, twoN_param, time_param)

  for (ln in lines) {
    tok <- strsplit(trimws(ln), "\\s+")[[1]]

    if (tok[[1]] %in% c("time", "twoN", "mixFrac", "param")) {
      # Use PARAM_RE to handle all subtype variants uniformly, including the
      # bounded-free form `{type} free [lo, hi] name = value` (per parse.c:211-235).
      m     <- regexec(PARAM_RE, ln, perl = TRUE)
      parts <- regmatches(ln, m)[[1]]
      if (length(parts) == 0) {
        rlang::abort(
          sprintf("Malformed %s declaration: %s", tok[[1]], ln),
          class = "legofit_lgo_unsupported"
        )
      }
      # parts indices (1 = full match, 2..7 = groups 1..6):
      p_type    <- parts[[2]]   # "time" | "twoN" | "mixFrac" | "param"
      p_subtype <- parts[[3]]   # "fixed" | "free" | "constrained"
      lo_str    <- parts[[4]]   # lo bound (or "" if no brackets)
      hi_str    <- parts[[5]]   # hi bound (or "" if no brackets)
      p_name    <- parts[[6]]   # parameter name
      rhs       <- parts[[7]]   # value or constrained expression

      # Bounds validation: brackets only valid on `free`.
      has_bounds <- nzchar(lo_str)
      if (has_bounds) {
        if (p_subtype != "free") {
          rlang::abort(
            sprintf(
              "Bounds `[lo, hi]` are only valid on `free` declarations; got `%s` at: %s",
              p_subtype, ln),
            class = "legofit_lgo_unsupported"
          )
        }
        lo_val <- suppressWarnings(as.numeric(lo_str))
        hi_val <- suppressWarnings(as.numeric(hi_str))
        if (is.na(lo_val) || is.na(hi_val)) {
          rlang::abort(
            sprintf("Non-numeric bounds in declaration: %s", ln),
            class = "legofit_lgo_unsupported"
          )
        }
        param_bounds[[p_name]] <- c(lo = lo_val, hi = hi_val)
      }

      p_val <- if (p_subtype == "constrained") {
        # `constrained` RHS is an arithmetic expression over previously
        # declared parameters. Evaluate safely (only +, -, *, /, ^, and
        # parenthesization allowed; identifiers must already be in the
        # params parsed so far).
        env <- list2env(
          setNames(lapply(params_rows, function(r) r$value),
                   vapply(params_rows, function(r) r$name, character(1))),
          parent = emptyenv()
        )
        # Parse the RHS expression in two steps so that (a) a syntax error in
        # parse() yields a clear message, and (b) an empty RHS (bare `=` with
        # nothing after it) produces a targeted diagnostic rather than the
        # confusing "subscript out of bounds" from parse(text="")[[1]].
        parsed <- tryCatch(
          parse(text = rhs),
          error = function(e) rlang::abort(
            sprintf("Syntax error in constrained expression for '%s': %s", p_name, rhs),
            class = "legofit_lgo_unsupported", parent = e)
        )
        if (length(parsed) == 0L)
          rlang::abort(
            sprintf("Empty right-hand side in constrained declaration: '%s =' has no value.", p_name),
            class = "legofit_lgo_unsupported")
        safe_eval_arith(parsed[[1]], env)
      } else {
        suppressWarnings(as.numeric(rhs))
      }
      params_rows[[length(params_rows) + 1]] <- list(
        name = p_name, value = p_val, type = p_subtype
      )

    } else if (tok[[1]] == "segment") {
      # segment <name> [t=<param>] twoN=<param> [samples=<int>]
      # Capture per-segment metadata (samples, twoN token, t token) into
      # segments_meta so it can be attached as the graph's nodes tibble.
      # Tolerant of whitespace around `=` and of keyword casing, matching
      # the time/twoN/mixFrac parser above. `samples=0` is not a sample.
      seg_name <- tok[[2]]
      if (length(tok) > 2) {
        rest <- paste(tok[-(1:2)], collapse = " ")
        meta <- list()
        m_s <- regmatches(rest, regexec(
          "(?:^|\\s)samples\\s*=\\s*([0-9]+)", rest,
          ignore.case = TRUE, perl = TRUE))[[1]]
        if (length(m_s) == 2) {
          n <- as.integer(m_s[[2]])
          if (!is.na(n) && n > 0) meta$samples <- n
        }
        m_tn <- regmatches(rest, regexec(
          "(?:^|\\s)twoN\\s*=\\s*([^\\s]+)", rest,
          ignore.case = TRUE, perl = TRUE))[[1]]
        if (length(m_tn) == 2) meta$twoN_param <- m_tn[[2]]
        m_t <- regmatches(rest, regexec(
          "(?:^|\\s)t\\s*=\\s*([^\\s]+)", rest,
          ignore.case = TRUE, perl = TRUE))[[1]]
        if (length(m_t) == 2) meta$time_param <- m_t[[2]]
        if (length(meta) > 0) segments_meta[[seg_name]] <- meta
      }
      next

    } else if (tok[[1]] == "derive") {
      # derive <child> from <parent>
      if (length(tok) != 4 || tok[[3]] != "from") {
        rlang::abort(
          sprintf("Malformed derive statement: %s", ln),
          class = "legofit_lgo_unsupported"
        )
      }
      derive_rows[[length(derive_rows) + 1]] <- list(child = tok[[2]], parent = tok[[4]])

    } else if (tok[[1]] == "mix") {
      # mix <child> from <high> + <mixfrac_param> * <low>
      # tok: mix, child, from, high, +, mixfrac_param, *, low
      if (length(tok) != 8 || tok[[3]] != "from" || tok[[5]] != "+" || tok[[7]] != "*") {
        rlang::abort(
          sprintf("Malformed mix statement: %s", ln),
          class = "legofit_lgo_unsupported"
        )
      }
      mix_rows[[length(mix_rows) + 1]] <- list(
        child       = tok[[2]],
        parent_high = tok[[4]],
        mf_param    = tok[[6]],
        parent_low  = tok[[8]]
      )

    } else {
      rlang::abort(
        sprintf("Unsupported `.lgo` construct: %s", ln),
        class = "legofit_lgo_unsupported"
      )
    }
  }

  # Build params lookup: name -> value
  param_val <- if (length(params_rows) > 0) {
    setNames(
      vapply(params_rows, function(r) r$value, numeric(1)),
      vapply(params_rows, function(r) r$name,  character(1))
    )
  } else {
    numeric(0)
  }

  # name -> declaration kind ("fixed"/"free"/"constrained"). Lets the writer
  # re-emit a captured param with its original fixed/free type when round
  # tripping; `constrained` is emitted as free at its evaluated value.
  param_type <- if (length(params_rows) > 0) {
    setNames(
      vapply(params_rows, function(r) r$type, character(1)),
      vapply(params_rows, function(r) r$name, character(1))
    )
  } else {
    character(0)
  }

  # Narrow ghost detection: a segment is a ghost iff
  # (a) name matches <dest>_<parent> for some admix dest in this file,
  # (b) it appears as the child of `derive <ghost> from <real>`, AND
  # (c) it appears as a mix parent.
  # The old greedy rule treated ANY mix parent as a ghost, which
  # collapsed real third-party intermediate segments (rha20.lgo's d2,
  # s, n, a2, nd2, y2). The narrow rule preserves them while still
  # recognising our own writer's <dest>_<parent> ghosts.
  admix_dests_seen <- unique(vapply(mix_rows, `[[`, character(1), "child"))
  derive_children  <- if (length(derive_rows) > 0) {
    vapply(derive_rows, `[[`, character(1), "child")
  } else character(0)

  ghost_names <- character(0)
  if (length(admix_dests_seen) > 0) {
    for (mr in mix_rows) {
      for (cand in c(mr$parent_high, mr$parent_low)) {
        # Condition (a): cand is "<dest>_<suffix>" for some admix dest in
        # this file, with a non-empty suffix. Tested by fixed-string prefix
        # match (not a regex) so node names containing regex metacharacters
        # (e.g. `+`, `[`, `(`) are matched literally and cannot corrupt or
        # crash the parse.
        matches_a <- any(vapply(admix_dests_seen, function(d) {
          prefix <- paste0(d, "_")
          startsWith(cand, prefix) && nchar(cand) > nchar(prefix)
        }, logical(1)))
        if (!matches_a) next
        # Condition (b): cand is a derive child (real parent in derive)
        i <- which(derive_children == cand)
        if (length(i) != 1) next
        # Condition (c): inherent — cand appeared as a mix parent here
        ghost_names <- c(ghost_names, cand)
      }
    }
    ghost_names <- unique(ghost_names)
  }

  ghost_to_real    <- list()
  real_derive_rows <- list()
  for (dr in derive_rows) {
    if (dr$child %in% ghost_names) {
      ghost_to_real[[dr$child]] <- dr$parent
    } else {
      real_derive_rows[[length(real_derive_rows) + 1]] <- dr
    }
  }

  # Build edge tibble from real-derive + mix rows
  edge_rows <- list()
  for (dr in real_derive_rows) {
    # The .lgo format does not carry edge drift directly; weight is unknown
    # on round-trip. NA_real_ is the honest sentinel.
    edge_rows[[length(edge_rows) + 1]] <- list(
      from = dr$parent, to = dr$child, type = "normal", weight = NA_real_
    )
  }

  for (mr in mix_rows) {
    mf_val <- if (mr$mf_param %in% names(param_val)) param_val[[mr$mf_param]] else NA_real_
    # Resolve ghost names back to real parents (graph_to_lgo ghosts are named
    # <dest>_<real>; resolve via the derive that produced them).
    real_high <- if (mr$parent_high %in% names(ghost_to_real)) {
      ghost_to_real[[mr$parent_high]]
    } else {
      mr$parent_high
    }
    real_low <- if (mr$parent_low %in% names(ghost_to_real)) {
      ghost_to_real[[mr$parent_low]]
    } else {
      mr$parent_low
    }
    edge_rows[[length(edge_rows) + 1]] <- list(
      from   = real_high,
      to     = mr$child,
      type   = "admix",
      weight = if (!is.na(mf_val)) 1 - mf_val else NA_real_
    )
    edge_rows[[length(edge_rows) + 1]] <- list(
      from   = real_low,
      to     = mr$child,
      type   = "admix",
      weight = if (!is.na(mf_val)) mf_val else NA_real_
    )
  }

  if (length(edge_rows) == 0) {
    if (as == "igraph") {
      ig <- igraph::make_empty_graph(n = 0, directed = TRUE)
      if (length(param_bounds) > 0)
        ig <- igraph::set_graph_attr(ig, "param_bounds", param_bounds)
      return(ig)
    }
    result <- tibble::tibble(from = character(), to = character(),
                             type = character(), weight = numeric())
    if (length(param_bounds) > 0) attr(result, "param_bounds") <- param_bounds
    return(result)
  }

  result <- do.call(rbind, lapply(edge_rows, as.data.frame, stringsAsFactors = FALSE))
  result <- tibble::as_tibble(result)

  # Attach per-segment metadata as the graph's nodes tibble. Exclude
  # synthetic ghosts and any segment not present as a real node in the
  # final edge tibble.
  if (length(segments_meta) > 0) {
    real_nodes <- unique(c(result$from, result$to))
    keep <- setdiff(intersect(names(segments_meta), real_nodes), ghost_names)
    if (length(keep) > 0) {
      nt <- tibble::tibble(
        name = keep,
        samples = vapply(keep, function(n) {
          v <- segments_meta[[n]]$samples
          if (is.null(v)) NA_integer_ else as.integer(v)
        }, integer(1), USE.NAMES = FALSE),
        twoN_param = vapply(keep, function(n) {
          v <- segments_meta[[n]]$twoN_param
          if (is.null(v)) NA_character_ else v
        }, character(1), USE.NAMES = FALSE),
        twoN = vapply(keep, function(n) {
          tp <- segments_meta[[n]]$twoN_param
          if (is.null(tp)) return(NA_real_)
          if (tp %in% names(param_val)) return(unname(param_val[[tp]]))
          lit <- suppressWarnings(as.numeric(tp))   # inline literal twoN=10000
          if (!is.na(lit)) lit else NA_real_
        }, numeric(1), USE.NAMES = FALSE),
        time_param = vapply(keep, function(n) {
          v <- segments_meta[[n]]$time_param
          if (is.null(v)) NA_character_ else v
        }, character(1), USE.NAMES = FALSE),
        time = vapply(keep, function(n) {
          tp <- segments_meta[[n]]$time_param
          if (is.null(tp)) return(NA_real_)
          if (tp %in% names(param_val)) return(unname(param_val[[tp]]))
          lit <- suppressWarnings(as.numeric(tp))   # inline literal t=0.5
          if (!is.na(lit)) lit else NA_real_
        }, numeric(1), USE.NAMES = FALSE),
        admix_event_time = NA_real_
      )
      # Carry the param fixed/free/constrained map on the nodes tibble itself,
      # so it travels with `nodes` and is dropped wherever `nodes` is dropped
      # (e.g. topology comparisons that null attr(.,"nodes")). graph_to_lgo
      # reads it back to re-emit captured params with their original type.
      attr(nt, "lgo_param_types") <- param_type
      attr(result, "nodes") <- nt
    }
  }

  # `as = "igraph"` return path
  if (as == "igraph") {
    edge_mat <- as.matrix(result[, c("from", "to")])
    ig <- igraph::graph_from_edgelist(edge_mat)
    for (col in setdiff(names(result), c("from", "to"))) {
      igraph::edge_attr(ig, col) <- result[[col]]
    }
    if (length(param_bounds) > 0) {
      ig <- igraph::set_graph_attr(ig, "param_bounds", param_bounds)
    }
    return(ig)
  }

  # Only attach param_bounds when non-empty to avoid changing the class/
  # attribute footprint of files that have no bounded declarations — this
  # keeps round-trip equality tests (which use expect_equal and compare
  # attributes) clean for standard .lgo files.
  if (length(param_bounds) > 0) attr(result, "param_bounds") <- param_bounds
  result
}

# Classify LEGOFIT parameter names into families for use by
# read_legofit_output and read_legofit_bootstrap.
# Vectorized over `name`. Returns a character vector with values in:
#   "admix_time" | "time" | "twoN" | "twoN_one" | "twoN_shared" | "mixFrac" | "unknown"
# Pattern ordering matters: T_admix_* must be checked before T_* (otherwise
# T_admix_M would classify as time for node "admix_M").
#' @keywords internal
classify_legofit_param <- function(name) {
  result <- rep("unknown", length(name))
  result[grepl("^T_admix_(.+)$",  name, perl = TRUE)] <- "admix_time"
  result[grepl("^T_(.+)$",        name, perl = TRUE) &
         result == "unknown"]                          <- "time"
  result[grepl("^twoN_(.+)$",     name, perl = TRUE)] <- "twoN"
  result[name == "one"]                               <- "twoN_one"
  result[name == "shared"]                            <- "twoN_shared"
  result[grepl("^m_(.+)$",        name, perl = TRUE)] <- "mixFrac"
  result
}

#' @keywords internal
validate_via_roundtrip <- function(lgo_text, edges) {
  parsed <- tryCatch(
    read_lgo(text = lgo_text, as = "edges"),
    error = function(e) {
      rlang::abort(
        c("graph_to_lgo produced .lgo that read_lgo cannot parse.",
          "i" = "This is a bug in graph_to_lgo or read_lgo.",
          "x" = paste("Parser error:", conditionMessage(e))),
        class = "legofit_validation_failed",
        parent = e
      )
    }
  )

  # read_lgo emits type="normal" for non-admix edges; the input may use
  # "edge" (qpgraph convention) or "normal" (random_sim convention) —
  # treat them as equivalent for the topology check.
  norm_type <- function(t) ifelse(t == "edge", "normal", t)
  expected <- edges  %>% dplyr::select(from, to, type) %>%
              dplyr::mutate(type = norm_type(type)) %>%
              dplyr::arrange(from, to)
  got      <- parsed %>% dplyr::select(from, to, type) %>%
              dplyr::mutate(type = norm_type(type)) %>%
              dplyr::arrange(from, to)
  # This is a topology check (from/to/type rows only). `edges` may carry a
  # `nodes` attribute (propagated by dplyr verbs onto `expected` but not onto
  # `got`, which is parsed from text), and read_lgo may attach `param_bounds`;
  # compare values only so neither extraneous attribute registers as a spurious
  # topology mismatch.
  if (!isTRUE(all.equal(expected, got, check.attributes = FALSE))) {
    rlang::abort(
      c("graph_to_lgo's output round-trips to a different topology.",
        "i" = "Likely a bug in parameter naming or .lgo assembly."),
      class = "legofit_validation_failed"
    )
  }
  invisible(TRUE)
}

#' Export an admixtools admixture graph to LEGOFIT .lgo format
#'
#' Converts an admixtools admixture graph (as an edge tibble or igraph) into
#' a LEGOFIT-compatible `.lgo` file.
#'
#' @param graph An edge tibble (columns `from`, `to`, `type`, `weight`) or an
#'   igraph with `type` and `weight` edge attributes.
#' @param file Output file path. If `NULL` (default), returns the `.lgo` text
#'   invisibly without writing.
#' @param samples Samples per leaf. A scalar (applied to all leaves) or a named
#'   integer vector (per-leaf overrides). Default `1`.
#' @param dates_terminal Starting time value for terminal (leaf) nodes. Default
#'   `0`.
#' @param outpop Name of an outgroup leaf to strip from the graph before export.
#'   `NULL` (default) skips stripping.
#' @param twoN Population size parameterization. `NULL` (default) emits
#'   `twoN fixed one=1` (coalescent units). A scalar numeric emits a fixed
#'   shared `twoN`. A named numeric vector emits per-segment free `twoN`
#'   parameters.
#' @param time_handling How to handle node times. One of:
#'   - `"fix_admix"` (default): admixture-destination nodes get no `t=`
#'     declaration (their time is implicit). Requires drift values that are
#'     additively consistent along every root-to-leaf path.
#'   - `"init"`: all non-leaf nodes are free starting points; respects an
#'     explicit `time` column on the edge tibble. Same consistency
#'     requirement as `"fix_admix"`.
#'   - `"free"`: all non-leaf times are declared `free` with depth-based
#'     starting values. Bypasses the topology walk entirely, so it accepts
#'     edge tibbles whose drifts are not additively consistent (e.g.,
#'     outputs from [qpgraph()], which optimizes drift for f-stat fit
#'     rather than time-walk consistency).
#' @param drift_to_time Function converting per-edge drift to branch length.
#'   Default [default_drift_to_time()] (identity). Ignored when
#'   `time_handling = "free"` or when an explicit `time` column is present.
#' @param fix_times Logical (default `FALSE`). When `TRUE`, internal-node times
#'   are emitted as `time fixed` at their computed absolute values rather than
#'   `time free`. Use this with a free `twoN` (a named `twoN=` vector) when you
#'   need to recover **absolute** effective population sizes: with both times and
#'   `twoN` free, site-pattern data fit only the `Δt/twoN` ratios, so the absolute
#'   scale is unidentified (the data constrain the ratios but not the overall
#'   scale). Fixing the times removes that degeneracy. With
#'   `"fix_admix"` or `"init"` the fixed values are the real absolute times, so
#'   the recovered `twoN` is absolute. With `"free"` (the only mode that exports
#'   additively-inconsistent graphs, e.g. most admixture topologies) the fixed
#'   values are topological-depth placeholders: still consistent and still enough
#'   to make `twoN` identifiable in a self-consistent fit, but not biologically
#'   calibrated absolute times.
#' @param validate If `TRUE` (default), round-trips the output through
#'   [read_lgo()] to confirm topology is preserved.
#' @return The `.lgo` text, invisibly.
#' @examples
#' # A 3-leaf no-admixture graph with explicit node times.
#' g <- tibble::tribble(
#'   ~from,  ~to,  ~type,    ~weight,  ~time,
#'   "xyz",  "xy", "normal", NA_real_, 1.5,
#'   "xyz",  "z",  "normal", NA_real_, 2,
#'   "xy",   "x",  "normal", NA_real_, 0.5,
#'   "xy",   "y",  "normal", NA_real_, 0.5
#' )
#' cat(graph_to_lgo(g, time_handling = "init"))
#' @export
graph_to_lgo <- function(graph,
                         file = NULL,
                         samples = 1,
                         dates_terminal = 0,
                         outpop = NULL,
                         twoN = NULL,
                         time_handling = c("fix_admix", "init", "free"),
                         drift_to_time = default_drift_to_time,
                         fix_times = FALSE,
                         validate = TRUE) {

  time_handling <- match.arg(time_handling)
  samples_supplied <- !missing(samples)   # only warn on an explicit samples= override

  edges <- coerce_to_edge_tibble(graph)
  edges <- validate_edge_tibble(edges)

  if (!is.null(outpop)) edges <- strip_outgroup(edges, outpop)

  params <- generate_param_names(edges)
  times  <- compute_times(edges, time_handling, drift_to_time, dates_terminal,
                          fix_times = fix_times)
  twoN_decls <- resolve_twoN_decls(edges, twoN)

  leaves       <- setdiff(edges$to, edges$from)
  all_nodes    <- unique(c(edges$from, edges$to))
  samples_vec  <- resolve_scalar_or_named(samples, leaves, default = NA_integer_)
  full_samples <- setNames(rep(NA_integer_, length(all_nodes)), all_nodes)
  full_samples[leaves] <- as.integer(samples_vec[leaves])

  # Merge per-node samples from the nodes tibble (precedence: nodes
  # tibble wins over the leaves-only default). Internal sampled nodes
  # (e.g. a Neanderthal ancestral to moderns) thus get emitted and are
  # counted by the unfittable check below.
  nodes_tbl <- attr(edges, "nodes")
  if (!is.null(nodes_tbl) && nrow(nodes_tbl) > 0) {
    ns <- nodes_tbl$samples
    have <- !is.na(ns) & nodes_tbl$name %in% all_nodes
    if (any(have)) {
      # Warn only when an explicit samples= disagrees with the nodes tibble; the
      # samples= default (1) must not trip this on a plain graph_to_lgo() round
      # trip of a graph whose captured leaf samples differ from 1.
      arg_conflict <- intersect(nodes_tbl$name[have], names(full_samples))
      arg_conflict <- arg_conflict[
        !is.na(full_samples[arg_conflict]) &
        full_samples[arg_conflict] != ns[match(arg_conflict, nodes_tbl$name)]]
      if (samples_supplied && length(arg_conflict) > 0)
        rlang::warn(
          c("`samples=` argument conflicts with the nodes tibble; using the nodes tibble.",
            "i" = paste("Nodes:", paste(arg_conflict, collapse = ", "))),
          class = "admixtools_samples_arg_overridden")
      full_samples[nodes_tbl$name[have]] <- as.integer(ns[have])
    }
  }

  # Round-trip fidelity: when the graph carries a nodes tibble with captured
  # source params, re-emit those param names, values, and fixed/free types so
  # the .lgo references real, declared parameters (valid LEGOFIT) instead of
  # regenerated names. We override the per-node name maps and declarations the
  # formatters consume; segments and the Parameters block then stay in sync. A
  # constrained-time param is emitted as `free` at its evaluated value (the
  # constraint expression is not reproduced). An explicit `twoN=` still wins,
  # so twoN is only rewritten when the caller left it NULL.
  if (!is.null(nodes_tbl) && nrow(nodes_tbl) > 0) {
    ptypes <- attr(nodes_tbl, "lgo_param_types")
    is_free <- function(pname) {
      ty <- if (!is.null(ptypes) && !is.na(pname) && pname %in% names(ptypes))
              ptypes[[pname]] else "free"
      !identical(ty, "fixed")   # fixed stays fixed; free/constrained -> free
    }
    in_graph <- nodes_tbl$name %in% all_nodes

    has_t <- in_graph & !is.na(nodes_tbl$time_param) & !is.na(nodes_tbl$time)
    if (any(has_t)) {
      tn <- nodes_tbl$name[has_t]
      params$time[tn] <- nodes_tbl$time_param[has_t]
      times$value[tn] <- nodes_tbl$time[has_t]
      times$free[tn]  <- vapply(nodes_tbl$time_param[has_t], is_free,
                                logical(1), USE.NAMES = FALSE)
    }

    has_n <- in_graph & !is.na(nodes_tbl$twoN_param) & !is.na(nodes_tbl$twoN)
    if (is.null(twoN) && any(has_n)) {
      twoN_decls$per_segment[nodes_tbl$name[has_n]] <- nodes_tbl$twoN_param[has_n]
      cap <- list()
      for (i in which(has_n))
        cap[[nodes_tbl$twoN_param[i]]] <-
          list(value = nodes_tbl$twoN[i], free = is_free(nodes_tbl$twoN_param[i]))
      # One declaration per distinct referenced twoN param. Captured params use
      # their value+type; any segment still on the generated default uses the
      # `one`=1 (fixed) convention from resolve_twoN_decls(NULL).
      refs <- unique(twoN_decls$per_segment)
      twoN_decls$declarations <- paste(vapply(refs, function(p) {
        info <- if (!is.null(cap[[p]])) cap[[p]] else list(value = 1, free = FALSE)
        sprintf("twoN %s %s=%g", if (isTRUE(info$free)) "free" else "fixed",
                p, info$value)
      }, character(1), USE.NAMES = FALSE), collapse = "\n")
    } else if (!is.null(twoN) && any(has_n)) {
      # An explicit `twoN=` wins (it defaults to NULL, so a non-NULL value is a
      # deliberate override of the captured effective sizes). Warn where it
      # differs from a captured value so the discard of source sizes is visible.
      # Dispatch scalar-vs-named by length to match resolve_twoN_decls(), which
      # treats ANY length-1 numeric as a scalar (a length-1 named vector is not
      # a per-node map there); keying such a vector by name would yield NA.
      cn      <- nodes_tbl$name[has_n]
      cap_val <- nodes_tbl$twoN[has_n]
      arg_val <- if (length(twoN) == 1) rep(twoN[[1]], length(cn)) else twoN[cn]
      differ  <- abs(arg_val - cap_val) > 1e-6 * pmax(1, abs(cap_val))
      if (isTRUE(any(differ)))
        rlang::warn(
          c("`twoN=` argument overrides captured twoN values from the nodes tibble.",
            "i" = paste("Nodes:", paste(cn[differ], collapse = ", "))),
          class = "admixtools_nodes_twoN_overridden")
    }
  }

  lgo_text <- assemble_lgo(edges, params, times, twoN_decls, full_samples,
                           fix_times = fix_times)

  is_sampled <- !is.na(full_samples) & full_samples > 0
  n_sampled <- sum(is_sampled)
  if (n_sampled < 2) {
    sampled_names <- names(full_samples)[is_sampled]
    rlang::warn(
      c("Resulting .lgo has < 2 sampled segments; LEGOFIT cannot fit it (zero observable site patterns).",
        "i" = sprintf("Sampled segments (%d): %s", n_sampled,
                      if (length(sampled_names) > 0) paste(sampled_names, collapse = ", ") else "none"),
        "i" = "Pass `samples = c(nodeA = 1, nodeB = 1, ...)` or use a graph with more leaves."),
      class = "legofit_unfittable_lgo"
    )
  }

  if (validate) validate_via_roundtrip(lgo_text, edges)
  if (!is.null(file)) writeLines(lgo_text, file)

  invisible(lgo_text)
}

# ===========================================================================
# Phase B — read_legofit_output helpers + public function (Steps 6-7)
# ===========================================================================

# ---------------------------------------------------------------------------
# Step 6 helpers: extract_param_section, parse_param_lines,
#                 extract_convergence_status, inform_convergence_status
# ---------------------------------------------------------------------------

# Extract a named parameter block from LEGOFIT stdout.
# Returns list(fixed = character(0), free = character(0)) always.
# Walks lines looking for `header`, then captures Fixed:/Free:
# subsections until a blank line, a `#` comment, or a letter-at-col-1.
#' @keywords internal
extract_param_section <- function(lines, header) {
  # Header matching is whitespace-tolerant.
  hdr_re <- paste0("^\\s*", gsub(" ", "\\\\s+", header), "\\s*$")
  hdr_idx <- which(grepl(hdr_re, lines, perl = TRUE))
  if (length(hdr_idx) == 0) return(list(fixed = character(0), free = character(0)))
  hdr_idx <- hdr_idx[[1]]

  fixed_lines <- character(0)
  free_lines  <- character(0)
  current     <- "free"   # default bucket

  for (i in seq_len(max(0L, length(lines) - hdr_idx)) + hdr_idx) {
    ln <- lines[[i]]
    # Termination: blank line
    if (!nzchar(trimws(ln))) break
    # Termination: comment / BranchLen header
    if (grepl("^\\s*#", ln)) break
    # Subsection toggle
    if (grepl("^Fixed:\\s*$", ln)) { current <- "fixed"; next }
    if (grepl("^Free:\\s*$",  ln)) { current <- "free";  next }
    # Termination: letter at col 1 that is not Fixed/Free (next section)
    if (grepl("^[A-Za-z]", ln) && !grepl("^(Fixed|Free):", ln)) break
    # Parameter row: leading whitespace + contains `=`
    if (grepl("^\\s+.*=", ln)) {
      if (current == "fixed") {
        fixed_lines <- c(fixed_lines, ln)
      } else {
        free_lines <- c(free_lines, ln)
      }
    }
  }
  list(fixed = fixed_lines, free = free_lines)
}

# Parse `name = value` lines into a tibble.
# Robust to whitespace around `=`.
#' @keywords internal
parse_param_lines <- function(lines) {
  if (length(lines) == 0) return(tibble::tibble(name = character(), value = numeric()))
  bad <- character(0)
  rows <- lapply(lines, function(ln) {
    ln    <- trimws(ln)
    eq    <- regexpr("=", ln, fixed = TRUE)
    if (eq == -1) { bad <<- c(bad, ln); return(NULL) }
    nm  <- trimws(substr(ln, 1, eq - 1))
    val <- suppressWarnings(as.numeric(trimws(substr(ln, eq + 1, nchar(ln)))))
    list(name = nm, value = val)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(bad) > 0) {
    rlang::inform(
      c("Some parameter lines could not be parsed (no `=` found).",
        "i" = paste("Offending:", paste(bad, collapse = "; "))),
      class = "legofit_invalid_input"
    )
  }
  if (length(rows) == 0) return(tibble::tibble(name = character(), value = numeric()))
  tibble::tibble(
    name  = vapply(rows, `[[`, character(1), "name"),
    value = vapply(rows, `[[`, numeric(1),   "value")
  )
}

# Extract DiffEv convergence status line from legofit output.
# Returns list(status, cost, spread) with NAs if not found.
#' @keywords internal
extract_convergence_status <- function(lines) {
  re <- paste0(
    "^DiffEv\\s+(reached_goal|finished_iterations)",
    "\\.\\s+cost=(\\S+)\\s+spread=(\\S+)\\s*$"
  )
  idx <- grep(re, lines, perl = TRUE)
  if (length(idx) == 0) {
    return(list(status = NA_character_, cost = NA_real_, spread = NA_real_))
  }
  m <- regmatches(lines[[idx[[1]]]], regexec(re, lines[[idx[[1]]]], perl = TRUE))[[1]]
  list(
    status = m[[2]],
    cost   = suppressWarnings(as.numeric(m[[3]])),
    spread = suppressWarnings(as.numeric(m[[4]]))
  )
}

# Emit convergence-aware informs.
# `path` is incorporated into .frequency_id so each file gets an independent
# once-per-session dedup key; without it a multi-file workflow would suppress
# all but the first convergence warn.
#' @keywords internal
inform_convergence_status <- function(convergence, path = NULL) {
  freq_suffix <- if (!is.null(path)) path else ""
  if (is.null(convergence) || is.na(convergence$status)) {
    rlang::inform(
      "Could not detect convergence status; fit may be incomplete.",
      class = "legofit_fit_status_unknown",
      .frequency = "once",
      .frequency_id = paste0("legofit_fit_status_unknown:", freq_suffix)
    )
    return(invisible())
  }
  if (convergence$status == "finished_iterations") {
    rlang::inform(
      c("LEGOFIT hit max iterations without reaching the convergence goal.",
        "i" = sprintf(
          "Final cost=%g, spread=%g. Consider re-running with more iterations or a tighter tolerance.",
          convergence$cost, convergence$spread)),
      class = "legofit_fit_incomplete",
      .frequency = "once",
      .frequency_id = paste0("legofit_fit_incomplete:", freq_suffix)
    )
  }
  invisible()
}

# ---------------------------------------------------------------------------
# Step 7 helpers: drop_ghost_params, pivot_params_to_edges, surface_mismatches
# ---------------------------------------------------------------------------

# Drop ghost-segment params from the param table.
# Ghost segments produce T_<ghost_name> params only if LEGOFIT
# bleeds segment naming into the param namespace — which standard LEGOFIT does
# NOT do (ghost segments reference T_admix_<dest>, not their own T_<ghost>).
# This function is therefore a no-op in practice; retained for future-proofing.
#' @keywords internal
drop_ghost_params <- function(params, edges) {
  params  # no-op: LEGOFIT does not emit ghost-segment-named params
}

# Pivot fitted params onto the edge tibble.
# Attaches time, admix_event_time, twoN, admix_prop, admix_event columns.
# Root and internal nodes with no incoming edge get their fitted times
# in attr(result, "node_times").
#' @keywords internal
pivot_params_to_edges <- function(params, edges) {
  all_nodes    <- unique(c(edges$from, edges$to))
  admix_dests  <- unique(edges$to[edges$type == "admix"])

  # Build fast lookup from param name to value
  pv <- setNames(params$value, params$name)

  # Per-edge time: T_<edges$to>
  edges$time <- vapply(edges$to, function(n) {
    nm <- paste0("T_", n)
    if (nm %in% names(pv)) pv[[nm]] else NA_real_
  }, numeric(1))

  # admix_event_time: T_admix_<edges$to> (only for admix-typed edges)
  edges$admix_event_time <- vapply(seq_len(nrow(edges)), function(i) {
    if (edges$type[[i]] != "admix") return(NA_real_)
    nm <- paste0("T_admix_", edges$to[[i]])
    if (nm %in% names(pv)) pv[[nm]] else NA_real_
  }, numeric(1))

  # twoN per edge: try twoN_<to> first; if all missing and `one`/`shared`
  # sentinel present, fill accordingly.
  has_one    <- "one"    %in% params$name
  has_shared <- "shared" %in% params$name
  shared_val <- if (has_shared) pv[["shared"]] else NA_real_

  edges$twoN <- vapply(edges$to, function(n) {
    nm <- paste0("twoN_", n)
    if (nm %in% names(pv)) return(pv[[nm]])
    if (has_shared) return(shared_val)
    if (has_one)    return(1)
    NA_real_
  }, numeric(1))

  # admix_prop: m_<to> on admix edges; weight on low-weight parent.
  # LEGOFIT emits m = low weight (same convention as graph_to_lgo).
  edges$admix_prop <- vapply(seq_len(nrow(edges)), function(i) {
    if (edges$type[[i]] != "admix") return(NA_real_)
    nm <- paste0("m_", edges$to[[i]])
    if (nm %in% names(pv)) pv[[nm]] else NA_real_
  }, numeric(1))

  # admix_event: param name grouping key for the admix event
  edges$admix_event <- vapply(seq_len(nrow(edges)), function(i) {
    if (edges$type[[i]] != "admix") return(NA_character_)
    paste0("m_", edges$to[[i]])
  }, character(1))

  # node_times attribute: fitted T_<n> for ALL graph nodes. Edges only carry
  # the time for their child endpoint; this attribute additionally exposes the
  # root (and every other node) so callers can always retrieve any node's fitted
  # time without reconstructing it from the edge tibble.
  node_times_all <- vapply(all_nodes, function(n) {
    nm <- paste0("T_", n)
    if (nm %in% names(pv)) pv[[nm]] else NA_real_
  }, numeric(1))
  names(node_times_all) <- all_nodes
  attr(edges, "node_times") <- node_times_all

  tibble::as_tibble(edges)
}

# Compute and emit structural-mismatch inform.
# Rule chain applied in strict order (order matters).
#
# `sentinel_names`: character vector of sentinel param names ("one", "shared")
#   extracted from the Initial/Fixed block. Passed separately so that
#   include_fixed=FALSE (which omits fixed params from `params`) still lets
#   Rules 4-5 correctly detect the coalescent-unit or scalar-twoN mode.
#
# `path`: incorporated into .frequency_id so each file gets an independent
#   once-per-session dedup key.
#' @keywords internal
surface_mismatches <- function(params, edges, include_fixed = TRUE,
                               path = NULL, sentinel_names = character(0)) {
  all_nodes   <- unique(c(edges$from, edges$to))
  admix_dests <- unique(edges$to[edges$type == "admix"])

  # Expected param names for this graph.
  # Guard paste0() with length() > 0: in R 4.6, paste0("prefix", character(0))
  # returns "prefix" rather than character(0) (recycling behavior), so we
  # explicitly suppress the computation when the vector is empty.
  time_names       <- if (length(all_nodes)   > 0) paste0("T_",       all_nodes)   else character(0)
  twoN_names       <- if (length(all_nodes)   > 0) paste0("twoN_",    all_nodes)   else character(0)
  admix_time_names <- if (length(admix_dests) > 0) paste0("T_admix_", admix_dests) else character(0)
  mfrac_names      <- if (length(admix_dests) > 0) paste0("m_",       admix_dests) else character(0)
  expected_names <- c(time_names, admix_time_names, twoN_names, mfrac_names)

  seen_names <- params$name
  # Merge in any sentinels passed from the Initial/Fixed block so that
  # had_one / had_shared detection is correct even when include_fixed=FALSE.
  all_for_sentinels <- union(seen_names, sentinel_names)

  # Rule 1: HARD CHECK — mutual exclusion of twoN modes
  if (all(c("one", "shared") %in% all_for_sentinels)) {
    rlang::abort(
      c("File mixes coalescent-unit (one=1) and scalar (shared=N) twoN modes.",
        "x" = "Both `one` and `shared` are present; these are mutually exclusive in LEGOFIT."),
      class = "legofit_invalid_input"
    )
  }

  # Rule 2: Remove sentinels from both sets (they are LEGOFIT internals)
  had_one    <- "one"    %in% all_for_sentinels
  had_shared <- "shared" %in% all_for_sentinels
  seen_names     <- setdiff(seen_names,     c("one", "shared"))
  expected_names <- setdiff(expected_names, c("one", "shared"))

  # Rule 3: compute missing and extra
  missing <- setdiff(expected_names, seen_names)
  extra   <- setdiff(seen_names, expected_names)

  # Rule 4: twoN-sentinel suppression. In coalescent-unit mode (`one`) or
  # scalar mode (`shared`) LEGOFIT emits no per-node twoN params, so the whole
  # expected twoN set legitimately goes missing. Rule 1 already guarantees
  # `one` and `shared` are mutually exclusive, so a single combined branch
  # covers both modes.
  twoN_expected <- twoN_names   # already computed above, reuse
  if ((had_one || had_shared) && all(twoN_expected %in% missing)) {
    missing <- setdiff(missing, twoN_expected)
  }

  # Rule 5: leaf-time suppression (only when include_fixed = FALSE)
  if (!include_fixed) {
    leaves <- setdiff(edges$to, edges$from)
    # Guard against R 4.6 paste0() recycling: paste0("T_", character(0))
    # returns "T_" rather than character(0), consistent with the guards above.
    leaf_time_params <- if (length(leaves) > 0) paste0("T_", leaves) else character(0)
    missing <- setdiff(missing, leaf_time_params)
  }

  # Rule 6: emit inform if any mismatch remains
  if (length(missing) > 0 || length(extra) > 0) {
    parts <- character(0)
    if (length(missing) > 0) {
      sample_missing <- paste(head(missing, 5), collapse = ", ")
      if (length(missing) > 5) sample_missing <- paste0(sample_missing, ", ...")
      parts <- c(parts, sprintf("Missing (%d): %s", length(missing), sample_missing))
    }
    if (length(extra) > 0) {
      sample_extra <- paste(head(extra, 5), collapse = ", ")
      if (length(extra) > 5) sample_extra <- paste0(sample_extra, ", ...")
      parts <- c(parts, sprintf("Extra (%d): %s", length(extra), sample_extra))
    }
    freq_suffix <- if (!is.null(path)) path else ""
    rlang::inform(
      c("Parameter set does not exactly match graph-derived expectations.",
        setNames(parts, rep("i", length(parts)))),
      class = "legofit_param_mismatch",
      .frequency = "once",
      .frequency_id = paste0("legofit_param_mismatch:", freq_suffix)
    )
  }
  invisible(NULL)
}

# Structural parameter-identifiability classification.
#
# Returns a tibble (parameter, family, identifiability, reason). Classes:
#   "fixed"            non-free parameter; an input, not an estimate.
#   "structural_none"  data CANNOT constrain it (value is arbitrary): admix-event
#                      times, whose segment subtends a single lineage and so
#                      admits no coalescence, meaning its duration cannot move
#                      any site-pattern frequency.
#   "scale_degenerate" identified only up to a global scale: a model with BOTH
#                      free times and free twoN fits only the Delta_t/twoN ratios,
#                      not absolute values. Fires only in named-twoN mode; with
#                      fixed `one`/`shared` twoN it does not.
#   "weak_identified" / "weak_unconstrained"
#                      a deep tree split time whose coalescent depth puts it in
#                      the downward-biased regime (see the deep-split tier below).
#   "identifiable"     the data constrain it (mixFrac, shallow tree split times).
#
# IMPORTANT: this is structural, derived from topology, fitted depths, and the
# model's free/fixed declarations, never from bootstrap CI width: in the failing
# regime the CI is narrow around a biased estimate, so CI width is not evidence
# of identifiability.
#' @keywords internal
classify_identifiability <- function(params, edges = NULL, node_times = NULL) {
  n <- nrow(params)
  if (n == 0) {
    return(tibble::tibble(parameter = character(), family = character(),
                          identifiability = character(), reason = character()))
  }
  free <- if (!is.null(params$free)) params$free else rep(TRUE, n)
  fam  <- params$family

  # Global scale degeneracy: at least one free time AND one free twoN.
  has_free_time <- any(free & fam == "time")
  has_free_twoN <- any(free & fam == "twoN")
  scale_degenerate <- has_free_time && has_free_twoN

  cls <- character(n); rsn <- character(n)
  for (i in seq_len(n)) {
    if (!isTRUE(free[i])) {
      cls[i] <- "fixed"
      rsn[i] <- "fixed parameter (an input, not an estimate)"
    } else if (fam[i] == "admix_time") {
      cls[i] <- "structural_none"
      rsn[i] <- "admix-event time: single-lineage segment admits no coalescence; fitted value is arbitrary"
    } else if (scale_degenerate && fam[i] %in% c("time", "twoN")) {
      cls[i] <- "scale_degenerate"
      rsn[i] <- "absolute scale unidentified: model has free times and free twoN; only Delta_t/twoN ratios are fitted"
    } else {
      cls[i] <- "identifiable"
      rsn[i] <- "data-constrained"
    }
  }

  # Deep-split weak-identification tier. A tree split time is recovered with a
  # downward bias that grows with its coalescent depth D = t/twoN. Past a few
  # coalescent units, coalescence within the branch saturates and the split time
  # stops moving the site-pattern spectrum (the single-lineage segment handled by
  # `structural_none` above is the limiting case). The onset depth depends on
  # tree shape: a split with a directly attached sampled leaf child is anchored
  # by that leaf's coalescences and holds far better than one whose children are
  # all internal clades. This tier only runs on the graph path (needs topology +
  # fitted depths), and only downgrades parameters the structural tier left
  # "identifiable".
  #
  # Threshold provenance. The cutoffs below are calibrated from coalescent
  # simulation, not assumed: independent ground-truth datasets simulated with
  # msprime across a depth x tree-shape x leaf-count grid, fit with `legofit -1`
  # (singletons on, which carry the deep-branch-length signal), comparing fitted
  # against true split times. A matched-pair experiment (same leaves and depth,
  # the focal node's leaf child toggled on/off) confirmed the leaf child is the
  # causal discriminator, not a proxy for global tree balance. Mean signed
  # relative error of the focal split, by depth D:
  #   no leaf child:   ~ -6% (D=5), -19% (D=6), -42% (D=8)  => well<=4, weak<=6
  #   has a leaf child:~ -2% (D=5),  -3% (D=6),  -4% (D=8)  => well<=5, weak<=9
  # These are an advisory: every flagged estimate is biased downward, so treat it
  # as a lower bound, and never infer good identification from a narrow bootstrap
  # CI (the failing regime is biased low WITH a narrow CI). Method reference:
  # Rogers (2019) "Legofit: estimating population history from genetic data",
  # BMC Bioinformatics 20:526.
  #
  # Large-tree margin. At >= 8 leaves the focal split degrades earlier than the
  # 4-7 leaf grid implies (even the leaf-child variant), and the exact large-tree
  # knees are not yet pinned, so both cutoffs are tightened by one coalescent unit
  # as a conservative margin. This changes the CLASS (a borderline deep split
  # becomes weak_unconstrained), not merely the reason text, so a genuinely
  # unconstrained split on a large tree is not reported as weak_identified.
  if (!is.null(edges) && !is.null(node_times)) {
    pv     <- stats::setNames(params$value, params$name)
    leaves <- setdiff(edges$to, edges$from)
    nleaf  <- length(leaves)
    margin <- if (nleaf >= 8) 1 else 0     # conservative large-tree tightening
    twoN_of <- function(node) {
      nm <- paste0("twoN_", node)
      if (nm %in% names(pv)) return(pv[[nm]])
      if ("shared" %in% names(pv)) return(pv[["shared"]])
      1                                    # `one` (coalescent units) or default
    }
    for (i in seq_len(n)) {
      if (cls[i] != "identifiable" || fam[i] != "time") next
      node <- sub("^T_", "", params$name[i])
      t    <- if (node %in% names(node_times)) node_times[[node]] else NA_real_
      tn   <- twoN_of(node)
      if (!is.finite(t) || !is.finite(tn) || tn <= 0) next
      D    <- t / tn
      has_leaf <- any(edges$to[edges$from == node] %in% leaves)
      base <- if (has_leaf) c(well = 5, weak = 9) else c(well = 4, weak = 6)
      well <- base[["well"]] - margin
      weak <- base[["weak"]] - margin
      big  <- if (margin > 0) "; >=8-leaf conservative margin applied" else ""
      if (D > weak) {
        cls[i] <- "weak_unconstrained"
        rsn[i] <- sprintf(
          "deep split (coalescent depth %.1f, %s leaf child): recovery effectively unconstrained, biased low%s",
          D, if (has_leaf) "has" else "no", big)
      } else if (D > well) {
        cls[i] <- "weak_identified"
        rsn[i] <- sprintf(
          "deep split (coalescent depth %.1f, %s leaf child): weakly identified, estimate biased low; treat as a lower bound%s",
          D, if (has_leaf) "has" else "no", big)
      }
    }
  }

  tibble::tibble(parameter = params$name, family = fam,
                 identifiability = cls, reason = rsn)
}

# ---------------------------------------------------------------------------
# Step 7 public function: read_legofit_output
# ---------------------------------------------------------------------------

#' Read a LEGOFIT fitted-output file into an admixtools edge tibble
#'
#' Parses the stdout of a `legofit` run into an admixtools edge tibble with
#' fitted parameter values attached as additional columns. When `graph` is
#' supplied, parameters are mapped to edges by segment name; when `NULL`, a
#' raw parameter tibble is returned.
#'
#' @section Parameter naming convention:
#' Family classification (the `family` column, and the edge mapping when `graph`
#' is supplied) keys off the `[graph_to_lgo()]` naming convention: `T_<node>`
#' (split time), `T_admix_<dest>` (admix-event time), `twoN_<seg>` (population
#' size), `m_<dest>` (mixFrac), plus the `one`/`shared` twoN sentinels. A
#' `.legofit` produced from a hand-written `.lgo` with arbitrary parameter names
#' (e.g. `Tnd`, `Tmnd`) parses correctly but every parameter classifies as
#' `"unknown"` and cannot be mapped onto graph edges. To read a third-party fit
#' into an edge tibble, export the matching topology with [graph_to_lgo()] so the
#' names follow the convention.
#'
#' @param path Path to a `.legofit` file (legofit stdout, UTF-8).
#' @param graph An edge tibble, igraph, or `NULL`. When supplied, parameter
#'   names are mapped to admixtools edge endpoints per the segment-naming
#'   convention from [graph_to_lgo()]. When `NULL`, returns a raw parameter
#'   table.
#' @param include_fixed Logical (default `TRUE`). When `TRUE`, fixed parameter
#'   values from the **Initial parameter values** block are included in the
#'   result. When `FALSE`, only free-fitted values are returned.
#'
#' @section Fitted columns (when `graph` is supplied):
#' \describe{
#'   \item{`time`}{Fitted `T_<to>` for each edge's younger endpoint.}
#'   \item{`admix_event_time`}{Fitted `T_admix_<to>` for admix-typed edges.}
#'   \item{`twoN`}{Fitted effective population size for `<to>`.}
#'   \item{`admix_prop`}{Fitted mixFrac (`m_<to>`) for admix-typed edges.}
#'   \item{`admix_event`}{The mixFrac parameter name (grouping ID).}
#' }
#'
#' @section Attributes:
#' \describe{
#'   \item{`node_times`}{Named numeric vector of fitted `T_<n>` for ALL graph
#'     nodes, including the root and other nodes that have no incoming edge and
#'     therefore do not appear in the edge-tibble `time` column. Example:
#'     \preformatted{
#'       result <- read_legofit_output("foo.legofit", graph = g)
#'       attr(result, "node_times")
#'       #>  xyz    xy     x     y     z
#'       #> 2.00  0.50  0.00  0.00  0.00
#'     }
#'   }
#'   \item{`fit_convergence`}{List with `status`, `cost`, `spread` from the
#'     DiffEv summary line. `status` is `"reached_goal"` (converged),
#'     `"finished_iterations"` (hit iteration cap), or `NA` (unknown).}
#'   \item{`identifiability`}{Per-parameter tibble (`parameter`, `family`,
#'     `identifiability`, `reason`) flagging which fitted values the site-pattern
#'     data can actually constrain. Classes: `"identifiable"` (data-constrained,
#'     e.g. mixFracs and shallow split times); `"structural_none"` (the value is
#'     arbitrary, e.g. admix-event times, whose single-lineage segment admits no
#'     coalescence); `"scale_degenerate"` (identified only up to a global scale,
#'     when a model has both free times and free `twoN` so only `Δt/twoN` ratios
#'     are fitted); `"weak_identified"` and `"weak_unconstrained"` (a deep tree
#'     split time whose coalescent depth `t/twoN` puts it in the downward-biased
#'     regime, with a looser onset depth when the split has a directly attached
#'     sampled leaf child, so treat these as lower bounds); and `"fixed"` (an
#'     input, not an estimate). This is a
#'     **structural** judgement from topology, fitted depths, and the model's
#'     declarations; it is
#'     never inferred from bootstrap CI width, because a non-identifiable
#'     parameter can carry a narrow CI around a badly biased estimate.}
#' }
#'
#' @return A tibble (with `node_times` and `fit_convergence` attributes) when
#'   `graph` is supplied; a raw parameter tibble (with `fit_convergence`) when
#'   `graph = NULL`.
#' @examples
#' # Construct a minimal `.legofit` output and read it back.
#' tmp <- tempfile(fileext = ".legofit")
#' writeLines(c(
#'   "Initial parameter values",
#'   "Fixed:",
#'   "       one = 1",
#'   "Free:",
#'   "     T_R = 2",
#'   "DiffEv reached_goal. cost=1e-10 spread=1e-7",
#'   "Fitted parameter values",
#'   "Free:",
#'   "     T_R = 2"
#' ), tmp)
#' read_legofit_output(tmp)
#' unlink(tmp)
#' @export
read_legofit_output <- function(path, graph = NULL, include_fixed = TRUE) {

  lines <- readLines(path, encoding = "UTF-8")

  # 2. Extract the two parameter blocks.
  # Always parse the Initial block to detect coalescent-unit sentinels (one=1 /
  # shared=N). surface_mismatches Rules 4-5 need these even when
  # include_fixed=FALSE (the sentinels live in the Fixed sub-block, not Fitted).
  fitted      <- extract_param_section(lines, header = "Fitted parameter values")
  initial_all <- extract_param_section(lines, header = "Initial parameter values")

  # Sentinel names extracted unconditionally for surface_mismatches.
  sentinel_names_from_file <- intersect(
    parse_param_lines(initial_all$fixed)$name, c("one", "shared")
  )

  # 3. Parse: fitted block uses `free` only; initial block uses `fixed` only.
  free_fitted   <- parse_param_lines(fitted$free)
  fixed_initial <- if (include_fixed) parse_param_lines(initial_all$fixed) else NULL

  # 3a. Convergence status
  convergence <- extract_convergence_status(lines)

  if (nrow(free_fitted) == 0 &&
      (is.null(fixed_initial) || nrow(fixed_initial) == 0)) {
    rlang::abort(
      paste0("File does not contain a `Fitted parameter values` block with ",
             "Free declarations.\n  path: ", path),
      class = "legofit_invalid_input"
    )
  }

  params <- dplyr::bind_rows(
    dplyr::mutate(free_fitted,  source = "fitted_free",    free = TRUE),
    if (!is.null(fixed_initial) && nrow(fixed_initial) > 0) {
      dplyr::mutate(fixed_initial, source = "initial_fixed", free = FALSE)
    }
  )

  # 4. Classify each param
  params$family <- classify_legofit_param(params$name)

  # 5. Optionally align to graph
  if (!is.null(graph)) {
    edges <- coerce_to_edge_tibble(graph)
    validate_edge_tibble(edges)

    # Read-time T_admix collision check
    admix_dests  <- unique(edges$to[edges$type == "admix"])
    all_nodes    <- unique(c(edges$from, edges$to))
    for (pname in params$name[params$family == "admix_time"]) {
      x <- sub("^T_admix_(.+)$", "\\1", pname, perl = TRUE)
      if (paste0("admix_", x) %in% all_nodes && x %in% admix_dests) {
        rlang::abort(
          c("Ambiguous parameter name in fitted output.",
            "x" = sprintf(
              "`%s` could refer to the admix-event time OR the time of node `admix_%s`.",
              pname, x),
            "i" = "Supply a graph without nodes named `admix_<admix_dest>`."),
          class = "legofit_invalid_input"
        )
      }
    }

    params <- drop_ghost_params(params, edges)
  }

  # 6. Return raw param table when no graph supplied.
  # Drop the internal `source` and `free` columns — they are implementation
  # details used for sentinel detection and graph alignment, not part of the
  # public API. Documented @return columns: name, value, family.
  # Also fire the convergence inform here (same as the graph path) so callers
  # using graph=NULL still learn when the fit hit max iterations.
  if (is.null(graph)) {
    result <- tibble::as_tibble(
      params[order(params$family, params$name), c("name", "value", "family")]
    )
    attr(result, "fit_convergence") <- convergence
    attr(result, "identifiability") <- classify_identifiability(params)
    inform_convergence_status(convergence, path = path)
    return(result)
  }

  # 7. Pivot to edge tibble
  result <- pivot_params_to_edges(params, edges)

  # 8. Surface structural mismatches
  surface_mismatches(params, edges, include_fixed = include_fixed,
                     path = path, sentinel_names = sentinel_names_from_file)

  # 9. Attach convergence + identifiability metadata and fire incomplete-fit inform
  # On the graph path, pass topology + fitted depths so the deep-split tier
  # (F-1b) can flag weakly-identified tree split times, not just the structural
  # cases the params-only classification covers.
  attr(result, "fit_convergence") <- convergence
  attr(result, "identifiability") <-
    classify_identifiability(params, edges, attr(result, "node_times"))
  inform_convergence_status(convergence, path = path)

  result
}

# ===========================================================================
# Phase C — read_legofit_bootstrap helpers + public function (Steps 8-9)
# ===========================================================================

# ---------------------------------------------------------------------------
# Step 8 helper: parse_bootci_output
# ---------------------------------------------------------------------------

# Parse bootci.py output table format.
# Verified format (2026-05-20):
#   # bootci.py run at: ...
#   # input: ...
#   # confidence: 0.950
#          par             est             low            high
#         T_xy      0.50000000      0.46250000      0.53750000
#    ...
# Optional lbl column when bootci.py -l <label> was used.
#' @keywords internal
parse_bootci_output <- function(lines) {
  # 1. Strip leading `#` comment lines; locate column header
  non_comment_idx <- which(!grepl("^\\s*#", lines))
  if (length(non_comment_idx) == 0) {
    rlang::abort(
      "No column header found in bootci.py output (all lines are comments).",
      class = "legofit_invalid_input"
    )
  }
  hdr_idx  <- non_comment_idx[[1]]
  hdr_toks <- strsplit(trimws(lines[[hdr_idx]]), "\\s+")[[1]]

  # 2. Validate header: must be c("par", "est", "low", "high") + optional "lbl"
  required_hdr <- c("par", "est", "low", "high")
  if (!identical(hdr_toks[seq_along(required_hdr)], required_hdr)) {
    rlang::abort(
      c("Unexpected column header in bootci.py output.",
        "x" = paste("Got:", paste(hdr_toks, collapse = " ")),
        "i" = "Expected: par est low high [lbl]"),
      class = "legofit_invalid_input"
    )
  }
  has_lbl  <- length(hdr_toks) >= 5 && hdr_toks[[5]] == "lbl"
  n_cols   <- if (has_lbl) 5L else 4L

  # 3. Parse data rows: guard seq() against the descending-range hazard when
  # the header is the last line (seq(n+1, n) counts down in R, not empty).
  # Also strip blank lines and any trailing comment lines (#...) that bootci.py
  # or downstream pipelines may append after the data table.
  if (hdr_idx >= length(lines)) {
    data_lines <- character(0)
  } else {
    data_lines <- lines[seq(hdr_idx + 1L, length(lines))]
    data_lines <- data_lines[nzchar(trimws(data_lines))]
    data_lines <- data_lines[!grepl("^\\s*#", data_lines)]
  }

  if (length(data_lines) == 0) {
    empty <- tibble::tibble(parameter = character(), est = numeric(),
                            low = numeric(),  high = numeric())
    if (has_lbl) empty$lbl <- character()
    return(empty)
  }

  rows <- lapply(data_lines, function(ln) {
    toks <- strsplit(trimws(ln), "\\s+")[[1]]
    if (length(toks) != n_cols) {
      rlang::abort(
        c("Unexpected number of columns in bootci.py data row.",
          "x" = paste("Got", length(toks), "columns, expected", n_cols),
          "i" = paste("Row:", ln)),
        class = "legofit_invalid_input"
      )
    }
    list(
      parameter = toks[[1]],
      est  = suppressWarnings(as.numeric(toks[[2]])),
      low  = suppressWarnings(as.numeric(toks[[3]])),
      high = suppressWarnings(as.numeric(toks[[4]])),
      lbl  = if (has_lbl) toks[[5]] else NA_character_
    )
  })

  result <- tibble::tibble(
    parameter = vapply(rows, `[[`, character(1), "parameter"),
    est       = vapply(rows, `[[`, numeric(1),   "est"),
    low       = vapply(rows, `[[`, numeric(1),   "low"),
    high      = vapply(rows, `[[`, numeric(1),   "high")
  )
  if (has_lbl) result$lbl <- vapply(rows, `[[`, character(1), "lbl")
  result
}

# ---------------------------------------------------------------------------
# Step 9 public function: read_legofit_bootstrap
# ---------------------------------------------------------------------------

#' Read a bootci.py bootstrap output file into a tidy CI table
#'
#' Parses the output of LEGOFIT's `bootci.py` utility into a tidy tibble with
#' per-parameter bootstrap confidence intervals. The point estimate column
#' (`point_estimate`) comes from `bootci.py`'s `est` column, which reflects
#' the real-data fit value.
#'
#' The confidence level used by `bootci.py` is stored in a comment line in the
#' file (e.g., `# confidence: 0.950`). This reader does not expose a `level`
#' argument; the bounds in the file are returned as-is. To change the level,
#' re-run `bootci.py` with the desired `--confidence` flag.
#'
#' @param path Path to a `bootci.py` output file (UTF-8 encoded).
#' @param graph An edge tibble, igraph, or `NULL`. When supplied, ghost-segment
#'   parameters are dropped (no-op for standard LEGOFIT output) and parameter
#'   families are classified.
#'
#' @return A tibble with columns `parameter`, `family`, `point_estimate`,
#'   `lo`, `hi`. If `bootci.py -l <label>` was used, an additional `lbl`
#'   column appears. The confidence level is embedded in the file as a comment
#'   (not returned as a column).
#' @examples
#' # Construct a minimal `bootci.py` output and read it back.
#' tmp <- tempfile(fileext = ".bootci")
#' writeLines(c(
#'   "# bootci.py run",
#'   "# confidence: 0.950",
#'   "       par             est             low            high",
#'   "      T_R      2.00000000      1.90000000      2.10000000"
#' ), tmp)
#' read_legofit_bootstrap(tmp)
#' unlink(tmp)
#' @export
read_legofit_bootstrap <- function(path, graph = NULL) {
  lines <- readLines(path, encoding = "UTF-8")

  # 1. Parse the CI table
  cis <- parse_bootci_output(lines)
  # → tibble(parameter, est, low, high [, lbl])

  # 2. Optionally align to graph (drops ghost params — no-op currently).
  # drop_ghost_params expects a column named 'name'; cis uses 'parameter'.
  # Rename around the call so the column contract is satisfied when the
  # function is eventually made non-trivial.
  if (!is.null(graph)) {
    edges <- coerce_to_edge_tibble(graph)
    validate_edge_tibble(edges)
    cis <- dplyr::rename(cis, name = parameter)
    cis <- drop_ghost_params(cis, edges)
    cis <- dplyr::rename(cis, parameter = name)
  }

  # 3. Rename columns to public names, then classify, then
  # reorder so final column order matches @return docs: parameter, family,
  # point_estimate, lo, hi [, lbl].
  cis <- dplyr::rename(cis, point_estimate = est, lo = low, hi = high)
  cis$family <- classify_legofit_param(cis$parameter)
  cis <- dplyr::select(cis, parameter, family, point_estimate, lo, hi,
                        dplyr::everything())

  cis
}
