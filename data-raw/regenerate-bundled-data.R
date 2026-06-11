# Regenerate bundled datasets in data/
#
# This script documents how the .rda files shipped under data/ are constructed
# and, in particular, how to refresh igraph objects stored in those datasets
# when a newer igraph release deprecates the on-disk representation.
#
# When to run
#   1. After an igraph release that triggers
#        "This graph was created by an old(er) igraph version. Call
#         `igraph::upgrade_graph()` on it ..."
#      messages from any vignette or test that touches a bundled graph.
#   2. When intentionally rebuilding a dataset from scratch (rare; the
#      construction recipes below depend on external genotype inputs not
#      shipped with the package).
#
# Datasets containing igraph objects (the ones this script actively touches):
#   example_igraph              single igraph
#   example_opt                 tibble with $graph list-column of igraph objects
#
# Datasets without igraph objects (recipes documented; not regenerated here):
#   example_f2_blocks, example_f2sim1, example_anno, example_graph,
#   example_triples, example_qpgraph_ref_results
#
# Storage convention (matches what already ships in data/):
#   RDX2 stream, bzip2 compression. Equivalent to
#     usethis::use_data(x, overwrite = TRUE, compress = "bzip2", version = 2)
#   Using base save() here so this script has no dev-only dependency.
#
# Note for maintainers
#   .Rbuildignore is git-ignored in this repo, so this directory ships in
#   the source tarball by default. If you want to exclude it from
#   R CMD build, add `^data-raw$` to your local .Rbuildignore.

stopifnot(basename(getwd()) == "admixtools")

library(admixtools)
library(igraph)
library(dplyr)
library(purrr)


# ---------------------------------------------------------------------------
# 1. Upgrade igraph objects in-place
# ---------------------------------------------------------------------------
# Loading the existing .rda files into a fresh session and running
# igraph::upgrade_graph() on each igraph object normalizes the internal
# representation to whatever the currently installed igraph expects. Vertex
# count, edge count, edge list, and vertex/edge attributes are preserved.

load("data/example_igraph.rda")
load("data/example_opt.rda")

example_igraph <- igraph::upgrade_graph(example_igraph)
example_opt$graph <- purrr::map(example_opt$graph, igraph::upgrade_graph)


# ---------------------------------------------------------------------------
# 2. Re-save with the existing storage convention
# ---------------------------------------------------------------------------

save(example_igraph, file = "data/example_igraph.rda",
     compress = "bzip2", version = 2)

save(example_opt,    file = "data/example_opt.rda",
     compress = "bzip2", version = 2)


# ---------------------------------------------------------------------------
# Construction recipes for the other bundled datasets
# ---------------------------------------------------------------------------
# These are not re-run here because they require external genotype data and,
# in some cases, an msprime simulation. Documented for future maintainers who
# need to rebuild from scratch.
#
# example_f2_blocks
#   Output of extract_f2() over a 7-population subset of a real-data genotype
#   panel (Chimp, Altai, Vindija, Denisova, Mbuti, French, Sardinian-style
#   ancient/modern set used throughout the vignettes).
#
# example_graph
#   A two-column character matrix of edges (from, to) describing the topology
#   used in the qpgraph vignette example.
#
# example_igraph
#   igraph::graph_from_edgelist(as.matrix(example_graph), directed = TRUE).
#   See commit f6f2b57 for the original `igraph::upgrade_graph()` pass.
#
# example_opt
#   find_graphs_old(example_f2_blocks, ...) on the 7-pop dataset, retaining
#   the per-generation graph search trace.
#
# example_triples
#   summarize_triples() over the same 7-pop dataset used for example_f2_blocks.
#
# example_qpgraph_ref_results
#   qpgraph(example_f2_blocks, example_graph) - the reference fit used to
#   sanity-check regressions in qpgraph internals.
#
# example_f2sim1
#   3d f2 array computed from msprime simulation output for a 5-pop scenario.
#   The simulation script itself is not in the package tree.
#
# example_anno
#   Sample annotation table (iid, group, lat, lon, ...) for the populations
#   referenced by example_f2_blocks; used by plot_graph_map() examples.
