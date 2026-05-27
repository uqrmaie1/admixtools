# Canonical 5-node, 1-admixture graph for unit tests.
# See LLD §9.1 for the topology diagram.
make_minimal_graph <- function() {
  tibble::tribble(
    ~from, ~to,  ~type,   ~weight,
    "R",   "A",  "normal",  0.05,
    "R",   "B",  "normal",  0.05,
    "A",   "M",  "admix", 0.60,
    "B",   "M",  "admix", 0.40,
    "M",   "X",  "normal",  0.02
  )
}

# Same graph as an igraph (for tests covering the igraph input path).
make_minimal_igraph <- function() {
  edges <- make_minimal_graph()
  ig <- igraph::graph_from_edgelist(as.matrix(edges[, c("from", "to")]))
  igraph::edge_attr(ig, "type")   <- edges$type
  igraph::edge_attr(ig, "weight") <- edges$weight
  ig
}

# Topology from rha20.lgo (Rogers et al. 2020, Sci Adv,
# doi:10.1126/sciadv.aay5483). Bundled with LEGOFIT under ISC license
# at src/rha20.lgo, citation embedded in the file as Li:N-505-43-S88.
#
# 18 nodes, 13 normal edges, 4 admixture events (mN, mS, mXY, mSND).
# Internal segments that were `samples=1` in LEGOFIT (n, v, a, y2, s,
# etc.) lose their sample status here — admixtools' edge tibble does
# not model sampled internal nodes. The 3 admixtools leaves of this
# topology are x, y, d.
#
# Branch lengths (`time` column) are derived from the published
# LEGOFIT-fitted absolute times in rha20.lgo: for each edge u->v, the
# branch length is t(u) - t(v) where t(*) are the rha20.lgo values
# (Txynds, Txynd, Tnd, Tav, Txy, Td, Ta, Tv, TmN, plus the constrained
# variables TmXY=0.5*(Txy+Tav), TmS=0.5*(Td+Tnd), TmSND=0.5*(Tnd+Txynd)).
# These satisfy the topology_walk consistency check.
#
# Use with: time_handling = "init", dates_terminal = c(x = 0, y = 0,
#                                                      d = 3484.25)
# to reconstruct the published times exactly.
make_rogers2020_graph <- function() {
  # Mixfracs from rha20.lgo (free-variable initial values)
  mS   <- 0.0190121
  mN   <- 0.0192768
  mXY  <- 0.015536
  mSND <- 0.0343425
  # Constrained-variable values (computed from the free variables)
  Td    <- 3484.25
  Tnd   <- 25416.9
  Txynd <- 25920    # fixed
  Tav   <- 16307
  Txy   <- 774.856
  TmN   <- 1
  TmXY  <- 0.5 * (Txy + Tav)        # 8540.928
  TmS   <- 0.5 * (Td + Tnd)         # 14450.575
  TmSND <- 0.5 * (Tnd + Txynd)      # 25668.45
  tibble::tribble(
    ~from,    ~to,    ~type,    ~weight,  ~time,
    # Normal edges (from `derive` statements in rha20.lgo)
    "v",      "n",    "normal", NA_real_, 2511.5 - TmN,         # Tv - TmN
    "a",      "v",    "normal", NA_real_, 5307.02 - 2511.5,     # Ta - Tv
    "av",     "a2",   "normal", NA_real_, Tav - TmXY,
    "nd",     "av",   "normal", NA_real_, Tnd - Tav,
    "nd",     "d2",   "normal", NA_real_, Tnd - TmS,
    "xy",     "x",    "normal", NA_real_, Txy - 0,
    "xy",     "y2",   "normal", NA_real_, Txy - TmN,
    "xy2",    "xy",   "normal", NA_real_, TmXY - Txy,
    "xynd",   "xy2",  "normal", NA_real_, Txynd - TmXY,
    "xynd",   "nd2",  "normal", NA_real_, Txynd - TmSND,
    "s2",     "s",    "normal", NA_real_, TmSND - TmS,
    "xynds",  "s2",   "normal", NA_real_, 82008.2 - TmSND,      # Txynds - TmSND
    "xynds",  "xynd", "normal", NA_real_, 82008.2 - Txynd,      # Txynds - Txynd
    # Admix edges: weight = admixture proportion; time = branch length
    # to the admix event. Both parents of each admix child sit at the
    # admix event's height.
    "d2",     "d",    "admix",  1 - mS,   TmS - Td,
    "s",      "d",    "admix",  mS,       TmS - Td,
    "y2",     "y",    "admix",  1 - mN,   TmN - 0,
    "n",      "y",    "admix",  mN,       TmN - 0,
    "a2",     "a",    "admix",  1 - mXY,  TmXY - 5307.02,        # TmXY - Ta
    "xy2",    "a",    "admix",  mXY,      TmXY - 5307.02,
    "nd2",    "nd",   "admix",  1 - mSND, TmSND - Tnd,
    "s2",     "nd",   "admix",  mSND,     TmSND - Tnd
  )
}

# 7-node graph with a bifurcating root and an outgroup leaf.
# Topology:
#   Root → Outgroup        (the outgroup leaf)
#   Root → IngroupRoot     (the new root after stripping)
#     IngroupRoot → A
#     IngroupRoot → B
#     A -.admix.-> M
#     B -.admix.-> M
#     M → X
make_minimal_graph_with_outgroup <- function() {
  tibble::tribble(
    ~from,         ~to,            ~type,    ~weight,
    "Root",        "Outgroup",     "normal",   0.10,
    "Root",        "IngroupRoot",  "normal",   0.10,
    "IngroupRoot", "A",            "normal",   0.05,
    "IngroupRoot", "B",            "normal",   0.05,
    "A",           "M",            "admix",  0.60,
    "B",           "M",            "admix",  0.40,
    "M",           "X",            "normal",   0.02
  )
}


# Singular-case fixture builder for qpadm rcond / loadings tests.
#
# Adds a clone of `src` to a 3d f2_blocks array as a new dimname `dst`.
# The cloned dim is row/col-for-row identical to src's row/col, which makes
# the unfudged f4_var rank-deficient (perfectly collinear). Combined with
# `fudge = 1e-12` on qpadm() this drives `f4_var_rcond` below the 1e-8
# auto-bar and triggers the loadings tibble. Originally lived in
# test-qpadm-singular-threshold.R; promoted here so test-qpadm_sweep.R can
# exercise the same singular path without duplicating the helper.
.clone_pop_in_f2_blocks = function(f2, src = "Mbuti.DG", dst = "Mbuti_clone") {
  if(length(dim(f2)) != 3L)
    stop("f2 must be a 3D array (n_pops x n_pops x n_blocks); got ndim = ",
         length(dim(f2)))
  nam = dimnames(f2)[[1]]
  if(!(src %in% nam)) stop("source pop not in f2_blocks: ", src)
  if(dst %in% nam)    stop("destination pop already exists: ", dst)
  n   = length(nam)
  nb  = dim(f2)[3]
  new_nam = c(nam, dst)
  out = array(NA_real_, dim = c(n + 1, n + 1, nb),
              dimnames = list(new_nam, new_nam, dimnames(f2)[[3]]))
  out[1:n, 1:n, ] = f2
  src_idx = match(src, nam)
  out[1:n,   n+1, ] = f2[1:n, src_idx, ]
  out[n+1,   1:n, ] = f2[src_idx, 1:n, ]
  out[n+1,   n+1, ] = 0
  out
}
