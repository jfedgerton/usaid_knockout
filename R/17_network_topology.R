# ==========================================================
# 17_network_topology.R
# Build network-level and node-level panels
#
# CHANGES:
#   1. Node-level vars for node DiD: usaid_funded, sole_usaid_funded,
#      n_funders, n_non_usaid_funders, usaid_exposure_share, distance_to_usaid
#   2. Bilateral donor partner dependency metrics
#   3. sole_funder_rate S(d) per Eq. 1
#   4. Density retained as descriptive only
#   5. Reads from harmonized networks (entity-resolved names)
#   6. Cascade simulations removed (DVs: fragmentation, weighted_frag,
#      node-level measures only)
# ==========================================================

library(igraph)
library(RColorBrewer)
library(dplyr)
library(fixest)
library(ggplot2)
library(purrr)

if (!grepl("Ecosystem_of_Aid", getwd())) setwd("Ecosystem_of_Aid")

# ----------------------------------------------------------
# Constants
# ----------------------------------------------------------
USAID_REGEX <- "United States Agency for International Development|\\bUSAID\\b|US AID"

BILATERAL_DONORS <- list(
  usaid  = "United States Agency for International Development|\\bUSAID\\b|US AID",
  dfid   = "Department for International Development|\\bDFID\\b|\\bFCDO\\b|Foreign.*Commonwealth.*Development",
  giz    = "Deutsche Gesellschaft|\\bGIZ\\b|\\bBMZ\\b|German Federal Ministry|Federal Ministry for Economic Cooperation.*Germany",
  jica   = "Japan(ese)? International Cooperation|\\bJICA\\b",
  afd    = "Agence Fran.aise de D.veloppement|\\bAFD\\b|French Development Agency",
  sida   = "Swedish International Development|\\bSida\\b|\\bSIDA\\b",
  cida   = "Canadian International Development|\\bCIDA\\b|Global Affairs Canada",
  norad  = "Norwegian Agency for Development|\\bNorad\\b|\\bNORAD\\b"
)

# ----------------------------------------------------------
# Helpers
# ----------------------------------------------------------
jaccard_sets <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

all_nodes_jaccard_lag1 <- function(g_lag, g, mode = c("out","in","all")) {
  mode <- match.arg(mode)
  common <- intersect(V(g_lag)$name, V(g)$name)
  out <- setNames(rep(NA_real_, length(common)), common)
  for (v in common) {
    n_lag <- neighbors(g_lag, v, mode = mode)$name
    n_cur <- neighbors(g,     v, mode = mode)$name
    out[v] <- jaccard_sets(n_lag, n_cur)
  }
  out
}

get_donor_nodes <- function(g, regex) {
  if (!inherits(g, "igraph") || vcount(g) == 0) return(character(0))
  nm <- V(g)$name
  if (is.null(nm)) return(character(0))
  nm <- as.character(nm)
  nm[grepl(regex, nm, ignore.case = TRUE)]
}

get_usaid_nodes <- function(g, regex = USAID_REGEX) {
  get_donor_nodes(g, regex)
}

calc_fragmentation <- function(g) {
  if (!inherits(g, "igraph") || vcount(g) == 0) return(1)
  cc <- components(g, mode = "weak")
  1 - (max(cc$csize) / vcount(g))
}

calc_weighted_fragmentation <- function(g) {
  if (!inherits(g, "igraph") || vcount(g) == 0) return(1)
  cc <- components(g, mode = "weak")
  s <- cc$csize; n <- vcount(g)
  1 - sum((s / n)^2)
}

calc_modularity_louvain <- function(g) {
  if (!inherits(g, "igraph") || vcount(g) < 2 || ecount(g) == 0) return(NA_real_)
  gu <- as.undirected(g, mode = "collapse")
  cl <- cluster_louvain(gu)
  modularity(cl)
}

partner_dependency_metrics <- function(g, donor_nodes) {
  donor_nodes <- intersect(donor_nodes, V(g)$name)
  if (length(donor_nodes) == 0) {
    return(list(n_partners = 0L, share_no_alt_funder = NA_real_,
                sole_funder_rate = NA_real_, mean_partner_total_degree = NA_real_,
                mean_partner_in_degree = NA_real_, mean_partner_out_degree = NA_real_))
  }
  partners <- unique(unlist(lapply(donor_nodes, function(d) V(g)[neighbors(g, d, mode = "out")]$name)))
  partners <- setdiff(partners, donor_nodes)
  if (length(partners) == 0) {
    return(list(n_partners = 0L, share_no_alt_funder = NA_real_,
                sole_funder_rate = NA_real_, mean_partner_total_degree = NA_real_,
                mean_partner_in_degree = NA_real_, mean_partner_out_degree = NA_real_))
  }
  in_deg  <- degree(g, mode = "in");  names(in_deg)  <- V(g)$name
  out_deg <- degree(g, mode = "out"); names(out_deg) <- V(g)$name
  tot_deg <- degree(g, mode = "all"); names(tot_deg) <- V(g)$name
  alt_n <- vapply(partners, function(p) {
    funders <- V(g)[neighbors(g, p, mode = "in")]$name
    length(setdiff(funders, donor_nodes))
  }, integer(1))
  list(n_partners = length(partners), share_no_alt_funder = mean(alt_n == 0),
       sole_funder_rate = mean(alt_n == 0),
       mean_partner_total_degree = mean(tot_deg[partners]),
       mean_partner_in_degree = mean(in_deg[partners]),
       mean_partner_out_degree = mean(out_deg[partners]))
}

# ----------------------------------------------------------
# Load networks
# ----------------------------------------------------------
# Read from harmonized networks if available, otherwise fall back to clean_nets
input_dir <- if (dir.exists("Data/clean_nets_harmonized")) "Data/clean_nets_harmonized" else "Data/clean_nets"
cat(sprintf("Reading networks from: %s\n", input_dir))
all_nets <- list.files(path = input_dir, full.names = TRUE, pattern = "\\.(rds|Rds)$")
all_nets <- all_nets[grepl("crs", all_nets)]

network_attributes_list <- list()
node_attributes_list <- list()

for (z in seq_along(all_nets)) {
  network_data <- readRDS(all_nets[z])
  networks <- network_data$networks
  index_name <- tolower(gsub("_(.*)", "", names(networks)[1]))

  for (i in seq_along(networks)) {
    g <- networks[[i]]
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    if (is.null(V(g)$name)) V(g)$name <- as.character(seq_len(vcount(g)))
    yr <- as.numeric(gsub("[^0-9]", "", names(networks))[i])

    is_usaid_vec <- as.integer(grepl(USAID_REGEX, V(g)$name, ignore.case = TRUE))
    usaid_nodes <- get_usaid_nodes(g)
    in_deg  <- degree(g, mode = "in")
    out_deg <- degree(g, mode = "out")
    tot_deg <- degree(g, mode = "all")

    # Lagged Jaccard
    jaccard_df <- NULL
    if (i != 1) {
      g_lag <- networks[[i - 1]]
      g_lag <- igraph::simplify(g_lag, remove.multiple = TRUE, remove.loops = TRUE)
      if (is.null(V(g_lag)$name)) V(g_lag)$name <- as.character(seq_len(vcount(g_lag)))
      j_out <- all_nodes_jaccard_lag1(g_lag, g, mode = "out")
      j_in  <- all_nodes_jaccard_lag1(g_lag, g, mode = "in")
      jaccard_df <- data.frame(node_name = names(j_in), j_in = as.numeric(j_in),
                                j_out = as.numeric(j_out[names(j_in)]), stringsAsFactors = FALSE)
    }

    # ==========================================================
    # NETWORK-LEVEL METRICS
    # ==========================================================
    network_attributes <- data.frame(
      network_name = paste0(all_nets[z], "_network_", i),
      state = index_name, year = yr,
      num_nodes = vcount(g), num_edges = ecount(g),
      density = edge_density(g, loops = FALSE),
      diameter = if (vcount(g) > 0) diameter(g, directed=TRUE, unconnected=TRUE) else NA_real_,
      flow_hierarchy = if (ecount(g) > 0) 1 - reciprocity(g, mode="ratio") else NA_real_,
      avg_in_degree = mean(in_deg), avg_out_degree = mean(out_deg),
      clustering_coef = if (ecount(g) > 0) transitivity(g, type="average") else NA_real_,
      reciprocity = if (ecount(g) > 0) reciprocity(g) else NA_real_,
      mean_distance = tryCatch(mean_distance(g, directed=TRUE, unconnected=TRUE), error=function(e) NA_real_),
      largest_component_size = if (vcount(g) > 0) max(components(g, mode="strong")$csize) else 0,
      indegree_centralization = centr_degree(g, mode="in")$centralization,
      outdegree_centralization = centr_degree(g, mode="out")$centralization,
      betweenness_centralization = centr_betw(g, directed=TRUE)$centralization,
      num_components = components(g, mode="strong")$no,
      edge_connectivity = edge_connectivity(g),
      vertex_connectivity = vertex_connectivity(g),
      stringsAsFactors = FALSE)

    if (vcount(g) <= 1) { network_attributes$network_efficiency <- NA_real_
    } else {
      dist_mat <- distances(g, mode="out"); diag(dist_mat) <- Inf; fin <- is.finite(dist_mat)
      network_attributes$network_efficiency <- if (any(fin)) mean(1/dist_mat[fin]) else NA_real_
    }
    k_rec <- in_deg[in_deg >= 1]
    network_attributes$funding_redundancy <- if (length(k_rec) > 0) 1 - mean(1/k_rec) else NA_real_
    if (sum(out_deg) > 0) { p_out <- out_deg/sum(out_deg); network_attributes$donor_concentration <- sum(p_out^2)
    } else { network_attributes$donor_concentration <- NA_real_ }

    network_attributes$fragmentation <- calc_fragmentation(g)
    network_attributes$weighted_fragmentation <- calc_weighted_fragmentation(g)
    network_attributes$modularity <- calc_modularity_louvain(g)

    dep <- partner_dependency_metrics(g, usaid_nodes)
    network_attributes$usaid_present <- as.integer(length(usaid_nodes) > 0)
    network_attributes$usaid_n_nodes <- length(usaid_nodes)
    network_attributes$usaid_n_partners <- dep$n_partners
    network_attributes$usaid_sole_funder_rate <- dep$sole_funder_rate
    network_attributes$usaid_share_partners_no_alt_funder <- dep$share_no_alt_funder
    network_attributes$usaid_mean_partner_total_degree <- dep$mean_partner_total_degree

    for (donor_key in names(BILATERAL_DONORS)) {
      d_nodes <- get_donor_nodes(g, BILATERAL_DONORS[[donor_key]])
      d_dep   <- partner_dependency_metrics(g, d_nodes)
      network_attributes[[paste0(donor_key,"_present")]]            <- as.integer(length(d_nodes) > 0)
      network_attributes[[paste0(donor_key,"_n_partners")]]         <- d_dep$n_partners
      network_attributes[[paste0(donor_key,"_sole_funder_rate")]]   <- d_dep$sole_funder_rate
      network_attributes[[paste0(donor_key,"_mean_partner_degree")]] <- d_dep$mean_partner_total_degree
    }

    network_attributes_list[[length(network_attributes_list) + 1]] <- network_attributes

    # ==========================================================
    # NODE-LEVEL ATTRIBUTES
    # ==========================================================
    hub_scores <- hub_score(g)$vector
    auth_scores <- authority_score(g)$vector
    eigen_centrality_safe <- tryCatch(eigen_centrality(g, directed=TRUE)$vector, error=function(e) rep(0, vcount(g)))
    local_clustering_safe <- transitivity(g, type = "local")
    local_clustering_safe[is.na(local_clustering_safe)] <- 0
    eccentricity_safe <- tryCatch(eccentricity(g, mode="out"), error=function(e) rep(0, vcount(g)))

    # tied_to_usaid (existing): any-mode neighbor of USAID
    if (length(usaid_nodes) > 0) {
      usaid_neighbors <- unique(unlist(lapply(usaid_nodes, function(x)
        c(neighbors(g, x, mode="out")$name, neighbors(g, x, mode="in")$name))))
      tied_to_usaid <- as.integer(V(g)$name %in% usaid_neighbors)
    } else {
      tied_to_usaid <- rep(0L, vcount(g))
    }

    # NEW: usaid_funded -- receives direct funding from USAID (out-neighbor)
    usaid_funded_vec <- rep(0L, vcount(g))
    if (length(usaid_nodes) > 0) {
      usaid_out_partners <- unique(unlist(lapply(usaid_nodes, function(x)
        V(g)[neighbors(g, x, mode = "out")]$name)))
      usaid_out_partners <- setdiff(usaid_out_partners, usaid_nodes)
      usaid_funded_vec <- as.integer(V(g)$name %in% usaid_out_partners)
    }

    # NEW: n_funders and n_non_usaid_funders
    n_funders_vec <- degree(g, mode = "in")
    n_non_usaid_funders_vec <- rep(NA_integer_, vcount(g))
    for (vi in seq_len(vcount(g))) {
      in_nbrs <- V(g)[neighbors(g, V(g)$name[vi], mode = "in")]$name
      n_non_usaid_funders_vec[vi] <- sum(!in_nbrs %in% usaid_nodes)
    }

    # NEW: sole_usaid_funded -- funded by USAID only, no other funders
    sole_usaid_funded_vec <- as.integer(usaid_funded_vec == 1L & n_non_usaid_funders_vec == 0L)

    # NEW: usaid_exposure_share -- fraction of neighbors that are USAID or USAID-funded
    usaid_connected_set <- unique(c(usaid_nodes, V(g)$name[usaid_funded_vec == 1L]))
    usaid_exposure_share_vec <- rep(NA_real_, vcount(g))
    for (vi in seq_len(vcount(g))) {
      all_nbrs <- V(g)[neighbors(g, V(g)$name[vi], mode = "all")]$name
      usaid_exposure_share_vec[vi] <- if (length(all_nbrs) == 0) 0 else mean(all_nbrs %in% usaid_connected_set)
    }

    # NEW: distance_to_usaid -- shortest undirected path
    distance_to_usaid_vec <- rep(NA_real_, vcount(g))
    if (length(usaid_nodes) > 0) {
      g_undir <- as.undirected(g, mode = "collapse")
      for (vi in seq_len(vcount(g))) {
        v_name <- V(g)$name[vi]
        if (v_name %in% usaid_nodes) {
          distance_to_usaid_vec[vi] <- 0
        } else {
          dists_to_usaid <- distances(g_undir, v = v_name, to = usaid_nodes)
          min_d <- min(dists_to_usaid, na.rm = TRUE)
          distance_to_usaid_vec[vi] <- if (is.finite(min_d)) min_d else NA_real_
        }
      }
    }

    # Assemble
    node_attributes <- data.frame(
      network_name = rep(paste0(all_nets[z], "_network_", i), vcount(g)),
      node_name = V(g)$name, state = rep(index_name, vcount(g)), year = yr,
      is_usaid = is_usaid_vec, tied_to_usaid = tied_to_usaid,
      usaid_funded = usaid_funded_vec, sole_usaid_funded = sole_usaid_funded_vec,
      in_degree = in_deg, out_degree = out_deg, total_degree = tot_deg,
      n_funders = n_funders_vec, n_non_usaid_funders = n_non_usaid_funders_vec,
      usaid_exposure_share = usaid_exposure_share_vec,
      distance_to_usaid = distance_to_usaid_vec,
      betweenness = betweenness(g, directed=TRUE, normalized=TRUE),
      closeness = closeness(g, mode="out", normalized=TRUE),
      pagerank = page_rank(g, directed=TRUE)$vector,
      authority_score = auth_scores, hub_score = hub_scores,
      eigenvector_centrality = eigen_centrality_safe,
      local_clustering = local_clustering_safe,
      eccentricity = eccentricity_safe,
      flow_balance = out_deg - in_deg,
      stringsAsFactors = FALSE)

    if (!is.null(jaccard_df)) {
      node_attributes <- left_join(node_attributes, jaccard_df, by = "node_name")
    } else { node_attributes$j_in <- NA_real_; node_attributes$j_out <- NA_real_ }
    node_attributes$rewire_in  <- ifelse(is.na(node_attributes$j_in),  NA_real_, 1 - node_attributes$j_in)
    node_attributes$rewire_out <- ifelse(is.na(node_attributes$j_out), NA_real_, 1 - node_attributes$j_out)

    node_attributes_list[[length(node_attributes_list) + 1]] <- node_attributes
  }
}

# ----------------------------------------------------------
# Combine and save
# ----------------------------------------------------------
network_attributes_df <- bind_rows(network_attributes_list)
node_attributes_df <- bind_rows(node_attributes_list)

saveRDS(network_attributes_df, file = "Data/clean_nets/network_attributes_directed.rds")
saveRDS(node_attributes_df, file = "Data/clean_nets/node_attributes_directed.rds")
