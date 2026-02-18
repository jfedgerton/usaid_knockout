# ==========================================================
# Build network + node panels, then run Ebola DiD
# with node-level exposure controls
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

# ----------------------------------------------------------
# Helpers
# ----------------------------------------------------------
jaccard_sets <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)  # stable: empty-empty
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

get_usaid_nodes <- function(g, regex = USAID_REGEX) {
  if (!inherits(g, "igraph") || vcount(g) == 0) return(character(0))
  nm <- V(g)$name
  if (is.null(nm)) return(character(0))
  nm <- as.character(nm)
  nm[grepl(regex, nm, ignore.case = TRUE)]
}

calc_fragmentation <- function(g) {
  if (!inherits(g, "igraph") || vcount(g) == 0) return(1)
  cc <- components(g, mode = "weak")
  1 - (max(cc$csize) / vcount(g))
}

calc_weighted_fragmentation <- function(g) {
  if (!inherits(g, "igraph") || vcount(g) == 0) return(1)
  cc <- components(g, mode = "weak")
  s <- cc$csize
  n <- vcount(g)
  1 - sum((s / n)^2)
}

calc_modularity_louvain <- function(g) {
  if (!inherits(g, "igraph") || vcount(g) < 2 || ecount(g) == 0) return(NA_real_)
  gu <- as_undirected(g, mode = "collapse")
  cl <- cluster_louvain(gu)
  modularity(cl)
}

partner_dependency_metrics <- function(g, donor_nodes) {
  donor_nodes <- intersect(donor_nodes, V(g)$name)
  if (length(donor_nodes) == 0) {
    return(list(n_partners = 0L,
                share_no_alt_funder = NA_real_,
                mean_partner_total_degree = NA_real_,
                mean_partner_in_degree = NA_real_,
                mean_partner_out_degree = NA_real_))
  }
  
  partners <- unique(unlist(lapply(donor_nodes, function(d) {
    V(g)[neighbors(g, d, mode = "out")]$name
  })))
  partners <- setdiff(partners, donor_nodes)
  
  if (length(partners) == 0) {
    return(list(n_partners = 0L,
                share_no_alt_funder = NA_real_,
                mean_partner_total_degree = NA_real_,
                mean_partner_in_degree = NA_real_,
                mean_partner_out_degree = NA_real_))
  }
  
  in_deg  <- degree(g, mode = "in");  names(in_deg)  <- V(g)$name
  out_deg <- degree(g, mode = "out"); names(out_deg) <- V(g)$name
  tot_deg <- degree(g, mode = "all"); names(tot_deg) <- V(g)$name
  
  alt_n <- vapply(partners, function(p) {
    funders <- V(g)[neighbors(g, p, mode = "in")]$name
    funders <- setdiff(funders, donor_nodes)
    length(funders)
  }, integer(1))
  
  list(
    n_partners = length(partners),
    share_no_alt_funder = mean(alt_n == 0),
    mean_partner_total_degree = mean(tot_deg[partners]),
    mean_partner_in_degree = mean(in_deg[partners]),
    mean_partner_out_degree = mean(out_deg[partners])
  )
}

cascade_drop_unfunded <- function(g, removed_seed) {
  removed <- unique(intersect(removed_seed, V(g)$name))
  if (length(removed) == 0) return(character(0))
  
  repeat {
    remaining <- setdiff(V(g)$name, removed)
    if (length(remaining) == 0) break
    
    new_fail <- character(0)
    for (v in remaining) {
      funders <- V(g)[neighbors(g, v, mode = "in")]$name
      if (length(funders) > 0 && all(funders %in% removed)) {
        new_fail <- c(new_fail, v)
      }
    }
    new_fail <- setdiff(unique(new_fail), removed)
    if (length(new_fail) == 0) break
    removed <- c(removed, new_fail)
  }
  
  unique(removed)
}

simulate_donor_withdrawal <- function(g, donor_nodes, scenario = c("direct","drop_unfunded")) {
  scenario <- match.arg(scenario)
  donor_nodes <- intersect(donor_nodes, V(g)$name)
  
  frag_pre  <- calc_fragmentation(g)
  wfrag_pre <- calc_weighted_fragmentation(g)
  cc_pre <- components(g, mode = "weak")
  comps_pre <- cc_pre$no
  giant_pre <- max(cc_pre$csize) / vcount(g)
  
  if (length(donor_nodes) == 0) {
    return(list(frag_delta = NA_real_,
                wfrag_delta = NA_real_,
                component_increase = NA_integer_,
                giant_share_delta = NA_real_,
                n_removed_total = NA_integer_))
  }
  
  removed <- donor_nodes
  if (scenario == "drop_unfunded") removed <- cascade_drop_unfunded(g, donor_nodes)
  
  g2 <- delete_vertices(g, removed)
  
  if (vcount(g2) == 0) {
    frag_post <- 1; wfrag_post <- 1; comps_post <- 0; giant_post <- 0
  } else {
    frag_post  <- calc_fragmentation(g2)
    wfrag_post <- calc_weighted_fragmentation(g2)
    cc_post <- components(g2, mode = "weak")
    comps_post <- cc_post$no
    giant_post <- max(cc_post$csize) / vcount(g2)
  }
  
  list(
    frag_delta = frag_post - frag_pre,
    wfrag_delta = wfrag_post - wfrag_pre,
    component_increase = comps_post - comps_pre,
    giant_share_delta = giant_post - giant_pre,
    n_removed_total = length(removed)
  )
}

# ----------------------------------------------------------
# Load networks
# ----------------------------------------------------------
all_nets <- list.files(path = "Data/clean_nets", full.names = TRUE, pattern = "Rds")
all_nets <- all_nets[grepl("crs", all_nets)]

network_attributes_list <- list()
node_attributes_list <- list()

for (z in seq_along(all_nets)) {
  
  network_data <- readRDS(all_nets[z])
  networks <- network_data$networks
  
  # state abbrev inferred from first network name in this file
  index_name <- gsub("_(.*)", "", names(networks)[1])
  
  for (i in seq_along(networks)) {
    
    g <- networks[[i]]
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    
    if (is.null(V(g)$name)) V(g)$name <- as.character(seq_len(vcount(g)))
    
    yr <- as.numeric(gsub("[^0-9]", "", names(networks))[i])
    
    # USAID flags
    is_usaid_vec <- as.integer(grepl(USAID_REGEX, V(g)$name, ignore.case = TRUE))
    usaid_nodes <- get_usaid_nodes(g)
    
    # Degrees
    in_deg  <- degree(g, mode = "in")
    out_deg <- degree(g, mode = "out")
    tot_deg <- degree(g, mode = "all")
    
    # Lagged Jaccard (only if prior network exists)
    jaccard_df <- NULL
    if (i != 1) {
      g_lag <- networks[[i - 1]]
      g_lag <- igraph::simplify(g_lag, remove.multiple = TRUE, remove.loops = TRUE)
      if (is.null(V(g_lag)$name)) V(g_lag)$name <- as.character(seq_len(vcount(g_lag)))
      
      j_out <- all_nodes_jaccard_lag1(g_lag, g, mode = "out")
      j_in  <- all_nodes_jaccard_lag1(g_lag, g, mode = "in")
      
      jaccard_df <- data.frame(
        node_name = names(j_in),
        j_in = as.numeric(j_in),
        j_out = as.numeric(j_out[names(j_in)]),
        stringsAsFactors = FALSE
      )
    }
    
    # --------------------------
    # Network-level metrics
    # --------------------------
    network_attributes <- data.frame(
      network_name = paste0(all_nets[z], "_network_", i),
      state = index_name,
      year = yr,
      
      num_nodes = vcount(g),
      num_edges = ecount(g),
      density = edge_density(g, loops = FALSE),
      diameter = if (vcount(g) > 0) diameter(g, directed = TRUE, unconnected = TRUE) else NA_real_,
      flow_hierarchy = if (ecount(g) > 0) 1 - reciprocity(g, mode = "ratio") else NA_real_,
      avg_in_degree = mean(in_deg),
      avg_out_degree = mean(out_deg),
      clustering_coef = if (ecount(g) > 0) transitivity(g, type = "average") else NA_real_,
      reciprocity = if (ecount(g) > 0) reciprocity(g) else NA_real_,
      mean_distance = tryCatch(mean_distance(g, directed = TRUE, unconnected = TRUE), error = function(e) NA_real_),
      largest_component_size = if (vcount(g) > 0) max(components(g, mode = "strong")$csize) else 0,
      indegree_centralization = centr_degree(g, mode = "in")$centralization,
      outdegree_centralization = centr_degree(g, mode = "out")$centralization,
      betweenness_centralization = centr_betw(g, directed = TRUE)$centralization,
      num_components = components(g, mode = "strong")$no,
      edge_connectivity = edge_connectivity(g),
      vertex_connectivity = vertex_connectivity(g),
      stringsAsFactors = FALSE
    )
    
    # Efficiency
    if (vcount(g) <= 1) {
      network_attributes$network_efficiency <- NA_real_
    } else {
      dist_mat <- distances(g, mode = "out")
      diag(dist_mat) <- Inf
      fin <- is.finite(dist_mat)
      network_attributes$network_efficiency <- if (any(fin)) mean(1 / dist_mat[fin]) else NA_real_
    }
    
    # Funding redundancy (system-wide)
    k_rec <- in_deg[in_deg >= 1]
    network_attributes$funding_redundancy <- if (length(k_rec) > 0) 1 - mean(1 / k_rec) else NA_real_
    
    # Donor concentration
    if (sum(out_deg) > 0) {
      p_out <- out_deg / sum(out_deg)
      network_attributes$donor_concentration <- sum(p_out^2)
    } else {
      network_attributes$donor_concentration <- NA_real_
    }
    
    # Simulation-aligned DVs + mechanism
    network_attributes$fragmentation <- calc_fragmentation(g)
    network_attributes$weighted_fragmentation <- calc_weighted_fragmentation(g)
    network_attributes$modularity <- calc_modularity_louvain(g)
    
    dep <- partner_dependency_metrics(g, usaid_nodes)
    network_attributes$usaid_present <- as.integer(length(usaid_nodes) > 0)
    network_attributes$usaid_n_nodes <- length(usaid_nodes)
    network_attributes$usaid_n_partners <- dep$n_partners
    network_attributes$usaid_share_partners_no_alt_funder <- dep$share_no_alt_funder
    network_attributes$usaid_mean_partner_total_degree <- dep$mean_partner_total_degree
    
    sim_direct <- simulate_donor_withdrawal(g, usaid_nodes, "direct")
    sim_drop   <- simulate_donor_withdrawal(g, usaid_nodes, "drop_unfunded")
    
    network_attributes$usaid_frag_delta_direct <- sim_direct$frag_delta
    network_attributes$usaid_frag_delta_drop_unfunded <- sim_drop$frag_delta
    
    network_attributes_list[[length(network_attributes_list) + 1]] <- network_attributes
    
    # --------------------------
    # Node-level attributes
    # --------------------------
    hits <- hits_scores(g)
    eigen_centrality_safe <- tryCatch(
      eigen_centrality(g, directed = TRUE)$vector,
      error = function(e) rep(0, vcount(g))
    )
    
    local_clustering_safe <- transitivity(g, type = "local")
    local_clustering_safe[is.na(local_clustering_safe)] <- 0
    
    eccentricity_safe <- tryCatch(eccentricity(g, mode = "out"), error = function(e) rep(0, vcount(g)))
    
    # tied_to_usaid: neighbor of USAID
    if (length(usaid_nodes) > 0) {
      usaid_neighbors <- unique(unlist(lapply(usaid_nodes, function(x) {
        c(neighbors(g, x, mode = "out")$name, neighbors(g, x, mode = "in")$name)
      })))
      tied_to_usaid <- as.integer(V(g)$name %in% usaid_neighbors)
    } else {
      tied_to_usaid <- rep(0L, vcount(g))
    }
    
    node_attributes <- data.frame(
      network_name = rep(paste0(all_nets[z], "_network_", i), vcount(g)),
      node_name = V(g)$name,
      state = rep(index_name, vcount(g)),
      year = yr,
      
      is_usaid = is_usaid_vec,
      tied_to_usaid = tied_to_usaid,
      
      in_degree = in_deg,
      out_degree = out_deg,
      total_degree = tot_deg,
      
      betweenness = betweenness(g, directed = TRUE, normalized = TRUE),
      closeness = closeness(g, mode = "out", normalized = TRUE),
      pagerank = page_rank(g, directed = TRUE)$vector,
      authority_score = hits$authority,
      hub_score = hits$hub,
      eigenvector_centrality = eigen_centrality_safe,
      local_clustering = local_clustering_safe,
      eccentricity = eccentricity_safe,
      flow_balance = out_deg - in_deg,
      stringsAsFactors = FALSE
    )
    
    # Join Jaccard if available; else set missing
    if (!is.null(jaccard_df)) {
      node_attributes <- left_join(node_attributes, jaccard_df, by = "node_name")
    } else {
      node_attributes$j_in <- NA_real_
      node_attributes$j_out <- NA_real_
    }
    
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
