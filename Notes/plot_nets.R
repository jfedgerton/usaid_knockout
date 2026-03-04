library("igraph")
library("RColorBrewer")
library("plyr")
library("dplyr")
library("fixest")

if (!grepl("Ecosystem_of_Aid", getwd()))
{
  setwd("Ecosystem_of_Aid")
}


get_initials <- function(name) {
  sapply(strsplit(name, " "), function(x) paste(substr(x, 1, 1), collapse = ""))
}

assign_colors_by_orgtype <- function(orgtypes) {
  unique_orgs <- unique(orgtypes)
  palette <- c(brewer.pal(min(length(unique_orgs), 8), "Dark2"), rep("black", 10))
  colors <- setNames(palette, unique_orgs)
  return(colors[orgtypes])
}

beautify_network_plot <- function(g, plot_title = NULL) {
  layout <- layout_with_fr(g)
  V(g)$label <- get_initials(V(g)$name)
  V(g)$size <- 5 + degree(g, mode = "out") / max(degree(g, mode = "out")) * 5
  V(g)$label.cex <- 0.7
  V(g)$label.color <- "black"
  V(g)$color <- assign_colors_by_orgtype(V(g)$orgtype)
  E(g)$color <- rgb(0.2, 0.2, 0.2, 0.4)
  E(g)$width <- 0.5
  plot(g, layout = layout,
       vertex.frame.color = "white",
       main = plot_title,
       vertex.label.dist = 0.2)
}

all_nets <- list.files(path = "Data/clean_nets", full.names = T, pattern = "Rds")
all_nets <- all_nets[grepl("crs", all_nets)]
network_attributes_list <- list()
node_attributes_list <- list()

for (z in seq_along(all_nets)) {
  network_data <- readRDS(all_nets[z])
  networks <- network_data$networks
  
  for (i in seq_along(networks)) {
    g <- networks[[i]]
    
    # Network-level metrics for directed network
    network_attributes <- data.frame(
      network_name = paste0(all_nets[z], "_network_", i),
      num_nodes = vcount(g),
      num_edges = ecount(g),
      density = edge_density(g, loops = FALSE),
      diameter = diameter(g, directed = TRUE, unconnected = TRUE),
      flow_hierarchy = 1 - reciprocity(g, mode = "ratio"),
      avg_in_degree = mean(degree(g, mode = "in")),
      avg_out_degree = mean(degree(g, mode = "out")),
      clustering_coef = transitivity(g, type = "average"),
      reciprocity = reciprocity(g),
      mean_distance = mean_distance(g, directed = TRUE),
      largest_component_size = max(components(g, mode = "strong")$csize),
      num_components = components(g, mode = "strong")$no,
      degree_centralization = centr_degree(g, mode = "out")$centralization,
      betweenness_centralization = centr_betw(g, directed = TRUE)$centralization,
      edge_connectivity = edge_connectivity(g),
      vertex_connectivity = vertex_connectivity(g),
      year = as.numeric(gsub("[^0-9]", "", names(networks))[i])
    )
    
    network_attributes_list[[length(network_attributes_list) + 1]] <- network_attributes
    
    hits <- hits_scores(g)
    
    suppressWarnings({eigen_centrality_safe <- tryCatch(
      eigen_centrality(g, directed = TRUE, scale = TRUE)$vector,
      error = function(e) rep(0, vcount(g))
    )})
    
    local_clustering_safe <- transitivity(g, type = "local")
    local_clustering_safe[is.na(local_clustering_safe)] <- 0
    
    eccentricity_safe <- tryCatch(
      eccentricity(g, mode = "out"),
      error = function(e) rep(0, vcount(g))
    )
    
    state_attr <- if("state" %in% vertex_attr_names(g)) V(g)$state else rep(NA, vcount(g))
    orgtype_attr <- if("orgtype" %in% vertex_attr_names(g)) V(g)$orgtype else rep(NA, vcount(g))
    
    node_attributes <- data.frame(
      network_name = rep(paste0(all_nets[z], "_network_", i), vcount(g)),
      node_name = V(g)$name,
      orgtype = orgtype_attr,
      state = state_attr, 
      in_degree = degree(g, mode = "in"),
      out_degree = degree(g, mode = "out"),
      total_degree = degree(g, mode = "all"),
      betweenness = betweenness(g, directed = TRUE, normalized = TRUE),
      closeness = closeness(g, mode = "out", normalized = TRUE),
      pagerank = page_rank(g, directed = TRUE)$vector,
      authority_score = hits$authority,
      hub_score = hits$hub,
      eigenvector_centrality = eigen_centrality_safe,
      local_clustering = local_clustering_safe,
      eccentricity = eccentricity_safe,
      flow_balance <- degree(g, mode = "out") - degree(g, mode = "in"),
      year = as.numeric(gsub("[^0-9]", "", names(networks))[i])
    )
    
    node_attributes_list[[length(node_attributes_list) + 1]] <- node_attributes
    
    beautify_network_plot(g, plot_title = paste0(toupper(gsub("_result_|\\.Rds", " ", all_nets[z])), " ", gsub("[^0-9]", "", names(networks))[i]))
    Sys.sleep(1)
  }
}

# Combine data frames
network_attributes_df <- do.call(rbind, network_attributes_list)
node_attributes_df <- do.call(rbind, node_attributes_list)
network_attributes_df <- network_attributes_df %>% 
  mutate(state = gsub("Data/clean_nets/", "", network_name),
         state = gsub("_result_crs.Rds_network_", "", state),
         state = gsub("\\d", "", state))

node_attributes_df <- node_attributes_df %>% 
  mutate(state = gsub("Data/clean_nets/", "", network_name),
         state = gsub("_result_crs.Rds_network_", "", state),
         state = gsub("\\d", "", state))
# Save results
saveRDS(network_attributes_df, file = "Data/clean_nets/network_attributes_directed.rds")
saveRDS(node_attributes_df, file = "Data/clean_nets/node_attributes_directed.rds")

network_attributes_df <- readRDS("Data/clean_nets/network_attributes_directed.rds")
node_attributes_df    <- readRDS("Data/clean_nets/node_attributes_directed.rds")

str(network_attributes_df)
str(node_attributes_df)
unique(network_attributes_df$state)
unique(node_attributes_df$year)

# Define treatment states for Haiti and Ebola crises
haiti_treated_states <- c("hti", "car", "mali", "sdn", "som")  # based on earlier examples
haiti_treatment_year <- 2010


ebola_treated_states <- c("lbn", "lib", "mali")
ebola_treatment_year <- 2014

# Label network_attributes_df
network_attributes_df$haiti_treated <- ifelse(network_attributes_df$state %in% haiti_treated_states & network_attributes_df$year >= haiti_treatment_year, 1, 0)
network_attributes_df$ebola_treated <- ifelse(network_attributes_df$state %in% ebola_treated_states & network_attributes_df$year >= ebola_treatment_year, 1, 0)
network_attributes_df$post_haiti <- ifelse(network_attributes_df$year >= haiti_treatment_year, 1, 0)
network_attributes_df$post_ebola <- ifelse(network_attributes_df$year >= ebola_treatment_year, 1, 0)

# Label node_attributes_df
node_attributes_df$haiti_treated <- ifelse(node_attributes_df$state %in% haiti_treated_states & node_attributes_df$year >= haiti_treatment_year, 1, 0)
node_attributes_df$ebola_treated <- ifelse(node_attributes_df$state %in% ebola_treated_states & node_attributes_df$year >= ebola_treatment_year, 1, 0)
node_attributes_df$post_haiti <- ifelse(node_attributes_df$year >= haiti_treatment_year, 1, 0)
node_attributes_df$post_ebola <- ifelse(node_attributes_df$year >= ebola_treatment_year, 1, 0)

# Verify new columns
head(network_attributes_df)
head(node_attributes_df)

# DiD for Haiti Earthquake (assuming your dataset is named 'df')
net_model_haiti <- feols(degree_centralization ~ haiti_treated:post_haiti | state + year, 
                     data = 
                       subset(network_attributes_df, year %in% (haiti_treatment_year - 3):(haiti_treatment_year + 2))
                     , cluster = ~state + year)
summary(net_model_haiti)

# DiD for Ebola Crisis
net_model_ebola <- feols(degree_centralization ~ ebola_treated:post_ebola | state + year, 
                     data = 
                       subset(network_attributes_df, year %in% (ebola_treatment_year - 3):(ebola_treatment_year + 2))
                     , cluster = ~state + year)
summary(net_model_ebola)

# DiD for Haiti Earthquake (assuming your dataset is named 'df')
node_model_haiti <- feols(out_degree ~ haiti_treated:post_haiti | state + year, 
                         data = 
                           subset(node_attributes_df, year %in% (haiti_treatment_year - 3):(haiti_treatment_year + 2))
                         , cluster = ~state + year)
summary(node_model_haiti)

# DiD for Ebola Crisis
node_model_ebola <- feols(betweenness ~ ebola_treated:post_ebola | state + year, 
                     data = 
                       subset(node_attributes_df, year %in% (ebola_treatment_year - 3):(ebola_treatment_year + 2))
                     , cluster = ~state + year)
summary(node_model_ebola)

