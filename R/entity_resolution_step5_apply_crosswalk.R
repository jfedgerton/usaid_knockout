###############################################################################
# Entity Resolution Step 5: Apply crosswalk to networks
#
# Input:  Data/entity_resolution/name_crosswalk.csv
#         Data/clean_nets/*crs*.Rds (original network files)
# Output: Data/clean_nets_harmonized/*crs*.Rds (harmonized network files)
#
# For each network:
#   1. Rename V(g)$name according to the crosswalk
#   2. If renaming creates duplicate nodes (multi-edges), merge them:
#      - Collapse duplicate nodes into one
#      - For parallel edges: keep max weight if weighted, sum otherwise
#   3. Re-save with identical list structure
#   4. Log changes per file
###############################################################################

library(igraph)

cat("=== Entity Resolution Step 5: Apply crosswalk to networks ===\n")

# --------------------------------------------------------------------------
# 1. Read crosswalk
# --------------------------------------------------------------------------

crosswalk_path <- "Data/entity_resolution/name_crosswalk.csv"

if (!file.exists(crosswalk_path)) {
  stop("ERROR: ", crosswalk_path, " not found.\n",
       "You need to complete manual review (Step 4) and produce the crosswalk first.\n",
       "The crosswalk should have columns: raw_name, canonical_name, match_type")
}

crosswalk <- read.csv(crosswalk_path, stringsAsFactors = FALSE)
cat(sprintf("Loaded crosswalk: %d mappings\n", nrow(crosswalk)))
cat(sprintf("  match_type breakdown:\n"))
for (mt in sort(unique(crosswalk$match_type))) {
  cat(sprintf("    %s: %d\n", mt, sum(crosswalk$match_type == mt)))
}

# Build lookup: raw_name -> canonical_name
name_map <- setNames(crosswalk$canonical_name, crosswalk$raw_name)

# --------------------------------------------------------------------------
# 2. Helper: apply name mapping and merge duplicate nodes
# --------------------------------------------------------------------------

harmonize_network <- function(g, name_map) {
  # Get current names
  old_names <- V(g)$name
  if (is.null(old_names)) return(list(graph = g, n_renamed = 0, n_merged = 0))

  old_names <- as.character(old_names)

  # Apply mapping (keep original if not in crosswalk)
  new_names <- ifelse(
    old_names %in% names(name_map),
    name_map[old_names],
    old_names
  )

  n_renamed <- sum(old_names != new_names)

  # Check for duplicates after renaming
  dup_names <- unique(new_names[duplicated(new_names)])
  n_merged <- length(dup_names)

  if (n_merged == 0) {
    # Simple case: just rename
    V(g)$name <- new_names
    return(list(graph = g, n_renamed = n_renamed, n_merged = 0))
  }

  # Complex case: need to merge nodes with same new name
  # Strategy: rebuild the edge list with new names, then collapse

  el <- as_edgelist(g)
  el_char <- matrix(as.character(el), ncol = 2)

  # Get edge weights if they exist
  has_weight <- "weight" %in% edge_attr_names(g)
  weights <- if (has_weight) E(g)$weight else rep(1, ecount(g))

  # Map old names to new names in edge list
  map_name <- function(nm) {
    idx <- match(nm, old_names)
    new_names[idx]
  }

  new_from <- map_name(el_char[, 1])
  new_to <- map_name(el_char[, 2])

  # Remove self-loops created by merging
  keep <- new_from != new_to
  new_from <- new_from[keep]
  new_to <- new_to[keep]
  weights <- weights[keep]

  # Aggregate parallel edges (take max weight)
  edge_key <- paste(new_from, new_to, sep = "|||")
  unique_edges <- unique(edge_key)

  agg_from <- character(length(unique_edges))
  agg_to <- character(length(unique_edges))
  agg_weight <- numeric(length(unique_edges))

  for (i in seq_along(unique_edges)) {
    idx <- which(edge_key == unique_edges[i])
    parts <- strsplit(unique_edges[i], "|||", fixed = TRUE)[[1]]
    agg_from[i] <- parts[1]
    agg_to[i] <- parts[2]
    agg_weight[i] <- max(weights[idx])  # Keep max weight
  }

  # Build new graph
  all_nodes <- sort(unique(c(agg_from, agg_to)))

  # Include isolated nodes (nodes that lost all edges due to self-loop removal)
  all_new_names <- sort(unique(new_names))
  all_nodes <- sort(unique(c(all_nodes, all_new_names)))

  if (length(agg_from) > 0) {
    new_el <- cbind(agg_from, agg_to)
    g_new <- graph_from_edgelist(new_el, directed = is_directed(g))
    # Add any isolated vertices
    missing <- setdiff(all_nodes, V(g_new)$name)
    if (length(missing) > 0) {
      g_new <- add_vertices(g_new, length(missing), name = missing)
    }
    if (has_weight) {
      E(g_new)$weight <- agg_weight
    }
  } else {
    # All edges became self-loops — graph is just isolated nodes
    g_new <- make_empty_graph(n = length(all_nodes), directed = is_directed(g))
    V(g_new)$name <- all_nodes
  }

  return(list(graph = g_new, n_renamed = n_renamed, n_merged = n_merged))
}

# --------------------------------------------------------------------------
# 3. Process all CRS network files
# --------------------------------------------------------------------------

output_dir <- "Data/clean_nets_harmonized"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

all_nets <- list.files(
  path = "Data/clean_nets",
  full.names = TRUE,
  pattern = "\\.(rds|Rds)$"
)
all_nets <- all_nets[grepl("crs", all_nets, ignore.case = TRUE)]

cat(sprintf("\nProcessing %d CRS network files\n", length(all_nets)))

total_renamed <- 0
total_merged <- 0
total_networks <- 0

log_entries <- list()

for (f in all_nets) {
  cat(sprintf("\n  Processing: %s\n", basename(f)))
  flush.console()

  network_data <- readRDS(f)
  networks <- network_data$networks

  file_renamed <- 0
  file_merged <- 0

  for (net_name in names(networks)) {
    g <- networks[[net_name]]
    result <- harmonize_network(g, name_map)

    networks[[net_name]] <- result$graph
    file_renamed <- file_renamed + result$n_renamed
    file_merged <- file_merged + result$n_merged
    total_networks <- total_networks + 1

    if (result$n_renamed > 0 || result$n_merged > 0) {
      cat(sprintf("    %s: %d renamed, %d merged nodes\n",
                  net_name, result$n_renamed, result$n_merged))
    }
  }

  total_renamed <- total_renamed + file_renamed
  total_merged <- total_merged + file_merged

  log_entries[[basename(f)]] <- data.frame(
    file = basename(f),
    n_networks = length(networks),
    n_renamed = file_renamed,
    n_merged = file_merged,
    stringsAsFactors = FALSE
  )

  # Save harmonized networks with same structure
  network_data$networks <- networks
  out_path <- file.path(output_dir, basename(f))
  saveRDS(network_data, out_path)
  cat(sprintf("    Saved to: %s\n", out_path))
}

# --------------------------------------------------------------------------
# 4. Write log
# --------------------------------------------------------------------------

log_df <- do.call(rbind, log_entries)
write.csv(log_df, "Data/entity_resolution/harmonization_log.csv", row.names = FALSE)

# --------------------------------------------------------------------------
# 5. Summary
# --------------------------------------------------------------------------

cat("\n=== Harmonization Summary ===\n")
cat(sprintf("Files processed:       %d\n", length(all_nets)))
cat(sprintf("Networks processed:    %d\n", total_networks))
cat(sprintf("Total names renamed:   %d\n", total_renamed))
cat(sprintf("Total nodes merged:    %d\n", total_merged))
cat(sprintf("\nOutput directory: %s\n", output_dir))
cat(sprintf("Log: Data/entity_resolution/harmonization_log.csv\n"))

cat("\nPer-file breakdown:\n")
for (i in seq_len(nrow(log_df))) {
  cat(sprintf("  %-35s  networks: %2d  renamed: %4d  merged: %2d\n",
              log_df$file[i], log_df$n_networks[i],
              log_df$n_renamed[i], log_df$n_merged[i]))
}

cat("\n=== Step 5 complete ===\n")
cat("Harmonized networks saved to Data/clean_nets_harmonized/\n")
cat("Update 17_network_topology.R to read from Data/clean_nets_harmonized/ instead of Data/clean_nets/\n")
