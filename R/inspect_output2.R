###############################################################################
# Inspect node-level panel + diagnostics
###############################################################################

cat("========================================================\n")
cat("=== NODE-LEVEL PANEL ===\n")
cat("========================================================\n\n")

n <- readRDS("Data/clean_nets/node_attributes_directed.rds")
cat("Class:", class(n), "\n")
cat("Dimensions:", nrow(n), "rows x", ncol(n), "cols\n\n")

cat("Column names:\n")
print(names(n))

cat("\nColumn types:\n")
print(sapply(n, class))

cat("\n--- Unique states ---\n")
if ("state" %in% names(n)) {
  cat("N unique states:", length(unique(n$state)), "\n")
  print(sort(unique(n$state)))
}

cat("\n--- Year range ---\n")
if ("year" %in% names(n)) {
  cat("Years:", min(n$year), "to", max(n$year), "\n")
  print(table(n$year))
}

cat("\n--- Unique nodes ---\n")
name_col <- if ("node_name" %in% names(n)) "node_name" else "name"
cat("Name column:", name_col, "\n")
cat("N unique node names:", length(unique(n[[name_col]])), "\n")

cat("\n--- Rows per state ---\n")
if ("state" %in% names(n)) {
  print(table(n$state))
}

cat("\n--- NA counts per column ---\n")
na_counts <- colSums(is.na(n))
print(na_counts[na_counts > 0])
if (all(na_counts == 0)) cat("No NAs in any column!\n")

cat("\n--- Summary stats for key numeric columns ---\n")
key_cols <- c("in_degree", "out_degree", "total_degree", "n_funders",
              "n_non_usaid_funders", "betweenness", "closeness",
              "pagerank", "eigenvector_centrality",
              "usaid_exposure_share", "distance_to_usaid",
              "rewire_in", "rewire_out")
key_cols <- key_cols[key_cols %in% names(n)]
if (length(key_cols) > 0) {
  print(summary(n[, key_cols]))
}

cat("\n--- USAID-related binary columns ---\n")
usaid_cols <- c("usaid_funded", "sole_usaid_funded")
usaid_cols <- usaid_cols[usaid_cols %in% names(n)]
for (col in usaid_cols) {
  cat(col, ":\n")
  print(table(n[[col]], useNA = "ifany"))
  cat("\n")
}

cat("\n--- First 5 rows ---\n")
print(head(n, 5))

cat("\n--- USAID-related node sample ---\n")
usaid_rows <- n[grepl("usaid|united states agency", n[[name_col]], ignore.case = TRUE), ]
cat("N USAID node-years:", nrow(usaid_rows), "\n")
if (nrow(usaid_rows) > 0) {
  cat("Unique USAID node names:\n")
  print(unique(usaid_rows[[name_col]]))
  cat("\nFirst 5 USAID rows:\n")
  print(head(usaid_rows, 5))
}

cat("\n\n========================================================\n")
cat("=== KEY DIAGNOSTICS ===\n")
cat("========================================================\n\n")

# Check network panel
d <- readRDS("Data/clean_nets/network_attributes_directed.rds")

cat("--- GIZ, JICA, AFD: All have 0 partners across all networks ---\n")
cat("GIZ n_partners: min=", min(d$giz_n_partners), " max=", max(d$giz_n_partners), "\n")
cat("JICA n_partners: min=", min(d$jica_n_partners), " max=", max(d$jica_n_partners), "\n")
cat("AFD n_partners: min=", min(d$afd_n_partners), " max=", max(d$afd_n_partners), "\n")

cat("\nGIZ present count:", sum(d$giz_present), "of", nrow(d), "\n")
cat("JICA present count:", sum(d$jica_present), "of", nrow(d), "\n")
cat("AFD present count:", sum(d$afd_present), "of", nrow(d), "\n")

cat("\n--- Edge/vertex connectivity always 0 ---\n")
cat("Edge connectivity: all 0?", all(d$edge_connectivity == 0), "\n")
cat("Vertex connectivity: all 0?", all(d$vertex_connectivity == 0), "\n")

cat("\n--- Reciprocity very low ---\n")
cat("Reciprocity: mean=", mean(d$reciprocity), " max=", max(d$reciprocity), "\n")

cat("\n--- Flow hierarchy near 1 ---\n")
cat("Flow hierarchy: mean=", mean(d$flow_hierarchy), " min=", min(d$flow_hierarchy), "\n")

cat("\n--- num_components matches num_nodes? ---\n")
cat("Cases where num_components == num_nodes - 1:",
    sum(abs(d$num_components - (d$num_nodes - 1)) <= 1), "of", nrow(d), "\n")
cat("Correlation:", cor(d$num_components, d$num_nodes), "\n")

cat("\n--- Fragmentation range check ---\n")
cat("Fragmentation: [", min(d$fragmentation), ",", max(d$fragmentation), "]\n")
cat("Weighted fragmentation: [", min(d$weighted_fragmentation), ",", max(d$weighted_fragmentation), "]\n")

cat("\n--- Distance to USAID distribution ---\n")
if ("distance_to_usaid" %in% names(n)) {
  cat("Summary:\n")
  print(summary(n$distance_to_usaid))
  cat("Inf count:", sum(is.infinite(n$distance_to_usaid)), "\n")
  cat("NA count:", sum(is.na(n$distance_to_usaid)), "\n")
  cat("Non-NA, non-Inf count:", sum(!is.na(n$distance_to_usaid) & !is.infinite(n$distance_to_usaid)), "\n")
}

cat("\n=== INSPECTION COMPLETE ===\n")
