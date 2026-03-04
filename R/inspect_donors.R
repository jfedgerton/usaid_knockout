###############################################################################
# Investigate GIZ/JICA/AFD 0-partner issue
###############################################################################

library(igraph)

# Load one harmonized network to check
nets <- readRDS("Data/clean_nets_harmonized/hti_result_crs.Rds")
g <- nets[[1]]
node_names <- V(g)$name

cat("=== Haiti 2005 network ===\n")
cat("Nodes:", vcount(g), "\n")
cat("Edges:", ecount(g), "\n\n")

# Search for GIZ-related names
cat("--- GIZ-related node names across ALL networks ---\n")
all_names <- c()
for (fname in list.files("Data/clean_nets_harmonized", full.names = TRUE)) {
  nets <- readRDS(fname)
  for (g in nets) {
    all_names <- c(all_names, V(g)$name)
  }
}
all_names <- unique(all_names)

cat("\nGIZ matches:\n")
giz_names <- all_names[grepl("giz|gesellschaft fur internationale zusammenarbeit|german|deutsche gesellschaft", all_names, ignore.case = TRUE)]
print(giz_names)

cat("\nJICA matches:\n")
jica_names <- all_names[grepl("jica|japan international cooperation|japan", all_names, ignore.case = TRUE)]
print(jica_names)

cat("\nAFD matches:\n")
afd_names <- all_names[grepl("afd|agence francaise|french development", all_names, ignore.case = TRUE)]
print(afd_names)

cat("\nFrance-related:\n")
france_names <- all_names[grepl("france|french|francais", all_names, ignore.case = TRUE)]
print(france_names)

cat("\n\n=== Now checking what the script 17 regex actually matches ===\n")

# The regex patterns from script 17
giz_pattern <- "deutsche gesellschaft|\\bgiz\\b|german.*(international|develop)"
jica_pattern <- "japan international cooperation|\\bjica\\b"
afd_pattern <- "agence francaise de developpement|\\bafd\\b"

cat("\nGIZ regex matches:\n")
giz_match <- all_names[grepl(giz_pattern, all_names, ignore.case = TRUE)]
print(giz_match)

cat("\nJICA regex matches:\n")
jica_match <- all_names[grepl(jica_pattern, all_names, ignore.case = TRUE)]
print(jica_match)

cat("\nAFD regex matches:\n")
afd_match <- all_names[grepl(afd_pattern, all_names, ignore.case = TRUE)]
print(afd_match)

# Check a specific network where GIZ is present
cat("\n\n=== Checking a GIZ node in a specific network ===\n")
nets <- readRDS("Data/clean_nets_harmonized/hti_result_crs.Rds")
for (i in seq_along(nets)) {
  g <- nets[[i]]
  giz_nodes <- V(g)$name[grepl(giz_pattern, V(g)$name, ignore.case = TRUE)]
  if (length(giz_nodes) > 0) {
    cat(sprintf("\nNetwork %d (year %d):\n", i, 2004 + i))
    for (gn in giz_nodes) {
      vid <- which(V(g)$name == gn)
      out_e <- incident(g, vid, mode = "out")
      in_e <- incident(g, vid, mode = "in")
      cat(sprintf("  Node: %s\n", gn))
      cat(sprintf("  Out-degree: %d, In-degree: %d\n", length(out_e), length(in_e)))
      if (length(out_e) > 0) {
        targets <- head(V(g)$name[head_of(g, out_e)], 5)
        cat("  First 5 out-neighbors:", paste(targets, collapse = "; "), "\n")
      }
    }
    break  # just show first match
  }
}
