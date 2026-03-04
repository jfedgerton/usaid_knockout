###############################################################################
# Investigate GIZ/JICA/AFD 0-partner issue
###############################################################################

library(igraph)

# Check structure
nets <- readRDS("Data/clean_nets_harmonized/hti_result_crs.Rds")
cat("Top-level names:", names(nets), "\n")
cat("networks class:", class(nets$networks), "\n")
cat("networks length:", length(nets$networks), "\n")
cat("networks names:", head(names(nets$networks), 5), "\n")
cat("First network class:", class(nets$networks[[1]]), "\n")

g <- nets$networks[[1]]
if (is_igraph(g)) {
  cat("\n=== First network (Haiti 2005) ===\n")
  cat("Nodes:", vcount(g), "\n")
  cat("Edges:", ecount(g), "\n")
  cat("Sample node names:\n")
  print(head(V(g)$name, 10))
}

# Now search all node names across harmonized files
cat("\n\n=== Searching all harmonized networks for donor names ===\n")
all_names <- c()
files <- list.files("Data/clean_nets_harmonized", full.names = TRUE, pattern = "\\.Rds$")
for (fname in files) {
  dat <- readRDS(fname)
  for (g in dat$networks) {
    if (is_igraph(g)) {
      all_names <- c(all_names, V(g)$name)
    }
  }
}
all_names <- unique(all_names)
cat("Total unique names across all harmonized networks:", length(all_names), "\n\n")

# GIZ
giz_pattern <- "deutsche gesellschaft|\\bgiz\\b|german.*(international|develop)"
cat("GIZ regex matches:\n")
giz_match <- all_names[grepl(giz_pattern, all_names, ignore.case = TRUE)]
if (length(giz_match) > 0) print(giz_match) else cat("  NONE\n")

# JICA
jica_pattern <- "japan international cooperation|\\bjica\\b"
cat("\nJICA regex matches:\n")
jica_match <- all_names[grepl(jica_pattern, all_names, ignore.case = TRUE)]
if (length(jica_match) > 0) print(jica_match) else cat("  NONE\n")

# AFD
afd_pattern <- "agence francaise de developpement|\\bafd\\b"
cat("\nAFD regex matches:\n")
afd_match <- all_names[grepl(afd_pattern, all_names, ignore.case = TRUE)]
if (length(afd_match) > 0) print(afd_match) else cat("  NONE\n")

# Broader searches
cat("\n--- Broader searches ---\n")

cat("\n'japan' matches:\n")
japan_names <- all_names[grepl("japan", all_names, ignore.case = TRUE)]
if (length(japan_names) > 0) print(japan_names) else cat("  NONE\n")

cat("\n'france' or 'french' or 'francais' matches:\n")
france_names <- all_names[grepl("france|french|francais", all_names, ignore.case = TRUE)]
if (length(france_names) > 0) print(france_names) else cat("  NONE\n")

cat("\n'german' matches:\n")
german_names <- all_names[grepl("german", all_names, ignore.case = TRUE)]
if (length(german_names) > 0) print(german_names) else cat("  NONE\n")

cat("\n'giz' matches (case insensitive):\n")
giz_names <- all_names[grepl("giz", all_names, ignore.case = TRUE)]
if (length(giz_names) > 0) print(giz_names) else cat("  NONE\n")

# Now check: where GIZ is present, does it have out-degree > 0?
cat("\n\n=== Checking GIZ node behavior in networks ===\n")
giz_has_partners <- 0
giz_no_partners <- 0

for (fname in files) {
  dat <- readRDS(fname)
  state_name <- gsub("Data/clean_nets_harmonized/|_result_crs\\.Rds", "", fname)
  for (i in seq_along(dat$networks)) {
    g <- dat$networks[[i]]
    if (!is_igraph(g)) next
    giz_nodes <- which(grepl(giz_pattern, V(g)$name, ignore.case = TRUE))
    if (length(giz_nodes) > 0) {
      for (vid in giz_nodes) {
        outdeg <- degree(g, vid, mode = "out")
        indeg <- degree(g, vid, mode = "in")
        if (outdeg > 0) {
          giz_has_partners <- giz_has_partners + 1
        } else {
          giz_no_partners <- giz_no_partners + 1
        }
      }
    }
  }
}

cat("GIZ with out-degree > 0:", giz_has_partners, "\n")
cat("GIZ with out-degree = 0:", giz_no_partners, "\n")

# Check USAID for comparison
cat("\n=== USAID detection check ===\n")
usaid_pattern <- "united states agency for international development|\\busaid\\b"
usaid_found <- 0
for (fname in files[1:3]) {
  dat <- readRDS(fname)
  for (i in 1:min(3, length(dat$networks))) {
    g <- dat$networks[[i]]
    if (!is_igraph(g)) next
    usaid_nodes <- which(grepl(usaid_pattern, V(g)$name, ignore.case = TRUE))
    if (length(usaid_nodes) > 0) {
      for (vid in usaid_nodes) {
        outdeg <- degree(g, vid, mode = "out")
        cat(sprintf("  %s year %d: USAID out-degree = %d\n",
                    basename(fname), 2004+i, outdeg))
        usaid_found <- usaid_found + 1
      }
    }
  }
}
cat("USAID found:", usaid_found, "times in first 3 files\n")
