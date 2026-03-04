library(igraph)

# Load all unique names
files <- list.files("Data/clean_nets_harmonized", full.names = TRUE, pattern = "\\.Rds$")
all_names <- c()
for (fname in files) {
  dat <- readRDS(fname)
  for (g in dat$networks) {
    if (is_igraph(g)) all_names <- c(all_names, V(g)$name)
  }
}
all_names <- unique(all_names)

# Test updated patterns
cat("=== Updated regex matches ===\n\n")

giz_pat <- "Deutsche Gesellschaft|\\bGIZ\\b|\\bBMZ\\b|German Federal Ministry|Federal Ministry for Economic Cooperation.*Germany"
cat("GIZ matches:\n")
print(all_names[grepl(giz_pat, all_names, ignore.case = TRUE)])

jica_pat <- "Japan(?:ese)? International Cooperation|\\bJICA\\b"
cat("\nJICA matches:\n")
print(all_names[grepl(jica_pat, all_names, ignore.case = TRUE, perl = TRUE)])

afd_pat <- "Agence Fran.aise de D.veloppement|\\bAFD\\b|French Development Agency"
cat("\nAFD matches:\n")
print(all_names[grepl(afd_pat, all_names, ignore.case = TRUE)])

# Check out-degree for matched nodes
cat("\n\n=== Out-degree check for matched donors ===\n")
for (fname in files) {
  dat <- readRDS(fname)
  state <- gsub("_result_crs\\.Rds", "", basename(fname))
  for (i in seq_along(dat$networks)) {
    g <- dat$networks[[i]]
    if (!is_igraph(g)) next

    # GIZ check (including BMZ)
    giz_nodes <- which(grepl(giz_pat, V(g)$name, ignore.case = TRUE))
    for (vid in giz_nodes) {
      od <- degree(g, vid, mode = "out")
      if (od > 0) {
        cat(sprintf("  GIZ/BMZ %s year %d: %s out-deg=%d\n", state, 2004+i, V(g)$name[vid], od))
      }
    }

    # JICA check
    jica_nodes <- which(grepl(jica_pat, V(g)$name, ignore.case = TRUE, perl = TRUE))
    for (vid in jica_nodes) {
      od <- degree(g, vid, mode = "out")
      if (od > 0) {
        cat(sprintf("  JICA %s year %d: %s out-deg=%d\n", state, 2004+i, V(g)$name[vid], od))
      }
    }

    # AFD check
    afd_nodes <- which(grepl(afd_pat, V(g)$name, ignore.case = TRUE))
    for (vid in afd_nodes) {
      od <- degree(g, vid, mode = "out")
      if (od > 0) {
        cat(sprintf("  AFD %s year %d: %s out-deg=%d\n", state, 2004+i, V(g)$name[vid], od))
      }
    }
  }
}

cat("\n=== BMZ (Federal Ministry for Economic Cooperation) stats ===\n")
bmz_present <- 0
bmz_has_partners <- 0
for (fname in files) {
  dat <- readRDS(fname)
  for (g in dat$networks) {
    if (!is_igraph(g)) next
    bmz_nodes <- which(grepl(giz_pat, V(g)$name, ignore.case = TRUE))
    if (length(bmz_nodes) > 0) {
      bmz_present <- bmz_present + 1
      total_out <- sum(sapply(bmz_nodes, function(v) degree(g, v, mode = "out")))
      if (total_out > 0) bmz_has_partners <- bmz_has_partners + 1
    }
  }
}
cat("Present:", bmz_present, "\n")
cat("Has partners (out-degree > 0):", bmz_has_partners, "\n")
