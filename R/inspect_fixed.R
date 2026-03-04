###############################################################################
# Verify the fixed donor metrics
###############################################################################

d <- readRDS("Data/clean_nets/network_attributes_directed.rds")
cat("=== NETWORK PANEL: Fixed donor metrics ===\n\n")
cat("Dimensions:", nrow(d), "x", ncol(d), "\n\n")

donors <- c("usaid", "dfid", "giz", "jica", "afd", "sida", "cida", "norad")

for (donor in donors) {
  pres_col <- paste0(donor, "_present")
  part_col <- paste0(donor, "_n_partners")
  sfr_col <- paste0(donor, "_sole_funder_rate")
  mpd_col <- paste0(donor, "_mean_partner_degree")

  present <- sum(d[[pres_col]])
  has_partners <- sum(d[[part_col]] > 0, na.rm = TRUE)
  max_partners <- max(d[[part_col]], na.rm = TRUE)

  cat(sprintf("%-6s: present=%3d/%d, has_partners=%3d, max_partners=%3d",
              toupper(donor), present, nrow(d), has_partners, max_partners))

  if (has_partners > 0) {
    sfr_mean <- mean(d[[sfr_col]], na.rm = TRUE)
    mpd_mean <- mean(d[[mpd_col]], na.rm = TRUE)
    cat(sprintf(", mean_sole_funder_rate=%.3f, mean_partner_degree=%.1f", sfr_mean, mpd_mean))
  }
  cat("\n")
}

cat("\n\n=== NODE PANEL check ===\n")
n <- readRDS("Data/clean_nets/node_attributes_directed.rds")
cat("Dimensions:", nrow(n), "x", ncol(n), "\n")
cat("Unique nodes:", length(unique(n$node_name)), "\n")
cat("States:", length(unique(n$state)), "\n")

# Quick sanity checks
cat("\n--- Key variable summaries ---\n")
cat("usaid_funded: ", sum(n$usaid_funded == 1), "yes /", sum(n$usaid_funded == 0), "no\n")
cat("sole_usaid_funded: ", sum(n$sole_usaid_funded == 1), "yes /", sum(n$sole_usaid_funded == 0), "no\n")
cat("distance_to_usaid NAs:", sum(is.na(n$distance_to_usaid)), "of", nrow(n), "\n")

cat("\nFragmentation: [", min(d$fragmentation), ",", max(d$fragmentation), "]\n")
cat("Weighted frag: [", min(d$weighted_fragmentation), ",", max(d$weighted_fragmentation), "]\n")
