###############################################################################
# THOROUGH MODEL OUTPUT DIAGNOSTIC
# Checks: structural integrity, value ranges, logical consistency,
#         cross-panel alignment, donor detection, edge cases, NA patterns
###############################################################################

library(igraph)
library(dplyr)

d <- readRDS("Data/clean_nets/network_attributes_directed.rds")
n <- readRDS("Data/clean_nets/node_attributes_directed.rds")

issues <- list()
add_issue <- function(severity, category, msg) {
  issues[[length(issues) + 1]] <<- list(severity = severity, category = category, msg = msg)
}

cat("###############################################################################\n")
cat("# PART 1: NETWORK-LEVEL PANEL CHECKS\n")
cat("###############################################################################\n\n")

# --- 1.1 Dimensions & Coverage ---
cat("=== 1.1 Dimensions & Coverage ===\n")
cat("Rows:", nrow(d), " Cols:", ncol(d), "\n")
n_states <- length(unique(d$state))
n_years <- length(unique(d$year))
cat("States:", n_states, " Years:", n_years, " Expected rows:", n_states * n_years, "\n")
if (nrow(d) != n_states * n_years) {
  add_issue("ERROR", "coverage", sprintf("Expected %d rows (states*years) but got %d", n_states * n_years, nrow(d)))
}

# Check balanced panel
state_year_counts <- table(d$state)
unbalanced <- names(state_year_counts[state_year_counts != n_years])
if (length(unbalanced) > 0) {
  add_issue("ERROR", "coverage", sprintf("Unbalanced panel! States with != %d years: %s", n_years, paste(unbalanced, collapse = ", ")))
} else {
  cat("Panel is balanced: every state has", n_years, "years\n")
}

# Check for duplicates
dupes <- d %>% group_by(state, year) %>% filter(n() > 1)
if (nrow(dupes) > 0) {
  add_issue("ERROR", "coverage", sprintf("Duplicate state-year rows found: %d", nrow(dupes)))
} else {
  cat("No duplicate state-year rows\n")
}

# --- 1.2 Value Range Checks ---
cat("\n=== 1.2 Value Range Checks ===\n")

# Density should be [0, 1]
if (any(d$density < 0 | d$density > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "density outside [0,1]")
} else cat("density [0,1]: OK\n")

# Fragmentation [0, 1]
if (any(d$fragmentation < 0 | d$fragmentation > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "fragmentation outside [0,1]")
} else cat("fragmentation [0,1]: OK\n")

if (any(d$weighted_fragmentation < 0 | d$weighted_fragmentation > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "weighted_fragmentation outside [0,1]")
} else cat("weighted_fragmentation [0,1]: OK\n")

# Modularity [0, 1] for non-overlapping communities
if (any(d$modularity < -0.5 | d$modularity > 1, na.rm = TRUE)) {
  add_issue("WARNING", "range", sprintf("modularity outside expected range: min=%.3f max=%.3f", min(d$modularity, na.rm=T), max(d$modularity, na.rm=T)))
} else cat("modularity [-0.5,1]: OK (min=%.3f, max=%.3f)\n", min(d$modularity, na.rm=T), max(d$modularity, na.rm=T))

# Reciprocity [0, 1]
if (any(d$reciprocity < 0 | d$reciprocity > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "reciprocity outside [0,1]")
} else cat(sprintf("reciprocity [0,1]: OK (max=%.4f)\n", max(d$reciprocity, na.rm = TRUE)))

# Centralization [0, 1]
for (col in c("indegree_centralization", "outdegree_centralization", "betweenness_centralization")) {
  if (any(d[[col]] < 0 | d[[col]] > 1, na.rm = TRUE)) {
    add_issue("ERROR", "range", sprintf("%s outside [0,1]", col))
  } else cat(sprintf("%s [0,1]: OK\n", col))
}

# flow_hierarchy [0, 1]
if (any(d$flow_hierarchy < 0 | d$flow_hierarchy > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "flow_hierarchy outside [0,1]")
} else cat(sprintf("flow_hierarchy [0,1]: OK (min=%.4f)\n", min(d$flow_hierarchy, na.rm = TRUE)))

# clustering_coef [0, 1]
if (any(d$clustering_coef < 0 | d$clustering_coef > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "clustering_coef outside [0,1]")
} else cat("clustering_coef [0,1]: OK\n")

# donor_concentration [0, 1] (HHI)
if (any(d$donor_concentration < 0 | d$donor_concentration > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "donor_concentration outside [0,1]")
} else cat("donor_concentration [0,1]: OK\n")

# funding_redundancy [0, 1]
if (any(d$funding_redundancy < 0 | d$funding_redundancy > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "funding_redundancy outside [0,1]")
} else cat("funding_redundancy [0,1]: OK\n")

# network_efficiency [0, 1]
if (any(d$network_efficiency < 0 | d$network_efficiency > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "network_efficiency outside [0,1]")
} else cat("network_efficiency [0,1]: OK\n")

# sole_funder_rate [0,1] for all donors
donors <- c("usaid", "dfid", "giz", "jica", "afd", "sida", "cida", "norad")
for (donor in donors) {
  sfr_col <- paste0(donor, "_sole_funder_rate")
  vals <- d[[sfr_col]]
  if (any(vals < 0 | vals > 1, na.rm = TRUE)) {
    add_issue("ERROR", "range", sprintf("%s outside [0,1]", sfr_col))
  }
}
cat("All donor sole_funder_rates [0,1]: OK\n")

# --- 1.3 Logical Consistency ---
cat("\n=== 1.3 Logical Consistency ===\n")

# avg_in_degree should equal avg_out_degree (directed graph property)
if (!all(abs(d$avg_in_degree - d$avg_out_degree) < 1e-10)) {
  add_issue("ERROR", "logic", "avg_in_degree != avg_out_degree (should be equal for directed graphs)")
} else cat("avg_in_degree == avg_out_degree: OK\n")

# num_edges / num_nodes should approximate avg_degree
expected_avg_deg <- d$num_edges / d$num_nodes
if (any(abs(d$avg_in_degree - expected_avg_deg) > 0.01)) {
  add_issue("WARNING", "logic", "avg_in_degree doesn't match num_edges/num_nodes")
} else cat("avg_in_degree ≈ num_edges/num_nodes: OK\n")

# num_edges >= num_nodes - 1 for connected, but we have disconnected nets
# At minimum num_edges >= 1
if (any(d$num_edges < 1)) {
  add_issue("WARNING", "logic", sprintf("Networks with 0 edges: %d", sum(d$num_edges < 1)))
} else cat("All networks have >= 1 edge: OK\n")

# density should equal num_edges / (num_nodes * (num_nodes - 1)) for directed
expected_density <- d$num_edges / (d$num_nodes * (d$num_nodes - 1))
if (any(abs(d$density - expected_density) > 0.001, na.rm = TRUE)) {
  mismatches <- which(abs(d$density - expected_density) > 0.001)
  add_issue("WARNING", "logic", sprintf("density doesn't match e/(n*(n-1)) in %d rows", length(mismatches)))
} else cat("density = e/(n*(n-1)): OK\n")

# usaid_n_partners should be 0 when usaid_present = 0
if (any(d$usaid_present == 0 & d$usaid_n_partners > 0)) {
  add_issue("ERROR", "logic", "usaid_n_partners > 0 when usaid_present = 0")
} else cat("usaid_present=0 => n_partners=0: OK\n")

# Check for each bilateral donor
for (donor in donors) {
  pres_col <- paste0(donor, "_present")
  part_col <- paste0(donor, "_n_partners")
  if (any(d[[pres_col]] == 0 & d[[part_col]] > 0, na.rm = TRUE)) {
    add_issue("ERROR", "logic", sprintf("%s has partners when not present", donor))
  }
}
cat("All donors: present=0 => n_partners=0: OK\n")

# largest_component_size should be <= num_nodes
if (any(d$largest_component_size > d$num_nodes)) {
  add_issue("ERROR", "logic", "largest_component_size > num_nodes")
} else cat("largest_component_size <= num_nodes: OK\n")

# edge_connectivity & vertex_connectivity always 0 — is this expected?
cat(sprintf("\nedge_connectivity: all zeros? %s\n", all(d$edge_connectivity == 0)))
cat(sprintf("vertex_connectivity: all zeros? %s\n", all(d$vertex_connectivity == 0)))
if (all(d$edge_connectivity == 0) && all(d$vertex_connectivity == 0)) {
  add_issue("NOTE", "logic", "edge_connectivity and vertex_connectivity are always 0 (expected for disconnected directed networks)")
}

# --- 1.4 Suspicious Patterns ---
cat("\n=== 1.4 Suspicious Patterns ===\n")

# Any constant columns?
for (col in names(d)) {
  if (is.numeric(d[[col]]) && !all(is.na(d[[col]]))) {
    if (sd(d[[col]], na.rm = TRUE) == 0) {
      add_issue("WARNING", "suspicious", sprintf("Column %s has zero variance (constant = %s)", col, unique(d[[col]][!is.na(d[[col]])])[1]))
    }
  }
}

# Very small networks
tiny <- d[d$num_nodes <= 5, ]
if (nrow(tiny) > 0) {
  cat(sprintf("Tiny networks (<=5 nodes): %d\n", nrow(tiny)))
  for (r in seq_len(nrow(tiny))) {
    cat(sprintf("  %s %d: %d nodes, %d edges\n", tiny$state[r], tiny$year[r], tiny$num_nodes[r], tiny$num_edges[r]))
  }
  add_issue("NOTE", "suspicious", sprintf("%d networks have <=5 nodes — may produce unstable metrics", nrow(tiny)))
}

# Check if usaid_sole_funder_rate and usaid_share_partners_no_alt_funder are identical
if (all(d$usaid_sole_funder_rate == d$usaid_share_partners_no_alt_funder |
        (is.na(d$usaid_sole_funder_rate) & is.na(d$usaid_share_partners_no_alt_funder)))) {
  add_issue("NOTE", "redundancy", "usaid_sole_funder_rate and usaid_share_partners_no_alt_funder are identical columns")
  cat("usaid_sole_funder_rate == usaid_share_partners_no_alt_funder: IDENTICAL (redundant)\n")
} else {
  cat("usaid_sole_funder_rate vs usaid_share_partners_no_alt_funder: DIFFERENT\n")
}

# --- 1.5 NA Pattern Analysis ---
cat("\n=== 1.5 NA Pattern Analysis ===\n")
na_counts <- colSums(is.na(d))
na_cols <- na_counts[na_counts > 0]
if (length(na_cols) > 0) {
  cat("Columns with NAs:\n")
  for (col_name in names(na_cols)) {
    pct <- round(100 * na_cols[col_name] / nrow(d), 1)
    cat(sprintf("  %-45s %3d NAs (%5.1f%%)\n", col_name, na_cols[col_name], pct))
  }
} else {
  cat("No NAs in any column\n")
}

# Check if NAs are structurally sensible (donor NAs = when donor not present)
for (donor in donors) {
  sfr_col <- paste0(donor, "_sole_funder_rate")
  part_col <- paste0(donor, "_n_partners")
  na_sfr <- is.na(d[[sfr_col]])
  zero_partners <- d[[part_col]] == 0
  if (!all(na_sfr == zero_partners)) {
    mismatch_n <- sum(na_sfr != zero_partners)
    add_issue("WARNING", "na_pattern", sprintf("%s: %d rows where NA pattern doesn't match zero-partners", sfr_col, mismatch_n))
  }
}
cat("Donor NA patterns (sole_funder_rate NA iff n_partners=0): checked\n")

# --- 1.6 Temporal Trends ---
cat("\n=== 1.6 Temporal Trends (any sudden jumps?) ===\n")
for (st in unique(d$state)) {
  sub <- d[d$state == st, ]
  sub <- sub[order(sub$year), ]
  for (col in c("num_nodes", "num_edges", "fragmentation")) {
    vals <- sub[[col]]
    if (length(vals) > 1) {
      pct_changes <- diff(vals) / pmax(vals[-length(vals)], 1) * 100
      big_jumps <- which(abs(pct_changes) > 200)
      if (length(big_jumps) > 0) {
        for (j in big_jumps) {
          add_issue("WARNING", "temporal", sprintf("%s %s: %d->%d jump %.0f%% (year %d->%d)",
                                                     st, col, vals[j], vals[j+1], pct_changes[j],
                                                     sub$year[j], sub$year[j+1]))
        }
      }
    }
  }
}
cat("Temporal jump check (>200% change year-over-year) complete\n")


cat("\n\n###############################################################################\n")
cat("# PART 2: NODE-LEVEL PANEL CHECKS\n")
cat("###############################################################################\n\n")

# --- 2.1 Dimensions ---
cat("=== 2.1 Dimensions ===\n")
cat("Rows:", nrow(n), " Cols:", ncol(n), "\n")
cat("Unique nodes:", length(unique(n$node_name)), "\n")
cat("Unique states:", length(unique(n$state)), "\n")

# --- 2.2 Value Range Checks ---
cat("\n=== 2.2 Value Range Checks ===\n")

# Degrees should be >= 0
for (col in c("in_degree", "out_degree", "total_degree")) {
  if (any(n[[col]] < 0, na.rm = TRUE)) {
    add_issue("ERROR", "range", sprintf("%s has negative values", col))
  }
}
cat("Degrees >= 0: OK\n")

# total_degree should equal in_degree + out_degree
if (any(abs(n$total_degree - (n$in_degree + n$out_degree)) > 0, na.rm = TRUE)) {
  bad_n <- sum(abs(n$total_degree - (n$in_degree + n$out_degree)) > 0, na.rm = TRUE)
  add_issue("ERROR", "logic", sprintf("total_degree != in_degree + out_degree in %d rows", bad_n))
} else cat("total_degree = in_degree + out_degree: OK\n")

# n_funders should equal in_degree (since funders are in-neighbors)
if (any(abs(n$n_funders - n$in_degree) > 0, na.rm = TRUE)) {
  bad_n <- sum(abs(n$n_funders - n$in_degree) > 0, na.rm = TRUE)
  add_issue("WARNING", "logic", sprintf("n_funders != in_degree in %d rows", bad_n))
} else cat("n_funders == in_degree: OK\n")

# n_non_usaid_funders <= n_funders
if (any(n$n_non_usaid_funders > n$n_funders, na.rm = TRUE)) {
  add_issue("ERROR", "logic", "n_non_usaid_funders > n_funders")
} else cat("n_non_usaid_funders <= n_funders: OK\n")

# usaid_exposure_share [0, 1]
if (any(n$usaid_exposure_share < 0 | n$usaid_exposure_share > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "usaid_exposure_share outside [0,1]")
} else cat("usaid_exposure_share [0,1]: OK\n")

# distance_to_usaid >= 0
if (any(n$distance_to_usaid < 0, na.rm = TRUE)) {
  add_issue("ERROR", "range", "distance_to_usaid has negative values")
} else cat("distance_to_usaid >= 0: OK\n")

# betweenness [0, 1] (normalized)
if (any(n$betweenness < 0 | n$betweenness > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "betweenness outside [0,1]")
} else cat("betweenness [0,1]: OK\n")

# closeness [0, 1] (normalized)
if (any(n$closeness < 0 | n$closeness > 1, na.rm = TRUE)) {
  add_issue("ERROR", "range", "closeness outside [0,1]")
} else cat("closeness [0,1]: OK\n")

# pagerank should sum to ~1 per network
cat("\n=== 2.3 PageRank Sum Check (should ≈ 1 per network) ===\n")
pr_sums <- n %>% group_by(network_name) %>% summarize(pr_sum = sum(pagerank), .groups = "drop")
pr_bad <- pr_sums[abs(pr_sums$pr_sum - 1) > 0.01, ]
if (nrow(pr_bad) > 0) {
  add_issue("ERROR", "logic", sprintf("PageRank doesn't sum to ~1 in %d networks (range: %.3f to %.3f)",
                                        nrow(pr_bad), min(pr_bad$pr_sum), max(pr_bad$pr_sum)))
  cat(sprintf("BAD: %d networks with PageRank sum != 1\n", nrow(pr_bad)))
} else {
  cat(sprintf("PageRank sums: min=%.4f max=%.4f — OK\n", min(pr_sums$pr_sum), max(pr_sums$pr_sum)))
}

# authority_score and hub_score range
cat("\n=== 2.4 Authority/Hub Score Checks ===\n")
for (col in c("authority_score", "hub_score")) {
  rng <- range(n[[col]], na.rm = TRUE)
  cat(sprintf("%s range: [%.6f, %.6f]\n", col, rng[1], rng[2]))
  if (rng[1] < -0.01) add_issue("ERROR", "range", sprintf("%s has negative values", col))
}

# eigenvector_centrality range
rng_eig <- range(n$eigenvector_centrality, na.rm = TRUE)
cat(sprintf("eigenvector_centrality range: [%.6f, %.6f]\n", rng_eig[1], rng_eig[2]))

# rewire_in and rewire_out [0, 1]
for (col in c("rewire_in", "rewire_out")) {
  if (any(n[[col]] < 0 | n[[col]] > 1, na.rm = TRUE)) {
    add_issue("ERROR", "range", sprintf("%s outside [0,1]", col))
  } else cat(sprintf("%s [0,1]: OK\n", col))
}

# --- 2.5 USAID Variable Consistency ---
cat("\n=== 2.5 USAID Variable Consistency ===\n")

# is_usaid should be binary
if (!all(n$is_usaid %in% c(0L, 1L))) {
  add_issue("ERROR", "logic", "is_usaid has non-binary values")
} else cat("is_usaid binary: OK\n")

# usaid_funded should be binary
if (!all(n$usaid_funded %in% c(0L, 1L))) {
  add_issue("ERROR", "logic", "usaid_funded has non-binary values")
} else cat("usaid_funded binary: OK\n")

# sole_usaid_funded should be binary
if (!all(n$sole_usaid_funded %in% c(0L, 1L))) {
  add_issue("ERROR", "logic", "sole_usaid_funded has non-binary values")
} else cat("sole_usaid_funded binary: OK\n")

# sole_usaid_funded=1 implies usaid_funded=1
if (any(n$sole_usaid_funded == 1 & n$usaid_funded == 0)) {
  add_issue("ERROR", "logic", "sole_usaid_funded=1 but usaid_funded=0")
} else cat("sole_usaid_funded=1 => usaid_funded=1: OK\n")

# usaid_funded=1 implies n_funders >= 1 (should receive from USAID so in-degree >= 1)
uf_but_no_funders <- sum(n$usaid_funded == 1 & n$n_funders == 0)
if (uf_but_no_funders > 0) {
  add_issue("ERROR", "logic", sprintf("usaid_funded=1 but n_funders=0 in %d rows", uf_but_no_funders))
} else cat("usaid_funded=1 => n_funders>=1: OK\n")

# sole_usaid_funded=1 implies n_non_usaid_funders=0
if (any(n$sole_usaid_funded == 1 & n$n_non_usaid_funders > 0)) {
  add_issue("ERROR", "logic", "sole_usaid_funded=1 but n_non_usaid_funders > 0")
} else cat("sole_usaid_funded=1 => n_non_usaid_funders=0: OK\n")

# is_usaid nodes should have distance_to_usaid = 0
usaid_rows <- n[n$is_usaid == 1, ]
if (any(usaid_rows$distance_to_usaid != 0, na.rm = TRUE)) {
  add_issue("ERROR", "logic", "is_usaid=1 but distance_to_usaid != 0")
} else cat("is_usaid=1 => distance_to_usaid=0: OK\n")

# is_usaid should NOT be usaid_funded (you don't fund yourself)
usaid_self_funded <- sum(n$is_usaid == 1 & n$usaid_funded == 1)
if (usaid_self_funded > 0) {
  add_issue("WARNING", "logic", sprintf("is_usaid=1 AND usaid_funded=1 in %d rows (USAID funding itself?)", usaid_self_funded))
} else cat("USAID not self-funded: OK\n")

# usaid_exposure_share should be 1.0 for is_usaid nodes
usaid_nodes_exp <- n[n$is_usaid == 1, "usaid_exposure_share"]
cat(sprintf("USAID node exposure_share: min=%.3f, max=%.3f, mean=%.3f\n",
            min(usaid_nodes_exp, na.rm=T), max(usaid_nodes_exp, na.rm=T), mean(usaid_nodes_exp, na.rm=T)))

# --- 2.6 Closeness NA Analysis ---
cat("\n=== 2.6 Closeness NA Analysis ===\n")
closeness_na <- sum(is.na(n$closeness))
closeness_pct <- round(100 * closeness_na / nrow(n), 1)
cat(sprintf("closeness NAs: %d (%.1f%% of rows)\n", closeness_na, closeness_pct))
if (closeness_pct > 50) {
  add_issue("WARNING", "na_pattern", sprintf("closeness has %.1f%% NAs — may limit its analytical usefulness", closeness_pct))
}

# Is closeness NA correlated with being isolated?
if (closeness_na > 0) {
  closeness_na_rows <- n[is.na(n$closeness), ]
  cat(sprintf("  Mean out_degree of closeness-NA rows: %.3f\n", mean(closeness_na_rows$out_degree)))
  cat(sprintf("  Mean out_degree of closeness-non-NA rows: %.3f\n", mean(n$out_degree[!is.na(n$closeness)])))
}


cat("\n\n###############################################################################\n")
cat("# PART 3: CROSS-PANEL CONSISTENCY\n")
cat("###############################################################################\n\n")

# --- 3.1 num_nodes vs actual node count ---
cat("=== 3.1 num_nodes vs actual node rows ===\n")
node_counts <- n %>% group_by(network_name) %>% summarize(actual_nodes = n(), .groups = "drop")
merged <- merge(d, node_counts, by = "network_name", all.x = TRUE)

mismatches <- merged[merged$num_nodes != merged$actual_nodes, ]
if (nrow(mismatches) > 0) {
  add_issue("ERROR", "cross_panel", sprintf("num_nodes doesn't match node panel rows in %d networks", nrow(mismatches)))
  cat("MISMATCHES:\n")
  print(mismatches[, c("state", "year", "num_nodes", "actual_nodes")])
} else {
  cat("num_nodes matches node panel row count in all networks: OK\n")
}

# --- 3.2 degree sums should match ---
cat("\n=== 3.2 Degree sums vs num_edges ===\n")
edge_check <- n %>% group_by(network_name) %>%
  summarize(sum_in = sum(in_degree), sum_out = sum(out_degree), .groups = "drop")
edge_merged <- merge(d[, c("network_name", "num_edges", "state", "year")], edge_check, by = "network_name")

# Sum of in-degrees = sum of out-degrees = num_edges
bad_edges <- edge_merged[edge_merged$sum_in != edge_merged$num_edges |
                          edge_merged$sum_out != edge_merged$num_edges, ]
if (nrow(bad_edges) > 0) {
  add_issue("ERROR", "cross_panel", sprintf("Sum of degrees != num_edges in %d networks", nrow(bad_edges)))
  cat("First 5 mismatches:\n")
  print(head(bad_edges))
} else {
  cat("Sum(in_degree) = Sum(out_degree) = num_edges in all networks: OK\n")
}

# --- 3.3 USAID presence consistency ---
cat("\n=== 3.3 USAID presence: network vs node panel ===\n")
usaid_node_check <- n %>% group_by(network_name) %>%
  summarize(any_usaid = max(is_usaid), n_usaid = sum(is_usaid), .groups = "drop")
usaid_merged <- merge(d[, c("network_name", "usaid_present", "usaid_n_nodes", "state", "year")],
                       usaid_node_check, by = "network_name")

pres_mismatch <- usaid_merged[usaid_merged$usaid_present != usaid_merged$any_usaid, ]
if (nrow(pres_mismatch) > 0) {
  add_issue("ERROR", "cross_panel", sprintf("usaid_present mismatch between panels in %d networks", nrow(pres_mismatch)))
} else cat("usaid_present matches node panel: OK\n")

nodes_mismatch <- usaid_merged[usaid_merged$usaid_n_nodes != usaid_merged$n_usaid, ]
if (nrow(nodes_mismatch) > 0) {
  add_issue("ERROR", "cross_panel", sprintf("usaid_n_nodes mismatch in %d networks", nrow(nodes_mismatch)))
} else cat("usaid_n_nodes matches node panel: OK\n")

# --- 3.4 usaid_n_partners vs usaid_funded count ---
cat("\n=== 3.4 usaid_n_partners vs usaid_funded count in node panel ===\n")
funded_check <- n %>% group_by(network_name) %>%
  summarize(n_funded = sum(usaid_funded), .groups = "drop")
funded_merged <- merge(d[, c("network_name", "usaid_n_partners", "state", "year")],
                        funded_check, by = "network_name")
partners_mismatch <- funded_merged[funded_merged$usaid_n_partners != funded_merged$n_funded, ]
if (nrow(partners_mismatch) > 0) {
  add_issue("ERROR", "cross_panel", sprintf("usaid_n_partners != count(usaid_funded=1) in %d networks", nrow(partners_mismatch)))
  cat("First 5 mismatches:\n")
  print(head(partners_mismatch[, c("state", "year", "usaid_n_partners", "n_funded")]))
} else {
  cat("usaid_n_partners matches usaid_funded count: OK\n")
}


cat("\n\n###############################################################################\n")
cat("# PART 4: DONOR DETECTION THOROUGHNESS\n")
cat("###############################################################################\n\n")

# Check all unique node names for potential donor misses
all_node_names <- unique(n$node_name)
cat("Total unique node names:", length(all_node_names), "\n\n")

# Look for potential USAID variants missed
usaid_pat <- "United States Agency for International Development|\\bUSAID\\b|US AID"
potential_us <- all_node_names[grepl("united states|\\bus\\b|america", all_node_names, ignore.case = TRUE)]
actual_usaid <- all_node_names[grepl(usaid_pat, all_node_names, ignore.case = TRUE)]
missed_us <- setdiff(potential_us, actual_usaid)
if (length(missed_us) > 0) {
  cat("Potential USAID-related names NOT matched by regex:\n")
  for (nm in missed_us) cat(sprintf("  %s\n", nm))
}

# Look for potential DFID/FCDO misses
dfid_pat <- "Department for International Development|\\bDFID\\b|\\bFCDO\\b|Foreign.*Commonwealth.*Development"
potential_uk <- all_node_names[grepl("uk|britain|british|united kingdom|england", all_node_names, ignore.case = TRUE)]
actual_dfid <- all_node_names[grepl(dfid_pat, all_node_names, ignore.case = TRUE)]
missed_uk <- setdiff(potential_uk, actual_dfid)
if (length(missed_uk) > 0) {
  cat("\nPotential UK/DFID-related names NOT matched by regex:\n")
  for (nm in head(missed_uk, 20)) cat(sprintf("  %s\n", nm))
}

# Check SIDA
sida_pat <- "Swedish International Development|\\bSida\\b|\\bSIDA\\b"
potential_swe <- all_node_names[grepl("swed|\\bsida\\b", all_node_names, ignore.case = TRUE)]
actual_sida <- all_node_names[grepl(sida_pat, all_node_names, ignore.case = TRUE)]
missed_swe <- setdiff(potential_swe, actual_sida)
if (length(missed_swe) > 0) {
  cat("\nPotential SIDA-related names NOT matched by regex:\n")
  for (nm in missed_swe) cat(sprintf("  %s\n", nm))
}

# Check CIDA
cida_pat <- "Canadian International Development|\\bCIDA\\b|Global Affairs Canada"
potential_can <- all_node_names[grepl("canad|\\bcida\\b", all_node_names, ignore.case = TRUE)]
actual_cida <- all_node_names[grepl(cida_pat, all_node_names, ignore.case = TRUE)]
missed_can <- setdiff(potential_can, actual_cida)
if (length(missed_can) > 0) {
  cat("\nPotential CIDA-related names NOT matched by regex:\n")
  for (nm in missed_can) cat(sprintf("  %s\n", nm))
}

# Check NORAD
norad_pat <- "Norwegian Agency for Development|\\bNorad\\b|\\bNORAD\\b"
potential_nor <- all_node_names[grepl("norw|\\bnorad\\b", all_node_names, ignore.case = TRUE)]
actual_norad <- all_node_names[grepl(norad_pat, all_node_names, ignore.case = TRUE)]
missed_nor <- setdiff(potential_nor, actual_norad)
if (length(missed_nor) > 0) {
  cat("\nPotential NORAD-related names NOT matched by regex:\n")
  for (nm in missed_nor) cat(sprintf("  %s\n", nm))
}


cat("\n\n###############################################################################\n")
cat("# PART 5: ENTITY RESOLUTION VERIFICATION\n")
cat("###############################################################################\n\n")

# Check for near-duplicate node names that survived entity resolution
cat("=== 5.1 Checking for remaining near-duplicates ===\n")
# Quick check: exact duplicates ignoring case
lower_names <- tolower(all_node_names)
case_dupes <- all_node_names[duplicated(lower_names) | duplicated(lower_names, fromLast = TRUE)]
if (length(case_dupes) > 0) {
  add_issue("WARNING", "entity_resolution", sprintf("Case-variant duplicates found: %s", paste(head(case_dupes, 10), collapse = "; ")))
  cat("Case-variant duplicates:\n")
  for (nm in case_dupes) cat(sprintf("  %s\n", nm))
} else {
  cat("No case-variant duplicates: OK\n")
}

# Check for duplicate node names within same state-year
cat("\n=== 5.2 Duplicate node names within same network ===\n")
intra_dupes <- n %>% group_by(network_name, node_name) %>% filter(n() > 1)
if (nrow(intra_dupes) > 0) {
  add_issue("ERROR", "entity_resolution", sprintf("Duplicate node_name within same network in %d rows", nrow(intra_dupes)))
  cat("Duplicate rows:\n")
  print(head(intra_dupes, 10))
} else {
  cat("No duplicate node_names within any network: OK\n")
}


cat("\n\n###############################################################################\n")
cat("# PART 6: JACCARD TURNOVER CHECKS\n")
cat("###############################################################################\n\n")

# Jaccard should only be NA for year 1
cat("=== 6.1 Jaccard NA pattern ===\n")
year1_rows <- n[n$year == min(n$year), ]
non_year1_rows <- n[n$year > min(n$year), ]

cat(sprintf("Year 1 (2005) rows: %d, all j_in NA? %s, all j_out NA? %s\n",
            nrow(year1_rows), all(is.na(year1_rows$j_in)), all(is.na(year1_rows$j_out))))

# For non-year-1, NAs should only be for new nodes (not in previous year)
non_y1_na <- sum(is.na(non_year1_rows$j_in))
cat(sprintf("Non-year-1 j_in NAs: %d of %d (%.1f%%) — these are nodes new to that state-year\n",
            non_y1_na, nrow(non_year1_rows), 100 * non_y1_na / nrow(non_year1_rows)))

# rewire = 1 - jaccard, so it should be NA when j is NA
rewire_na_match <- all(is.na(n$rewire_in) == is.na(n$j_in)) && all(is.na(n$rewire_out) == is.na(n$j_out))
if (!rewire_na_match) {
  add_issue("ERROR", "logic", "rewire NA pattern doesn't match Jaccard NA pattern")
} else cat("rewire NA pattern matches Jaccard: OK\n")

# Check rewire = 1 - jaccard where non-NA
non_na_j <- !is.na(n$j_in) & !is.na(n$rewire_in)
if (any(abs(n$rewire_in[non_na_j] - (1 - n$j_in[non_na_j])) > 1e-10)) {
  add_issue("ERROR", "logic", "rewire_in != 1 - j_in")
} else cat("rewire_in = 1 - j_in: OK\n")


cat("\n\n###############################################################################\n")
cat("# PART 7: ADDITIONAL EDGE CASES\n")
cat("###############################################################################\n\n")

# Check for Inf values in any numeric column
cat("=== 7.1 Infinite values ===\n")
for (col in names(n)) {
  if (is.numeric(n[[col]])) {
    n_inf <- sum(is.infinite(n[[col]]))
    if (n_inf > 0) {
      add_issue("ERROR", "edge_case", sprintf("Column %s has %d Inf values", col, n_inf))
    }
  }
}
for (col in names(d)) {
  if (is.numeric(d[[col]])) {
    n_inf <- sum(is.infinite(d[[col]]))
    if (n_inf > 0) {
      add_issue("ERROR", "edge_case", sprintf("Network panel column %s has %d Inf values", col, n_inf))
    }
  }
}
cat("No Inf values in either panel: checked\n")

# Check for NaN values
cat("\n=== 7.2 NaN values ===\n")
for (col in names(n)) {
  if (is.numeric(n[[col]])) {
    n_nan <- sum(is.nan(n[[col]]))
    if (n_nan > 0) {
      add_issue("ERROR", "edge_case", sprintf("Node panel column %s has %d NaN values", col, n_nan))
    }
  }
}
for (col in names(d)) {
  if (is.numeric(d[[col]])) {
    n_nan <- sum(is.nan(d[[col]]))
    if (n_nan > 0) {
      add_issue("ERROR", "edge_case", sprintf("Network panel column %s has %d NaN values", col, n_nan))
    }
  }
}
cat("NaN check: done\n")

# Check that USAID_REGEX and BILATERAL_DONORS usaid pattern are consistent
cat("\n=== 7.3 USAID regex consistency ===\n")
usaid_main <- all_node_names[grepl("United States Agency for International Development|\\bUSAID\\b|US AID", all_node_names, ignore.case = TRUE)]
cat(sprintf("USAID main regex matches %d unique names: %s\n", length(usaid_main), paste(usaid_main, collapse = "; ")))

# Check that usaid_mean_partner_degree == usaid_mean_partner_total_degree
cat("\n=== 7.4 Redundant USAID columns ===\n")
eq_check <- all(d$usaid_mean_partner_degree == d$usaid_mean_partner_total_degree |
                 (is.na(d$usaid_mean_partner_degree) & is.na(d$usaid_mean_partner_total_degree)))
cat(sprintf("usaid_mean_partner_degree == usaid_mean_partner_total_degree? %s\n", eq_check))
if (eq_check) {
  add_issue("NOTE", "redundancy", "usaid_mean_partner_degree and usaid_mean_partner_total_degree are identical")
}

# Check for extremely small networks affecting centrality measures
cat("\n=== 7.5 Small network centrality stability ===\n")
small_nets <- d[d$num_nodes <= 10, ]
if (nrow(small_nets) > 0) {
  cat(sprintf("Networks with <= 10 nodes: %d\n", nrow(small_nets)))
  cat(sprintf("  Max centralization in small nets: in=%.3f, out=%.3f, betw=%.4f\n",
              max(small_nets$indegree_centralization, na.rm = TRUE),
              max(small_nets$outdegree_centralization, na.rm = TRUE),
              max(small_nets$betweenness_centralization, na.rm = TRUE)))
}


cat("\n\n###############################################################################\n")
cat("# SUMMARY OF ALL ISSUES\n")
cat("###############################################################################\n\n")

if (length(issues) == 0) {
  cat("NO ISSUES FOUND!\n")
} else {
  errors <- Filter(function(x) x$severity == "ERROR", issues)
  warnings <- Filter(function(x) x$severity == "WARNING", issues)
  notes <- Filter(function(x) x$severity == "NOTE", issues)

  cat(sprintf("Total issues: %d (ERRORS: %d, WARNINGS: %d, NOTES: %d)\n\n",
              length(issues), length(errors), length(warnings), length(notes)))

  if (length(errors) > 0) {
    cat("=== ERRORS (must fix) ===\n")
    for (e in errors) cat(sprintf("  [%s] %s\n", e$category, e$msg))
    cat("\n")
  }

  if (length(warnings) > 0) {
    cat("=== WARNINGS (should investigate) ===\n")
    for (w in warnings) cat(sprintf("  [%s] %s\n", w$category, w$msg))
    cat("\n")
  }

  if (length(notes) > 0) {
    cat("=== NOTES (informational) ===\n")
    for (n_item in notes) cat(sprintf("  [%s] %s\n", n_item$category, n_item$msg))
    cat("\n")
  }
}

cat("\n=== DIAGNOSTIC COMPLETE ===\n")
