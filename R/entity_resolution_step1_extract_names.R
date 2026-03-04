###############################################################################
# Entity Resolution Step 1: Extract all unique organization names from networks
#
# Input:  Data/clean_nets/*crs*.Rds (igraph network objects)
# Output: Data/entity_resolution/raw_names.csv
#         Data/entity_resolution/raw_names_summary.csv
#
# Loops over all CRS network files, extracts V(g)$name from every network,
# and produces:
#   raw_names.csv        — one row per (raw_name, state, year) with n_appearances
#   raw_names_summary.csv — deduplicated: raw_name, n_states, n_years, first_year, last_year
###############################################################################

library(igraph)

cat("=== Entity Resolution Step 1: Extracting organization names ===\n")

# --------------------------------------------------------------------------
# 1. List all CRS network files
# --------------------------------------------------------------------------

all_nets <- list.files(
  path = "Data/clean_nets",
  full.names = TRUE,
  pattern = "\\.(rds|Rds)$"
)
all_nets <- all_nets[grepl("crs", all_nets, ignore.case = TRUE)]
cat(sprintf("Found %d CRS network files\n", length(all_nets)))

# --------------------------------------------------------------------------
# 2. Extract names from every network in every file
# --------------------------------------------------------------------------

records <- list()
record_idx <- 0

for (f in all_nets) {
  cat(sprintf("  Processing: %s\n", basename(f)))
  flush.console()

  network_data <- readRDS(f)
  networks <- network_data$networks

  for (net_name in names(networks)) {
    g <- networks[[net_name]]
    nm <- V(g)$name

    if (is.null(nm)) {
      cat(sprintf("    WARNING: %s has no vertex names, skipping\n", net_name))
      next
    }

    nm <- as.character(nm)

    # Parse state and year from network name (e.g., "hti_2005")
    parts <- strsplit(net_name, "_")[[1]]
    state <- tolower(parts[1])
    year <- as.integer(parts[length(parts)])

    for (name in nm) {
      record_idx <- record_idx + 1
      records[[record_idx]] <- data.frame(
        raw_name = name,
        state = state,
        year = year,
        stringsAsFactors = FALSE
      )
    }
  }
}

cat(sprintf("Extracted %d total name-state-year records\n", record_idx))

# --------------------------------------------------------------------------
# 3. Build raw_names.csv: one row per (raw_name, state, year) with count
# --------------------------------------------------------------------------

all_records <- do.call(rbind, records)

# Count how many networks each exact string appears in (per state-year)
raw_names <- aggregate(
  year ~ raw_name + state + year,
  data = all_records,
  FUN = length
)
# The aggregate above doesn't make sense for counting — let's do it properly
# We want: for each (raw_name, state, year), how many times does it appear?
# Since each org appears at most once per network, the count per (name, state, year)
# is always 1. The useful count is across all networks.

raw_names <- unique(all_records[, c("raw_name", "state", "year")])

# n_appearances: how many (state, year) pairs this raw_name appears in
name_counts <- as.data.frame(table(all_records$raw_name), stringsAsFactors = FALSE)
colnames(name_counts) <- c("raw_name", "n_appearances")

raw_names <- merge(raw_names, name_counts, by = "raw_name", all.x = TRUE)
raw_names <- raw_names[order(raw_names$raw_name, raw_names$state, raw_names$year), ]

cat(sprintf("raw_names.csv: %d rows, %d unique names\n",
            nrow(raw_names), length(unique(raw_names$raw_name))))

write.csv(raw_names, "Data/entity_resolution/raw_names.csv", row.names = FALSE)

# --------------------------------------------------------------------------
# 4. Build raw_names_summary.csv: deduplicated summary per raw_name
# --------------------------------------------------------------------------

summary_df <- data.frame(
  raw_name = character(),
  n_states = integer(),
  n_years = integer(),
  first_year = integer(),
  last_year = integer(),
  n_appearances = integer(),
  stringsAsFactors = FALSE
)

unique_names <- sort(unique(all_records$raw_name))

for (i in seq_along(unique_names)) {
  nm <- unique_names[i]
  subset <- all_records[all_records$raw_name == nm, ]
  summary_df[i, ] <- list(
    raw_name = nm,
    n_states = length(unique(subset$state)),
    n_years = length(unique(subset$year)),
    first_year = min(subset$year),
    last_year = max(subset$year),
    n_appearances = nrow(subset)
  )
}

summary_df <- summary_df[order(-summary_df$n_appearances), ]

write.csv(summary_df, "Data/entity_resolution/raw_names_summary.csv", row.names = FALSE)

# --------------------------------------------------------------------------
# 5. Report basic stats
# --------------------------------------------------------------------------

cat("\n=== Summary Statistics ===\n")
cat(sprintf("Total unique raw names:        %d\n", length(unique_names)))
cat(sprintf("Total (name, state, year):     %d\n", nrow(raw_names)))
cat(sprintf("Total CRS files processed:     %d\n", length(all_nets)))
cat(sprintf("Names appearing in 1 network:  %d\n",
            sum(summary_df$n_appearances == 1)))
cat(sprintf("Names appearing in 2+ states:  %d\n",
            sum(summary_df$n_states >= 2)))
cat(sprintf("Names appearing in 5+ years:   %d\n",
            sum(summary_df$n_years >= 5)))
cat(sprintf("Names appearing in 10+ years:  %d\n",
            sum(summary_df$n_years >= 10)))

cat("\nTop 20 most frequent names:\n")
top20 <- head(summary_df, 20)
for (i in seq_len(nrow(top20))) {
  cat(sprintf("  [%3d appearances, %2d states] %s\n",
              top20$n_appearances[i], top20$n_states[i], top20$raw_name[i]))
}

cat("\nNames with accented characters:\n")
accented <- unique_names[grepl("[^\\x00-\\x7F]", unique_names, perl = TRUE)]
cat(sprintf("  %d names contain non-ASCII characters\n", length(accented)))
if (length(accented) > 0) {
  cat("  First 20:\n")
  for (nm in head(accented, 20)) {
    cat(sprintf("    %s\n", nm))
  }
}

cat("\n=== Step 1 complete ===\n")
cat("Output files:\n")
cat("  Data/entity_resolution/raw_names.csv\n")
cat("  Data/entity_resolution/raw_names_summary.csv\n")
