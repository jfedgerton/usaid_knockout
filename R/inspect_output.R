###############################################################################
# Inspect model output files
###############################################################################

cat("=== NETWORK-LEVEL PANEL ===\n")
d <- readRDS("Data/clean_nets/network_attributes_directed.rds")
cat("Class:", class(d), "\n")
cat("Dimensions:", nrow(d), "rows x", ncol(d), "cols\n\n")

cat("Column names:\n")
print(names(d))

cat("\nColumn types:\n")
print(sapply(d, class))

cat("\n--- Summary stats for numeric columns ---\n")
nums <- d[, sapply(d, is.numeric)]
print(summary(nums))

cat("\n--- Unique states ---\n")
if ("state" %in% names(d)) {
  cat("N unique states:", length(unique(d$state)), "\n")
  print(sort(unique(d$state)))
}

cat("\n--- Year range ---\n")
if ("year" %in% names(d)) {
  cat("Years:", min(d$year), "to", max(d$year), "\n")
  print(table(d$year))
}

cat("\n--- Rows per state ---\n")
if ("state" %in% names(d)) {
  print(table(d$state))
}

cat("\n--- NA counts per column ---\n")
na_counts <- colSums(is.na(d))
print(na_counts[na_counts > 0])
if (all(na_counts == 0)) cat("No NAs in any column!\n")

cat("\n--- First 10 rows ---\n")
print(head(d, 10))

cat("\n\n========================================================\n")
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
if ("node_name" %in% names(n)) {
  cat("N unique node names:", length(unique(n$node_name)), "\n")
} else if ("name" %in% names(n)) {
  cat("N unique node names:", length(unique(n$name)), "\n")
}

cat("\n--- NA counts per column ---\n")
na_counts_n <- colSums(is.na(n))
print(na_counts_n[na_counts_n > 0])
if (all(na_counts_n == 0)) cat("No NAs in any column!\n")

cat("\n--- Summary stats for numeric columns ---\n")
nums_n <- n[, sapply(n, is.numeric)]
print(summary(nums_n))

cat("\n--- First 10 rows ---\n")
print(head(n, 10))

cat("\n--- Sample of USAID-related nodes ---\n")
name_col <- if ("node_name" %in% names(n)) "node_name" else "name"
usaid_rows <- n[grepl("usaid|united states agency", n[[name_col]], ignore.case = TRUE), ]
cat("N USAID rows:", nrow(usaid_rows), "\n")
if (nrow(usaid_rows) > 0) {
  cat("Sample:\n")
  print(head(usaid_rows, 10))
}

cat("\n=== INSPECTION COMPLETE ===\n")
