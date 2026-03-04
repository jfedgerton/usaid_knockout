###############################################################################
# Investigate GIZ/JICA/AFD 0-partner issue - check file structure first
###############################################################################

library(igraph)

# Check structure of harmonized file
nets <- readRDS("Data/clean_nets_harmonized/hti_result_crs.Rds")
cat("Class:", class(nets), "\n")
cat("Length:", length(nets), "\n")
cat("Names:", head(names(nets), 5), "\n")

# Check what's inside
cat("\nFirst element class:", class(nets[[1]]), "\n")

# If it's a list of lists, check deeper
if (!is.igraph(nets[[1]])) {
  cat("First element structure:\n")
  str(nets[[1]], max.level = 1)

  # Try named access
  if (is.list(nets[[1]])) {
    cat("\nFirst sub-element class:", class(nets[[1]][[1]]), "\n")
  }
}

# Try getting igraph object
get_graph <- function(x) {
  if (is.igraph(x)) return(x)
  if (is.list(x)) {
    for (item in x) {
      if (is.igraph(item)) return(item)
    }
  }
  return(NULL)
}

g <- get_graph(nets[[1]])
if (is.null(g)) {
  g <- get_graph(nets)
}

if (!is.null(g)) {
  cat("\nFound igraph object!\n")
  cat("Nodes:", vcount(g), "\n")
  cat("Edges:", ecount(g), "\n")
} else {
  cat("\nCould not find igraph object. Trying raw file...\n")
  nets_raw <- readRDS("Data/clean_nets/hti_result_crs.Rds")
  cat("Raw class:", class(nets_raw), "\n")
  cat("Raw length:", length(nets_raw), "\n")
  cat("Raw names:", head(names(nets_raw), 5), "\n")
  cat("First element class:", class(nets_raw[[1]]), "\n")

  if (is.igraph(nets_raw[[1]])) {
    g <- nets_raw[[1]]
    cat("Nodes:", vcount(g), "\n")
    cat("Sample names:", head(V(g)$name, 5), "\n")
  }
}
