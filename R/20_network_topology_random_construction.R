# ==========================================================
# 20_network_topology_random_construction.R
# PURPOSE:
#   Rebuild network- AND node-level panels after introducing
#   stochastic perturbations to the edge list *before* computing
#   metrics. Tests robustness of both network-level and
#   node-level DiD results to measurement error in network
#   construction.
#
# CHANGES FROM ORIGINAL:
#   1. Node-level attributes computed per perturbation
#   2. Node-level DiD analysis added (balanced panel, spillover)
#   3. Density removed from dv_focus (descriptive only)
#   4. Node DVs aligned with revised 17/18
# ==========================================================

library(igraph)
library(dplyr)
library(fixest)
library(readr)
library(ggplot2)

# ----------------------------------------------------------
# PROJECT ROOT
# ----------------------------------------------------------
PROJECT_ROOT <- "Ecosystem_of_Aid"
if (!grepl(PROJECT_ROOT, getwd())) setwd(PROJECT_ROOT)

# ----------------------------------------------------------
# Constants
# ----------------------------------------------------------
USAID_REGEX <- "United States Agency for International Development|\\bUSAID\\b|US AID"

# ----------------------------------------------------------
# Random construction settings
# ----------------------------------------------------------
out_dir <- "Data/results/robustness/haiti"

do_edge_perturb <- TRUE

perturb_drop_p    <- 0.05
perturb_rewire_p  <- 0.03
perturb_add_p     <- 0.01
perturb_jitter_sd <- 0.30

edge_perturb_scenario <- "combo"


# ----------------------------------------------------------
# Helpers (same as 17_network_topology.R)
# ----------------------------------------------------------

get_usaid_nodes <- function(g, regex = USAID_REGEX) {
  if (!inherits(g, "igraph") || vcount(g) == 0) return(character(0))
  nm <- V(g)$name
  if (is.null(nm)) return(character(0))
  nm <- as.character(nm)
  nm[grepl(regex, nm, ignore.case = TRUE)]
}

calc_fragmentation <- function(g) {
  if (!inherits(g, "igraph") || vcount(g) == 0) return(1)
  cc <- components(g, mode = "weak")
  1 - (max(cc$csize) / vcount(g))
}

calc_weighted_fragmentation <- function(g) {
  if (!inherits(g, "igraph") || vcount(g) == 0) return(1)
  cc <- components(g, mode = "weak")
  s <- cc$csize
  n <- vcount(g)
  1 - sum((s / n)^2)
}

calc_modularity_louvain <- function(g) {
  if (!inherits(g, "igraph") || vcount(g) < 2 || ecount(g) == 0) return(NA_real_)
  gu <- as_undirected(g, mode = "collapse")
  cl <- cluster_louvain(gu)
  modularity(cl)
}

partner_dependency_metrics <- function(g, donor_nodes) {
  donor_nodes <- intersect(donor_nodes, V(g)$name)
  if (length(donor_nodes) == 0) {
    return(list(n_partners = 0L,
                share_no_alt_funder = NA_real_,
                mean_partner_total_degree = NA_real_,
                mean_partner_in_degree = NA_real_,
                mean_partner_out_degree = NA_real_))
  }

  partners <- unique(unlist(lapply(donor_nodes, function(d) {
    V(g)[neighbors(g, d, mode = "out")]$name
  })))
  partners <- setdiff(partners, donor_nodes)

  if (length(partners) == 0) {
    return(list(n_partners = 0L,
                share_no_alt_funder = NA_real_,
                mean_partner_total_degree = NA_real_,
                mean_partner_in_degree = NA_real_,
                mean_partner_out_degree = NA_real_))
  }

  in_deg  <- degree(g, mode = "in");  names(in_deg)  <- V(g)$name
  out_deg <- degree(g, mode = "out"); names(out_deg) <- V(g)$name
  tot_deg <- degree(g, mode = "all"); names(tot_deg) <- V(g)$name

  alt_n <- vapply(partners, function(p) {
    funders <- V(g)[neighbors(g, p, mode = "in")]$name
    funders <- setdiff(funders, donor_nodes)
    length(funders)
  }, integer(1))

  list(
    n_partners = length(partners),
    share_no_alt_funder = mean(alt_n == 0),
    mean_partner_total_degree = mean(tot_deg[partners]),
    mean_partner_in_degree = mean(in_deg[partners]),
    mean_partner_out_degree = mean(out_deg[partners])
  )
}

cascade_drop_unfunded <- function(g, removed_seed) {
  removed <- unique(intersect(removed_seed, V(g)$name))
  if (length(removed) == 0) return(character(0))

  repeat {
    remaining <- setdiff(V(g)$name, removed)
    if (length(remaining) == 0) break

    new_fail <- character(0)
    for (v in remaining) {
      funders <- V(g)[neighbors(g, v, mode = "in")]$name
      if (length(funders) > 0 && all(funders %in% removed)) {
        new_fail <- c(new_fail, v)
      }
    }
    new_fail <- setdiff(unique(new_fail), removed)
    if (length(new_fail) == 0) break
    removed <- c(removed, new_fail)
  }

  unique(removed)
}

simulate_donor_withdrawal <- function(g, donor_nodes, scenario = c("direct","drop_unfunded")) {
  scenario <- match.arg(scenario)
  donor_nodes <- intersect(donor_nodes, V(g)$name)

  frag_pre  <- calc_fragmentation(g)
  wfrag_pre <- calc_weighted_fragmentation(g)
  cc_pre <- components(g, mode = "weak")
  comps_pre <- cc_pre$no
  giant_pre <- max(cc_pre$csize) / vcount(g)

  if (length(donor_nodes) == 0) {
    return(list(frag_delta = NA_real_,
                wfrag_delta = NA_real_,
                component_increase = NA_integer_,
                giant_share_delta = NA_real_,
                n_removed_total = NA_integer_))
  }

  removed <- donor_nodes
  if (scenario == "drop_unfunded") removed <- cascade_drop_unfunded(g, donor_nodes)

  g2 <- delete_vertices(g, removed)

  if (vcount(g2) == 0) {
    frag_post <- 1; wfrag_post <- 1; comps_post <- 0; giant_post <- 0
  } else {
    frag_post  <- calc_fragmentation(g2)
    wfrag_post <- calc_weighted_fragmentation(g2)
    cc_post <- components(g2, mode = "weak")
    comps_post <- cc_post$no
    giant_post <- max(cc_post$csize) / vcount(g2)
  }

  list(
    frag_delta = frag_post - frag_pre,
    wfrag_delta = wfrag_post - wfrag_pre,
    component_increase = comps_post - comps_pre,
    giant_share_delta = giant_post - giant_pre,
    n_removed_total = length(removed)
  )
}

# ----------------------------------------------------------
# Load networks
# ----------------------------------------------------------
dir.create("Data/results/robustness/permute", showWarnings = FALSE, recursive = TRUE)
all_nets <- list.files(path = "Data/clean_nets", full.names = TRUE, pattern = "Rds")
all_nets <- all_nets[grepl("crs", all_nets)]


# ##############################################################
# PART 1: PERTURB NETWORKS AND COMPUTE PANELS
# ##############################################################

for (tag in 1:100) {

  seed_now <- 123 + tag
  set.seed(seed_now)

  network_attributes_list <- list()
  node_attributes_list <- list()

  for (z in seq_along(all_nets)) {

    network_data <- readRDS(all_nets[z])
    networks <- network_data$networks

    index_name <- gsub("_(.*)", "", names(networks)[1])
    index_name <- tolower(index_name)

    for (i in seq_along(networks)) {

      g <- networks[[i]]
      g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

      if (is.null(V(g)$name)) V(g)$name <- as.character(seq_len(vcount(g)))

      yr_txt <- gsub("[^0-9]", "", names(networks)[i])
      yr <- suppressWarnings(as.numeric(yr_txt))
      if (!is.finite(yr)) next

      # ==========================================================
      # RANDOM CONSTRUCTION: perturb edge list BEFORE metrics
      # ==========================================================
      if (do_edge_perturb && ecount(g) > 0) {

        ed <- igraph::as_data_frame(g, what = "edges")

        if (!("from" %in% names(ed))) names(ed)[1] <- "from"
        if (!("to" %in% names(ed))) names(ed)[2] <- "to"
        if (!("weight" %in% names(ed))) ed$weight <- NA_real_

        nodes_all <- as.character(V(g)$name)

        # 1) DROP EDGES
        if (edge_perturb_scenario %in% c("drop_edges", "combo")) {
          keep <- stats::runif(nrow(ed)) > perturb_drop_p
          ed <- ed[keep, , drop = FALSE]
        }

        # 2) REWIRE TO-ENDPOINT
        if (nrow(ed) > 0 && edge_perturb_scenario %in% c("rewire_to", "combo")) {
          rw <- stats::runif(nrow(ed)) < perturb_rewire_p
          if (any(rw)) {
            ed$to[rw] <- sample(nodes_all, size = sum(rw), replace = TRUE)
            ed <- ed[ed$from != ed$to, , drop = FALSE]
          }
        }

        # 3) ADD RANDOM NOISE EDGES
        if (edge_perturb_scenario %in% c("add_noise_edges", "combo")) {
          n_add <- max(1, round(max(1, nrow(ed)) * perturb_add_p))
          add_from <- sample(nodes_all, size = n_add, replace = TRUE)
          add_to <- sample(nodes_all, size = n_add, replace = TRUE)

          add_df <- data.frame(from = add_from, to = add_to, weight = NA_real_)
          add_df <- add_df[add_df$from != add_df$to, , drop = FALSE]

          ed <- rbind(ed, add_df)
        }

        # 4) JITTER WEIGHTS (if present)
        if (nrow(ed) > 0 && edge_perturb_scenario %in% c("jitter_weights", "combo")) {
          if (!all(is.na(ed$weight))) {
            mult <- exp(stats::rnorm(nrow(ed), mean = 0, sd = perturb_jitter_sd))
            ed$weight <- ed$weight * mult
          }
        }

        # Rebuild graph
        if (nrow(ed) == 0) {
          g <- igraph::make_empty_graph(n = length(nodes_all), directed = TRUE)
          igraph::V(g)$name <- nodes_all
        } else {
          g_new <- igraph::graph_from_data_frame(ed[, c("from", "to", "weight")],
                                                  directed = TRUE, vertices = nodes_all)
          g_new <- igraph::simplify(g_new, remove.multiple = TRUE, remove.loops = TRUE)

          if ("weight" %in% igraph::edge_attr_names(g_new)) {
            w2 <- igraph::edge_attr(g_new, "weight")
            if (all(is.na(w2))) g_new <- igraph::delete_edge_attr(g_new, "weight")
          }

          g <- g_new
        }
      }

      # ==========================================================
      # NETWORK-LEVEL METRICS (on perturbed graph)
      # ==========================================================
      usaid_nodes <- get_usaid_nodes(g)
      pm_usaid <- partner_dependency_metrics(g, usaid_nodes)
      wd_usaid <- simulate_donor_withdrawal(g, usaid_nodes, scenario = "drop_unfunded")

      network_attributes <- data.frame(
        tag = tag,
        seed = seed_now,
        state = index_name,
        year = yr,

        num_nodes = vcount(g),
        num_edges = ecount(g),
        density = edge_density(g, loops = FALSE),

        fragmentation = calc_fragmentation(g),
        weighted_fragmentation = calc_weighted_fragmentation(g),
        modularity = calc_modularity_louvain(g),

        usaid_present = as.integer(length(usaid_nodes) > 0),

        usaid_n_partners = pm_usaid$n_partners,
        usaid_share_no_alt_funder = pm_usaid$share_no_alt_funder,
        usaid_mean_partner_total_degree = pm_usaid$mean_partner_total_degree,
        usaid_mean_partner_in_degree = pm_usaid$mean_partner_in_degree,
        usaid_mean_partner_out_degree = pm_usaid$mean_partner_out_degree,

        usaid_frag_delta_drop_unfunded = wd_usaid$frag_delta,
        usaid_wfrag_delta_drop_unfunded = wd_usaid$wfrag_delta,
        usaid_component_increase_drop_unfunded = wd_usaid$component_increase,
        usaid_giant_share_delta_drop_unfunded = wd_usaid$giant_share_delta,
        usaid_n_removed_total_drop_unfunded = wd_usaid$n_removed_total,

        rc_scenario = edge_perturb_scenario,
        rc_drop_p = perturb_drop_p,
        rc_rewire_p = perturb_rewire_p,
        rc_add_p = perturb_add_p,
        rc_jitter_sd = perturb_jitter_sd,

        stringsAsFactors = FALSE
      )

      network_attributes_list[[length(network_attributes_list) + 1]] <- network_attributes

      # ==========================================================
      # NODE-LEVEL ATTRIBUTES (on perturbed graph)
      # ==========================================================
      in_deg  <- degree(g, mode = "in")
      out_deg <- degree(g, mode = "out")
      tot_deg <- degree(g, mode = "all")

      is_usaid_vec <- as.integer(grepl(USAID_REGEX, V(g)$name, ignore.case = TRUE))

      # tied_to_usaid: any-mode neighbor of USAID
      tied_to_usaid_vec <- rep(0L, vcount(g))
      if (length(usaid_nodes) > 0) {
        usaid_neighbors <- unique(unlist(lapply(usaid_nodes, function(x)
          c(neighbors(g, x, mode = "out")$name, neighbors(g, x, mode = "in")$name))))
        tied_to_usaid_vec <- as.integer(V(g)$name %in% usaid_neighbors)
      }

      # usaid_funded: receives direct funding from USAID (out-neighbor)
      usaid_funded_vec <- rep(0L, vcount(g))
      if (length(usaid_nodes) > 0) {
        usaid_out_partners <- unique(unlist(lapply(usaid_nodes, function(x)
          V(g)[neighbors(g, x, mode = "out")]$name)))
        usaid_out_partners <- setdiff(usaid_out_partners, usaid_nodes)
        usaid_funded_vec <- as.integer(V(g)$name %in% usaid_out_partners)
      }

      # n_funders and n_non_usaid_funders
      n_funders_vec <- in_deg
      n_non_usaid_funders_vec <- rep(NA_integer_, vcount(g))
      for (vi in seq_len(vcount(g))) {
        in_nbrs <- V(g)[neighbors(g, V(g)$name[vi], mode = "in")]$name
        n_non_usaid_funders_vec[vi] <- sum(!in_nbrs %in% usaid_nodes)
      }

      # sole_usaid_funded
      sole_usaid_funded_vec <- as.integer(usaid_funded_vec == 1L & n_non_usaid_funders_vec == 0L)

      # usaid_exposure_share: fraction of neighbors that are USAID or USAID-funded
      usaid_connected_set <- unique(c(usaid_nodes, V(g)$name[usaid_funded_vec == 1L]))
      usaid_exposure_share_vec <- rep(0, vcount(g))
      for (vi in seq_len(vcount(g))) {
        all_nbrs <- V(g)[neighbors(g, V(g)$name[vi], mode = "all")]$name
        if (length(all_nbrs) > 0) {
          usaid_exposure_share_vec[vi] <- mean(all_nbrs %in% usaid_connected_set)
        }
      }

      node_attributes <- data.frame(
        tag = tag,
        state = index_name,
        year = yr,
        node_name = V(g)$name,

        is_usaid = is_usaid_vec,
        tied_to_usaid = tied_to_usaid_vec,
        usaid_funded = usaid_funded_vec,
        sole_usaid_funded = sole_usaid_funded_vec,

        in_degree = in_deg,
        out_degree = out_deg,
        total_degree = tot_deg,
        n_funders = n_funders_vec,
        n_non_usaid_funders = n_non_usaid_funders_vec,
        usaid_exposure_share = usaid_exposure_share_vec,

        stringsAsFactors = FALSE
      )

      node_attributes_list[[length(node_attributes_list) + 1]] <- node_attributes
    }
  }

  network_attributes_df <- dplyr::bind_rows(network_attributes_list)
  node_attributes_df    <- dplyr::bind_rows(node_attributes_list)

  saveRDS(
    network_attributes_df,
    file = file.path("Data/results/robustness/permute/",
                     paste0("network_attributes_directed_random_", tag, ".rds"))
  )
  saveRDS(
    node_attributes_df,
    file = file.path("Data/results/robustness/permute/",
                     paste0("node_attributes_directed_random_", tag, ".rds"))
  )

  message("Saved random-construction outputs | tag = ", tag, " | seed = ", seed_now)
}


# ##############################################################
# PART 2: NETWORK-LEVEL DiD ON PERTURBED DATA
# ##############################################################

treat_state <- "hti"
event_year  <- 2010
windows     <- c(2, 3, 4, 5)
include_full_sample <- TRUE
drop_states <- character(0)

perm_files <- list.files(
  path = "Data/results/robustness/permute/",
  pattern = "^network_attributes_directed_random_[0-9]+\\.rds$",
  full.names = TRUE
)

# DVs: density removed, aligned with manuscript
net_dv_focus <- c(
  "fragmentation",
  "weighted_fragmentation",
  "modularity",
  "usaid_frag_delta_drop_unfunded",
  "usaid_n_partners"
)

net_rows <- list()

for (f in perm_files) {

  message("\rReading: ", f)

  dat0 <- readRDS(f)

  if (!("state" %in% names(dat0))) stop("Missing column: state in file: ", f)
  if (!("year"  %in% names(dat0))) stop("Missing column: year in file: ", f)

  dat0$state   <- tolower(dat0$state)
  dat0$treated <- as.integer(dat0$state == treat_state)
  dat0$post    <- as.integer(dat0$year >= event_year)
  dat0$tp      <- dat0$treated * dat0$post

  if (!("tag" %in% names(dat0))) {
    tag_txt <- gsub("^.*network_attributes_directed_random_([0-9]+)\\.rds$", "\\1", f)
    dat0$tag <- as.integer(tag_txt)
  }
  if (!("seed" %in% names(dat0))) dat0$seed <- NA_integer_

  dvs_here <- net_dv_focus[net_dv_focus %in% names(dat0)]
  if (length(dvs_here) == 0) { message("  -> Skipping."); next }

  for (dv in dvs_here) {

    # Full sample
    if (include_full_sample) {
      dat <- dat0[!dat0$state %in% drop_states & !is.na(dat0[[dv]]), ]
      fit <- tryCatch(
        feols(as.formula(paste0(dv, " ~ tp | state + year")),
              data = dat, cluster = ~state),
        error = function(e) NULL
      )
      if (!is.null(fit) && ("tp" %in% names(coef(fit)))) {
        net_rows[[length(net_rows) + 1]] <- data.frame(
          tag = unique(dat0$tag)[1], seed = unique(dat0$seed)[1],
          level = "network", dep_var = dv, window = "all",
          estimate = as.numeric(coef(fit)["tp"]),
          std_error = as.numeric(se(fit)["tp"]),
          p_value = as.numeric(pvalue(fit)["tp"]),
          n = nobs(fit), n_units = length(unique(dat$state)),
          stringsAsFactors = FALSE
        )
      }
    }

    # Windowed
    for (w in windows) {
      yr_range <- (event_year - w):(event_year + w)
      dat <- dat0[!dat0$state %in% drop_states & dat0$year %in% yr_range & !is.na(dat0[[dv]]), ]

      fit <- tryCatch(
        feols(as.formula(paste0(dv, " ~ tp | state + year")),
              data = dat, cluster = ~state),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      if (!("tp" %in% names(coef(fit)))) next

      net_rows[[length(net_rows) + 1]] <- data.frame(
        tag = unique(dat0$tag)[1], seed = unique(dat0$seed)[1],
        level = "network", dep_var = dv, window = paste0("+-", w),
        estimate = as.numeric(coef(fit)["tp"]),
        std_error = as.numeric(se(fit)["tp"]),
        p_value = as.numeric(pvalue(fit)["tp"]),
        n = nobs(fit), n_units = length(unique(dat$state)),
        stringsAsFactors = FALSE
      )
    }
  }
}


# ##############################################################
# PART 3: NODE-LEVEL DiD ON PERTURBED DATA
# ##############################################################

node_perm_files <- list.files(
  path = "Data/results/robustness/permute/",
  pattern = "^node_attributes_directed_random_[0-9]+\\.rds$",
  full.names = TRUE
)

node_dv_focus <- c("in_degree", "sole_usaid_funded", "n_funders",
                   "n_non_usaid_funders", "total_degree")

node_rows <- list()

for (f in node_perm_files) {

  message("\rReading node file: ", f)

  nd0 <- readRDS(f)
  nd0$state <- tolower(nd0$state)

  if (!("tag" %in% names(nd0))) {
    tag_txt <- gsub("^.*node_attributes_directed_random_([0-9]+)\\.rds$", "\\1", f)
    nd0$tag <- as.integer(tag_txt)
  }

  # Construct balanced panel for this perturbation
  pre_orgs_p  <- unique(nd0[nd0$year < event_year, c("state", "node_name")])
  post_orgs_p <- unique(nd0[nd0$year >= event_year, c("state", "node_name")])
  bal_orgs_p  <- merge(pre_orgs_p, post_orgs_p, by = c("state", "node_name"))
  nd_bal <- merge(nd0, bal_orgs_p, by = c("state", "node_name"))

  nd_bal$treated <- as.integer(nd_bal$state == treat_state)
  nd_bal$post    <- as.integer(nd_bal$year >= event_year)
  nd_bal$tp      <- nd_bal$treated * nd_bal$post
  nd_bal$org_id  <- paste0(nd_bal$state, "||", nd_bal$node_name)

  if (nrow(nd_bal) < 20) { message("  -> Too few rows, skipping."); next }

  # Pre-treatment exposure
  nd_pre <- nd_bal[nd_bal$year < event_year, ]

  pre_funded_p <- aggregate(usaid_funded ~ state + node_name,
    data = nd_pre, FUN = max, na.rm = TRUE)
  names(pre_funded_p)[3] <- "pre_usaid_funded"

  pre_exp_p <- aggregate(usaid_exposure_share ~ state + node_name,
    data = nd_pre, FUN = mean, na.rm = TRUE)
  names(pre_exp_p)[3] <- "pre_usaid_exposure_share"

  nd_bal <- merge(nd_bal, pre_funded_p, by = c("state", "node_name"), all.x = TRUE)
  nd_bal <- merge(nd_bal, pre_exp_p,    by = c("state", "node_name"), all.x = TRUE)

  dvs_here <- node_dv_focus[node_dv_focus %in% names(nd_bal)]
  if (length(dvs_here) == 0) next

  current_tag <- unique(nd0$tag)[1]

  for (dv in dvs_here) {
    dat <- nd_bal[!is.na(nd_bal[[dv]]), ]

    # Baseline ATE (full sample)
    fit_base <- tryCatch(
      feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
            data = dat, cluster = ~org_id),
      error = function(e) NULL
    )
    if (!is.null(fit_base) && "tp" %in% names(coef(fit_base))) {
      node_rows[[length(node_rows) + 1]] <- data.frame(
        tag = current_tag, level = "node", spec = "baseline",
        dep_var = dv, window = "all",
        estimate = as.numeric(coef(fit_base)["tp"]),
        std_error = as.numeric(se(fit_base)["tp"]),
        p_value = as.numeric(pvalue(fit_base)["tp"]),
        n = nobs(fit_base), n_units = length(unique(dat$org_id)),
        stringsAsFactors = FALSE
      )
    }

    # Heterogeneous: direct USAID funding
    fit_het <- tryCatch(
      feols(as.formula(paste0(dv, " ~ tp + tp:pre_usaid_funded | org_id + year")),
            data = dat, cluster = ~org_id),
      error = function(e) NULL
    )
    if (!is.null(fit_het)) {
      for (cf in names(coef(fit_het))) {
        node_rows[[length(node_rows) + 1]] <- data.frame(
          tag = current_tag, level = "node", spec = "heterogeneous_direct",
          dep_var = dv, window = "all",
          estimate = as.numeric(coef(fit_het)[cf]),
          std_error = as.numeric(se(fit_het)[cf]),
          p_value = as.numeric(pvalue(fit_het)[cf]),
          n = nobs(fit_het), n_units = length(unique(dat$org_id)),
          stringsAsFactors = FALSE
        )
      }
    }

    # Spillover: exposure share
    fit_sp <- tryCatch(
      feols(as.formula(paste0(dv, " ~ tp + tp:pre_usaid_exposure_share | org_id + year")),
            data = dat, cluster = ~org_id),
      error = function(e) NULL
    )
    if (!is.null(fit_sp)) {
      for (cf in names(coef(fit_sp))) {
        node_rows[[length(node_rows) + 1]] <- data.frame(
          tag = current_tag, level = "node", spec = "spillover_exposure",
          dep_var = dv, window = "all",
          estimate = as.numeric(coef(fit_sp)[cf]),
          std_error = as.numeric(se(fit_sp)[cf]),
          p_value = as.numeric(pvalue(fit_sp)[cf]),
          n = nobs(fit_sp), n_units = length(unique(dat$org_id)),
          stringsAsFactors = FALSE
        )
      }
    }

    # Windowed baseline
    for (w in windows) {
      yr_range <- (event_year - w):(event_year + w)
      dat_w <- nd_bal[nd_bal$year %in% yr_range & !is.na(nd_bal[[dv]]), ]

      fit_w <- tryCatch(
        feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
              data = dat_w, cluster = ~org_id),
        error = function(e) NULL
      )
      if (is.null(fit_w)) next
      if (!("tp" %in% names(coef(fit_w)))) next

      node_rows[[length(node_rows) + 1]] <- data.frame(
        tag = current_tag, level = "node", spec = "baseline",
        dep_var = dv, window = paste0("+-", w),
        estimate = as.numeric(coef(fit_w)["tp"]),
        std_error = as.numeric(se(fit_w)["tp"]),
        p_value = as.numeric(pvalue(fit_w)["tp"]),
        n = nobs(fit_w), n_units = length(unique(dat_w$org_id)),
        stringsAsFactors = FALSE
      )
    }
  }
}


# ##############################################################
# PART 4: COMBINE AND PLOT
# ##############################################################

# --- Network-level results ---
net_tbl <- bind_rows(net_rows)
net_tbl$lb <- net_tbl$estimate - 1.96 * net_tbl$std_error
net_tbl$ub <- net_tbl$estimate + 1.96 * net_tbl$std_error
net_tbl <- net_tbl %>%
  arrange(dep_var, tag, factor(window, levels = c("+-2","+-3","+-4","+-5","all")))

write_csv(net_tbl, file.path(out_dir, "haiti_random_construction_twfe_results_long.csv"))
saveRDS(net_tbl, file.path(out_dir, "haiti_random_construction_twfe_results_long.rds"))

# --- Node-level results ---
node_tbl <- bind_rows(node_rows)
node_tbl$lb <- node_tbl$estimate - 1.96 * node_tbl$std_error
node_tbl$ub <- node_tbl$estimate + 1.96 * node_tbl$std_error

write_csv(node_tbl, file.path(out_dir, "haiti_random_construction_node_results_long.csv"))
saveRDS(node_tbl, file.path(out_dir, "haiti_random_construction_node_results_long.rds"))

# --- Network-level plots ---
net_plot_dvs <- c("fragmentation", "weighted_fragmentation",
                  "usaid_frag_delta_drop_unfunded", "usaid_n_partners")
net_plot_labels <- c("Fragmentation", "Weighted fragmentation",
                     "Drop unfunded", "USAID partners")

net_tbl_plot <- net_tbl %>%
  filter(dep_var %in% net_plot_dvs) %>%
  mutate(dep_var = factor(dep_var, levels = net_plot_dvs, labels = net_plot_labels))

dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

net_permute_est <- ggplot(net_tbl_plot, aes(x = estimate, fill = window)) +
  geom_density(alpha = 0.35) +
  facet_wrap(~dep_var, scales = "free") +
  theme_minimal() +
  labs(fill = "Time Window", title = "Network-level: Permutation coefficients",
       y = "Density", x = "Coefficient") +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "bottom")

p_values_net <- ggplot(net_tbl_plot, aes(x = p_value, y = window, fill = window)) +
  geom_boxplot() +
  facet_wrap(~dep_var, scales = "free") +
  theme_minimal() +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  labs(fill = "Time Window", title = "Network-level: p-values across perturbations",
       x = "p-values") +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "plots/p_value_network_measure_noise.png"),
       plot = p_values_net, width = 10, height = 5, dpi = 200)
ggsave(file.path(out_dir, "plots/estimate_network_measure_noise.png"),
       plot = net_permute_est, width = 10, height = 5, dpi = 200)

# --- Node-level plots ---
node_tbl_baseline <- node_tbl %>% filter(spec == "baseline")

if (nrow(node_tbl_baseline) > 0) {
  node_permute_est <- ggplot(node_tbl_baseline, aes(x = estimate, fill = window)) +
    geom_density(alpha = 0.35) +
    facet_wrap(~dep_var, scales = "free") +
    theme_minimal() +
    labs(fill = "Time Window", title = "Node-level: Permutation coefficients (baseline)",
         y = "Density", x = "Coefficient") +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = "bottom")

  p_values_node <- ggplot(node_tbl_baseline, aes(x = p_value, y = window, fill = window)) +
    geom_boxplot() +
    facet_wrap(~dep_var, scales = "free") +
    theme_minimal() +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    labs(fill = "Time Window", title = "Node-level: p-values across perturbations",
         x = "p-values") +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = "none")

  ggsave(file.path(out_dir, "plots/p_value_node_measure_noise.png"),
         plot = p_values_node, width = 10, height = 5, dpi = 200)
  ggsave(file.path(out_dir, "plots/estimate_node_measure_noise.png"),
         plot = node_permute_est, width = 10, height = 5, dpi = 200)
}

cat("\n=== Random construction robustness complete (network + node level) ===\n")
