# ==========================================================
# PURPOSE:
#   Rebuild network- and node-level panels after introducing
#   stochastic perturbations to the edge list *before* computing
#   network metrics.
# NOTES:
#   - This is a robustness/stress-test construction, not the
#     primary data pipeline.
#   - The perturbation is applied per state-year graph.
#
# ==========================================================

library(igraph)
library(dplyr)
library(fixest)
library(readr)
library(ggplot2)
# ----------------------------------------------------------
# PROJECT ROOT (edit if needed)
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

## Export path
out_dir <- "Data/results/robustness/haiti"

# Tag appended to output filenames


# Turn on/off perturbations
do_edge_perturb <- TRUE

# Perturbation knobs (used in the "combo" scenario below)
perturb_drop_p    <- 0.05   # Bernoulli drop per existing edge
perturb_rewire_p  <- 0.03   # Bernoulli rewire per remaining edge (rewire TO endpoint)
perturb_add_p     <- 0.01   # add ~ round(m * add_p) random edges where m = current edges
perturb_jitter_sd <- 0.30   # multiplicative lognormal jitter for weights (if present)

# One of: "drop_edges","rewire_to","add_noise_edges","jitter_weights","combo"
edge_perturb_scenario <- "combo"


# ----------------------------------------------------------
# Helpers (same as your baseline construction)
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
dir.create("Data/results/robustness/permute", showWarnings = F, recursive = T)
all_nets <- list.files(path = "Data/clean_nets", full.names = TRUE, pattern = "Rds")
all_nets <- all_nets[grepl("crs", all_nets)]



# ----------------------------------------------------------
# Main loop over simulations
# ----------------------------------------------------------
for (tag in 1:100) {
  
  seed_now <- 123 + tag
  set.seed(seed_now)
  
  network_attributes_list <- list()
  
  for (z in seq_along(all_nets)) {
    
    network_data <- readRDS(all_nets[z])
    networks <- network_data$networks
    
    # state abbrev inferred from first network name in this file
    index_name <- gsub("_(.*)", "", names(networks)[1])
    index_name <- tolower(index_name)
    
    for (i in seq_along(networks)) {
      
      g <- networks[[i]]
      g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
      
      if (is.null(V(g)$name)) V(g)$name <- as.character(seq_len(vcount(g)))
      
      # Parse year from network name safely
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
        
        # Rebuild graph (robust to ed becoming empty)
        if (nrow(ed) == 0) {
          g <- igraph::make_empty_graph(n = length(nodes_all), directed = TRUE)
          igraph::V(g)$name <- nodes_all
        } else {
          g_new <- igraph::graph_from_data_frame(ed[, c("from", "to", "weight")],
                                                 directed = TRUE, vertices = nodes_all)
          g_new <- igraph::simplify(g_new, remove.multiple = TRUE, remove.loops = TRUE)
          
          # Safe remove weight attr if all NA (avoid E(g)$weight <- NULL)
          if ("weight" %in% igraph::edge_attr_names(g_new)) {
            w2 <- igraph::edge_attr(g_new, "weight")
            if (all(is.na(w2))) g_new <- igraph::delete_edge_attr(g_new, "weight")
          }
          
          g <- g_new
        }
      }
      
      # --------------------------
      # Compute key metrics
      # (assumes these helper functions exist in your session/script)
      # --------------------------
      # REQUIRED: get_usaid_nodes(), calc_fragmentation(), calc_weighted_fragmentation(), calc_modularity_louvain()
      usaid_nodes <- get_usaid_nodes(g)
      
      # ---- USAID partner metrics ----
      pm_usaid <- partner_dependency_metrics(g, usaid_nodes)
      
      # ---- USAID withdrawal deltas (drop_unfunded) ----
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
        
        # --- NEW: USAID partner + dependency metrics ---
        usaid_n_partners = pm_usaid$n_partners,
        usaid_share_no_alt_funder = pm_usaid$share_no_alt_funder,
        usaid_mean_partner_total_degree = pm_usaid$mean_partner_total_degree,
        usaid_mean_partner_in_degree = pm_usaid$mean_partner_in_degree,
        usaid_mean_partner_out_degree = pm_usaid$mean_partner_out_degree,
        
        # --- NEW: USAID withdrawal deltas (drop_unfunded) ---
        usaid_frag_delta_drop_unfunded = wd_usaid$frag_delta,
        usaid_wfrag_delta_drop_unfunded = wd_usaid$wfrag_delta,
        usaid_component_increase_drop_unfunded = wd_usaid$component_increase,
        usaid_giant_share_delta_drop_unfunded = wd_usaid$giant_share_delta,
        usaid_n_removed_total_drop_unfunded = wd_usaid$n_removed_total,
        
        # Traceability
        rc_scenario = edge_perturb_scenario,
        rc_drop_p = perturb_drop_p,
        rc_rewire_p = perturb_rewire_p,
        rc_add_p = perturb_add_p,
        rc_jitter_sd = perturb_jitter_sd,
        
        stringsAsFactors = FALSE
      )
      
      network_attributes_list[[length(network_attributes_list) + 1]] <- network_attributes
    }
  }
  
  network_attributes_df <- dplyr::bind_rows(network_attributes_list)
  
  saveRDS(
    network_attributes_df,
    file = file.path("Data/results/robustness/permute/", paste0("network_attributes_directed_random_", tag, ".rds"))
  )
  
  message("Saved random-construction outputs | tag = ", tag, " | seed = ", seed_now)
}

treat_state <- "hti"
event_year  <- 2010

# ----------------------------------------------------------
# Windows (match Script 18)
# ----------------------------------------------------------
windows <- c(2, 3, 4, 5)
include_full_sample <- TRUE

# Optional: exclude states (leave empty unless you used this elsewhere)
drop_states <- character(0)

# ----------------------------------------------------------
# Find permutation files
# ----------------------------------------------------------
perm_files <- list.files(
  path = "Data/results/robustness/permute/",
  pattern = "^network_attributes_directed_random_[0-9]+\\.rds$",
  full.names = TRUE
)


# ----------------------------------------------------------
# Outcomes to estimate (edit as needed)
# Keep this conservative unless your random pipeline computed the full set.
# ----------------------------------------------------------
dv_focus <- c(
  "fragmentation",
  "weighted_fragmentation",
  "modularity",
  "usaid_frag_delta_drop_unfunded",
  "usaid_n_partners",
  "density"
)

# ----------------------------------------------------------
# Storage for results
# ----------------------------------------------------------
rows <- list()

# ----------------------------------------------------------
# Main loop over permutation tags/files
# ----------------------------------------------------------
for (f in perm_files) {
  
  message("\rReading: ", f)
  
  dat0 <- readRDS(f)
  
  # Basic sanity checks
  if (!("state" %in% names(dat0))) stop("Missing column: state in file: ", f)
  if (!("year"  %in% names(dat0))) stop("Missing column: year in file: ", f)
  
  dat0 <- dat0 %>%
    mutate(
      state = tolower(state),
      treated = as.integer(state == treat_state),
      post    = as.integer(year >= event_year),
      tp      = treated * post
    )
  
  # If tag/seed not present, infer tag from filename
  if (!("tag" %in% names(dat0))) {
    tag_txt <- gsub("^.*network_attributes_directed_random_([0-9]+)\\.rds$", "\\1", f)
    dat0$tag <- as.integer(tag_txt)
  }
  if (!("seed" %in% names(dat0))) {
    dat0$seed <- NA_integer_
  }
  
  # Keep only dv names that exist in this file
  dvs_here <- dv_focus[dv_focus %in% names(dat0)]
  if (length(dvs_here) == 0) {
    message("  -> Skipping (no dv_focus columns found in file).")
    next
  }
  
  # --------------------------------------------------------
  # Run models for each DV and window
  # --------------------------------------------------------
  for (dv in dvs_here) {
    
    # Full sample
    if (include_full_sample) {
      
      dat <- dat0 %>%
        filter(!(state %in% drop_states)) %>%
        filter(!is.na(.data[[dv]]))
      
      fit <- tryCatch(
        feols(as.formula(paste0(dv, " ~ tp | state + year")),
              data = dat, cluster = ~state),
        error = function(e) NULL
      )
      
      if (!is.null(fit) && ("tp" %in% names(coef(fit)))) {
        rows[[length(rows) + 1]] <- data.frame(
          tag = unique(dat0$tag)[1],
          seed = unique(dat0$seed)[1],
          dep_var = dv,
          window = "all",
          estimate = as.numeric(coef(fit)["tp"]),
          std_error = as.numeric(se(fit)["tp"]),
          p_value = as.numeric(pvalue(fit)["tp"]),
          n = nobs(fit),
          n_states = length(unique(dat$state)),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Windowed
    for (w in windows) {
      
      dat <- dat0 %>%
        filter(!(state %in% drop_states)) %>%
        filter(year %in% (event_year - w):(event_year + w)) %>%
        filter(!is.na(.data[[dv]]))
      
      fit <- tryCatch(
        feols(as.formula(paste0(dv, " ~ tp | state + year")),
              data = dat, cluster = ~state),
        error = function(e) NULL
      )
      
      if (is.null(fit)) next
      if (!("tp" %in% names(coef(fit)))) next
      
      rows[[length(rows) + 1]] <- data.frame(
        tag = unique(dat0$tag)[1],
        seed = unique(dat0$seed)[1],
        dep_var = dv,
        window = paste0("+-", w),
        estimate = as.numeric(coef(fit)["tp"]),
        std_error = as.numeric(se(fit)["tp"]),
        p_value = as.numeric(pvalue(fit)["tp"]),
        n = nobs(fit),
        n_states = length(unique(dat$state)),
        stringsAsFactors = FALSE
      )
    }
  }
}

# ----------------------------------------------------------
# Combine + add CIs + save
# ----------------------------------------------------------
tbl <- bind_rows(rows) %>%
  mutate(
    lb = estimate - 1.96 * std_error,
    ub = estimate + 1.96 * std_error
  ) %>%
  arrange(dep_var, tag, factor(window, levels = c("+-2","+-3","+-4","+-5","all")))


write_csv(tbl, file.path(out_dir, "haiti_random_construction_twfe_results_long.csv"))
saveRDS(tbl, file.path(out_dir, "haiti_random_construction_twfe_results_long.rds"))
tbl <- tbl %>% 
  filter(dep_var %in% c(
    "fragmentation",
    "weighted_fragmentation",
    "usaid_frag_delta_drop_unfunded",
    "usaid_n_partners"
  )) %>% 
  mutate(dep_var = factor(
    dep_var, 
    levels = c(
      "fragmentation",
      "weighted_fragmentation",
      "usaid_frag_delta_drop_unfunded",
      "usaid_n_partners"
    ),
    labels = c(
      "Fragmentation",
      "Weighted_fragmentation",
      "Drop unfunded",
      "USAID partners"
    )
  )) 
net_permute_est <- ggplot(tbl, aes(x = estimate, fill = window)) + 
  geom_density(alpha = 0.35) + 
  facet_wrap(~dep_var, scales = "free") + 
  theme_minimal() + 
  labs(fill = "Time Window", title = "Permutation coefficients", y = "Density", x = "Coefficient") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "bottom") 


p_values <- ggplot(tbl, aes(x = p_value, y = window, fill = window)) + 
  geom_boxplot() + 
  facet_wrap(~dep_var, scales = "free") + 
  theme_minimal() + 
  geom_hline(yintercept = 0.05) + 
  labs(fill = "Time Window", title = "Window", x = "p-values") + 
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "none") 


ggsave(
  filename = file.path(out_dir, "plots/p_value_network_measure_noise.png"),
  plot = p_values, width = 10, height = 5, dpi = 200
)

ggsave(
  filename = file.path(out_dir, "plots/estimate_network_measure_noise.png"),
  plot = net_permute_est, width = 10, height = 5, dpi = 200
)