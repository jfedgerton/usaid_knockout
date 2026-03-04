# ==========================================================
# 18_haiti_earthquake.R
# Haiti DiD: Network-level + Node-level analysis
#
# DESIGN:
#   Network-level: TWFE DiD (state + year FE, cluster by state)
#     DVs: fragmentation, weighted_fragmentation, usaid_n_partners
#   Node-level: TWFE DiD (org + year FE, cluster by org)
#     DVs: in_degree, sole_usaid_funded, n_funders,
#           n_non_usaid_funders, total_degree, rewire_in, rewire_out
#     Spillover: pre-treatment USAID exposure share (Aronow & Samii 2017)
#     Balanced panel: restrict to orgs present both pre and post
#     Entry/exit modeled as separate outcome
# ==========================================================

library(dplyr)
library(fixest)
library(readr)
library(ggplot2)

if (!grepl("Ecosystem_of_Aid", getwd())) setwd("Ecosystem_of_Aid")

out_dir <- "Data/results/haiti_main_effects"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "node_level"), showWarnings = FALSE, recursive = TRUE)

treat_state <- "hti"
event_year  <- 2010
windows     <- c(2, 3, 4, 5)


# ##############################################################
# PART 1: NETWORK-LEVEL DiD
# ##############################################################

net <- readRDS("Data/clean_nets/network_attributes_directed.rds")
net$state   <- tolower(net$state)
net$treated <- as.integer(net$state == treat_state)
net$post    <- as.integer(net$year >= event_year)
net$tp      <- net$treated * net$post

net_dvs <- c("fragmentation", "weighted_fragmentation", "usaid_n_partners")
net_dvs <- net_dvs[net_dvs %in% names(net)]

cat("Network-level DVs:", paste(net_dvs, collapse = ", "), "\n")
cat("States:", length(unique(net$state)), "| Years:", length(unique(net$year)), "\n")

# --- 1a: Full-sample estimates ---
net_rows <- list()

for (dv in net_dvs) {
  dat <- net[!is.na(net[[dv]]), ]
  fit <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp | state + year")),
          data = dat, cluster = ~state),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  if (!("tp" %in% names(coef(fit)))) next

  net_rows[[length(net_rows) + 1]] <- data.frame(
    level     = "network",
    dep_var   = dv,
    window    = "all",
    estimate  = as.numeric(coef(fit)["tp"]),
    std_error = as.numeric(se(fit)["tp"]),
    p_value   = as.numeric(pvalue(fit)["tp"]),
    n         = nobs(fit),
    n_states  = length(unique(dat$state)),
    stringsAsFactors = FALSE
  )
}

# --- 1b: Windowed estimates ---
for (dv in net_dvs) {
  for (w in windows) {
    yr_range <- (event_year - w):(event_year + w)
    dat <- net[net$year %in% yr_range & !is.na(net[[dv]]), ]

    fit <- tryCatch(
      feols(as.formula(paste0(dv, " ~ tp | state + year")),
            data = dat, cluster = ~state),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    if (!("tp" %in% names(coef(fit)))) next

    net_rows[[length(net_rows) + 1]] <- data.frame(
      level     = "network",
      dep_var   = dv,
      window    = paste0("+-", w),
      estimate  = as.numeric(coef(fit)["tp"]),
      std_error = as.numeric(se(fit)["tp"]),
      p_value   = as.numeric(pvalue(fit)["tp"]),
      n         = nobs(fit),
      n_states  = length(unique(dat$state)),
      stringsAsFactors = FALSE
    )
  }
}

net_tbl <- bind_rows(net_rows)
net_tbl$lb <- net_tbl$estimate - 1.96 * net_tbl$std_error
net_tbl$ub <- net_tbl$estimate + 1.96 * net_tbl$std_error
net_tbl <- net_tbl %>%
  arrange(dep_var, factor(window, levels = c("+-2","+-3","+-4","+-5","all")))

write_csv(net_tbl, file.path(out_dir, "haiti_network_level_results.csv"))
cat("\n=== Network-level results ===\n")
print(net_tbl)

# --- 1c: Coefficient plots ---
for (dv in unique(net_tbl$dep_var)) {
  dfp <- net_tbl[net_tbl$dep_var == dv, ]
  p <- ggplot(dfp, aes(x = window, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2) +
    theme_classic() +
    labs(title = paste0("Haiti network-level DiD: ", dv),
         x = "Window around 2010", y = "DiD estimate (tp)")
  ggsave(file.path(out_dir, "figures", paste0("network_", dv, ".png")),
         p, width = 7, height = 4)
}


# ##############################################################
# PART 2: NODE-LEVEL DiD
# ##############################################################

nodes <- readRDS("Data/clean_nets/node_attributes_directed.rds")
nodes$state <- tolower(nodes$state)

cat("\nNode panel dimensions:", nrow(nodes), "rows,",
    length(unique(nodes$node_name)), "unique orgs,",
    length(unique(nodes$state)), "states\n")

# ----------------------------------------------------------
# Step 2.1: Construct balanced panel
#   Restrict to orgs observed in at least one pre AND one post year.
#   This avoids composition bias from entry/exit.
# ----------------------------------------------------------
pre_orgs  <- unique(nodes[nodes$year < event_year, c("state", "node_name")])
post_orgs <- unique(nodes[nodes$year >= event_year, c("state", "node_name")])
balanced_orgs <- merge(pre_orgs, post_orgs, by = c("state", "node_name"))

nodes_bal <- merge(nodes, balanced_orgs, by = c("state", "node_name"))

cat("Balanced panel:", nrow(balanced_orgs), "unique state-org pairs\n")
cat("Balanced panel rows:", nrow(nodes_bal), "\n")

# ----------------------------------------------------------
# Step 2.2: Treatment indicators
# ----------------------------------------------------------
nodes_bal$treated <- as.integer(nodes_bal$state == treat_state)
nodes_bal$post    <- as.integer(nodes_bal$year >= event_year)
nodes_bal$tp      <- nodes_bal$treated * nodes_bal$post
nodes_bal$org_id  <- paste0(nodes_bal$state, "||", nodes_bal$node_name)

# ----------------------------------------------------------
# Step 2.3: Pre-treatment exposure variables (time-invariant)
#   Computed from pre-period averages/maxima to avoid endogeneity.
#   These serve as heterogeneous treatment / spillover moderators.
# ----------------------------------------------------------
nodes_pre <- nodes_bal[nodes_bal$year < event_year, ]

# For each org: was it ever directly funded by USAID pre-treatment?
pre_usaid_funded <- aggregate(
  usaid_funded ~ state + node_name,
  data = nodes_pre, FUN = max, na.rm = TRUE
)
names(pre_usaid_funded)[3] <- "pre_usaid_funded"

# Was it ever sole-USAID-funded pre-treatment?
pre_sole <- aggregate(
  sole_usaid_funded ~ state + node_name,
  data = nodes_pre, FUN = max, na.rm = TRUE
)
names(pre_sole)[3] <- "pre_sole_usaid_funded"

# Was it ever tied to USAID (any neighbor)?
pre_tied <- aggregate(
  tied_to_usaid ~ state + node_name,
  data = nodes_pre, FUN = max, na.rm = TRUE
)
names(pre_tied)[3] <- "pre_tied_to_usaid"

# Mean pre-treatment USAID exposure share (spillover measure)
pre_exposure <- aggregate(
  usaid_exposure_share ~ state + node_name,
  data = nodes_pre, FUN = mean, na.rm = TRUE
)
names(pre_exposure)[3] <- "pre_usaid_exposure_share"

# Mean pre-treatment in-degree (controls for baseline connectivity)
pre_indeg <- aggregate(
  in_degree ~ state + node_name,
  data = nodes_pre, FUN = mean, na.rm = TRUE
)
names(pre_indeg)[3] <- "pre_in_degree"

# Merge all pre-treatment variables
nodes_bal <- merge(nodes_bal, pre_usaid_funded, by = c("state", "node_name"), all.x = TRUE)
nodes_bal <- merge(nodes_bal, pre_sole,         by = c("state", "node_name"), all.x = TRUE)
nodes_bal <- merge(nodes_bal, pre_tied,         by = c("state", "node_name"), all.x = TRUE)
nodes_bal <- merge(nodes_bal, pre_exposure,     by = c("state", "node_name"), all.x = TRUE)
nodes_bal <- merge(nodes_bal, pre_indeg,        by = c("state", "node_name"), all.x = TRUE)

cat("Pre-treatment exposure vars merged.\n")
cat("  USAID-funded orgs (pre):", sum(nodes_bal$pre_usaid_funded == 1, na.rm = TRUE) / nrow(nodes_bal), "\n")

# ----------------------------------------------------------
# Step 2.4: Define node-level DVs
#   These directly reflect the cascade simulation mechanism:
#   - in_degree / n_funders: funder redundancy per org
#   - sole_usaid_funded: binary dependency indicator
#   - n_non_usaid_funders: alternative funder count
#   - total_degree: overall connectivity (peripherality)
#   - rewire_in / rewire_out: partnership turnover
# ----------------------------------------------------------
node_dvs <- c("in_degree", "n_funders", "sole_usaid_funded",
              "n_non_usaid_funders", "total_degree",
              "rewire_in", "rewire_out")
node_dvs <- node_dvs[node_dvs %in% names(nodes_bal)]

cat("Node-level DVs:", paste(node_dvs, collapse = ", "), "\n")

# ----------------------------------------------------------
# Step 2.5A: Pooled node-level DiD (all countries)
#   Spec: DV ~ tp | org_id + year, cluster = org_id
#   This is the cross-country DiD: Haiti orgs vs orgs in other
#   fragile states, before vs after 2010.
# ----------------------------------------------------------
cat("\n=== Pooled node-level DiD (all countries) ===\n")
node_rows <- list()

for (dv in node_dvs) {
  dat <- nodes_bal[!is.na(nodes_bal[[dv]]), ]
  n_orgs_dv <- length(unique(dat$org_id))

  # --- Baseline ATE ---
  fit_base <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
          data = dat, cluster = ~org_id),
    error = function(e) NULL
  )
  if (!is.null(fit_base) && "tp" %in% names(coef(fit_base))) {
    node_rows[[length(node_rows) + 1]] <- data.frame(
      spec      = "baseline",
      dep_var   = dv,
      term      = "tp",
      estimate  = as.numeric(coef(fit_base)["tp"]),
      std_error = as.numeric(se(fit_base)["tp"]),
      p_value   = as.numeric(pvalue(fit_base)["tp"]),
      n         = nobs(fit_base),
      n_orgs    = n_orgs_dv,
      stringsAsFactors = FALSE
    )
    cat("  ", dv, "baseline: tp =", round(coef(fit_base)["tp"], 4),
        " p =", round(pvalue(fit_base)["tp"], 4), "\n")
  }

  # --- Heterogeneous effect by direct USAID funding (pre-treatment) ---
  fit_het <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp + tp:pre_usaid_funded | org_id + year")),
          data = dat, cluster = ~org_id),
    error = function(e) NULL
  )
  if (!is.null(fit_het)) {
    for (cf in names(coef(fit_het))) {
      node_rows[[length(node_rows) + 1]] <- data.frame(
        spec      = "heterogeneous_direct",
        dep_var   = dv,
        term      = cf,
        estimate  = as.numeric(coef(fit_het)[cf]),
        std_error = as.numeric(se(fit_het)[cf]),
        p_value   = as.numeric(pvalue(fit_het)[cf]),
        n         = nobs(fit_het),
        n_orgs    = n_orgs_dv,
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Spillover: exposure share (Aronow & Samii) ---
  fit_sp <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp + tp:pre_usaid_exposure_share | org_id + year")),
          data = dat, cluster = ~org_id),
    error = function(e) NULL
  )
  if (!is.null(fit_sp)) {
    for (cf in names(coef(fit_sp))) {
      node_rows[[length(node_rows) + 1]] <- data.frame(
        spec      = "spillover_exposure",
        dep_var   = dv,
        term      = cf,
        estimate  = as.numeric(coef(fit_sp)[cf]),
        std_error = as.numeric(se(fit_sp)[cf]),
        p_value   = as.numeric(pvalue(fit_sp)[cf]),
        n         = nobs(fit_sp),
        n_orgs    = n_orgs_dv,
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Full model: direct + spillover ---
  fit_full <- tryCatch(
    feols(as.formula(paste0(
      dv, " ~ tp + tp:pre_usaid_funded + tp:pre_usaid_exposure_share | org_id + year")),
      data = dat, cluster = ~org_id),
    error = function(e) NULL
  )
  if (!is.null(fit_full)) {
    for (cf in names(coef(fit_full))) {
      node_rows[[length(node_rows) + 1]] <- data.frame(
        spec      = "full_direct_spillover",
        dep_var   = dv,
        term      = cf,
        estimate  = as.numeric(coef(fit_full)[cf]),
        std_error = as.numeric(se(fit_full)[cf]),
        p_value   = as.numeric(pvalue(fit_full)[cf]),
        n         = nobs(fit_full),
        n_orgs    = n_orgs_dv,
        stringsAsFactors = FALSE
      )
    }
  }
}

node_tbl <- bind_rows(node_rows)
write_csv(node_tbl, file.path(out_dir, "node_level", "node_did_pooled.csv"))

cat("\n=== Pooled node-level results written ===\n")


# ----------------------------------------------------------
# Step 2.5B: Haiti-only node-level DiD
#   Within Haiti, treatment = pre-treatment USAID funding status.
#   Question: did USAID-funded orgs change differently from
#   non-USAID-funded orgs in Haiti after the earthquake?
# ----------------------------------------------------------
cat("\n=== Haiti-only node-level DiD ===\n")
nodes_hti <- nodes_bal[nodes_bal$state == treat_state, ]

cat("Haiti balanced panel:", length(unique(nodes_hti$org_id)), "orgs,",
    nrow(nodes_hti), "rows\n")

haiti_rows <- list()

for (dv in node_dvs) {
  dat <- nodes_hti[!is.na(nodes_hti[[dv]]), ]
  if (nrow(dat) < 10) next
  n_orgs_dv <- length(unique(dat$org_id))

  # Direct effect: pre_usaid_funded x post
  dat$usaid_x_post <- dat$pre_usaid_funded * dat$post
  fit_direct <- tryCatch(
    feols(as.formula(paste0(dv, " ~ usaid_x_post | org_id + year")),
          data = dat, cluster = ~org_id),
    error = function(e) NULL
  )
  if (!is.null(fit_direct) && "usaid_x_post" %in% names(coef(fit_direct))) {
    haiti_rows[[length(haiti_rows) + 1]] <- data.frame(
      spec      = "haiti_direct",
      dep_var   = dv,
      term      = "usaid_x_post",
      estimate  = as.numeric(coef(fit_direct)["usaid_x_post"]),
      std_error = as.numeric(se(fit_direct)["usaid_x_post"]),
      p_value   = as.numeric(pvalue(fit_direct)["usaid_x_post"]),
      n         = nobs(fit_direct),
      n_orgs    = n_orgs_dv,
      stringsAsFactors = FALSE
    )
    cat("  ", dv, "haiti_direct:", round(coef(fit_direct)["usaid_x_post"], 4),
        " p =", round(pvalue(fit_direct)["usaid_x_post"], 4), "\n")
  }

  # Spillover effect: pre_usaid_exposure_share x post
  dat$exposure_x_post <- dat$pre_usaid_exposure_share * dat$post
  fit_spill <- tryCatch(
    feols(as.formula(paste0(dv, " ~ exposure_x_post | org_id + year")),
          data = dat, cluster = ~org_id),
    error = function(e) NULL
  )
  if (!is.null(fit_spill) && "exposure_x_post" %in% names(coef(fit_spill))) {
    haiti_rows[[length(haiti_rows) + 1]] <- data.frame(
      spec      = "haiti_spillover",
      dep_var   = dv,
      term      = "exposure_x_post",
      estimate  = as.numeric(coef(fit_spill)["exposure_x_post"]),
      std_error = as.numeric(se(fit_spill)["exposure_x_post"]),
      p_value   = as.numeric(pvalue(fit_spill)["exposure_x_post"]),
      n         = nobs(fit_spill),
      n_orgs    = n_orgs_dv,
      stringsAsFactors = FALSE
    )
  }

  # Joint: direct + spillover
  fit_joint <- tryCatch(
    feols(as.formula(paste0(dv, " ~ usaid_x_post + exposure_x_post | org_id + year")),
          data = dat, cluster = ~org_id),
    error = function(e) NULL
  )
  if (!is.null(fit_joint)) {
    for (cf in names(coef(fit_joint))) {
      haiti_rows[[length(haiti_rows) + 1]] <- data.frame(
        spec      = "haiti_joint",
        dep_var   = dv,
        term      = cf,
        estimate  = as.numeric(coef(fit_joint)[cf]),
        std_error = as.numeric(se(fit_joint)[cf]),
        p_value   = as.numeric(pvalue(fit_joint)[cf]),
        n         = nobs(fit_joint),
        n_orgs    = n_orgs_dv,
        stringsAsFactors = FALSE
      )
    }
  }
}

haiti_node_tbl <- bind_rows(haiti_rows)
write_csv(haiti_node_tbl, file.path(out_dir, "node_level", "node_did_haiti_only.csv"))


# ----------------------------------------------------------
# Step 2.5C: Windowed node-level DiD (pooled, baseline spec)
#   Matches the network-level window analysis.
# ----------------------------------------------------------
cat("\n=== Windowed node-level DiD ===\n")
node_win_rows <- list()

for (dv in node_dvs) {
  for (w in windows) {
    yr_range <- (event_year - w):(event_year + w)
    dat <- nodes_bal[nodes_bal$year %in% yr_range & !is.na(nodes_bal[[dv]]), ]
    n_orgs_dv <- length(unique(dat$org_id))

    fit <- tryCatch(
      feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
            data = dat, cluster = ~org_id),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    if (!("tp" %in% names(coef(fit)))) next

    node_win_rows[[length(node_win_rows) + 1]] <- data.frame(
      dep_var   = dv,
      window    = paste0("+-", w),
      estimate  = as.numeric(coef(fit)["tp"]),
      std_error = as.numeric(se(fit)["tp"]),
      p_value   = as.numeric(pvalue(fit)["tp"]),
      n         = nobs(fit),
      n_orgs    = n_orgs_dv,
      stringsAsFactors = FALSE
    )
  }
}

node_win_tbl <- bind_rows(node_win_rows)
node_win_tbl$lb <- node_win_tbl$estimate - 1.96 * node_win_tbl$std_error
node_win_tbl$ub <- node_win_tbl$estimate + 1.96 * node_win_tbl$std_error
write_csv(node_win_tbl, file.path(out_dir, "node_level", "node_did_windowed.csv"))


# ##############################################################
# PART 3: ENTRY/EXIT AS SEPARATE OUTCOME
# ##############################################################
# Entry/exit is NOT a node-level DV because orgs that enter post-
# treatment are not in the balanced panel. Instead, model the
# count of orgs per state-year as a network-level outcome.

cat("\n=== Entry/exit analysis ===\n")

# Use the FULL (unbalanced) node panel for this
entry_exit <- nodes %>%
  group_by(state, year) %>%
  summarise(
    n_orgs          = n_distinct(node_name),
    n_usaid_funded  = sum(usaid_funded, na.rm = TRUE),
    n_new_orgs      = NA_integer_,     # placeholder; computed below
    n_exiting_orgs  = NA_integer_,
    .groups = "drop"
  )

# Compute new entrants and exits by state
all_states <- unique(nodes$state)
entry_rows <- list()

for (st in all_states) {
  st_data <- nodes[nodes$state == st, ]
  st_years <- sort(unique(st_data$year))

  for (j in seq_along(st_years)) {
    yr_now  <- st_years[j]
    orgs_now <- unique(st_data$node_name[st_data$year == yr_now])

    if (j == 1) {
      n_new <- NA_integer_
      n_exit <- NA_integer_
    } else {
      yr_prev  <- st_years[j - 1]
      orgs_prev <- unique(st_data$node_name[st_data$year == yr_prev])
      n_new  <- length(setdiff(orgs_now, orgs_prev))
      n_exit <- length(setdiff(orgs_prev, orgs_now))
    }

    entry_rows[[length(entry_rows) + 1]] <- data.frame(
      state = st, year = yr_now,
      n_new_orgs = n_new, n_exiting_orgs = n_exit,
      stringsAsFactors = FALSE
    )
  }
}

entry_exit_detail <- bind_rows(entry_rows)
entry_exit <- merge(
  entry_exit[, c("state", "year", "n_orgs", "n_usaid_funded")],
  entry_exit_detail,
  by = c("state", "year"), all.x = TRUE
)

entry_exit$treated <- as.integer(entry_exit$state == treat_state)
entry_exit$post    <- as.integer(entry_exit$year >= event_year)
entry_exit$tp      <- entry_exit$treated * entry_exit$post

# DiD on entry/exit counts
entry_dvs <- c("n_orgs", "n_usaid_funded", "n_new_orgs", "n_exiting_orgs")
entry_results <- list()

for (dv in entry_dvs) {
  dat <- entry_exit[!is.na(entry_exit[[dv]]), ]
  fit <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp | state + year")),
          data = dat, cluster = ~state),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  if (!("tp" %in% names(coef(fit)))) next

  entry_results[[length(entry_results) + 1]] <- data.frame(
    dep_var   = dv,
    estimate  = as.numeric(coef(fit)["tp"]),
    std_error = as.numeric(se(fit)["tp"]),
    p_value   = as.numeric(pvalue(fit)["tp"]),
    n         = nobs(fit),
    stringsAsFactors = FALSE
  )
  cat("  ", dv, ": tp =", round(coef(fit)["tp"], 3),
      " p =", round(pvalue(fit)["tp"], 4), "\n")
}

entry_tbl <- bind_rows(entry_results)
write_csv(entry_tbl, file.path(out_dir, "node_level", "entry_exit_did.csv"))
write_csv(entry_exit, file.path(out_dir, "node_level", "entry_exit_by_state_year.csv"))


# ##############################################################
# PART 4: DIAGNOSTICS
# ##############################################################

cat("\n=== Diagnostics ===\n")

# Balanced panel summary by state
bal_summary <- nodes_bal %>%
  group_by(state) %>%
  summarise(
    n_orgs        = n_distinct(org_id),
    n_years       = n_distinct(year),
    n_rows        = n(),
    pct_usaid_funded_pre  = mean(usaid_funded[year < event_year], na.rm = TRUE),
    pct_sole_usaid_pre    = mean(sole_usaid_funded[year < event_year], na.rm = TRUE),
    mean_exposure_pre     = mean(usaid_exposure_share[year < event_year], na.rm = TRUE),
    .groups = "drop"
  )
write_csv(bal_summary, file.path(out_dir, "node_level", "balanced_panel_summary.csv"))
cat("Balanced panel summary:\n")
print(bal_summary)

# Pre-treatment balance check: are USAID-funded orgs different from
# non-USAID-funded orgs on observables before treatment?
nodes_pre_hti <- nodes_bal[nodes_bal$state == treat_state & nodes_bal$year < event_year, ]
if (nrow(nodes_pre_hti) > 10) {
  balance_vars <- c("in_degree", "total_degree", "betweenness", "pagerank")
  balance_vars <- balance_vars[balance_vars %in% names(nodes_pre_hti)]

  bal_rows <- list()
  for (bv in balance_vars) {
    funded_mean     <- mean(nodes_pre_hti[[bv]][nodes_pre_hti$usaid_funded == 1], na.rm = TRUE)
    not_funded_mean <- mean(nodes_pre_hti[[bv]][nodes_pre_hti$usaid_funded == 0], na.rm = TRUE)
    ttest <- tryCatch(
      t.test(nodes_pre_hti[[bv]] ~ nodes_pre_hti$usaid_funded),
      error = function(e) NULL
    )
    bal_rows[[length(bal_rows) + 1]] <- data.frame(
      variable       = bv,
      mean_funded    = funded_mean,
      mean_notfunded = not_funded_mean,
      diff           = funded_mean - not_funded_mean,
      p_value        = if (!is.null(ttest)) ttest$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  balance_tbl <- bind_rows(bal_rows)
  write_csv(balance_tbl, file.path(out_dir, "node_level", "pretreatment_balance_haiti.csv"))
  cat("\nPre-treatment balance (Haiti, funded vs not):\n")
  print(balance_tbl)
}


cat("\n=== Haiti DiD complete (network + node level) ===\n")
