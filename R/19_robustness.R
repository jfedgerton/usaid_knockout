# ==========================================================
# HAITI ROBUSTNESS SUITE (NETWORK TOPOLOGY DiD / EVENT STUDY)
# Focus: Haiti earthquake inflow case (HTI treated, 2010)
# Data: Data/clean_nets/network_attributes_directed.rds
#
# Outputs -> Data/results/robustness/haiti/
#
# Checks included:
#   A) Baseline pooled DiD (main donor pool + full donor pool)
#   B) Event study + joint test of leads
#   C) Placebo years (windowed)
#   D) Placebo treated units (treat each control as "treated")
#   E) Leave-one-out controls (influence)
#   F) Random-drop 1–5 controls (stress test) + plots
#   G) Window sensitivity grid (vary pre/post windows)
#   H) Wild cluster bootstrap p-values (if available)
#   I) Permutation/randomization p-values (via placebo-treated)
#   J) Synthetic control / synthDiD (if available) for headline DVs
# ==========================================================

# ----------------------------
# Packages
# ----------------------------
pkgs <- c("dplyr", "fixest", "readr", "ggplot2", "tidyr", "stringr", "fwildclusterboot", "synthdid")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(dplyr)
library(fixest)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)
library(fwildclusterboot)
library(synthdid)

# ----------------------------
# Paths
# ----------------------------
if (!grepl("Ecosystem_of_Aid", getwd())) setwd("Ecosystem_of_Aid")

out_dir <- "Data/results/robustness/haiti"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# User settings
# ----------------------------
treated_state <- "hti"
event_year <- 2010

# "Main" donor pool exclusion (your primary spec)
exclude_main <- c("som", "ssud")

# Sample years (you said you have 2005–2020)
min_year <- 2005
max_year <- 2020

# DVs to analyze (will auto-drop any missing columns)
dv_list <- c(
  "fragmentation",
  "weighted_fragmentation",
  "modularity",
  "density",
  "network_efficiency",
  "mean_distance",
  "donor_concentration",
  "indegree_centralization",
  "outdegree_centralization",
  "betweenness_centralization",
  "usaid_share_partners_no_alt_funder",
  "usaid_frag_delta_direct",
  "usaid_frag_delta_drop_unfunded",
  "usaid_n_partners",
  "usaid_n_nodes"
)

# For computationally heavy checks, restrict to headline outcomes:
dv_focus <- c("fragmentation", "weighted_fragmentation", "modularity", "density")

# Random-drop stress test settings
set.seed(1234)
random_drop_draws <- 500      # increase to 1000+ if you want
random_drop_k_vals <- 1:5     # drop 1–5 controls per draw

# Placebo-year test window (keeps design comparable given 2005–2020)
placebo_pre_w <- 3
placebo_post_w <- 3

# Window sensitivity grid (pre, post)
pre_grid <- c(2, 3, 4, 5)
post_grid <- c(2, 3, 5, 8, 10)

# Synthetic DiD outcomes (keep small)
synth_dvs <- c("fragmentation", "weighted_fragmentation", "modularity")

# ----------------------------
# Load data
# ----------------------------
net <- readRDS("Data/clean_nets/network_attributes_directed.rds")

# Basic sanity checks
stopifnot(all(c("state", "year") %in% names(net)))

net <- net %>%
  filter(year >= min_year, year <= max_year) %>%
  mutate(
    treated = as.integer(state == treated_state),
    post = as.integer(year >= event_year),
    tp = treated * post,
    rel_year = year - event_year
  )

# Keep only DVs that exist
dv_list <- dv_list[dv_list %in% names(net)]
dv_focus <- dv_focus[dv_focus %in% names(net)]
synth_dvs <- synth_dvs[synth_dvs %in% names(net)]

if (length(dv_list) == 0) stop("No DVs found in net. Check variable names.")

# Helper: summarize a single pooled DiD for one DV + one dataset
run_pooled_did_one <- function(dat, dv, tp_var = "tp") {
  dat2 <- dat %>% filter(!is.na(.data[[dv]]))
  fit <- tryCatch(
    feols(as.formula(paste0(dv, " ~ ", tp_var, " | state + year")),
          data = dat2, cluster = ~ state),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  b <- as.numeric(coef(fit)[tp_var])
  se <- as.numeric(se(fit)[tp_var])
  p <- as.numeric(pvalue(fit)[tp_var])
  
  data.frame(
    dep_var = dv,
    estimate = b,
    std_error = se,
    p_value = p,
    lb = b - 1.96 * se,
    ub = b + 1.96 * se,
    n = nrow(dat2),
    n_states = length(unique(dat2$state)),
    stringsAsFactors = FALSE
  )
}

# NOTE: You said you prefer no functions at top.
# The function above is small and only used to avoid copying the same 15 lines 50 times.
# If you want it inlined, tell me and I’ll rewrite it fully inline.

# ==========================================================
# A) Baseline pooled DiD (main donor pool + full donor pool)
# ==========================================================
net_main <- net %>% filter(!(state %in% exclude_main))
net_full <- net

baseline_main <- bind_rows(lapply(dv_list, function(dv) run_pooled_did_one(net_main, dv))) %>%
  mutate(spec = "baseline pooled DiD (exclude som+ssud)") %>%
  arrange(p_value)

baseline_full <- bind_rows(lapply(dv_list, function(dv) run_pooled_did_one(net_full, dv))) %>%
  mutate(spec = "baseline pooled DiD (ALL states)") %>%
  arrange(p_value)

write_csv(baseline_main, file.path(out_dir, "baseline_pooled_did_main.csv"))
write_csv(baseline_full, file.path(out_dir, "baseline_pooled_did_allstates.csv"))

# ==========================================================
# B) Event study + joint test of leads (pre-trends)
# ==========================================================
event_rows <- list()
lead_tests <- list()
for (dv in dv_list) {
  dat <- net_main %>% filter(!is.na(.data[[dv]]))
  
  fit_es <- tryCatch(
    feols(as.formula(paste0(dv, " ~ i(rel_year, treated, ref = -1) | state + year")),
          data = dat, cluster = ~ state),
    error = function(e) NULL
  )
  if (is.null(fit_es)) next
  
  ct <- as.data.frame(coeftable(fit_es))
  ct$term <- rownames(ct)
  rownames(ct) <- NULL
  
  ct <- ct %>% filter(str_detect(term, "^rel_year::"))
  
  rel_raw <- ct$term %>%
    str_replace("^rel_year::", "") %>%
    str_replace("(:treated|#treated)$", "")
  
  ct$rel_year <- suppressWarnings(as.integer(rel_raw))
  ct$dep_var <- dv
  
  event_rows[[length(event_rows) + 1]] <- ct
  
  # Joint test: leads (rel_year <= -2) = 0
  lead_terms <- ct %>% filter(!is.na(rel_year), rel_year <= -2) %>% pull(term)
  
  # keep only terms that exist in model coefficients
  lead_terms <- lead_terms[lead_terms %in% names(coef(fit_es))]
  
  lead_p <- NA_real_
  lead_stat <- NA_real_
  lead_df <- length(lead_terms)
  
  if (lead_df > 0) {
    V <- tryCatch(vcov(fit_es)[lead_terms, lead_terms, drop = FALSE], error = function(e) NULL)
    
    if (!is.null(V)) {
      vars <- diag(V)
      keep <- which(is.finite(vars) & vars > 1e-12)
      lead_terms2 <- lead_terms[keep]
      
      lead_df <- length(lead_terms2)
      
      if (lead_df > 0) {
        b2 <- coef(fit_es)[lead_terms2]
        V2 <- vcov(fit_es)[lead_terms2, lead_terms2, drop = FALSE]
        
        lead_stat <- as.numeric(t(b2) %*% solve(V2) %*% b2)
        lead_p <- 1 - pchisq(lead_stat, df = lead_df)
      }
    }
  }
  
  lead_tests[[length(lead_tests) + 1]] <- data.frame(
    dep_var = dv,
    n_leads = lead_df,
    wald_stat = lead_stat,
    p_value = lead_p,
    stringsAsFactors = FALSE
  )
}


eventstudy_coefs <- bind_rows(event_rows)
eventstudy_leadtests <- bind_rows(lead_tests)

write_csv(eventstudy_coefs, file.path(out_dir, "eventstudy_coefs.csv"))
write_csv(eventstudy_leadtests, file.path(out_dir, "eventstudy_leads_joint_test.csv"))

# ==========================================================
# C) Placebo years (windowed)
# ==========================================================
placebo_years <- setdiff(
  seq(min_year + placebo_pre_w, max_year - placebo_post_w),
  event_year
)

placebo_rows <- list()

for (dv in dv_focus) {
  # True estimate in the same window size around event_year
  true_window <- net_main %>%
    filter(year >= (event_year - placebo_pre_w), year <= (event_year + placebo_post_w))
  
  true_fit <- run_pooled_did_one(true_window, dv)
  if (is.null(true_fit)) next
  true_b <- true_fit$estimate
  
  for (y0 in placebo_years) {
    dat <- net_main %>%
      filter(year >= (y0 - placebo_pre_w), year <= (y0 + placebo_post_w)) %>%
      mutate(
        post_pl = as.integer(year >= y0),
        tp_pl = treated * post_pl
      )
    
    res <- run_pooled_did_one(dat, dv, tp_var = "tp_pl")
    if (is.null(res)) next
    
    placebo_rows[[length(placebo_rows) + 1]] <- data.frame(
      dep_var = dv,
      placebo_year = y0,
      estimate = res$estimate,
      std_error = res$std_error,
      p_value = res$p_value,
      true_estimate = true_b,
      stringsAsFactors = FALSE
    )
  }
}

placebo_years_df <- bind_rows(placebo_rows)
write_csv(placebo_years_df, file.path(out_dir, "placebo_years_windowed.csv"))

# Simple plots
if (nrow(placebo_years_df) > 0) {
  for (dv in unique(placebo_years_df$dep_var)) {
    pdat <- placebo_years_df %>% filter(dep_var == dv)
    
    g <- ggplot(pdat, aes(x = placebo_year, y = estimate)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point() +
      geom_vline(xintercept = event_year, linetype = "dotted") +
      labs(
        title = paste0("Placebo-year DiD (windowed) — ", dv),
        subtitle = paste0("Window: [year-", placebo_pre_w, ", year+", placebo_post_w, "], dotted line = true 2010"),
        x = "Placebo year",
        y = "Estimated Treated×Post"
      )
    
    ggsave(
      filename = file.path(out_dir, "plots", paste0("placebo_years_", dv, ".png")),
      plot = g, width = 10, height = 5, dpi = 200
    )
  }
}

# ==========================================================
# D) Placebo treated units (treat each control as treated)
# ==========================================================
states_all <- sort(unique(net_main$state))
controls <- setdiff(states_all, treated_state)

placebo_units_rows <- list()

for (dv in dv_focus) {
  # True effect
  true_res <- run_pooled_did_one(net_main, dv)
  if (is.null(true_res)) next
  true_b <- true_res$estimate
  
  for (s in controls) {
    dat <- net_main %>%
      mutate(
        treated_pl = as.integer(state == s),
        tp_pl = treated_pl * post
      )
    
    res <- run_pooled_did_one(dat, dv, tp_var = "tp_pl")
    if (is.null(res)) next
    
    placebo_units_rows[[length(placebo_units_rows) + 1]] <- data.frame(
      dep_var = dv,
      placebo_treated = s,
      estimate = res$estimate,
      std_error = res$std_error,
      p_value = res$p_value,
      true_estimate = true_b,
      stringsAsFactors = FALSE
    )
  }
}

placebo_treated_units_df <- bind_rows(placebo_units_rows)
write_csv(placebo_treated_units_df, file.path(out_dir, "placebo_treated_units.csv"))

# Randomization p-value from placebo-treated distribution
placebo_treated_summary <- placebo_treated_units_df %>%
  group_by(dep_var) %>%
  summarise(
    true_estimate = first(true_estimate),
    placebo_mean = mean(estimate, na.rm = TRUE),
    placebo_sd = sd(estimate, na.rm = TRUE),
    p_perm_two_sided = mean(abs(estimate) >= abs(first(true_estimate)), na.rm = TRUE),
    n_placebos = sum(!is.na(estimate)),
    .groups = "drop"
  )

write_csv(placebo_treated_summary, file.path(out_dir, "placebo_treated_units_summary.csv"))

# Plot placebo-treated distribution
if (nrow(placebo_treated_units_df) > 0) {
  for (dv in unique(placebo_treated_units_df$dep_var)) {
    pdat <- placebo_treated_units_df %>% filter(dep_var == dv)
    true_b <- unique(pdat$true_estimate)
    
    g <- ggplot(pdat, aes(x = estimate)) +
      geom_histogram(bins = 40) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = true_b, linetype = "dotted") +
      labs(
        title = paste0("Placebo treated-units distribution — ", dv),
        subtitle = "Histogram of placebo Treated×Post (each control treated in turn); dotted = true HTI estimate",
        x = "Estimated Treated×Post",
        y = "Count"
      )
    
    ggsave(
      filename = file.path(out_dir, "plots", paste0("placebo_treated_units_hist_", dv, ".png")),
      plot = g, width = 10, height = 5, dpi = 200
    )
  }
}

# ==========================================================
# E) Leave-one-out controls (influence)
# ==========================================================
loo_rows <- list()
controls_for_loo <- setdiff(sort(unique(net_main$state)), treated_state)

for (dv in dv_focus) {
  base_res <- run_pooled_did_one(net_main, dv)
  if (is.null(base_res)) next
  
  for (s in controls_for_loo) {
    dat <- net_main %>% filter(state != s)
    res <- run_pooled_did_one(dat, dv)
    if (is.null(res)) next
    
    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      dep_var = dv,
      dropped_state = s,
      estimate = res$estimate,
      std_error = res$std_error,
      base_estimate = base_res$estimate,
      stringsAsFactors = FALSE
    )
  }
}

loo_df <- bind_rows(loo_rows)
write_csv(loo_df, file.path(out_dir, "leave_one_out_controls.csv"))

# ==========================================================
# F) Random-drop 1–5 controls (stress test) + plots
# Start from FULL donor pool (including som/ssud) and drop 1–5 OTHER states.
# ==========================================================
states_full <- sort(unique(net_full$state))
candidate_drop <- setdiff(states_full, c(treated_state))  # "other states"

random_rows <- list()

for (dv in dv_focus) {
  base_full <- run_pooled_did_one(net_full, dv)
  if (is.null(base_full)) next
  
  for (b in 1:random_drop_draws) {
    k <- sample(random_drop_k_vals, 1)
    drop_states <- sample(candidate_drop, k)
    
    dat <- net_full %>% filter(!(state %in% drop_states))
    res <- run_pooled_did_one(dat, dv)
    if (is.null(res)) next
    
    random_rows[[length(random_rows) + 1]] <- data.frame(
      dep_var = dv,
      draw = b,
      k_dropped = k,
      estimate = res$estimate,
      std_error = res$std_error,
      base_estimate_full = base_full$estimate,
      stringsAsFactors = FALSE
    )
  }
}

random_drop_df <- bind_rows(random_rows)
write_csv(random_drop_df, file.path(out_dir, "random_drop_stress_test.csv"))

# Plot random-drop distributions
if (nrow(random_drop_df) > 0) {
  for (dv in unique(random_drop_df$dep_var)) {
    pdat <- random_drop_df %>% filter(dep_var == dv)
    base_b <- unique(pdat$base_estimate_full)
    
    g <- ggplot(pdat, aes(x = estimate)) +
      geom_histogram(bins = 40) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = base_b, linetype = "dotted") +
      labs(
        title = paste0("Random-drop stress test — ", dv),
        subtitle = paste0("Start from ALL states; drop 1–5 random OTHER states (excluding som/ssud). Dotted = baseline(all-states)"),
        x = "Estimated Treated×Post",
        y = "Count"
      )
    
    ggsave(
      filename = file.path(out_dir, "plots", paste0("random_drop_hist_", dv, ".png")),
      plot = g, width = 10, height = 5, dpi = 200
    )
  }
}

# ==========================================================
# G) Window sensitivity grid (vary pre/post around 2010)
# ==========================================================
win_rows <- list()

for (dv in dv_focus) {
  for (pre_w in pre_grid) {
    for (post_w in post_grid) {
      dat <- net_main %>%
        filter(year >= (event_year - pre_w), year <= (event_year + post_w))
      
      res <- run_pooled_did_one(dat, dv)
      if (is.null(res)) next
      
      win_rows[[length(win_rows) + 1]] <- data.frame(
        dep_var = dv,
        pre_w = pre_w,
        post_w = post_w,
        estimate = res$estimate,
        std_error = res$std_error,
        p_value = res$p_value,
        stringsAsFactors = FALSE
      )
    }
  }
}

window_sens_df <- bind_rows(win_rows)
write_csv(window_sens_df, file.path(out_dir, "window_sensitivity.csv"))

# Plot (heatmap-style) window sensitivity
if (nrow(window_sens_df) > 0) {
  for (dv in unique(window_sens_df$dep_var)) {
    pdat <- window_sens_df %>% filter(dep_var == dv)
    
    g <- ggplot(pdat, aes(x = post_w, y = pre_w, fill = estimate)) +
      geom_tile() +
      labs(
        title = paste0("Window sensitivity (pooled DiD) — ", dv),
        x = "Post window length",
        y = "Pre window length"
      )
    
    ggsave(
      filename = file.path(out_dir, "plots", paste0("window_sensitivity_heat_", dv, ".png")),
      plot = g, width = 8, height = 6, dpi = 200
    )
  }
}

# ==========================================================
# H) Wild cluster bootstrap p-values
# ==========================================================
wb_rows <- list()

B <- 999

suppressMessages({
  suppressWarnings({
    
    for (dv in dv_focus) {
      
      dat <- net_main %>% dplyr::filter(!is.na(.data[[dv]]))
      
      # baseline model (clustered)
      fit <- tryCatch(
        feols(as.formula(paste0(dv, " ~ tp | state + year")),
              data = dat, cluster = ~state),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      
      # bootstrap distribution of tp
      states <- unique(dat$state)
      G <- length(states)
      b <- rep(NA_real_, B)
      
      for (b_i in 1:B) {
        
        boot_states <- sample(states, size = G, replace = TRUE)
        state_counts <- table(boot_states)  # multiplicity per state
        
        # multiplicity weight per row
        w <- state_counts[as.character(dat$state)]
        w[is.na(w)] <- 0
        w <- as.numeric(w)
        
        fit_b <- tryCatch(
          feols(as.formula(paste0(dv, " ~ tp | state + year")),
                data = dat, weights = w),
          error = function(e) NULL
        )
        if (is.null(fit_b)) next
        
        b[b_i] <- unname(coef(fit_b)["tp"])
      }
      
      b <- b[is.finite(b)]
      if (length(b) < 200) next  # basic sanity
      
      beta_hat <- unname(coef(fit)["tp"])
      
      # percentile CI
      ci <- stats::quantile(b, probs = c(0.025, 0.975), na.rm = TRUE)
      
      # simple two-sided bootstrap p-value around 0
      b_bar <- mean(b, na.rm = TRUE)
      p_boot <- mean(abs(b - b_bar) >= abs(beta_hat - b_bar), na.rm = TRUE)
      
      wb_rows[[length(wb_rows) + 1]] <- data.frame(
        dep_var = dv,
        beta_hat = beta_hat,
        p_boot = p_boot,
        ci_lo = ci[[1]],
        ci_hi = ci[[2]],
        B_used = length(b),
        stringsAsFactors = FALSE
      )
    }
    
  })
})
wb_out <- do.call(rbind, wb_rows)
write_csv(wb_out, file.path(out_dir, "wild_cluster_bootstrap_pvalues.csv"))
 
# ==========================================================
# I) Permutation inference p-values
# (Already computed via placebo-treated distribution summary)
# ==========================================================
# You already have:
#   placebo_treated_units_summary.csv
# which includes p_perm_two_sided for each DV in dv_focus.

# ==========================================================
# J) Synthetic control / synthDiD 
# ==========================================================
synth_rows <- list()

for (dv in synth_dvs) {
  dat <- net_full %>%
    select(state, year, treated, post, tp, all_of(dv)) %>%
    filter(!is.na(.data[[dv]]))
    
    # Enforce balanced panel for this DV
  years_all <- sort(unique(dat$year))
  need_T <- length(years_all)
    
  keep_states <- dat %>%
    group_by(state) %>%
    summarise(nT = n_distinct(year), .groups = "drop") %>%
    filter(nT == need_T) %>%
    pull(state)
    
  datb <- dat %>% filter(state %in% keep_states)
    
  # Treatment indicator for synthdid needs to be 1 in treated periods
  datb <- datb %>%
    mutate(treat_on = as.integer(state == treated_state & year >= event_year))
    
  # Build panel matrices
  pm <- tryCatch(
    synthdid::panel.matrices(datb, unit = "state", time = "year",
                             outcome = dv, treatment = "treat_on"),
    error = function(e) NULL
  )
  if (is.null(pm)) next
  
  est <- tryCatch(
    synthdid::synthdid_estimate(pm$Y, pm$N0, pm$T0),
    error = function(e) NULL
  )
  if (is.null(est)) next
  
  # SE via vcov (jackknife if available)
  se <- NA_real_
  ci_lb <- NA_real_
  ci_ub <- NA_real_
  
  V <- tryCatch(
    vcov(est, method = "placebo", replications = 500),
    error = function(e) { message("vcov placebo failed: ", conditionMessage(e)); NULL }
  )
  
  if (!is.null(V)) {
    se <- sqrt(as.numeric(V))
    ci_lb <- as.numeric(est) - 1.96 * se
    ci_ub <- as.numeric(est) + 1.96 * se
  }
    
  synth_rows[[length(synth_rows) + 1]] <- data.frame(
    dep_var = dv,
    synthdid_estimate = as.numeric(est),
    std_error = se,
    lb = ci_lb,
    ub = ci_ub,
    n_units = nrow(pm$Y),
    n_years = ncol(pm$Y),
    stringsAsFactors = FALSE
  )
    
  # Save a standard synthdid plot
  # Pull weights
  # controls are the first N0 rows of pm$Y
  # ---- Manual synthdid plot (robust) ----
  
  years_plot <- sort(unique(datb$year))
  
  control_units <- rownames(pm$Y)[1:pm$N0]
  ctrl <- synthdid::synthdid_controls(est)
  
  w_sparse <- as.numeric(ctrl[, 1])
  names(w_sparse) <- rownames(ctrl)
  
  # normalize names for matching
  control_units_clean <- trimws(tolower(control_units))
  names(w_sparse) <- trimws(tolower(names(w_sparse)))
  
  omega <- rep(0, pm$N0)
  names(omega) <- control_units_clean
  
  omega[names(w_sparse)] <- w_sparse
  
  # diagnostics + rescale instead of hard stop
  missing_in_controls <- setdiff(names(w_sparse), names(omega))
  unmapped_weight <- if (length(missing_in_controls) > 0) sum(w_sparse[missing_in_controls]) else 0
  
  message("DV = ", dv,
          " | sum(w_sparse) = ", round(sum(w_sparse), 6),
          " | sum(omega) = ", round(sum(omega), 6),
          " | unmapped donors = ", length(missing_in_controls),
          " | unmapped weight = ", round(unmapped_weight, 6))
  
  if (sum(omega) == 0) {
    message("Omega is all zeros after mapping; skipping plot for dv = ", dv)
    next
  }
  if (abs(sum(omega) - 1) > 1e-6) {
    omega <- omega / sum(omega)
    message("Rescaled omega to sum to 1 (after mapping).")
  }
  
  y_treated <- as.numeric(pm$Y[pm$N0 + 1, ])
  Y_controls <- pm$Y[1:pm$N0, , drop = FALSE]
  y_synth <- as.numeric(crossprod(omega, Y_controls))
  
  treat_start_year <- years_plot[pm$T0 + 1]
  
  plot_df <- data.frame(
    year = years_plot,
    treated = y_treated,
    synthetic = y_synth
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = year)) +
    ggplot2::geom_line(ggplot2::aes(y = synthetic, color = "synthetic control")) +
    ggplot2::geom_line(ggplot2::aes(y = treated, color = "treated")) +
    ggplot2::geom_vline(xintercept = treat_start_year, linetype = "dashed") +
    ggplot2::labs(
      title = paste0("synthDiD — ", dv),
      x = NULL,
      y = dv,
      color = NULL
    ) +
    ggplot2::theme_minimal()
  
  png(file.path(out_dir, "plots", paste0("synthdid_plot_", dv, ".png")),
      width = 1600, height = 900, res = 150)
  print(p)
  dev.off()
}
  
synthdid_df <- bind_rows(synth_rows)
write_csv(synthdid_df, file.path(out_dir, "synthdid_summary.csv"))
