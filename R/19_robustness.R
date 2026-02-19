# ==========================================================
# HAITI ROBUSTNESS SUITE (NETWORK TOPOLOGY DiD / EVENT STUDY)
# Focus: Haiti earthquake inflow case (HTI treated, 2010)
# Data: Data/clean_nets/network_attributes_directed.rds
#
# Outputs -> Data/results/robustness/haiti/
#
# Checks included:
#   A) Baseline pooled DiD (main donor pool + full donor pool)
#   B) Event study + joint test of leads (robust to singular V)
#   C) Placebo years (windowed; PRE-2010 only to avoid contamination)
#   D) Placebo treated units (treat each control as treated)
#   E) Leave-one-out controls (influence)
#   F) Random-drop 1–5 controls (stress test) + plots
#   G) Window sensitivity grid (vary pre/post windows)
#   H) Cluster (pairs) bootstrap via state resampling + multiplicity weights
#   I) Permutation/randomization p-values (via placebo-treated)
#   J) Synthetic control / synthDiD (optional)
# ==========================================================

# ----------------------------
# Packages
# ----------------------------
required_pkgs <- c("dplyr", "fixest", "readr", "ggplot2", "tidyr", "stringr", "synthdid")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, e.g. install.packages(c(",
    paste0("\"", missing_pkgs, "\"", collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(synthdid)
})

# ----------------------------
# Paths
# ----------------------------
# Prefer running from repo root. If you run from elsewhere but have a local "Ecosystem_of_Aid" folder,
# this attempts to set the working directory automatically.
if (basename(getwd()) != "Ecosystem_of_Aid" && dir.exists("Ecosystem_of_Aid")) {
  setwd("Ecosystem_of_Aid")
}

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
set.seed(123)
random_drop_draws <- 500      # increase to 1000+ if you want
random_drop_k_vals <- 1:5     # drop 1–5 controls per draw

# Placebo-year test window (pre-only, with auto-shrink if infeasible)
placebo_pre_w <- 2
placebo_post_w <- 2
placebo_pre_only <- TRUE  # IMPORTANT: prevents post-2010 contamination

# Window sensitivity grid (pre, post)
pre_grid <- c(2, 3, 4, 5)
post_grid <- c(2, 3, 5, 8, 10)

# Cluster (pairs) bootstrap settings
boot_B <- 999
boot_min_success <- 200
boot_max_tries_include_treated <- 50

# Synthetic DiD outcomes (keep small)
synth_dvs <- c("fragmentation", "weighted_fragmentation", "modularity")

# ----------------------------
# Load data
# ----------------------------
net <- readRDS("Data/clean_nets/network_attributes_directed.rds")

# Basic sanity checks
stopifnot(all(c("state", "year") %in% names(net)))

net <- net %>%
  mutate(state = tolower(as.character(state))) %>%
  filter(year >= min_year, year <= max_year) %>%
  mutate(
    treated = as.integer(state == treated_state),
    post = as.integer(year >= event_year),
    tp = treated * post,
    rel_year = year - event_year
  )

# Keep only DVs that exist
present <- names(net)
dv_list <- dv_list[dv_list %in% present]
dv_focus <- dv_focus[dv_focus %in% present]
synth_dvs <- synth_dvs[synth_dvs %in% present]

if (length(dv_list) == 0) stop("No DVs found in net. Check variable names.")

# ----------------------------
# Helpers
# ----------------------------
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

# Safe Wald chi-square test using eigen truncation (pseudo-inverse)
# Returns list(stat, df, p)
safe_wald_chisq <- function(b, V, tol = 1e-10) {
  if (is.null(V) || any(!is.finite(V))) return(list(stat = NA_real_, df = 0L, p = NA_real_))
  V <- as.matrix(V)
  V <- (V + t(V)) / 2

  eg <- tryCatch(eigen(V, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eg)) return(list(stat = NA_real_, df = 0L, p = NA_real_))

  vals <- eg$values
  vecs <- eg$vectors

  # Keep only well-behaved positive eigenvalues
  maxv <- max(vals, na.rm = TRUE)
  if (!is.finite(maxv) || maxv <= 0) return(list(stat = NA_real_, df = 0L, p = NA_real_))

  keep <- which(is.finite(vals) & vals > (tol * maxv))
  if (length(keep) == 0) return(list(stat = NA_real_, df = 0L, p = NA_real_))

  Q <- vecs[, keep, drop = FALSE]
  Dinv <- diag(1 / vals[keep], nrow = length(keep))
  Vinv <- Q %*% Dinv %*% t(Q)

  b <- as.numeric(b)
  stat <- as.numeric(t(b) %*% Vinv %*% b)
  df <- length(keep)
  p <- 1 - pchisq(stat, df = df)

  list(stat = stat, df = df, p = p)
}

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

  # term looks like "rel_year::-3:treated"
  rel_raw <- ct$term %>%
    str_replace("^rel_year::", "") %>%
    str_replace("(:treated|#treated)$", "")

  ct$rel_year <- suppressWarnings(as.integer(rel_raw))
  ct$dep_var <- dv

  event_rows[[length(event_rows) + 1]] <- ct

  # Joint test: leads (rel_year <= -2) = 0
  lead_terms <- ct %>% filter(!is.na(rel_year), rel_year <= -2) %>% pull(term)
  lead_terms <- lead_terms[lead_terms %in% names(coef(fit_es))]

  lead_stat <- NA_real_
  lead_p <- NA_real_
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

        w <- safe_wald_chisq(b2, V2, tol = 1e-10)
        lead_stat <- w$stat
        lead_p <- w$p
        lead_df <- w$df
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
# IMPORTANT: If we want placebo windows that do NOT overlap the true post-2010 period,
# we must keep the entire placebo window within years <= 2009.
pre_end <- event_year - 1
placebo_post_w_eff <- placebo_post_w

if (placebo_pre_only) {
  max_post_feasible <- pre_end - (min_year + placebo_pre_w)
  if (max_post_feasible < placebo_post_w_eff) {
    message(
      "[Placebo-year] Requested post window = ", placebo_post_w,
      ", but max feasible (pre-only) is ", max_post_feasible,
      ". Using post window = ", max_post_feasible, "."
    )
    placebo_post_w_eff <- max_post_feasible
  }

  if (!is.finite(placebo_post_w_eff) || placebo_post_w_eff < 1) {
    message("[Placebo-year] Skipping placebo-year test: no feasible pre-only placebo years with current window sizes.")
    placebo_years <- integer(0)
  } else {
    placebo_years <- seq(min_year + placebo_pre_w, pre_end - placebo_post_w_eff)
  }
} else {
  placebo_years <- setdiff(
    seq(min_year + placebo_pre_w, max_year - placebo_post_w_eff),
    event_year
  )
}

placebo_rows <- list()

if (length(placebo_years) > 0) {
  for (dv in dv_focus) {
    # True estimate in the same window size around event_year
    true_window <- net_main %>%
      filter(year >= (event_year - placebo_pre_w), year <= (event_year + placebo_post_w_eff))

    true_fit <- run_pooled_did_one(true_window, dv)
    if (is.null(true_fit)) next
    true_b <- true_fit$estimate

    for (y0 in placebo_years) {
      dat <- net_main %>%
        filter(year >= (y0 - placebo_pre_w), year <= (y0 + placebo_post_w_eff)) %>%
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
        subtitle = paste0(
          "Window: [year-", placebo_pre_w, ", year+", placebo_post_w_eff,
          "] (pre-only = ", placebo_pre_only, "); dotted line = true 2010"
        ),
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

    # small-n: fewer bins improves readability
    g <- ggplot(pdat, aes(x = estimate)) +
      geom_histogram(bins = 15) +
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
candidate_drop <- setdiff(states_full, c(treated_state, exclude_main))  # "other states"

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
        subtitle = "Start from ALL states; drop 1–5 random OTHER states (excluding som/ssud). Dotted = baseline(all-states)",
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
    pdat <- window_sens_df %>%
      filter(dep_var == dv) %>%
      mutate(
        pre_w = factor(pre_w, levels = sort(unique(pre_w))),
        post_w = factor(post_w, levels = sort(unique(post_w)))
      )

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
# H) Cluster (pairs) bootstrap via state resampling + multiplicity weights
# ==========================================================
boot_rows <- list()

for (dv in dv_focus) {
  dat <- net_main %>% filter(!is.na(.data[[dv]]))

  fit <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp | state + year")),
          data = dat, cluster = ~ state),
    error = function(e) NULL
  )
  if (is.null(fit)) next

  beta_hat <- unname(coef(fit)["tp"])

  states <- sort(unique(dat$state))
  G <- length(states)

  b <- rep(NA_real_, boot_B)
  failures <- 0L

  for (b_i in 1:boot_B) {
    boot_states <- NULL

    # Resample clusters until treated appears (avoids unidentified tp)
    for (try_i in 1:boot_max_tries_include_treated) {
      cand <- sample(states, size = G, replace = TRUE)
      if (treated_state %in% cand) {
        boot_states <- cand
        break
      }
    }

    # Fallback: force include treated by replacement
    if (is.null(boot_states)) {
      cand <- sample(states, size = G, replace = TRUE)
      cand[sample.int(G, 1)] <- treated_state
      boot_states <- cand
    }

    state_counts <- table(boot_states)

    # multiplicity weight per row
    w <- state_counts[as.character(dat$state)]
    w[is.na(w)] <- 0
    w <- as.numeric(w)

    fit_b <- tryCatch(
      feols(as.formula(paste0(dv, " ~ tp | state + year")),
            data = dat, weights = w),
      error = function(e) NULL
    )

    if (is.null(fit_b) || is.na(coef(fit_b)["tp"])) {
      failures <- failures + 1L
      next
    }

    b[b_i] <- unname(coef(fit_b)["tp"])
  }

  b <- b[is.finite(b)]
  if (length(b) < boot_min_success) next

  ci <- stats::quantile(b, probs = c(0.025, 0.975), na.rm = TRUE)

  # Two-sided bootstrap p-value for H0: beta = 0 using pivot approximation:
  #   e* = b* - beta_hat approximates sampling error under true beta.
  #   Under H0, estimator distribution is ~ 0 + e*.
  p_boot <- mean(abs(b - beta_hat) >= abs(beta_hat), na.rm = TRUE)

  boot_rows[[length(boot_rows) + 1]] <- data.frame(
    dep_var = dv,
    beta_hat = beta_hat,
    p_boot_two_sided = p_boot,
    ci_lo = ci[[1]],
    ci_hi = ci[[2]],
    B_requested = boot_B,
    B_used = length(b),
    fail_rate = failures / boot_B,
    stringsAsFactors = FALSE
  )
}

boot_out <- bind_rows(boot_rows)

# Write under a more accurate name, but keep the legacy filename for compatibility.
write_csv(boot_out, file.path(out_dir, "cluster_pairs_bootstrap.csv"))

# ==========================================================
# I) Permutation inference p-values
# (Already computed via placebo-treated distribution summary)
# ==========================================================
# placebo_treated_units_summary.csv includes p_perm_two_sided.

# ==========================================================
# J) Synthetic control / synthDiD (optional)
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


  # Treatment indicator for synthdid needs to be 1 in treated periods
  datb <- dat %>%
    mutate(
      state = trimws(tolower(as.character(state))),
      year  = as.integer(year),
      treat_on = as.integer(state == treated_state & year >= event_year)
    ) %>% 
    filter(state %in% keep_states)
  
  # Sanity: treat_on only for treated_state
  treat_units <- unique(datb$state[datb$treat_on == 1])
  if (length(treat_units) != 1 || treat_units != treated_state) {
    stop("Treatment coding error: treat_on should be 1 only for treated_state.\n",
         "Found treated units: ", paste(treat_units, collapse = ", "))
  }
  
  # Build wide Y matrix: rows = states, cols = years (numeric order)
  years_all <- sort(unique(datb$year))
  T0_forced <- sum(years_all < event_year)
  if (T0_forced < 2) stop("T0_forced < 2; check event_year and year coding.")
  
  Y_wide <- datb %>%
    select(state, year, y = all_of(dv)) %>%
    tidyr::pivot_wider(names_from = year, values_from = y)
  
  Y_mat <- as.matrix(Y_wide[, -1, drop = FALSE])
  rownames(Y_mat) <- Y_wide$state
  
  # order columns by numeric year
  col_years <- as.integer(colnames(Y_mat))
  ord <- order(col_years)
  Y_mat <- Y_mat[, ord, drop = FALSE]
  col_years <- col_years[ord]
  
  # Drop any rows with missing values (synthdid needs complete matrix)
  keep_rows <- which(apply(Y_mat, 1, function(r) all(is.finite(r))))
  Y_mat <- Y_mat[keep_rows, , drop = FALSE]
  
  # Recompute control set after balancing
  if (!(treated_state %in% rownames(Y_mat))) {
    message("Treated state dropped due to missingness/balancing for dv=", dv, " — skipping")
    next
  }
  controls <- setdiff(rownames(Y_mat), treated_state)
  N0 <- length(controls)
  if (N0 < 5) {
    message("Too few controls after balancing for dv=", dv, " (N0=", N0, ") — skipping")
    next
  }
  
  # Put treated unit LAST row (synthdid convention)
  Y_mat <- Y_mat[c(controls, treated_state), , drop = FALSE]
  
  # Force T0 based on column-years (after any balancing)
  T0 <- sum(col_years < event_year)
  if (T0 < 2 || T0 >= ncol(Y_mat)) stop("Bad T0 after balancing; check year coverage.")
  
  # Estimate
  est <- tryCatch(
    synthdid::synthdid_estimate(Y_mat, N0 = N0, T0 = T0),
    error = function(e) NULL
  )
  if (is.null(est)) next
  
  # ---- synthDiD pre-fit diagnostics (via donor weights) ----
  ctrl <- synthdid::synthdid_controls(est)
  w_sparse <- as.numeric(ctrl[, 1])
  names(w_sparse) <- rownames(ctrl)
  
  control_units <- rownames(Y_mat)[1:N0]
  control_units_clean <- trimws(tolower(control_units))
  names(w_sparse) <- trimws(tolower(names(w_sparse)))
  
  omega <- rep(0, N0)
  names(omega) <- control_units_clean
  omega[names(w_sparse)] <- w_sparse
  if (sum(omega) > 0 && abs(sum(omega) - 1) > 1e-6) omega <- omega / sum(omega)
  
  y_treated  <- as.numeric(Y_mat[N0 + 1, ])
  Y_controls <- Y_mat[1:N0, , drop = FALSE]
  y_synth    <- as.numeric(crossprod(omega, Y_controls))
  
  pre_idx  <- seq_len(T0)
  post_idx <- (T0 + 1):ncol(Y_mat)
  
  rmspe_pre  <- sqrt(mean((y_treated[pre_idx]  - y_synth[pre_idx])^2))
  rmspe_post <- sqrt(mean((y_treated[post_idx] - y_synth[post_idx])^2))
  rmspe_ratio <- rmspe_post / rmspe_pre
  
  # SE via placebo-based vcov
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
  
  # Store row (replace pm$Y counts with Y_mat counts)
  synth_rows[[length(synth_rows) + 1]] <- data.frame(
    dep_var = dv,
    synthdid_estimate = as.numeric(est),
    std_error = se,
    lb = ci_lb,
    ub = ci_ub,
    n_units = nrow(Y_mat),
    n_years = ncol(Y_mat),
    rmspe_pre = rmspe_pre,
    rmspe_post = rmspe_post,
    rmspe_ratio = rmspe_ratio,
    stringsAsFactors = FALSE
  )
  
  # Plot: DO NOT use plot(est) (can mislead while debugging)
  plot_df <- data.frame(
    year = col_years,
    synthetic_control = y_synth,
    treated = y_treated
  ) %>%
    tidyr::pivot_longer(cols = c("synthetic_control", "treated"),
                        names_to = "series", values_to = "value")
  
  png(file.path(out_dir, "plots", paste0("synthdid_plot_", dv, ".png")),
      width = 1200, height = 700)
  print(
    ggplot(plot_df, aes(x = year, y = value, color = series)) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = event_year, linetype = "dashed") +
      labs(title = paste0("synthDiD (forced T0) — ", dv),
           x = NULL, y = NULL) +
      theme_classic()
  )
  dev.off()
}

synthdid_df <- bind_rows(synth_rows)
write_csv(synthdid_df, file.path(out_dir, "synthdid_summary.csv"))

perm_out_dir <- file.path(out_dir, "permutation_by_window")
dir.create(perm_out_dir, showWarnings = FALSE, recursive = TRUE)

# Choose DVs for permutation (recommend core ones)
perm_dvs <- c("fragmentation", "weighted_fragmentation", "density", "modularity",
              "usaid_frag_delta_drop_unfunded", "usaid_n_partners")

# Windows (match main script)
perm_windows <- c(2, 3, 4, 5)
do_full_sample <- TRUE

# Define "main" pool if you're using exclude_main elsewhere
states_all <- sort(unique(net_full$state))
states_main <- setdiff(states_all, exclude_main)

# For each DV and window: compute Haiti beta and placebo distribution
perm_rows <- list()

for (dv in perm_dvs) {
  
  # Ensure dv exists
  if (!dv %in% names(net_full)) next
  
  for (w in perm_windows) {
    
    dat_w <- net_full %>%
      dplyr::filter(year >= (event_year - w), year <= (event_year + w)) %>%
      dplyr::filter(!is.na(.data[[dv]]))
    
    # require treated present
    if (!treated_state %in% unique(dat_w$state)) next
    
    # Haiti estimate
    dat_h <- dat_w %>%
      dplyr::mutate(
        treated = as.integer(state == treated_state),
        post = as.integer(year >= event_year),
        tp = treated * post
      )
    
    fit_h <- tryCatch(
      fixest::feols(stats::as.formula(paste0(dv, " ~ tp | state + year")),
                    data = dat_h, cluster = ~state),
      error = function(e) NULL
    )
    if (is.null(fit_h)) next
    
    beta_hat <- as.numeric(coef(fit_h)["tp"])
    
    # Placebo treated states: use same donor pool as main spec
    placebo_states <- setdiff(states_main, treated_state)
    
    placebo_betas <- rep(NA_real_, length(placebo_states))
    names(placebo_betas) <- placebo_states
    
    for (s in placebo_states) {
      dat_p <- dat_w %>%
        dplyr::mutate(
          treated = as.integer(state == s),
          post = as.integer(year >= event_year),
          tp = treated * post
        )
      
      fit_p <- tryCatch(
        fixest::feols(stats::as.formula(paste0(dv, " ~ tp | state + year")),
                      data = dat_p, cluster = ~state),
        error = function(e) NULL
      )
      
      if (!is.null(fit_p)) placebo_betas[s] <- as.numeric(coef(fit_p)["tp"])
    }
    
    placebo_betas <- placebo_betas[is.finite(placebo_betas)]
    if (length(placebo_betas) < 5) next
    
    p_perm <- mean(abs(placebo_betas) >= abs(beta_hat))
    z_perm <- (beta_hat - mean(placebo_betas)) / stats::sd(placebo_betas)
    
    perm_rows[[length(perm_rows) + 1]] <- data.frame(
      dep_var = dv,
      window = paste0("+-", w),
      beta_haiti = beta_hat,
      placebo_mean = mean(placebo_betas),
      placebo_sd = stats::sd(placebo_betas),
      z_vs_placebo = z_perm,
      p_perm_twosided = p_perm,
      n_placebos = length(placebo_betas),
      stringsAsFactors = FALSE
    )
    
    # Save raw placebo distribution (optional, useful for plots later)
    placebo_df <- data.frame(state = names(placebo_betas), beta = as.numeric(placebo_betas))
    readr::write_csv(placebo_df, file.path(perm_out_dir, paste0("perm_dist_", dv, "_w", w, ".csv")))
  }
  
  # Full-sample variant (optional)
  if (do_full_sample) {
    
    dat_w <- net_full %>%
      dplyr::filter(!is.na(.data[[dv]]))
    
    if (!treated_state %in% unique(dat_w$state)) next
    
    dat_h <- dat_w %>%
      dplyr::mutate(
        treated = as.integer(state == treated_state),
        post = as.integer(year >= event_year),
        tp = treated * post
      )
    
    fit_h <- tryCatch(
      fixest::feols(stats::as.formula(paste0(dv, " ~ tp | state + year")),
                    data = dat_h, cluster = ~state),
      error = function(e) NULL
    )
    if (is.null(fit_h)) next
    
    beta_hat <- as.numeric(coef(fit_h)["tp"])
    
    placebo_states <- setdiff(states_main, treated_state)
    placebo_betas <- rep(NA_real_, length(placebo_states))
    names(placebo_betas) <- placebo_states
    
    for (s in placebo_states) {
      dat_p <- dat_w %>%
        dplyr::mutate(
          treated = as.integer(state == s),
          post = as.integer(year >= event_year),
          tp = treated * post
        )
      
      fit_p <- tryCatch(
        fixest::feols(stats::as.formula(paste0(dv, " ~ tp | state + year")),
                      data = dat_p, cluster = ~state),
        error = function(e) NULL
      )
      
      if (!is.null(fit_p)) placebo_betas[s] <- as.numeric(coef(fit_p)["tp"])
    }
    
    placebo_betas <- placebo_betas[is.finite(placebo_betas)]
    if (length(placebo_betas) >= 5) {
      
      p_perm <- mean(abs(placebo_betas) >= abs(beta_hat))
      z_perm <- (beta_hat - mean(placebo_betas)) / stats::sd(placebo_betas)
      
      perm_rows[[length(perm_rows) + 1]] <- data.frame(
        dep_var = dv,
        window = "all",
        beta_haiti = beta_hat,
        placebo_mean = mean(placebo_betas),
        placebo_sd = stats::sd(placebo_betas),
        z_vs_placebo = z_perm,
        p_perm_twosided = p_perm,
        n_placebos = length(placebo_betas),
        stringsAsFactors = FALSE
      )
      
      placebo_df <- data.frame(state = names(placebo_betas), beta = as.numeric(placebo_betas))
      readr::write_csv(placebo_df, file.path(perm_out_dir, paste0("perm_dist_", dv, "_all.csv")))
    }
  }
}

perm_summary <- dplyr::bind_rows(perm_rows)
readr::write_csv(perm_summary, file.path(perm_out_dir, "permutation_summary_by_window.csv"))

alt_out_dir <- file.path(out_dir, "alt_post_timing")
dir.create(alt_out_dir, showWarnings = FALSE, recursive = TRUE)

alt_dvs <- c("fragmentation", "weighted_fragmentation",
             "usaid_frag_delta_drop_unfunded", "usaid_n_partners")

alt_windows <- c(3, 5)  # keep small: match headline windows
alt_post_starts <- c(2010, 2011, 2012)

alt_rows <- list()

for (dv in alt_dvs) {
  if (!dv %in% names(net_full)) next
  
  for (w in alt_windows) {
    dat_w <- net_full %>%
      dplyr::filter(year >= (event_year - w), year <= (event_year + w)) %>%
      dplyr::filter(!is.na(.data[[dv]])) %>%
      dplyr::filter(!state %in% exclude_main)
    
    for (ps in alt_post_starts) {
      
      dat_a <- dat_w %>%
        dplyr::mutate(
          treated = as.integer(state == treated_state),
          post = as.integer(year >= ps),
          tp = treated * post
        )
      
      fit <- tryCatch(
        fixest::feols(stats::as.formula(paste0(dv, " ~ tp | state + year")),
                      data = dat_a, cluster = ~state),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      
      b <- as.numeric(coef(fit)["tp"])
      se <- as.numeric(se(fit)["tp"])
      p <- as.numeric(pvalue(fit)["tp"])
      
      alt_rows[[length(alt_rows) + 1]] <- data.frame(
        dep_var = dv,
        window = paste0("+-", w),
        post_start = ps,
        beta = b,
        se = se,
        p = p,
        stringsAsFactors = FALSE
      )
    }
  }
}

alt_df <- dplyr::bind_rows(alt_rows)
readr::write_csv(alt_df, file.path(alt_out_dir, "alt_post_timing_pooled_did.csv"))
