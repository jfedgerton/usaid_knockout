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
  # "density",   # REMOVED: not a manuscript DV (descriptive only)
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
dv_focus <- c("fragmentation", "weighted_fragmentation", "usaid_n_partners")

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
perm_dvs <- c("fragmentation", "weighted_fragmentation", "modularity",
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


# ##############################################################
# ==============================================================
# NODE-LEVEL ROBUSTNESS CHECKS
# ==============================================================
# Parallels the network-level checks above for the node-level
# DiD design described in 18_haiti_earthquake.R.
#
# Checks:
#   K) Node-level event study (pre-trends)
#   L) Node-level window sensitivity
#   M) Full vs balanced panel comparison
#   N) Node-level leave-one-out controls
#   O) Alternative exposure definitions
#   P) Node-level cluster bootstrap (org-level)
# ##############################################################

node_out_dir <- file.path(out_dir, "node_level")
dir.create(node_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(node_out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# Load node data and construct balanced panel
# ----------------------------------------------------------
nodes <- readRDS("Data/clean_nets/node_attributes_directed.rds")
nodes$state <- tolower(nodes$state)

# Balanced panel: orgs present in at least one pre AND one post year
pre_orgs_r  <- unique(nodes[nodes$year < event_year, c("state", "node_name")])
post_orgs_r <- unique(nodes[nodes$year >= event_year, c("state", "node_name")])
balanced_orgs_r <- merge(pre_orgs_r, post_orgs_r, by = c("state", "node_name"))
nodes_bal <- merge(nodes, balanced_orgs_r, by = c("state", "node_name"))

nodes_bal$treated <- as.integer(nodes_bal$state == treated_state)
nodes_bal$post    <- as.integer(nodes_bal$year >= event_year)
nodes_bal$tp      <- nodes_bal$treated * nodes_bal$post
nodes_bal$org_id  <- paste0(nodes_bal$state, "||", nodes_bal$node_name)
nodes_bal$rel_year <- nodes_bal$year - event_year

cat("\nNode balanced panel:", nrow(balanced_orgs_r), "org-state pairs,",
    nrow(nodes_bal), "total rows\n")

# Pre-treatment exposure variables
nodes_pre_r <- nodes_bal[nodes_bal$year < event_year, ]

pre_usaid_funded_r <- aggregate(usaid_funded ~ state + node_name,
  data = nodes_pre_r, FUN = max, na.rm = TRUE)
names(pre_usaid_funded_r)[3] <- "pre_usaid_funded"

pre_exposure_r <- aggregate(usaid_exposure_share ~ state + node_name,
  data = nodes_pre_r, FUN = mean, na.rm = TRUE)
names(pre_exposure_r)[3] <- "pre_usaid_exposure_share"

pre_tied_r <- aggregate(tied_to_usaid ~ state + node_name,
  data = nodes_pre_r, FUN = max, na.rm = TRUE)
names(pre_tied_r)[3] <- "pre_tied_to_usaid"

nodes_bal <- merge(nodes_bal, pre_usaid_funded_r, by = c("state", "node_name"), all.x = TRUE)
nodes_bal <- merge(nodes_bal, pre_exposure_r,     by = c("state", "node_name"), all.x = TRUE)
nodes_bal <- merge(nodes_bal, pre_tied_r,         by = c("state", "node_name"), all.x = TRUE)

# Exclude states from main node pool
nodes_bal_main <- nodes_bal[!nodes_bal$state %in% exclude_main, ]

# Node DVs
node_dv_list <- c("in_degree", "n_funders", "sole_usaid_funded",
                   "n_non_usaid_funders", "total_degree",
                   "rewire_in", "rewire_out")
node_dv_list <- node_dv_list[node_dv_list %in% names(nodes_bal)]

# Headline DVs for expensive checks
node_dv_focus <- c("in_degree", "sole_usaid_funded", "n_funders")
node_dv_focus <- node_dv_focus[node_dv_focus %in% names(nodes_bal)]

cat("Node DVs:", paste(node_dv_list, collapse = ", "), "\n")
cat("Node DVs (focus):", paste(node_dv_focus, collapse = ", "), "\n")


# ==========================================================
# K) Node-level event study (pre-trends)
# ==========================================================
cat("\n=== K) Node-level event study ===\n")
node_es_rows <- list()
node_lead_tests <- list()

for (dv in node_dv_list) {
  dat <- nodes_bal_main[!is.na(nodes_bal_main[[dv]]), ]

  fit_es <- tryCatch(
    feols(as.formula(paste0(dv, " ~ i(rel_year, treated, ref = -1) | org_id + year")),
          data = dat, cluster = ~org_id),
    error = function(e) NULL
  )
  if (is.null(fit_es)) next

  ct <- as.data.frame(coeftable(fit_es))
  ct$term <- rownames(ct)
  rownames(ct) <- NULL

  ct <- ct[grepl("^rel_year::", ct$term), ]
  rel_raw <- ct$term
  rel_raw <- sub("^rel_year::", "", rel_raw)
  rel_raw <- sub("(:treated|#treated)$", "", rel_raw)
  ct$rel_year <- suppressWarnings(as.integer(rel_raw))
  ct$dep_var <- dv

  node_es_rows[[length(node_es_rows) + 1]] <- ct

  # Joint test of leads (rel_year <= -2)
  lead_terms <- ct$term[!is.na(ct$rel_year) & ct$rel_year <= -2]
  lead_terms <- lead_terms[lead_terms %in% names(coef(fit_es))]
  lead_df <- length(lead_terms)

  lead_stat <- NA_real_
  lead_p <- NA_real_

  if (lead_df > 0) {
    V <- tryCatch(vcov(fit_es)[lead_terms, lead_terms, drop = FALSE],
                  error = function(e) NULL)
    if (!is.null(V)) {
      w <- safe_wald_chisq(coef(fit_es)[lead_terms], V, tol = 1e-10)
      lead_stat <- w$stat
      lead_p <- w$p
      lead_df <- w$df
    }
  }

  node_lead_tests[[length(node_lead_tests) + 1]] <- data.frame(
    dep_var = dv, n_leads = lead_df,
    wald_stat = lead_stat, p_value = lead_p,
    stringsAsFactors = FALSE
  )

  cat("  ", dv, ": lead test p =", round(lead_p, 4), "\n")
}

node_es_coefs <- bind_rows(node_es_rows)
node_es_leads <- bind_rows(node_lead_tests)

write_csv(node_es_coefs, file.path(node_out_dir, "node_eventstudy_coefs.csv"))
write_csv(node_es_leads, file.path(node_out_dir, "node_eventstudy_leads_joint_test.csv"))

# Event study plots
for (dv in unique(node_es_coefs$dep_var)) {
  pdat <- node_es_coefs[node_es_coefs$dep_var == dv & !is.na(node_es_coefs$rel_year), ]
  if (nrow(pdat) == 0) next

  names(pdat)[names(pdat) == "Estimate"] <- "est"
  names(pdat)[names(pdat) == "Std. Error"] <- "se_val"

  if (!all(c("est", "se_val") %in% names(pdat))) next

  pdat$lb_es <- pdat$est - 1.96 * pdat$se_val
  pdat$ub_es <- pdat$est + 1.96 * pdat$se_val

  p <- ggplot(pdat, aes(x = rel_year, y = est)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey50") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lb_es, ymax = ub_es), width = 0.3) +
    theme_classic() +
    labs(title = paste0("Node-level event study: ", dv),
         x = "Years relative to 2010", y = "Coefficient")
  ggsave(file.path(node_out_dir, "plots", paste0("node_es_", dv, ".png")),
         p, width = 8, height = 5)
}


# ==========================================================
# L) Node-level window sensitivity
# ==========================================================
cat("\n=== L) Node-level window sensitivity ===\n")
node_win_rows <- list()

for (dv in node_dv_focus) {
  for (pre_w in pre_grid) {
    for (post_w in post_grid) {
      yr_range <- (event_year - pre_w):(event_year + post_w)
      dat <- nodes_bal_main[nodes_bal_main$year %in% yr_range & !is.na(nodes_bal_main[[dv]]), ]

      fit <- tryCatch(
        feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
              data = dat, cluster = ~org_id),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      if (!("tp" %in% names(coef(fit)))) next

      node_win_rows[[length(node_win_rows) + 1]] <- data.frame(
        dep_var   = dv,
        pre_w     = pre_w,
        post_w    = post_w,
        estimate  = as.numeric(coef(fit)["tp"]),
        std_error = as.numeric(se(fit)["tp"]),
        p_value   = as.numeric(pvalue(fit)["tp"]),
        n         = nobs(fit),
        stringsAsFactors = FALSE
      )
    }
  }
}

node_win_sens <- bind_rows(node_win_rows)
write_csv(node_win_sens, file.path(node_out_dir, "node_window_sensitivity.csv"))

# Heatmap plots
for (dv in unique(node_win_sens$dep_var)) {
  pdat <- node_win_sens[node_win_sens$dep_var == dv, ]
  pdat$pre_w  <- factor(pdat$pre_w,  levels = sort(unique(pdat$pre_w)))
  pdat$post_w <- factor(pdat$post_w, levels = sort(unique(pdat$post_w)))

  g <- ggplot(pdat, aes(x = post_w, y = pre_w, fill = estimate)) +
    geom_tile() +
    labs(title = paste0("Node-level window sensitivity: ", dv),
         x = "Post window", y = "Pre window")
  ggsave(file.path(node_out_dir, "plots", paste0("node_win_sens_", dv, ".png")),
         g, width = 8, height = 6)
}


# ==========================================================
# M) Full (unbalanced) vs balanced panel comparison
# ==========================================================
cat("\n=== M) Full vs balanced panel comparison ===\n")
nodes_full <- nodes[!nodes$state %in% exclude_main, ]
nodes_full$treated <- as.integer(nodes_full$state == treated_state)
nodes_full$post    <- as.integer(nodes_full$year >= event_year)
nodes_full$tp      <- nodes_full$treated * nodes_full$post
nodes_full$org_id  <- paste0(nodes_full$state, "||", nodes_full$node_name)

panel_comp_rows <- list()

for (dv in node_dv_focus) {
  # Balanced panel estimate
  dat_b <- nodes_bal_main[!is.na(nodes_bal_main[[dv]]), ]
  fit_b <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
          data = dat_b, cluster = ~org_id),
    error = function(e) NULL
  )

  # Full (unbalanced) panel estimate
  dat_f <- nodes_full[!is.na(nodes_full[[dv]]), ]
  fit_f <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
          data = dat_f, cluster = ~org_id),
    error = function(e) NULL
  )

  if (!is.null(fit_b) && "tp" %in% names(coef(fit_b))) {
    panel_comp_rows[[length(panel_comp_rows) + 1]] <- data.frame(
      dep_var = dv, panel = "balanced",
      estimate  = as.numeric(coef(fit_b)["tp"]),
      std_error = as.numeric(se(fit_b)["tp"]),
      p_value   = as.numeric(pvalue(fit_b)["tp"]),
      n = nobs(fit_b), n_orgs = length(unique(dat_b$org_id)),
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(fit_f) && "tp" %in% names(coef(fit_f))) {
    panel_comp_rows[[length(panel_comp_rows) + 1]] <- data.frame(
      dep_var = dv, panel = "full_unbalanced",
      estimate  = as.numeric(coef(fit_f)["tp"]),
      std_error = as.numeric(se(fit_f)["tp"]),
      p_value   = as.numeric(pvalue(fit_f)["tp"]),
      n = nobs(fit_f), n_orgs = length(unique(dat_f$org_id)),
      stringsAsFactors = FALSE
    )
  }
}

panel_comp <- bind_rows(panel_comp_rows)
write_csv(panel_comp, file.path(node_out_dir, "node_balanced_vs_full_panel.csv"))
cat("Panel comparison:\n")
print(panel_comp)


# ==========================================================
# N) Node-level leave-one-out controls
# ==========================================================
cat("\n=== N) Node-level leave-one-out ===\n")
node_loo_rows <- list()
node_controls <- setdiff(sort(unique(nodes_bal_main$state)), treated_state)

for (dv in node_dv_focus) {
  # Baseline
  dat_base <- nodes_bal_main[!is.na(nodes_bal_main[[dv]]), ]
  fit_base <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
          data = dat_base, cluster = ~org_id),
    error = function(e) NULL
  )
  if (is.null(fit_base)) next
  base_est <- as.numeric(coef(fit_base)["tp"])

  for (s in node_controls) {
    dat_loo <- dat_base[dat_base$state != s, ]
    fit_loo <- tryCatch(
      feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
            data = dat_loo, cluster = ~org_id),
      error = function(e) NULL
    )
    if (is.null(fit_loo)) next
    if (!("tp" %in% names(coef(fit_loo)))) next

    node_loo_rows[[length(node_loo_rows) + 1]] <- data.frame(
      dep_var       = dv,
      dropped_state = s,
      estimate      = as.numeric(coef(fit_loo)["tp"]),
      std_error     = as.numeric(se(fit_loo)["tp"]),
      base_estimate = base_est,
      stringsAsFactors = FALSE
    )
  }
}

node_loo_df <- bind_rows(node_loo_rows)
write_csv(node_loo_df, file.path(node_out_dir, "node_leave_one_out.csv"))


# ==========================================================
# O) Alternative exposure definitions
# ==========================================================
# Test whether results are sensitive to how we define USAID
# exposure at the node level. Compare:
#   (1) pre_usaid_funded (direct: USAID out-neighbor)
#   (2) pre_tied_to_usaid (any-mode neighbor of USAID)
#   (3) pre_usaid_exposure_share (continuous: share of neighbors connected)
cat("\n=== O) Alternative exposure definitions ===\n")
alt_exp_rows <- list()

exposure_vars <- c("pre_usaid_funded", "pre_tied_to_usaid", "pre_usaid_exposure_share")
exposure_vars <- exposure_vars[exposure_vars %in% names(nodes_bal_main)]

for (dv in node_dv_focus) {
  dat <- nodes_bal_main[!is.na(nodes_bal_main[[dv]]), ]

  for (ev in exposure_vars) {
    # Interaction with tp
    fml <- as.formula(paste0(dv, " ~ tp + tp:", ev, " | org_id + year"))
    fit <- tryCatch(
      feols(fml, data = dat, cluster = ~org_id),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    for (cf in names(coef(fit))) {
      alt_exp_rows[[length(alt_exp_rows) + 1]] <- data.frame(
        dep_var       = dv,
        exposure_var  = ev,
        term          = cf,
        estimate      = as.numeric(coef(fit)[cf]),
        std_error     = as.numeric(se(fit)[cf]),
        p_value       = as.numeric(pvalue(fit)[cf]),
        n             = nobs(fit),
        stringsAsFactors = FALSE
      )
    }
  }
}

alt_exp_df <- bind_rows(alt_exp_rows)
write_csv(alt_exp_df, file.path(node_out_dir, "node_alternative_exposure_defs.csv"))


# ==========================================================
# P) Node-level cluster bootstrap (org-level resampling)
# ==========================================================
# Parallels section H but at the node level.
# Resamples organizations (clusters) within each state, preserving
# the state-level treatment structure.
cat("\n=== P) Node-level cluster bootstrap ===\n")

node_boot_B <- 999
node_boot_min_success <- 200

node_boot_rows <- list()

for (dv in node_dv_focus) {
  dat <- nodes_bal_main[!is.na(nodes_bal_main[[dv]]), ]

  fit <- tryCatch(
    feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
          data = dat, cluster = ~org_id),
    error = function(e) NULL
  )
  if (is.null(fit)) next

  beta_hat <- as.numeric(coef(fit)["tp"])

  # Get list of orgs by state for stratified resampling
  orgs_by_state <- split(unique(dat$org_id), dat$state[match(unique(dat$org_id), dat$org_id)])

  b_vals <- rep(NA_real_, node_boot_B)
  failures <- 0L

  for (b_i in seq_len(node_boot_B)) {
    # Stratified resample: within each state, resample orgs with replacement
    boot_orgs <- character(0)
    for (st in names(orgs_by_state)) {
      st_orgs <- orgs_by_state[[st]]
      boot_orgs <- c(boot_orgs, sample(st_orgs, length(st_orgs), replace = TRUE))
    }

    # Build multiplicity weights
    org_counts <- table(boot_orgs)
    w <- as.numeric(org_counts[dat$org_id])
    w[is.na(w)] <- 0

    fit_b <- tryCatch(
      feols(as.formula(paste0(dv, " ~ tp | org_id + year")),
            data = dat, weights = w),
      error = function(e) NULL
    )

    if (is.null(fit_b) || is.na(coef(fit_b)["tp"])) {
      failures <- failures + 1L
      next
    }

    b_vals[b_i] <- as.numeric(coef(fit_b)["tp"])
  }

  b_vals <- b_vals[is.finite(b_vals)]
  if (length(b_vals) < node_boot_min_success) next

  ci <- quantile(b_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  p_boot <- mean(abs(b_vals - beta_hat) >= abs(beta_hat), na.rm = TRUE)

  node_boot_rows[[length(node_boot_rows) + 1]] <- data.frame(
    dep_var          = dv,
    beta_hat         = beta_hat,
    p_boot_two_sided = p_boot,
    ci_lo            = ci[[1]],
    ci_hi            = ci[[2]],
    B_requested      = node_boot_B,
    B_used           = length(b_vals),
    fail_rate        = failures / node_boot_B,
    stringsAsFactors = FALSE
  )

  cat("  ", dv, ": beta =", round(beta_hat, 4),
      " boot p =", round(p_boot, 4),
      " CI = [", round(ci[[1]], 4), ",", round(ci[[2]], 4), "]\n")
}

node_boot_out <- bind_rows(node_boot_rows)
write_csv(node_boot_out, file.path(node_out_dir, "node_cluster_bootstrap.csv"))


cat("\n=== All robustness checks complete (network + node level) ===\n")
