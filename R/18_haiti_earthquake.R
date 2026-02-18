# ==========================================================
# Haiti main effects (simple)
# - Runs pooled TWFE DiD for HTI with multiple windows
# - Exports ONE table: estimate(tp), SE, CI, p for each DV x window
# - Optional: one clean plot per DV (effect vs window)
# ==========================================================

library(dplyr)
library(fixest)
library(readr)
library(ggplot2)

if (!grepl("Ecosystem_of_Aid", getwd())) setwd("Ecosystem_of_Aid")

out_dir <- "Data/results/haiti_main_effects"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# Load data
# ----------------------------------------------------------
net <- readRDS("Data/clean_nets/network_attributes_directed.rds")

treat_state <- "hti"
event_year  <- 2010

net <- net %>%
  mutate(
    state = tolower(state),
    treated = as.integer(state == treat_state),
    post    = as.integer(year >= event_year),
    tp      = treated * post
  )

# ----------------------------------------------------------
# Outcomes (edit to your headline set)
# ----------------------------------------------------------
dv_focus <- c(
  "fragmentation",
  "weighted_fragmentation",
  "modularity",
  "usaid_frag_delta_drop_unfunded",
  "usaid_n_partners"
)
dv_focus <- dv_focus[dv_focus %in% names(net)]
if (length(dv_focus) == 0) stop("No dv_focus found in data. Check names.")

# ----------------------------------------------------------
# Windows
#   - Use NA / "all" for full-sample pooled estimate (optional)
# ----------------------------------------------------------
windows <- c(2, 3, 4, 5)
include_full_sample <- TRUE

# Optional exclusions (only if you want them in main effects; otherwise set empty)
drop_states <- character(0)

# ----------------------------------------------------------
# Estimate
# ----------------------------------------------------------
rows <- list()

for (dv in dv_focus) {
  
  # Full sample (optional reference)
  if (include_full_sample) {
    
    dat <- net %>%
      filter(!(state %in% drop_states)) %>%
      filter(!is.na(.data[[dv]]))
    
    fit <- tryCatch(
      feols(as.formula(paste0(dv, " ~ tp | state + year")),
            data = dat, cluster = ~state),
      error = function(e) NULL
    )
    
    if (!is.null(fit) && ("tp" %in% names(coef(fit)))) {
      rows[[length(rows) + 1]] <- data.frame(
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
  
  # Windowed estimates
  for (w in windows) {
    
    dat <- net %>%
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

tbl <- bind_rows(rows) %>%
  mutate(
    lb = estimate - 1.96 * std_error,
    ub = estimate + 1.96 * std_error
  ) %>%
  arrange(dep_var, factor(window, levels = c("+-2","+-3", "+-4", "+-5","all")))

write_csv(tbl, file.path(out_dir, "haiti_main_effects_windows.csv"))
print(tbl)


for (dv in unique(tbl$dep_var)) {
  
  dfp <- tbl %>%
    filter(dep_var == dv) %>%
    mutate(
      window_num = dplyr::case_when(
        window == "+-2" ~ 2,
        window == "+-3" ~ 3,
        window == "+-4" ~ 4,
        window == "+-5" ~ 5,
        window == "all" ~ 99,
        TRUE ~ NA_real_
      )
    )
  
  p <- ggplot(dfp, aes(x = window, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point() +
    geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2) +
    theme_classic() +
    labs(title = paste0("Haiti main effect across windows: ", dv),
         x = "Window around 2010",
         y = "DiD estimate (tp)")
  
  ggsave(file.path(out_dir, "figures", paste0("main_effect_windows_", dv, ".png")),
         p, width = 7, height = 4)
}

