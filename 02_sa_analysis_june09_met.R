# PROPENSITY SCORE MODELS ###############

library(tidyverse)
library(haven)
library(rlang)

## 1. SETUP & DIRECTORIES -------------------------------------------------------
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/Sensitivities/"
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
clust_vars <- c("studyid_rev", "episode", "time0")

# Default variables
default_out_vars <- c("sens_discont", "time_discont", 
                      "sens_od_acute", "time_od_acute", 
                      "sens_all_acute", "time_all_acute")

baselinevars <- c("sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
                  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
                  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
                  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss", "prior_episode")

continuous_vars <- c("ddose_prev", "year", "ret_length", "prior_episode")

## 2. DATA PREPARATION ----------------------------------------------------------
cat("Calculating cumulative OAT experience from longitudinal data...\n")

# A. Calculate eli_cum_oat from the comprehensive longitudinal file
dat_cum_oat <- read_sas(file.path(DATA_DIR, "TT_analytic_alls_m.sas7bdat")) |>
  arrange(studyid_rev, episode, myepisode, time0) |>
  group_by(studyid_rev, episode) |>
  mutate(
    gap_start = if_else(myepisode == 1, replace_na(gap_betwn_epi, 0), 0),
    cum_gap_days_epi = first(cumsum(gap_start)),
    cum_gap_days = cum_gap_days_epi[1],
    time_on_OAT = as.numeric(time0 - f_oat_dt),
    cumulative_OAT = time_on_OAT - cum_gap_days,
    eli_cum_oat = if_else(cumulative_OAT > 365, 1, 0)
  ) |>
  ungroup() |>
  # Keep only the merging keys and the new variable to avoid duplicating columns
  select(studyid_rev, episode, time0, eli_cum_oat) 


cat("Loading and merging with primary analytic cohort...\n")

# B. Load main analytic dataset and join the new eli_cum_oat variable
raw_first <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat")) |>
  # Left join ensures we only keep the rows already validated for the analytic cohort
  left_join(dat_cum_oat, by = c("studyid_rev", "episode", "time0")) |>
  mutate(
    # Collapse reduction categories
    trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456),
    prior_episode = as.numeric(episode - 1),
    
    # Setup conditions for easy filtering in the SA loop
    high_ddose = if_else(ddose_prev >= 80, 1, 0),
    post_fentanyl = if_else(time0 >= as.Date("2012-04-01"), 1, 0),
    high_fentanyl_area = if_else(clnt_lha >= 200 & clnt_lha < 400 & time0 >= as.Date("2016-04-01"), 1, 0)
  ) |>
  rename(trt = trt_s456) |>
  filter(trt %in% TRT_LEVELS) |>
  mutate(trt = factor(trt, levels = TRT_LEVELS))

# C. Split by trial
dat_tr2_base <- raw_first |> filter(trial_gap %in% 3:4)
dat_tr3_base <- raw_first |> filter(trial_cat == "5-14")

rm(raw_first); gc()

## 3. DEFINE SENSITIVITY SCENARIOS ----------------------------------------------
sa_scenarios <- list(
  SA02_Induction      = list(filter_expr = expr(induction == 1), out_vars = default_out_vars),
  SA03_Dose_LT80      = list(filter_expr = expr(high_ddose == 0), out_vars = default_out_vars),
  SA04_Dose_GE80      = list(filter_expr = expr(high_ddose == 1), out_vars = default_out_vars),
  
  # Aligned with primary analysis strata
  SA05_Ret_LE30       = list(filter_expr = expr(ret_length <= 30), out_vars = default_out_vars),
  SA06_Ret_GT30       = list(filter_expr = expr(ret_length > 30), out_vars = default_out_vars),
  
  SA09_HighFentArea   = list(filter_expr = expr(high_fentanyl_area == 1), out_vars = default_out_vars),
  SA10_PostFentEra    = list(filter_expr = expr(post_fentanyl == 1), out_vars = default_out_vars),
  
  # 1-year cumulative OAT
  SA11_CumOAT_GE1yr   = list(filter_expr = expr(eli_cum_oat == 1), out_vars = default_out_vars),
  
  SA13_30dWashout     = list(filter_expr = expr(gap30_eli == 1), out_vars = default_out_vars),
  
  # Outcome Swaps
  SA12_NoTruncation   = list(
    filter_expr = expr(TRUE), # No filtering, keep full cohort
    out_vars = c("s_discont", "t_discont", "s_od_acute", "t_od_acute", "s_all_acute", "t_all_acute")
  ),
  SA16_Composite_OD   = list(
    filter_expr = expr(TRUE), # No filtering, keep full cohort
    out_vars = c("sens_discont", "time_discont", "sens_od_comp_outcome", "time_od_comp_outcome")
  )
)

## 4. EXECUTION ENGINE ----------------------------------------------------------
run_sa_model <- function(base_dat, trial_name, sa_name, sa_params) {
  cat(sprintf("\nRunning %s for %s...\n", sa_name, trial_name))
  
  # 1. Apply specific filtering
  analytic_dat <- base_dat |> filter(!!sa_params$filter_expr)
  
  # Special pre-processing for SA16 composite outcome
  if (sa_name == "SA16_Composite_OD") {
    analytic_dat <- analytic_dat |> mutate(
      time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
      sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE)
    )
  }
  
  # Positivity & sample size safety check
  n_obs <- nrow(analytic_dat)
  if(n_obs < 500) {
    warning(sprintf("Skipping %s for %s: Insufficient N (%d) - Positivity violation likely.", sa_name, trial_name, n_obs))
    return(NULL)
  }
  
  # 2. Fit PS Model
  # Uses the robust ps_model_wt function from your 00_functions.R
  ps_fit <- ps_model_wt(
    dat = analytic_dat,
    exposure = "trt",
    baselinevars = baselinevars,
    continuous_vars = continuous_vars,
    out_vars = sa_params$out_vars,
    clust_vars = clust_vars
  )
  
  # 3. Save Output
  out_file <- paste0(RES_DIR, "ps_wt_", trial_name, "_", sa_name, ".rds")
  saveRDS(ps_fit, out_file)
  cat(sprintf("  ✓ Saved %d observations to %s\n", n_obs, out_file))
}


## 5. RUN ALL SCENARIOS LOOP ----------------------------------------------------
# Execute sequentially across all definitions and both trials
for (sa_name in names(sa_scenarios)) {
  params <- sa_scenarios[[sa_name]]
  
  # Run Trial 2 (3-4 days)
  run_sa_model(dat_tr2_base, "tr2_34d", sa_name, params)
  
  # Run Trial 3 (5-14 days)
  run_sa_model(dat_tr3_base, "tr3_514d", sa_name, params)
}




# OUTCOME MODELS ##########

library(tidyverse)
library(survival)
library(data.table)

predict_cuminc_rmst_fast <- function(Dat_expanded, trt_var, covs, horizons,
                                     cluster_var, samples = 1000, conf_int = TRUE,
                                     ref_group = "no reversion", wt_type = "wt_trunc",
                                     rmst_tau = 365) {
  
  # 1. Fit pooled logistic regression
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  model <- glm(plr_formula, data = Dat_expanded,
               family = quasibinomial(), weights = Dat_expanded[[wt_type]])
  
  V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                        sandwich = TRUE, fix = FALSE)
  coefs <- coef(model)
  
  # 2. Reference group & Factor Level Preservation
  is_trt_factor <- is.factor(Dat_expanded[[trt_var]])
  if (is_trt_factor) {
    trt_levels <- levels(Dat_expanded[[trt_var]])
  } else {
    trt_levels <- as.character(unique(Dat_expanded[[trt_var]]))
  }
  
  if (is.null(ref_group)) ref_group <- trt_levels[1]
  
  # 3. Pre-compute Design Matrices using data.table for extreme speed
  dt <- as.data.table(Dat_expanded)
  setorderv(dt, c(cluster_var, "stop")) 
  
  # Build a list of model matrices (one for each treatment level) ONCE.
  mm_list <- lapply(trt_levels, function(lvl) {
    dt_temp <- copy(dt)
    
    # Explicitly preserve factor levels so we don't trigger the contrasts error
    if (is_trt_factor) {
      dt_temp[[trt_var]] <- factor(lvl, levels = trt_levels)
    } else {
      dt_temp[[trt_var]] <- lvl
    }
    
    # Use xlev = model$xlevels to strictly enforce the model's known factor levels
    model.matrix(delete.response(terms(model)), data = dt_temp, xlev = model$xlevels)
  })
  names(mm_list) <- trt_levels
  
  # 4. Vectorized Risk Function (Runs strictly on pre-computed matrices)
  fast_risk_fun <- function(coef_draw) {
    res_list <- vector("list", length(trt_levels))
    
    for (i in seq_along(trt_levels)) {
      lvl <- trt_levels[i]
      
      # Matrix multiplication for all observations at once
      dt[, p_event := plogis(mm_list[[lvl]] %*% coef_draw)]
      
      # Fast cumulative product of survival probabilities grouped by cluster
      dt[, cum_surv := cumprod(1 - p_event), by = cluster_var]
      
      # Extract mean risk at exact specified horizons
      risks_h <- numeric(length(horizons))
      for (j in seq_along(horizons)) {
        h <- horizons[j]
        risk_val <- 1 - dt[stop <= h, .(last_surv = cum_surv[.N]), by = cluster_var][, mean(last_surv)]
        risks_h[j] <- risk_val
      }
      
      res_list[[i]] <- data.table(trt = lvl, horizon = horizons, risk = risks_h)
    }
    rbindlist(res_list)
  }
  
  # 5. Point Estimates
  point <- fast_risk_fun(coefs)
  
  # Return early if no CI requested
  if (!conf_int) {
    res <- list(risks = as_tibble(point))
    
    if (!is.null(rmst_tau)) {
      if (max(horizons) < rmst_tau) stop("max(horizons) must be >= rmst_tau to compute RMST.")
      
      point_rmst <- point %>%
        filter(horizon <= rmst_tau) %>%
        group_by(trt) %>%
        arrange(horizon) %>%
        mutate(
          survival = 1 - risk,
          time_diff = horizon - dplyr::lag(horizon, default = 0),
          area = dplyr::lag(survival, default = 1) * time_diff
        ) %>%
        summarise(RMST = sum(area), .groups = "drop") %>%
        mutate(
          ref_RMST = RMST[trt == ref_group],
          RMST_diff = ifelse(trt == ref_group, 0, RMST - ref_RMST)
        ) %>% select(-ref_RMST)
      
      res$rmst <- point_rmst
    }
    return(res)
  }
  
  # 6. Sequential Simulation for Confidence Intervals (Safe Mode)
  sim_coefs <- mvtnorm::rmvnorm(n = samples, mean = coefs, sigma = V)
  
  # Standard lapply - completely avoids parallel processing instability
  sim_list <- lapply(seq_len(nrow(sim_coefs)), function(i) {
    fast_risk_fun(sim_coefs[i, ])
  })
  
  sim_df <- rbindlist(sim_list, idcol = "sim")
  
  # 7. Compute Final Statistics (Standard dplyr format)
  risk_ci <- sim_df %>%
    group_by(trt, horizon) %>%
    summarise(
      risk_lower = quantile(risk, 0.025),
      risk_upper = quantile(risk, 0.975),
      risk = mean(risk), 
      .groups = "drop"
    )
  
  sim_contrasts <- sim_df %>%
    group_by(sim, horizon) %>%
    mutate(ref_risk = risk[trt == ref_group]) %>%
    filter(trt != ref_group) %>%
    mutate(RD = risk - ref_risk, RR = risk / ref_risk)
  
  contrast_ci <- sim_contrasts %>%
    group_by(trt, horizon) %>%
    summarise(
      RD_lower = quantile(RD, 0.025),
      RD_upper = quantile(RD, 0.975),
      RD = mean(RD),
      RR_lower = quantile(RR, 0.025),
      RR_upper = quantile(RR, 0.975),
      RR = mean(RR),
      .groups = "drop"
    )
  
  res <- list(risks = risk_ci, contrasts = contrast_ci, ref_group = ref_group)
  
  # 8. Compute RMST distributions if requested
  if (!is.null(rmst_tau)) {
    if (max(horizons) < rmst_tau) stop("max(horizons) must be >= rmst_tau to compute RMST.")
    
    sim_rmst <- sim_df %>%
      filter(horizon <= rmst_tau) %>%
      group_by(sim, trt) %>%
      arrange(horizon) %>%
      mutate(
        survival = 1 - risk,
        time_diff = horizon - dplyr::lag(horizon, default = 0),
        area = dplyr::lag(survival, default = 1) * time_diff
      ) %>%
      summarise(RMST = sum(area), .groups = "drop")
    
    rmst_ci <- sim_rmst %>%
      group_by(trt) %>%
      summarise(
        RMST_lower = quantile(RMST, 0.025),
        RMST_upper = quantile(RMST, 0.975),
        RMST = mean(RMST),
        .groups = "drop"
      )
    
    rmst_contrasts <- sim_rmst %>%
      group_by(sim) %>%
      mutate(ref_RMST = RMST[trt == ref_group]) %>%
      filter(trt != ref_group) %>%
      mutate(RMST_diff = RMST - ref_RMST) %>%
      group_by(trt) %>%
      summarise(
        RMST_diff_lower = quantile(RMST_diff, 0.025),
        RMST_diff_upper = quantile(RMST_diff, 0.975),
        RMST_diff = mean(RMST_diff),
        .groups = "drop"
      )
    
    res$rmst <- rmst_ci
    res$rmst_contrasts <- rmst_contrasts
  }
  
  return(res)
}

## 1. SETUP & DIRECTORIES -------------------------------------------------------
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/Sensitivities/"
RISK_DIR  <- paste0(RES_DIR, "Risk/")
dir.create(RISK_DIR, showWarnings = FALSE, recursive = TRUE)

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(unique(c(HORIZONS_CURVE, c(30, 60, 90, 180))))

## 2. DEFINE SCENARIO OUTCOME MAPPINGS ------------------------------------------
# We map the exact time and event variables needed for each specific SA
# If an SA isn't explicitly listed here, it falls back to the "default" setup.
sa_outcome_map <- list(
  default = list(
    od  = list(time = "time_od_acute", event = "sens_od_acute"),
    dis = list(time = "time_discont", event = "sens_discont")
  ),
  SA12_NoTruncation = list(
    od  = list(time = "t_od_acute", event = "s_od_acute"),
    dis = list(time = "t_discont", event = "s_discont")
  ),
  SA16_Composite_OD = list(
    od  = list(time = "time_od_comp_outcome", event = "sens_od_comp_outcome"),
    dis = list(time = "time_discont", event = "sens_discont") # Discont doesn't change
  )
)

# List of all SA names generated in the previous step
sa_names <- c("SA02_Induction", "SA03_Dose_LT80", "SA04_Dose_GE80", 
              "SA05_Ret_LE30", "SA06_Ret_GT30", "SA09_HighFentArea", 
              "SA10_PostFentEra", "SA11_CumOAT_GE1yr", "SA13_30dWashout", 
              "SA12_NoTruncation", "SA16_Composite_OD")

## 3. OUTCOME EXECUTION ENGINE --------------------------------------------------
run_sa_outcomes <- function(trial_name, sa_name) {
  
  cat(sprintf("\n--- Running Outcomes for %s | %s ---\n", trial_name, sa_name))
  
  # A. Load the propensity score weights
  ps_file <- paste0(RES_DIR, "ps_wt_", trial_name, "_", sa_name, ".rds")
  if (!file.exists(ps_file)) {
    warning(sprintf("File missing, skipping: %s", ps_file))
    return(NULL)
  }
  
  dat_weights <- readRDS(ps_file)$data_weights
  
  # B. Determine specific covariates and outcomes
  # Trial 2 uses prior_episode; Trial 3 uses year
  adj_covs <- if (trial_name == "tr2_34d") c("ddose_prev", "prior_episode") else c("ddose_prev", "year")
  
  # Lookup outcome variables (use specific if exists, otherwise default)
  out_map <- if (sa_name %in% names(sa_outcome_map)) sa_outcome_map[[sa_name]] else sa_outcome_map[["default"]]
  
  results <- list(od = list(), dis = list())
  
  # C. OVERDOSE OUTCOME
  cat("  -> Processing Overdose...\n")
  long_dat_od <- expand_person_period(
    dat_weights, 
    time_var = out_map$od$time, 
    event_var = out_map$od$event, 
    interval_length = 14
  )
  
  # 1. Point Estimates (Curve)
  results$od$curve <- predict_cuminc_rmst_fast(
    Dat_expanded = long_dat_od, trt_var = "trt", covs = adj_covs,
    horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
    samples = 499, conf_int = FALSE, ref_group = "no reversion", wt_type = "wt_trunc"
  )
  
  # 2. Confidence Intervals (Simulation)
  results$od$ci <- predict_cuminc_rmst_fast(
    Dat_expanded = long_dat_od, trt_var = "trt", covs = adj_covs,
    horizons = HORIZONS_CI, cluster_var = "studyid_rev",
    samples = 499, conf_int = TRUE, ref_group = "no reversion", wt_type = "wt_trunc"
  )
  
  rm(long_dat_od); gc()
  
  # D. DISCONTINUATION OUTCOME
  cat("  -> Processing Discontinuation...\n")
  long_dat_dis <- expand_person_period(
    dat_weights, 
    time_var = out_map$dis$time, 
    event_var = out_map$dis$event, 
    interval_length = 14
  )
  
  # 1. Point Estimates (Curve)
  results$dis$curve <- predict_cuminc_rmst_fast(
    Dat_expanded = long_dat_dis, trt_var = "trt", covs = adj_covs,
    horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
    samples = 499, conf_int = FALSE, ref_group = "no reversion", wt_type = "wt_trunc"
  )
  
  # 2. Confidence Intervals (Simulation)
  results$dis$ci <- predict_cuminc_rmst_fast(
    Dat_expanded = long_dat_dis, trt_var = "trt", covs = adj_covs,
    horizons = HORIZONS_CI, cluster_var = "studyid_rev",
    samples = 499, conf_int = TRUE, ref_group = "no reversion", wt_type = "wt_trunc"
  )
  
  rm(long_dat_dis); gc()
  
  # E. SAVE FINAL OBJECT
  out_file <- paste0(RISK_DIR, "risk_", trial_name, "_", sa_name, ".rds")
  saveRDS(results, out_file)
  cat(sprintf("  ✓ Successfully saved risk estimates to %s\n", out_file))
}

## 4. RUN AUTOMATED LOOP --------------------------------------------------------
for (sa in sa_names) {
  # Run Trial 2
  run_sa_outcomes(trial_name = "tr2_34d", sa_name = sa)
  
  # Run Trial 3
  run_sa_outcomes(trial_name = "tr3_514d", sa_name = sa)
}



# MULTIPLE TIME ZERO ###########
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/Sensitivities/"
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
clust_vars <- c("studyid_rev", "episode", "time0")

# Default variables
default_out_vars <- c("sens_discont", "time_discont", 
                      "sens_od_acute", "time_od_acute", 
                      "sens_all_acute", "time_all_acute")

baselinevars <- c("sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
                  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
                  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
                  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

continuous_vars <- c("ddose_prev", "year", "ret_length")

## 1. DATA PREPARATION ----------------------------------------------------------
raw_first <- read_sas(file.path(DATA_DIR, "TT_analytic_Rprev_sens_nw_m.sas7bdat")) |>
  arrange(studyid_rev, episode, time0) |>
  mutate(
    trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456),
    prior_episode = as.numeric(episode - 1),
  ) |>
  rename(trt = trt_s456) |>
  filter(trt %in% TRT_LEVELS) |>
  mutate(trt = factor(trt, levels = TRT_LEVELS))|>
  group_by(studyid_rev, episode)|>
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()
  )|> 
  ungroup()|> 
  filter(trial_gap %in% 3:14)|>
  group_by(studyid_rev, episode)|>
  mutate(
    days_since_prev = as.numeric(time0 - lag(time0)),
    eligible30 = is.na(days_since_prev) | days_since_prev >= 30,
    eligible15 = is.na(days_since_prev) | days_since_prev >= 15
  )|>
  ungroup()


# C. Split by trial
dat_tr2_sa2 <- raw_first |> filter(eligible30)|>  filter(trial_gap %in% 3:4)
dat_tr3_sa2 <- raw_first |> filter(eligible30)|> filter(trial_gap %in% 5:14)



## 2. PS MODELS ########

ps_tr2_sa2 <- ps_model_wt(
  dat             = dat_tr2_sa2,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_gap", "prior_episode"),
  continuous_vars = c(continuous_vars, "prior_gap", "prior_episode"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_sa2, paste0(RES_DIR, "ps_weight_tr2_34d_sa2.rds"))



ps_tr3_sa2 <- ps_model_wt(
  dat             = dat_tr3_sa2,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_gap", "prior_episode"),
  continuous_vars = c(continuous_vars, "prior_gap", "prior_episode"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr3_sa2, paste0(RES_DIR, "ps_weight_tr3_514d_sa2.rds"))




## 3. OUTCOME MODELS ########
predict_cuminc_rmst <- function(Dat_expanded, trt_var, covs, horizons,
                                cluster_var, samples = 1000, conf_int = TRUE,
                                ref_group = "no reversion", wt_type="wt_trunc",
                                rmst_tau = 365) {
  
  # Formula for pooled logistic regression
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  # Fit GLM using quasibinomial to suppress fractional weight warnings
  model <- glm(plr_formula, data = Dat_expanded,
               family = quasibinomial(), weights = Dat_expanded[[wt_type]])
  
  # Cluster-robust covariance
  V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                        sandwich = TRUE, fix = FALSE)
  
  coefs <- coef(model)
  
  # Function to calculate risk at a given horizon
  risk_fun <- function(coef_draw, h) {
    nd <- Dat_expanded %>% filter(stop <= h)
    
    risks <- lapply(unique(Dat_expanded[[trt_var]]), function(trt_level) {
      nd_trt <- nd
      nd_trt[[trt_var]] <- trt_level
      
      mm_trt <- model.matrix(delete.response(terms(model)), nd_trt)
      p_trt <- plogis(mm_trt %*% coef_draw)
      nd_trt$p_event <- p_trt
      
      df_pred <- nd_trt %>%
        group_by(.data[[cluster_var]]) %>%
        summarise(surv = prod(1 - p_event), .groups = "drop") %>%
        mutate(risk = 1 - surv)
      
      tibble(trt = as.character(trt_level), horizon = h, risk = mean(df_pred$risk))
    })
    
    bind_rows(risks)
  }
  
  # Reference group check
  if(is.null(ref_group)) {
    ref_group <- as.character(unique(Dat_expanded[[trt_var]]))[1]
  }
  
  # Point estimates
  point <- bind_rows(lapply(horizons, function(h) risk_fun(coefs, h)))
  
  if(!conf_int) {
    res <- list(risks = point)
    
    # Point estimate for RMST if requested
    if(!is.null(rmst_tau)) {
      if(max(horizons) < rmst_tau) stop("max(horizons) must be >= rmst_tau to compute RMST.")
      
      point_rmst <- point %>%
        filter(horizon <= rmst_tau) %>%
        group_by(trt) %>%
        arrange(horizon) %>%
        mutate(
          survival = 1 - risk,
          time_diff = horizon - dplyr::lag(horizon, default = 0),
          area = dplyr::lag(survival, default = 1) * time_diff
        ) %>%
        summarise(RMST = sum(area), .groups = "drop") %>%
        mutate(
          ref_RMST = RMST[trt == ref_group],
          RMST_diff = ifelse(trt == ref_group, 0, RMST - ref_RMST)
        ) %>% dplyr::select(-ref_RMST) # <-- Explicitly namespaced to avoid MASS collision
      
      res$rmst <- point_rmst
    }
    return(res)
  }
  
  # Simulate coefficient draws for CI
  sim_coefs <- mvtnorm::rmvnorm(n = samples, mean = coefs, sigma = V)
  
  sim_list <- lapply(seq_len(nrow(sim_coefs)), function(i) {
    bind_rows(lapply(horizons, function(h) risk_fun(sim_coefs[i,], h)))
  })
  sim_df <- bind_rows(sim_list, .id = "sim")
  
  # Risk CI
  risk_ci <- sim_df %>%
    group_by(trt, horizon) %>%
    summarise(
      risk_lower = quantile(risk, 0.025),
      risk_upper = quantile(risk, 0.975),
      risk = mean(risk),  # point estimate from simulation
      .groups = "drop"
    )
  
  # Compute RD and RR in each simulation
  sim_contrasts <- sim_df %>%
    group_by(sim, horizon) %>%
    mutate(ref_risk = risk[trt == ref_group]) %>%
    filter(trt != ref_group) %>%
    mutate(
      RD = risk - ref_risk,
      RR = risk / ref_risk
    )
  
  # Summarise RD and RR
  contrast_ci <- sim_contrasts %>%
    group_by(trt, horizon) %>%
    summarise(
      RD_lower = quantile(RD, 0.025),
      RD_upper = quantile(RD, 0.975),
      RD = mean(RD),
      RR_lower = quantile(RR, 0.025),
      RR_upper = quantile(RR, 0.975),
      RR = mean(RR),
      .groups = "drop"
    )
  
  res <- list(
    risks = risk_ci,
    contrasts = contrast_ci,
    ref_group = ref_group
  )
  
  # Compute RMST distribution if requested
  if (!is.null(rmst_tau)) {
    if (max(horizons) < rmst_tau) stop("max(horizons) must be >= rmst_tau to compute RMST.")
    
    sim_rmst <- sim_df %>%
      filter(horizon <= rmst_tau) %>%
      group_by(sim, trt) %>%
      arrange(horizon) %>%
      mutate(
        survival = 1 - risk,
        time_diff = horizon - dplyr::lag(horizon, default = 0),
        area = dplyr::lag(survival, default = 1) * time_diff
      ) %>%
      summarise(RMST = sum(area), .groups = "drop")
    
    rmst_ci <- sim_rmst %>%
      group_by(trt) %>%
      summarise(
        RMST_lower = quantile(RMST, 0.025),
        RMST_upper = quantile(RMST, 0.975),
        RMST = mean(RMST),
        .groups = "drop"
      )
    
    rmst_contrasts <- sim_rmst %>%
      group_by(sim) %>%
      mutate(ref_RMST = RMST[trt == ref_group]) %>%
      filter(trt != ref_group) %>%
      mutate(RMST_diff = RMST - ref_RMST) %>%
      group_by(trt) %>%
      summarise(
        RMST_diff_lower = quantile(RMST_diff, 0.025),
        RMST_diff_upper = quantile(RMST_diff, 0.975),
        RMST_diff = mean(RMST_diff),
        .groups = "drop"
      )
    
    res$rmst <- rmst_ci
    res$rmst_contrasts <- rmst_contrasts
  }
  
  return(res)
}


HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))

risk <- list(
  tr2 = list(od = list(), dis = list()),
  tr3 = list(od = list(), dis = list())
)

### MET interruption of 3-4 days #############

ps_fit <- readRDS(paste0(RES_DIR, "ps_weight_tr2_34d_sa2.rds"))

# Overdose
long_dat_od <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)


risk$tr2$od$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat_od, 
  trt_var = "trt",
  covs = c("ddose_prev", "prior_episode"),
  horizons = HORIZONS_CURVE, 
  cluster_var = "studyid_rev",
  samples = 199, 
  conf_int = FALSE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc", 
  rmst_tau = 365
)


risk$tr2$od$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat_od, 
  trt_var = "trt",
  covs = c("ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI, 
  cluster_var = "studyid_rev",
  samples = 199, 
  conf_int = TRUE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc", 
  rmst_tau = 365
)

# Discontinuation
long_dat_dis <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr2$dis$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat_dis, 
  trt_var = "trt",
  covs = c("ddose_prev", "prior_episode"),
  horizons = HORIZONS_CURVE, 
  cluster_var = "studyid_rev",
  samples = 199, 
  conf_int = FALSE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc", 
  rmst_tau = 365
)


risk$tr2$dis$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat_dis, 
  trt_var = "trt",
  covs = c("ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI, 
  cluster_var = "studyid_rev",
  samples = 199, 
  conf_int = TRUE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc", 
  rmst_tau = 365
)

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_34d_sa2.rds"))





### MET interruption of 5-14 days #############

ps_fit <- readRDS(paste0(RES_DIR, "ps_weight_tr3_514d_sa2.rds"))

# Overdose
long_dat_od <- expand_person_period(ps_fit$data_weights,
                                    time_var  = "time_od_acute",
                                    event_var = "sens_od_acute",
                                    interval_length = 14)


risk$tr3$od$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat_od, 
  trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, 
  cluster_var = "studyid_rev",
  samples = 199, 
  conf_int = FALSE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc", 
  rmst_tau = 365
)


risk$tr3$od$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat_od, 
  trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, 
  cluster_var = "studyid_rev",
  samples = 199, 
  conf_int = TRUE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc", 
  rmst_tau = 365
)

# Discontinuation
long_dat_dis <- expand_person_period(ps_fit$data_weights,
                                     time_var  = "time_discont",
                                     event_var = "sens_discont",
                                     interval_length = 14)

risk$tr3$dis$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat_dis, 
  trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, 
  cluster_var = "studyid_rev",
  samples = 199, 
  conf_int = FALSE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc", 
  rmst_tau = 365
)


risk$tr3$dis$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat_dis, 
  trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, 
  cluster_var = "studyid_rev",
  samples = 199, 
  conf_int = TRUE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc", 
  rmst_tau = 365
)

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_514d_sa2.rds"))


# FOUR SA #############

library(tidyverse)
library(haven)
library(survival)
library(sandwich)
library(data.table)

## 1. SETUP & DIRECTORIES -------------------------------------------------------
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"
SA_DIR    <- paste0(RES_DIR, "Sensitivities/")
RISK_DIR  <- paste0(SA_DIR, "Risk/")

dir.create(SA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RISK_DIR, showWarnings = FALSE, recursive = TRUE)

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(unique(c(HORIZONS_CURVE, c(30, 60, 90, 180))))

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
clust_vars <- c("studyid_rev", "episode", "time0")

# Default variables
default_out_vars <- c("sens_discont", "time_discont", 
                      "sens_od_acute", "time_od_acute", 
                      "sens_all_acute", "time_all_acute")

baselinevars <- c("sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
                  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
                  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
                  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

continuous_vars <- c("ddose_prev", "year", "ret_length")

# CORE OUTCOME FUNCTION
predict_cuminc_rmst_fast <- function(Dat_expanded, trt_var, covs, horizons,
                                     cluster_var, samples = 1000, conf_int = TRUE,
                                     ref_group = "no reversion", wt_type = "wt_trunc",
                                     rmst_tau = 365) {
  
  # 1. Fit pooled logistic regression
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  model <- glm(plr_formula, data = Dat_expanded,
               family = quasibinomial(), weights = Dat_expanded[[wt_type]])
  
  V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                        sandwich = TRUE, fix = FALSE)
  coefs <- coef(model)
  
  # 2. Reference group & Factor Level Preservation
  is_trt_factor <- is.factor(Dat_expanded[[trt_var]])
  if (is_trt_factor) {
    trt_levels <- levels(Dat_expanded[[trt_var]])
  } else {
    trt_levels <- as.character(unique(Dat_expanded[[trt_var]]))
  }
  
  if (is.null(ref_group)) ref_group <- trt_levels[1]
  
  # 3. Pre-compute Design Matrices using data.table for extreme speed
  dt <- as.data.table(Dat_expanded)
  setorderv(dt, c(cluster_var, "stop")) 
  
  # Build a list of model matrices (one for each treatment level) ONCE.
  mm_list <- lapply(trt_levels, function(lvl) {
    dt_temp <- copy(dt)
    if (is_trt_factor) {
      dt_temp[[trt_var]] <- factor(lvl, levels = trt_levels)
    } else {
      dt_temp[[trt_var]] <- lvl
    }
    model.matrix(delete.response(terms(model)), data = dt_temp, xlev = model$xlevels)
  })
  names(mm_list) <- trt_levels
  
  # 4. Vectorized Risk Function (Runs strictly on pre-computed matrices)
  fast_risk_fun <- function(coef_draw) {
    res_list <- vector("list", length(trt_levels))
    for (i in seq_along(trt_levels)) {
      lvl <- trt_levels[i]
      dt[, p_event := plogis(mm_list[[lvl]] %*% coef_draw)]
      dt[, cum_surv := cumprod(1 - p_event), by = cluster_var]
      
      risks_h <- numeric(length(horizons))
      for (j in seq_along(horizons)) {
        h <- horizons[j]
        risk_val <- 1 - dt[stop <= h, .(last_surv = cum_surv[.N]), by = cluster_var][, mean(last_surv)]
        risks_h[j] <- risk_val
      }
      res_list[[i]] <- data.table(trt = lvl, horizon = horizons, risk = risks_h)
    }
    rbindlist(res_list)
  }
  
  # 5. Point Estimates
  point <- fast_risk_fun(coefs)
  
  if (!conf_int) {
    res <- list(risks = as_tibble(point))
    if (!is.null(rmst_tau)) {
      if (max(horizons) < rmst_tau) stop("max(horizons) must be >= rmst_tau to compute RMST.")
      point_rmst <- point %>%
        filter(horizon <= rmst_tau) %>%
        group_by(trt) %>%
        arrange(horizon) %>%
        mutate(
          survival = 1 - risk,
          time_diff = horizon - dplyr::lag(horizon, default = 0),
          area = dplyr::lag(survival, default = 1) * time_diff
        ) %>%
        summarise(RMST = sum(area), .groups = "drop") %>%
        mutate(
          ref_RMST = RMST[trt == ref_group],
          RMST_diff = ifelse(trt == ref_group, 0, RMST - ref_RMST)
        ) %>% select(-ref_RMST)
      res$rmst <- point_rmst
    }
    return(res)
  }
  
  # 6. Sequential Simulation for Confidence Intervals
  sim_coefs <- mvtnorm::rmvnorm(n = samples, mean = coefs, sigma = V)
  sim_list <- lapply(seq_len(nrow(sim_coefs)), function(i) {
    fast_risk_fun(sim_coefs[i, ])
  })
  sim_df <- rbindlist(sim_list, idcol = "sim")
  
  # 7. Compute Final Statistics
  risk_ci <- sim_df %>%
    group_by(trt, horizon) %>%
    summarise(
      risk_lower = quantile(risk, 0.025),
      risk_upper = quantile(risk, 0.975),
      risk = mean(risk), 
      .groups = "drop"
    )
  
  sim_contrasts <- sim_df %>%
    group_by(sim, horizon) %>%
    mutate(ref_risk = risk[trt == ref_group]) %>%
    filter(trt != ref_group) %>%
    mutate(RD = risk - ref_risk, RR = risk / ref_risk)
  
  contrast_ci <- sim_contrasts %>%
    group_by(trt, horizon) %>%
    summarise(
      RD_lower = quantile(RD, 0.025),
      RD_upper = quantile(RD, 0.975),
      RD = mean(RD),
      RR_lower = quantile(RR, 0.025),
      RR_upper = quantile(RR, 0.975),
      RR = mean(RR),
      .groups = "drop"
    )
  
  res <- list(risks = risk_ci, contrasts = contrast_ci, ref_group = ref_group)
  
  # 8. Compute RMST distributions
  if (!is.null(rmst_tau)) {
    if (max(horizons) < rmst_tau) stop("max(horizons) must be >= rmst_tau to compute RMST.")
    sim_rmst <- sim_df %>%
      filter(horizon <= rmst_tau) %>%
      group_by(sim, trt) %>%
      arrange(horizon) %>%
      mutate(
        survival = 1 - risk,
        time_diff = horizon - dplyr::lag(horizon, default = 0),
        area = dplyr::lag(survival, default = 1) * time_diff
      ) %>%
      summarise(RMST = sum(area), .groups = "drop")
    
    rmst_ci <- sim_rmst %>%
      group_by(trt) %>%
      summarise(RMST_lower = quantile(RMST, 0.025), RMST_upper = quantile(RMST, 0.975), RMST = mean(RMST), .groups = "drop")
    
    rmst_contrasts <- sim_rmst %>%
      group_by(sim) %>%
      mutate(ref_RMST = RMST[trt == ref_group]) %>%
      filter(trt != ref_group) %>%
      mutate(RMST_diff = RMST - ref_RMST) %>%
      group_by(trt) %>%
      summarise(RMST_diff_lower = quantile(RMST_diff, 0.025), RMST_diff_upper = quantile(RMST_diff, 0.975), RMST_diff = mean(RMST_diff), .groups = "drop")
    
    res$rmst <- rmst_ci
    res$rmst_contrasts <- rmst_contrasts
  }
  return(res)
}

## 2. AUTOMATION WRAPPERS -------------------------------------------------------
run_ps_safely <- function(dat, ps_out_file, exposure_var = "trt", base_vars, cont_vars, out_vars, c_vars = c("studyid_rev", "episode", "time0")) {
  if (file.exists(ps_out_file)) {
    cat(sprintf("[SKIPPED] PS Model already exists: %s\n", basename(ps_out_file)))
    return(readRDS(ps_out_file))
  }
  cat(sprintf("[RUNNING] Fitting PS Model: %s\n", basename(ps_out_file)))
  ps_fit <- ps_model_wt(
    dat = dat, exposure = exposure_var, baselinevars = base_vars,
    continuous_vars = cont_vars, out_vars = out_vars, clust_vars = c_vars
  )
  saveRDS(ps_fit, ps_out_file)
  return(ps_fit)
}

run_outcome_safely <- function(ps_file, risk_out_file, time_var, event_var, covs, samples = 299) {
  if (file.exists(risk_out_file)) {
    cat(sprintf("  [SKIPPED] Outcome already exists: %s\n", basename(risk_out_file)))
    return(invisible())
  }
  
  cat(sprintf("  [RUNNING] Generating Outcome: %s\n", basename(risk_out_file)))
  ps_fit <- readRDS(ps_file)
  
  long_dat <- expand_person_period(
    ps_fit$data_weights, 
    time_var = time_var, 
    event_var = event_var, 
    interval_length = 14
  )
  
  res <- list()
  res$curve <- predict_cuminc_rmst_fast(
    Dat_expanded = long_dat, trt_var = "trt", covs = covs,
    horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
    samples = samples, conf_int = FALSE, ref_group = "no reversion", wt_type = "wt_trunc", rmst_tau = 365
  )
  
  res$ci <- predict_cuminc_rmst_fast(
    Dat_expanded = long_dat, trt_var = "trt", covs = covs,
    horizons = HORIZONS_CI, cluster_var = "studyid_rev",
    samples = samples, conf_int = TRUE, ref_group = "no reversion", wt_type = "wt_trunc", rmst_tau = 365
  )
  
  final_obj <- list()
  trial_key <- ifelse(grepl("tr2", risk_out_file), "tr2", "tr3")
  outcome_key <- ifelse(grepl("discont", risk_out_file), "dis", "od")
  
  final_obj[[trial_key]][[outcome_key]] <- res
  saveRDS(final_obj[[trial_key]], risk_out_file)
  
  rm(long_dat, ps_fit, res, final_obj); gc()
}

## 3. ANALYSIS 1: ALL-CAUSE ACUTE CARE ------------------------------------------
cat("\n======================================================\n")
cat("=== ANALYSIS 1: ALL-CAUSE ACUTE CARE               ===\n")
cat("======================================================\n")

ps_file_tr2_main <- paste0(RES_DIR, "PScore/ps_weight_tr2_34d_all.rds")
ps_file_tr3_main <- paste0(RES_DIR, "PScore/ps_weight_tr3_all.rds")

# Trial 2
run_outcome_safely(
  ps_file = ps_file_tr2_main,
  risk_out_file = paste0(RISK_DIR, "risk_tr2_34d_all_acute.rds"),
  time_var = "time_all_acute", event_var = "sens_all_acute",
  covs = c("ddose_prev", "prior_episode"), samples = 299
)

# Trial 3
run_outcome_safely(
  ps_file = ps_file_tr3_main,
  risk_out_file = paste0(RISK_DIR, "risk_tr3_514d_all_acute.rds"),
  time_var = "time_all_acute", event_var = "sens_all_acute",
  covs = c("ddose_prev", "year"), samples = 299
)


## 4. ANALYSIS 2: 30-DAY GAP DISCONTINUATION ------------------------------------
cat("\n======================================================\n")
cat("=== ANALYSIS 2: 30-DAY GAP DISCONTINUATION         ===\n")
cat("======================================================\n")

# File targets
ps_tr2_dis30_file <- paste0(SA_DIR, "ps_weight_tr2_34d_discont30.rds")
ps_tr3_dis30_file <- paste0(SA_DIR, "ps_weight_tr3_514d_discont30.rds")

# Check and build dedicated PS models if missing
if (!file.exists(ps_tr2_dis30_file) || !file.exists(ps_tr3_dis30_file)) {
  cat("Loading data to generate Discont-30 PS Models...\n")
  
  clean_trt_standard <- function(dat) {
    dat |>
      mutate(
        trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456),
        trt_s456 = factor(trt_s456, levels = TRT_LEVELS),
        prior_episode = as.numeric(episode - 1)
      ) |>
      rename(trt = trt_s456) |>
      filter(trt %in% TRT_LEVELS) |>
      group_by(studyid_rev, episode) |>
      mutate(prior_gap = cumsum(gap_days), prior_gap_num = row_number()) %>%
      ungroup()
  }
  
  raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))
  
  dat_tr2_dis30 <- raw_all |> filter(trial_gap %in% 3:4) |> clean_trt_standard() 
  dat_tr3_dis30 <- raw_all |> filter(trial_cat == "5-14") |> clean_trt_standard()
  rm(raw_all); gc()
  
  # Ensure the specific discont30 outcome vars are included
  out_vars_dis30 <- c("sens_discont30", "time_discont30", "sens_od_acute", "time_od_acute", "sens_all_acute", "time_all_acute")
  
  run_ps_safely(dat_tr2_dis30, ps_tr2_dis30_file, "trt", c(baselinevars, "prior_episode"), c(continuous_vars, "prior_episode"), out_vars_dis30)
  run_ps_safely(dat_tr3_dis30, ps_tr3_dis30_file, "trt", c(baselinevars, "prior_episode"), c(continuous_vars, "prior_episode"), out_vars_dis30)
  rm(dat_tr2_dis30, dat_tr3_dis30); gc()
}

# Run Outcomes
# Trial 2
run_outcome_safely(
  ps_file = ps_tr2_dis30_file,
  risk_out_file = paste0(RISK_DIR, "risk_tr2_34d_discont30.rds"),
  time_var = "time_discont30", event_var = "sens_discont30",
  covs = c("ddose_prev", "prior_episode"), samples = 299
)

# Trial 3 
run_outcome_safely(
  ps_file = ps_tr3_dis30_file,
  risk_out_file = paste0(RISK_DIR, "risk_tr3_514d_discont30.rds"),
  time_var = "time_discont30", event_var = "sens_discont30",
  covs = c("ddose_prev", "year"), samples = 299
)


## 5. ANALYSIS 3: 5-CATEGORY EXPOSURE -------------------------------------------
cat("\n======================================================\n")
cat("=== ANALYSIS 3: 5-CATEGORY DOSE EXPOSURE           ===\n")
cat("======================================================\n")

ps_tr2_5cat_file <- paste0(SA_DIR, "ps_weight_tr2_34d_5cat.rds")
ps_tr3_5cat_file <- paste0(SA_DIR, "ps_weight_tr3_514d_5cat.rds")

if (!file.exists(ps_tr2_5cat_file) || !file.exists(ps_tr3_5cat_file)) {
  cat("Loading data to generate 5-Category PS Models...\n")
  
  TRT_LEVELS_5 <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase <50%", "increase >=50%")
  
  raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))
  
  clean_trt_5cat <- function(dat) {
    dat |>
      mutate(
        trt_s456 = case_when(
          ddose_rev == 0 ~ "no reversion",
          ddose_rev > 0 & ddose_rev < 50  ~ "1-<50%",
          (!is.na(ddose) & ddose < ddose_prev & ddose <= 40) | ddose_rev >= 50 ~ ">=50%/<=40mg",
          ddose_rev < 0 & ddose_rev > -50  ~ "increase <50%",
          ddose_rev <= -50  ~ "increase >=50%", 
          TRUE ~ NA_character_
        ),
        trt_s456 = factor(trt_s456, levels = TRT_LEVELS_5),
        prior_episode = as.numeric(episode - 1)
      ) |>
      rename(trt = trt_s456) |>
      filter(trt %in% TRT_LEVELS_5) |>
      group_by(studyid_rev, episode) |>
      mutate(prior_gap = cumsum(gap_days), prior_gap_num = row_number()) %>%
      ungroup()
  }
  
  dat_tr2_5cat <- raw_all |> filter(trial_gap %in% 3:4) |> clean_trt_5cat() 
  dat_tr3_5cat <- raw_all |> filter(trial_cat == "5-14") |> clean_trt_5cat()
  rm(raw_all); gc()
  
  run_ps_safely(dat_tr2_5cat, ps_tr2_5cat_file, "trt", c(baselinevars, "prior_episode"), c(continuous_vars, "prior_episode"), default_out_vars)
  run_ps_safely(dat_tr3_5cat, ps_tr3_5cat_file, "trt", c(baselinevars, "prior_episode"), c(continuous_vars, "prior_episode"), default_out_vars)
  rm(dat_tr2_5cat, dat_tr3_5cat); gc()
}

# TR2 - Overdose
run_outcome_safely(
  ps_file = ps_tr2_5cat_file, risk_out_file = paste0(RISK_DIR, "risk_tr2_34d_5cat_od.rds"),
  time_var = "time_od_acute", event_var = "sens_od_acute", covs = c("ddose_prev", "prior_episode"), samples = 299
)
# TR2 - Discontinuation
run_outcome_safely(
  ps_file = ps_tr2_5cat_file, risk_out_file = paste0(RISK_DIR, "risk_tr2_34d_5cat_dis.rds"), 
  time_var = "time_discont", event_var = "sens_discont", covs = c("ddose_prev", "prior_episode"), samples = 299
)

# TR3 - Overdose
run_outcome_safely(
  ps_file = ps_tr3_5cat_file, risk_out_file = paste0(RISK_DIR, "risk_tr3_514d_5cat_od.rds"),
  time_var = "time_od_acute", event_var = "sens_od_acute", covs = c("ddose_prev", "year"), samples = 299
)
# TR3 - Discontinuation
run_outcome_safely(
  ps_file = ps_tr3_5cat_file, risk_out_file = paste0(RISK_DIR, "risk_tr3_514d_5cat_dis.rds"),
  time_var = "time_discont", event_var = "sens_discont", covs = c("ddose_prev", "year"), samples = 299
)

cat("\n=== PIPELINE COMPLETE ===\n")





## 6. ANALYSIS 4: 4 DAYS GAP ------------------------------------
cat("=== ANALYSIS 4: 4 DAYS GAP ===\n")
# ── Shared settings ────────────────────────────────────────────────────────────
TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")

out_vars   <- c("sens_discont", "time_discont",
                "sens_od_acute", "time_od_acute",
                "sens_all_acute", "time_all_acute")

clust_vars <- c("studyid_rev", "episode", "time0")

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss"
)
continuous_vars <- c("ddose_prev", "year", "ret_length")

clean_trt_all <- function(dat) {
  dat |>
    mutate(
      trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)
    ) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS),
           prior_episode = as.numeric(episode - 1)) |>
    group_by(studyid_rev, episode) |>
    mutate(
      prior_gap     = cumsum(gap_days),
      prior_gap_num = row_number()
    ) %>%
    ungroup()
  
}



# First-interruption-in-any-episode 
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 4) |>
  clean_trt_all() 

ps_tr2_all <- ps_model_wt(
  dat             = dat_tr2_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_episode"),
  continuous_vars = c(continuous_vars, "prior_episode"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_all, paste0(RES_DIR, "Sensitivities/", "ps_weight_tr2_4d_all.rds"))


ps_fit <- readRDS(paste0(RES_DIR, "Sensitivities/", "ps_weight_tr2_4d_all.rds"))



HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(unique(c(HORIZONS_CURVE, c(30, 60, 90, 180))))


results <- list(od = list(), dis = list())

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

results$od$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 299, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

results$od$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 299, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)

results$dis$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"),  
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 299, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

results$dis$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 299, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(results, paste0(RES_DIR, "Sensitivities/Risk/", "risk_tr2_4d_rmst.rds"))



# DISCONT 30 DAYS ########


library(tidyverse)
library(haven)
library(survival)
library(sandwich)
library(data.table)

## 1. SETUP & DIRECTORIES -------------------------------------------------------
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"
SA_DIR    <- paste0(RES_DIR, "Sensitivities/")
RISK_DIR  <- paste0(SA_DIR, "Risk/")

dir.create(SA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RISK_DIR, showWarnings = FALSE, recursive = TRUE)

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(unique(c(HORIZONS_CURVE, c(30, 60, 90, 180))))

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
clust_vars <- c("studyid_rev", "episode", "time0")

baselinevars <- c("sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
                  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
                  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
                  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

continuous_vars <- c("ddose_prev", "year", "ret_length")

## 2. AUTOMATION WRAPPERS (Loaded for context) ----------------------------------
run_ps_safely <- function(dat, ps_out_file, exposure_var = "trt", base_vars, cont_vars, out_vars, c_vars = c("studyid_rev", "episode", "time0")) {
  if (file.exists(ps_out_file)) {
    cat(sprintf("[SKIPPED] PS Model already exists: %s\n", basename(ps_out_file)))
    return(readRDS(ps_out_file))
  }
  cat(sprintf("[RUNNING] Fitting PS Model: %s\n", basename(ps_out_file)))
  ps_fit <- ps_model_wt(dat = dat, exposure = exposure_var, baselinevars = base_vars, continuous_vars = cont_vars, out_vars = out_vars, clust_vars = c_vars)
  saveRDS(ps_fit, ps_out_file)
  return(ps_fit)
}


run_outcome_safely <- function(ps_file, risk_out_file, time_var, event_var, covs, samples = 299) {
  if (file.exists(risk_out_file)) {
    cat(sprintf("  [SKIPPED] Outcome already exists: %s\n", basename(risk_out_file)))
    return(invisible())
  }
  cat(sprintf("  [RUNNING] Generating Outcome: %s\n", basename(risk_out_file)))
  ps_fit <- readRDS(ps_file)
  
  long_dat <- expand_person_period(ps_fit$data_weights, time_var = time_var, event_var = event_var, interval_length = 14)
  res <- list()
  
  res$curve <- predict_cuminc_rmst_fast(Dat_expanded = long_dat, trt_var = "trt", covs = covs, horizons = HORIZONS_CURVE, cluster_var = "studyid_rev", samples = samples, conf_int = FALSE, ref_group = "no reversion", wt_type = "wt_trunc", rmst_tau = 365)
  res$ci <- predict_cuminc_rmst_fast(Dat_expanded = long_dat, trt_var = "trt", covs = covs, horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = samples, conf_int = TRUE, ref_group = "no reversion", wt_type = "wt_trunc", rmst_tau = 365)
  
  final_obj <- list()
  trial_key <- ifelse(grepl("tr2", risk_out_file), "tr2", "tr3")
  final_obj[[trial_key]][["dis"]] <- res # Hardcoded to "dis" since this is exclusively a discont analysis
  saveRDS(final_obj[[trial_key]], risk_out_file)
  rm(long_dat, ps_fit, res, final_obj); gc()
}

## 3. ANALYSIS: 30-DAY GAP DISCONTINUATION --------------------------------------
cat("\n======================================================\n")
cat("=== ANALYSIS: 30-DAY GAP DISCONTINUATION           ===\n")
cat("======================================================\n")

# File targets
ps_tr2_dis30_file <- paste0(SA_DIR, "ps_weight_tr2_34d_discont30.rds")
ps_tr3_dis30_file <- paste0(SA_DIR, "ps_weight_tr3_514d_discont30.rds")

# We only load and clean data if the PS models are missing
if (!file.exists(ps_tr2_dis30_file) || !file.exists(ps_tr3_dis30_file)) {
  cat("Loading data to generate Discont-30 PS Models...\n")
  
  # Standard 4-category cleaning
  clean_trt_standard <- function(dat) {
    dat |>
      mutate(
        trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456),
        trt_s456 = factor(trt_s456, levels = TRT_LEVELS),
        prior_episode = as.numeric(episode - 1)
      ) |>
      rename(trt = trt_s456) |>
      filter(trt %in% TRT_LEVELS) |>
      group_by(studyid_rev, episode) |>
      mutate(prior_gap = cumsum(gap_days), prior_gap_num = row_number()) %>%
      ungroup()
  }
  
  raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))
  
  dat_tr2_dis30 <- raw_all |> filter(trial_gap %in% 3:4) |> clean_trt_standard() 
  dat_tr3_dis30 <- raw_all |> filter(trial_cat == "5-14") |> clean_trt_standard()
  rm(raw_all); gc()
  
  # Crucial: Define the out_vars specifically including time_discont30
  out_vars_dis30 <- c("sens_discont30", "time_discont30", "sens_od_acute", "time_od_acute", "sens_all_acute", "time_all_acute")
  
  run_ps_safely(dat_tr2_dis30, ps_tr2_dis30_file, "trt", c(baselinevars, "prior_episode"), c(continuous_vars, "prior_episode"), out_vars_dis30)
  run_ps_safely(dat_tr3_dis30, ps_tr3_dis30_file, "trt", c(baselinevars, "prior_episode"), c(continuous_vars, "prior_episode"), out_vars_dis30)
  rm(dat_tr2_dis30, dat_tr3_dis30); gc()
}

# Run Outcomes using the dedicated PS models
cat("\nRunning Outcomes...\n")

# Trial 2
run_outcome_safely(
  ps_file = ps_tr2_dis30_file,
  risk_out_file = paste0(RISK_DIR, "risk_tr2_34d_discont30.rds"),
  time_var = "time_discont30", event_var = "sens_discont30",
  covs = c("ddose_prev", "prior_episode"), samples = 299
)

# Trial 3 
run_outcome_safely(
  ps_file = ps_tr3_dis30_file,
  risk_out_file = paste0(RISK_DIR, "risk_tr3_514d_discont30.rds"),
  time_var = "time_discont30", event_var = "sens_discont30",
  covs = c("ddose_prev", "year"), samples = 299
)

cat("\n=== 30-DAY DISCONTINUATION PIPELINE COMPLETE ===\n")