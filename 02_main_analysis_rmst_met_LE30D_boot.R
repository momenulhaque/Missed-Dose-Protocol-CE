library(tidyverse)
library(haven)
library(nnet)
library(splines)
library(sandwich)
library(mvtnorm)
library(tableone)
library(survey)
library(cobalt)
library(patchwork)
library(cowplot)

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")   # ps_model_wt(), predict_cuminc(), expand_person_period()


# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/LE30D/"
dir.create(paste0(RES_DIR, "PScore"),   showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Risk"),     showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Tables"),   showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Figures"),  showWarnings = FALSE)


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

# Human-readable variable labels for Table 1 and love plots
VAR_LABELS <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
  home_grp_b    = "Unstable housing in past 5 years",
  incar_grp_bl  = "Incarcerated in the past year",
  drug_acute    = "Drug-related acute care visits in the past month",
  cci_gt0       = "Charlson Comorbidity Index (>0)",
  copd_asthma   = "Asthma or COPD in the past year",
  sed_psy_med   = "Dispensations of psychiatric or sedative medication in the past month",
  po            = "Dispensations opioid other than OAT in the past month",
  sud           = "Other substance use disorder",
  aud           = "Alcohol use disorder",
  smd           = "Serious mental disorder",
  hcv           = "Hepatitis C",
  tob           = "Tobacco use disorder (past year)",
  udt           = "Urine drug test in the past month",
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_gap     = "Total OAT interruptions since initiation",
  prior_episode = "Number of prior episode"
)

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

# 1. LOAD AND FILTER RAW SAS DATA #########

raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat")) %>% 
  mutate(t_death = rev_dod - time0,
         s_death_od = ifelse(!is.na(t_death) & t_death == time_od_acute, 1, 0),
         s_death_dis = ifelse(!is.na(t_death) & t_death == time_discont, 1, 0))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all() |>
  filter(ret_length <= 30)

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all() |>
  filter(ret_length <= 30)



bootstrap_ipw_surv_pipeline <- function(
    dat_baseline, 
    dat_weighted,         # NEW: Pass pre-computed weights for exact point estimates
    ps_args,              
    exp_time_var,         
    exp_event_var,        
    exp_interval,         
    trt_var,              
    out_covs,             
    horizons,             
    ci_horizons = 30,     
    cluster_var = "studyid_rev",
    wt_type = "wt_trunc",
    samples = 500,
    ref_group = "no reversion"
) {
  
  cat("\n--- 1. Computing Point Estimates (Using Provided Weighted Data) ---\n")
  
  # Use the exact weights from the original PS fit to lock in the point estimates
  dat_long_orig <- expand_person_period(
    dat_weighted, 
    time_var = exp_time_var, event_var = exp_event_var, interval_length = exp_interval
  )
  
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", out_covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  model_orig <- glm(plr_formula, data = dat_long_orig, family = quasibinomial(), weights = dat_long_orig[[wt_type]])
  
  is_trt_factor <- is.factor(dat_long_orig[[trt_var]])
  trt_levels <- if (is_trt_factor) levels(dat_long_orig[[trt_var]]) else as.character(unique(dat_long_orig[[trt_var]]))
  
  dt_orig <- as.data.table(dat_long_orig)
  setorderv(dt_orig, c(cluster_var, "stop"))
  
  mm_list_orig <- lapply(trt_levels, function(lvl) {
    dt_temp <- copy(dt_orig)
    dt_temp[[trt_var]] <- if (is_trt_factor) factor(lvl, levels = trt_levels) else lvl
    model.matrix(delete.response(terms(model_orig)), data = dt_temp, xlev = model_orig$xlevels)
  })
  names(mm_list_orig) <- trt_levels
  
  # Helper to calculate risks over a given set of horizons
  calc_risk <- function(dt, mm_list, coefs, id_var, calc_horizons) {
    res_list <- vector("list", length(trt_levels))
    for (i in seq_along(trt_levels)) {
      lvl <- trt_levels[i]
      dt[, p_event := plogis(mm_list[[lvl]] %*% coefs)]
      dt[, cum_surv := cumprod(1 - p_event), by = id_var]
      
      risks_h <- numeric(length(calc_horizons))
      for (j in seq_along(calc_horizons)) {
        h <- calc_horizons[j]
        risks_h[j] <- 1 - dt[stop <= h, .(last_surv = cum_surv[.N]), by = id_var][, mean(last_surv)]
      }
      res_list[[i]] <- data.table(trt = lvl, horizon = calc_horizons, risk = risks_h)
    }
    rbindlist(res_list)
  }
  
  # Calculate point estimates for ALL requested horizons
  point_est <- calc_risk(dt_orig, mm_list_orig, coef(model_orig), cluster_var, horizons)
  
  point_contrasts <- point_est %>%
    group_by(horizon) %>%
    mutate(ref_risk = risk[trt == ref_group]) %>%
    filter(trt != ref_group) %>%
    mutate(RD = risk - ref_risk, RR = risk / ref_risk) %>%
    select(trt, horizon, RD, RR) %>%
    ungroup()
  
  if (samples == 0) return(list(risks = as_tibble(point_est), contrasts = point_contrasts, ref_group = ref_group)) 
  
  cat("\n--- 2. Non-Parametric Bootstrap (Including PS Weighting) ---\n")
  
  dat_long_raw <- expand_person_period(
    dat_baseline, 
    time_var = exp_time_var, event_var = exp_event_var, interval_length = exp_interval
  )
  
  clusters <- unique(dat_baseline[[cluster_var]])
  n_clusters <- length(clusters)
  sim_list <- vector("list", samples)
  
  pb <- txtProgressBar(min = 0, max = samples, style = 3)
  for (i in seq_len(samples)) {
    boot_model <- tryCatch({
      boot_clusters <- sample(clusters, size = n_clusters, replace = TRUE)
      boot_idx <- data.table(orig_id = boot_clusters, boot_id = seq_len(n_clusters))
      
      boot_base <- merge(boot_idx, as.data.table(dat_baseline), by.x = "orig_id", by.y = cluster_var, allow.cartesian = TRUE)
      boot_base[[cluster_var]] <- boot_base$boot_id 
      
      boot_ps_res <- do.call(ps_model_wt, c(list(dat = as.data.frame(boot_base)), ps_args))
      boot_weights <- as.data.table(boot_ps_res$data_weights)[, c(cluster_var, "episode", wt_type), with = FALSE]
      
      boot_long <- merge(boot_idx, as.data.table(dat_long_raw), by.x = "orig_id", by.y = cluster_var, allow.cartesian = TRUE)
      boot_long[[cluster_var]] <- boot_long$boot_id
      boot_long <- merge(boot_long, boot_weights, by = c(cluster_var, "episode"))
      setorderv(boot_long, c(cluster_var, "stop"))
      
      suppressWarnings(glm(plr_formula, data = boot_long, family = quasibinomial(), weights = boot_long[[wt_type]]))
    }, error = function(e) NULL)
    
    if (!is.null(boot_model) && boot_model$converged) {
      mm_list_boot <- lapply(trt_levels, function(lvl) {
        dt_temp <- copy(boot_long)
        dt_temp[[trt_var]] <- if (is_trt_factor) factor(lvl, levels = trt_levels) else lvl
        model.matrix(delete.response(terms(boot_model)), data = dt_temp, xlev = boot_model$xlevels)
      })
      names(mm_list_boot) <- trt_levels
      
      sim_list[[i]] <- calc_risk(boot_long, mm_list_boot, coef(boot_model), cluster_var, ci_horizons)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  sim_list <- sim_list[!sapply(sim_list, is.null)]
  sim_df <- rbindlist(sim_list, idcol = "sim")
  
  cat(sprintf("\nSuccessfully completed %d/%d iterations.\n", length(sim_list), samples))
  cat("\n--- 3. Structuring Final Output ---\n")
  
  risk_ci <- sim_df %>%
    group_by(trt, horizon) %>%
    summarise(
      risk_lower = quantile(risk, 0.025, na.rm = TRUE),
      risk_upper = quantile(risk, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) 
  
  risk_out <- point_est %>% 
    left_join(risk_ci, by = c("trt", "horizon"))
  
  contrast_ci <- sim_df %>%
    group_by(sim, horizon) %>%
    mutate(ref_risk = risk[trt == ref_group]) %>%
    filter(trt != ref_group) %>%
    mutate(RD = risk - ref_risk, RR = risk / ref_risk) %>%
    group_by(trt, horizon) %>%
    summarise(
      RD_lower = quantile(RD, 0.025, na.rm = TRUE),
      RD_upper = quantile(RD, 0.975, na.rm = TRUE),
      RR_lower = quantile(RR, 0.025, na.rm = TRUE),
      RR_upper = quantile(RR, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  contrast_out <- point_contrasts %>% 
    left_join(contrast_ci, by = c("trt", "horizon"))
  
  return(list(risks = risk_out, contrasts = contrast_out, ref_group = ref_group))
}




HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))

# ── INITIALIZE RESULTS STRUCTURE ──

# ── INITIALIZE RESULTS STRUCTURE ──

risk <- list(
  tr2 = list(od = list(), dis = list()),
  tr3 = list(od = list(), dis = list())
)

# ── Trial 2 ──
cat("\nComputing Trial 2 Risks...\n")
ps_fit_tr2 <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))

risk$tr2$od$ci <- bootstrap_ipw_surv_pipeline(
  dat_baseline = dat_tr2_all, 
  dat_weighted = ps_fit_tr2$data_weights,  # <-- Locks in the exact point estimates
  
  ps_args = list(
    exposure        = "trt",
    baselinevars    =  c(baselinevars, "prior_episode"),
    continuous_vars = c(continuous_vars, "prior_episode"),
    out_vars        = out_vars,
    clust_vars      = clust_vars 
  ),
  
  exp_time_var    = "time_od_acute",
  exp_event_var   = "sens_od_acute",
  exp_interval    = 14,
  
  trt_var         = "trt",
  out_covs        = c("year", "ddose_prev", "prior_episode"),
  
  horizons        = HORIZONS_CI, 
  ci_horizons     = HORIZONS_CI,             
  cluster_var     = "studyid_rev",
  samples         = 999,         
  ref_group       = "no reversion"
)

risk$tr2$dis$ci <- bootstrap_ipw_surv_pipeline(
  dat_baseline = dat_tr2_all, 
  dat_weighted = ps_fit_tr2$data_weights,  # <-- Locks in the exact point estimates
  ps_args = list(
    exposure        = "trt",
    baselinevars    =  c(baselinevars, "prior_episode"),
    continuous_vars = c(continuous_vars, "prior_episode"),
    out_vars        = out_vars,
    clust_vars      = clust_vars 
  ),
  exp_time_var    = "time_discont",
  exp_event_var   = "sens_discont",
  exp_interval    = 14,
  trt_var         = "trt",
  out_covs        = c("year", "ddose_prev", "prior_episode"),
  horizons        = HORIZONS_CI, 
  ci_horizons     = HORIZONS_CI,
  cluster_var     = "studyid_rev",
  samples         = 999,
  ref_group       = "no reversion"
)

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_boot.rds"))

# ── Trial 3 ──
cat("\nComputing Trial 3 Risks...\n")
ps_fit_tr3 <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

risk$tr3$od$ci <- bootstrap_ipw_surv_pipeline(
  dat_baseline = dat_tr3_all, 
  dat_weighted = ps_fit_tr3$data_weights,  # <-- Locks in the exact point estimates
  ps_args = list(
    exposure        = "trt",
    baselinevars    =  c(baselinevars, "prior_episode"),
    continuous_vars = c(continuous_vars, "prior_episode"),
    out_vars        = out_vars,
    clust_vars      = clust_vars 
  ),
  exp_time_var    = "time_od_acute",
  exp_event_var   = "sens_od_acute",
  exp_interval    = 14,
  trt_var         = "trt",
  out_covs        = c("year", "ddose_prev", "prior_episode"),
  horizons        = HORIZONS_CI, 
  ci_horizons     = HORIZONS_CI,
  cluster_var     = "studyid_rev",
  samples         = 999,
  ref_group       = "no reversion"
)

risk$tr3$dis$ci <- bootstrap_ipw_surv_pipeline(
  dat_baseline = dat_tr3_all, 
  dat_weighted = ps_fit_tr3$data_weights,  # <-- Locks in the exact point estimates
  ps_args = list(
    exposure        = "trt",
    baselinevars    =  c(baselinevars, "prior_episode"),
    continuous_vars = c(continuous_vars, "prior_episode"),
    out_vars        = out_vars,
    clust_vars      = clust_vars 
  ),
  exp_time_var    = "time_discont",
  exp_event_var   = "sens_discont",
  exp_interval    = 14,
  trt_var         = "trt",
  out_covs        = c("year", "ddose_prev", "prior_episode"),
  horizons        = HORIZONS_CI, 
  ci_horizons     = HORIZONS_CI,
  cluster_var     = "studyid_rev",
  samples         = 999,
  ref_group       = "no reversion"
)

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_boot.rds"))

