library(tidyverse)
library(haven)
library(nnet)
library(splines)
library(sandwich)
library(mvtnorm)
library(survival)

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")   # ps_model_wt(), predict_cuminc(), expand_person_period()

# ── Paths & Settings ──────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/BNX/Sensitivities/"
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
clust_vars <- c("studyid_rev", "episode", "time0")


clean_trt_all <- function(dat) {
  dat |>
    mutate(
      trt_s5 = if_else(trt_s5 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s5)
    ) |>
    rename(trt = trt_s5) |>
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



default_out_vars <- c("sens_discont", "time_discont", 
                      "sens_od_acute", "time_od_acute", 
                      "sens_all_acute", "time_all_acute")

# Add the new Per-Protocol variables to our out_vars list so they aren't dropped
pp_out_vars <- c(
  default_out_vars, 
  "time_od_acute_pp", "sens_od_acute_pp", "is_censored_od",
  "time_discont_pp", "sens_discont_pp", "is_censored_dis"
)

baselinevars <- c("sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
                  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
                  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
                  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

continuous_vars <- c("ddose_prev", "year", "ret_length")


# ── 1. DATA PREPARATION ───────────────────────────────────────────────────────
raw_first <- read_sas(file.path(DATA_DIR, "tt_analytic_rprev_sens_nw_b.sas7bdat")) |>
  arrange(studyid_rev, episode, time0) |>
  mutate(
    trt_s5 = if_else(trt_s5 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s5)
  ) |>
  rename(trt = trt_s5) |>
  filter(trt %in% TRT_LEVELS) |>
  mutate(trt = factor(trt, levels = TRT_LEVELS)) |>
  group_by(studyid_rev, episode) |>
  mutate(
    prior_gap = cumsum(gap_days) # Add prior interruptions per reviewer request
  ) |> 
  ungroup()

# Step 1A: Isolate ALL gaps >= 3 days to serve as our registry of "next" gaps for censoring
all_gaps <- raw_first %>% 
  filter(trial_gap >= 3) %>% 
  select(studyid_rev, episode, time0, trial_gap) %>% 
  rename(next_time0 = time0, next_trial_gap = trial_gap)


# First-interruption-in-any-episode 
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_b.sas7bdat"))

dat_tr4_sa2 <- raw_all |>
  filter(trial_gap %in% 3:5) |>
  clean_trt_all() 

dat_tr5_sa2 <- raw_all |>
  filter(trial_cat == "6-14") |>
  clean_trt_all()

# ── 2. APPLY PER-PROTOCOL CENSORING ───────────────────────────────────────────
apply_pp_censoring <- function(dat_index) {
  
  # Find the exact next gap >= 3 days for each index event
  censor_dates <- dat_index %>%
    left_join(all_gaps, by = c("studyid_rev", "episode"), relationship = "many-to-many") %>%
    filter(next_time0 > time0) %>%          # Must happen AFTER the index reinitiation
    group_by(studyid_rev, episode, time0) %>%
    arrange(next_time0) %>%
    slice(1) %>%                            # Grab only the very next gap
    ungroup() %>%
    select(studyid_rev, episode, time0, next_time0, next_trial_gap)
  
  # Merge back to the trial data and create Per-Protocol outcomes
  dat_pp <- dat_index %>%
    left_join(censor_dates, by = c("studyid_rev", "episode", "time0")) %>%
    mutate(
      # The gap started 'next_trial_gap' days prior to the next reinitiation
      censor_date = next_time0 - next_trial_gap,
      time_to_censor = as.numeric(censor_date - time0),
      
      # If no subsequent gap exists, set censoring to max follow-up (365 days)
      time_to_censor = replace_na(time_to_censor, 365),
      
      # --- OVERDOSE ---
      time_od_acute_pp = pmin(time_od_acute, time_to_censor),
      sens_od_acute_pp = ifelse(time_od_acute <= time_to_censor, sens_od_acute, 0),
      is_censored_od   = ifelse(time_to_censor < time_od_acute & time_to_censor < 365, 1, 0),
      
      # --- DISCONTINUATION ---
      time_discont_pp = pmin(time_discont, time_to_censor),
      sens_discont_pp = ifelse(time_discont <= time_to_censor, sens_discont, 0),
      is_censored_dis = ifelse(time_to_censor < time_discont & time_to_censor < 365, 1, 0)
    )
  
  return(dat_pp)
}

# Apply to both trials
dat_tr4_pp <- apply_pp_censoring(dat_tr4_sa2)
dat_tr5_pp <- apply_pp_censoring(dat_tr5_sa2)
dat_tr5_pp %>% count(trt)

# ── 3. IPCW FUNCTION (Fixed Matrix Bug) ───────────────────────────────────────
add_ipcw <- function(dat, time_var, censor_var, covariates) {
  # Denominator model: Probability of remaining uncensored given TRT + Covariates
  f_den <- as.formula(paste("Surv(", time_var, ",", censor_var, ") ~ trt +", paste(covariates, collapse = " + ")))
  fit_den <- coxph(f_den, data = dat)
  surv_den <- exp(-predict(fit_den, type = "expected")) # Direct vector calculation
  
  # Numerator model: Probability of remaining uncensored given TRT only (stabilizer)
  f_num <- as.formula(paste("Surv(", time_var, ",", censor_var, ") ~ trt"))
  fit_num <- coxph(f_num, data = dat)
  surv_num <- exp(-predict(fit_num, type = "expected"))
  
  # Protect against extreme weights
  surv_den <- pmax(surv_den, 0.01)
  
  return(surv_num / surv_den)
}


# BNX 3-5 days ###############

ps_tr4_pp <- ps_model_wt(
  dat = dat_tr4_pp, exposure = "trt", out_vars = pp_out_vars, 
  clust_vars = clust_vars, baselinevars = c(baselinevars, "prior_episode"), 
  continuous_vars = c(continuous_vars, "prior_episode")
)
dat_tr4_weighted <- ps_tr4_pp$data_weights

# B. Calculate IPCW and Combined Weights for Overdose
dat_tr4_weighted$ipcw_od <- add_ipcw(
  dat = dat_tr4_weighted, time_var = "time_od_acute_pp", 
  censor_var = "is_censored_od", covariates = c(baselinevars, "prior_episode")
)

dat_tr4_weighted$wt_pp_od_raw <- dat_tr4_weighted$wt_stab * dat_tr4_weighted$ipcw_od
t_low_od <- quantile(dat_tr4_weighted$wt_pp_od_raw, 0.01, na.rm = TRUE)
t_high_od <- quantile(dat_tr4_weighted$wt_pp_od_raw, 0.99, na.rm = TRUE)
dat_tr4_weighted$wt_pp_od <- pmin(pmax(dat_tr4_weighted$wt_pp_od_raw, t_low_od), t_high_od)

# C. Calculate IPCW and Combined Weights for Discontinuation
dat_tr4_weighted$ipcw_dis <- add_ipcw(
  dat = dat_tr4_weighted, time_var = "time_discont_pp", 
  censor_var = "is_censored_dis", covariates = c(baselinevars, "prior_episode")
)

dat_tr4_weighted$wt_pp_dis_raw <- dat_tr4_weighted$wt_stab * dat_tr4_weighted$ipcw_dis
t_low_dis <- quantile(dat_tr4_weighted$wt_pp_dis_raw, 0.01, na.rm = TRUE)
t_high_dis <- quantile(dat_tr4_weighted$wt_pp_dis_raw, 0.99, na.rm = TRUE)
dat_tr4_weighted$wt_pp_dis <- pmin(pmax(dat_tr4_weighted$wt_pp_dis_raw, t_low_dis), t_high_dis)


# D. Outcome Models
HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))
risk_pp <- list(tr4 = list(od = list(), dis = list()), tr5 = list(od = list(), dis = list()))

# Overdose (Using prior_gap instead of prior_episode)
long_dat_tr4_od <- expand_person_period(dat_tr4_weighted, time_var = "time_od_acute_pp", event_var = "sens_od_acute_pp", interval_length = 14)

risk_pp$tr4$od$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat_tr4_od, trt_var = "trt", covs = "year", 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE, 
  ref_group = "no reversion", wt_type = "wt_pp_od"
)

# Discontinuation (Using prior_gap)
long_dat_tr4_dis <- expand_person_period(dat_tr4_weighted, time_var = "time_discont_pp", event_var = "sens_discont_pp", interval_length = 14)
risk_pp$tr4$dis$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat_tr4_dis, trt_var = "trt", covs = "year", 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE, 
  ref_group = "no reversion", wt_type = "wt_pp_dis"
)

saveRDS(risk_pp$tr4, paste0(RES_DIR, "risk/risk_tr4_per_protocol.rds"))


# BNX 6-14 days ###############

# A. Calculate Baseline IPTW
ps_tr5_pp <- ps_model_wt(
  dat = dat_tr5_pp, exposure = "trt", out_vars = pp_out_vars, 
  clust_vars = clust_vars, baselinevars = c(baselinevars, "prior_episode"), 
  continuous_vars = c(continuous_vars, "prior_episode")
)
dat_tr5_weighted <- ps_tr5_pp$data_weights

# B. Calculate IPCW and Combined Weights for Overdose
dat_tr5_weighted$ipcw_od <- add_ipcw(
  dat = dat_tr5_weighted, time_var = "time_od_acute_pp", 
  censor_var = "is_censored_od", covariates = c(baselinevars, "prior_episode")
)
dat_tr5_weighted$wt_pp_od_raw <- dat_tr5_weighted$wt_stab * dat_tr5_weighted$ipcw_od
t_low_od_3 <- quantile(dat_tr5_weighted$wt_pp_od_raw, 0.01, na.rm = TRUE)
t_high_od_3 <- quantile(dat_tr5_weighted$wt_pp_od_raw, 0.99, na.rm = TRUE)
dat_tr5_weighted$wt_pp_od <- pmin(pmax(dat_tr5_weighted$wt_pp_od_raw, t_low_od_3), t_high_od_3)

# C. Calculate IPCW and Combined Weights for Discontinuation
dat_tr5_weighted$ipcw_dis <- add_ipcw(
  dat = dat_tr5_weighted, time_var = "time_discont_pp", 
  censor_var = "is_censored_dis", covariates = c(baselinevars, "prior_episode")
)
dat_tr5_weighted$wt_pp_dis_raw <- dat_tr5_weighted$wt_stab * dat_tr5_weighted$ipcw_dis
t_low_dis_3 <- quantile(dat_tr5_weighted$wt_pp_dis_raw, 0.01, na.rm = TRUE)
t_high_dis_3 <- quantile(dat_tr5_weighted$wt_pp_dis_raw, 0.99, na.rm = TRUE)
dat_tr5_weighted$wt_pp_dis <- pmin(pmax(dat_tr5_weighted$wt_pp_dis_raw, t_low_dis_3), t_high_dis_3)

# D. Outcome Models
# Overdose (Using prior_gap)
long_dat_tr5_od <- expand_person_period(dat_tr5_weighted, time_var = "time_od_acute_pp", event_var = "sens_od_acute_pp", interval_length = 14)
risk_pp$tr5$od$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat_tr5_od, trt_var = "trt", covs = NULL, 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE, 
  ref_group = "no reversion", wt_type = "wt_pp_od"
)

# Discontinuation (Using prior_gap)
long_dat_tr5_dis <- expand_person_period(dat_tr5_weighted, time_var = "time_discont_pp", event_var = "sens_discont_pp", interval_length = 14)
risk_pp$tr5$dis$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat_tr5_dis, trt_var = "trt", covs = NULL, 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE, 
  ref_group = "no reversion", wt_type = "wt_pp_dis"
)

saveRDS(risk_pp$tr5, paste0(RES_DIR, "risk/risk_tr5_per_protocol.rds"))

cat("\n✓ Per-Protocol Analysis Complete. Results saved.\n")




res = readRDS(paste0(RES_DIR, "risk/risk_tr5_per_protocol.rds"))


res$od$ci

