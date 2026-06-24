# Load Packages ##################

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
library(autoCovariateSelection)

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")   # ps_model_wt(), predict_cuminc(), expand_person_period()

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/BNX/"
dir.create(paste0(RES_DIR, "PScore"),   showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Risk"),     showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Tables"),   showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Figures"),  showWarnings = FALSE)


# ── Shared settings ────────────────────────────────────────────────────────────
TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")

out_vars   <- c("sens_discont", "time_discont",
                "sens_od_acute", "time_od_acute",
                "sens_all_acute", "time_all_acute")

clust_vars <- c("studyid_rev", "episode", "time0")

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss", "prior_episode"
)
continuous_vars <- c("ddose_prev", "year", "ret_length", "prior_episode")

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
  prior_episode = "Total prior episode"
)

# Shared pre-processing: collapse all reduction categories into one

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
      prior_gap_num = row_number(),
      unique_id = paste0(studyid_rev, "_", time0)
    ) %>%
    ungroup()
  
}

# 1. LOAD AND FILTER RAW SAS DATA #########
# First-interruption-in-any-episode 
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_b.sas7bdat"))

dat_tr4_all <- raw_all |>
  filter(trial_gap %in% 3:5) |>
  clean_trt_all() 

dat_tr5_all <- raw_all |>
  filter(trial_cat == "6-14") |>
  clean_trt_all()


rm(raw_all); gc()



# hdPS analysis  ###########################
## PS fitting ########
### bnx 3-5 days #########################

dat3_all <- get_dat_all_prev_b(dat1=dat_tr4_all)
dim(dat3_all)
# empirical covariates
vars.empirical <- setdiff(names(dat3_all), names(dat_tr4_all))
length(vars.empirical)
# 84
dat3_all <- dat3_all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical_le200 <- select_top200_empirical_covariates(dat.all=dat3_all,
                                                           time_var = "time_discont",
                                                           event_var = "sens_discont",
                                                           max_vars = 200, nfolds = 5,
                                                           vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat3_all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical_le200),
                        continuous_vars = c("age_st", "ddose_prev", "year", 
                                            "ret_length", "prior_episode"))



saveRDS(ps_fit, paste0(RES_DIR, "PScore/", "ps_weight_tr4_hdPS.rds"))


### met 5-14 days #########################

dat5_all <- get_dat_all_prev_b(dat1=dat_tr5_all)
dim(dat5_all)
# empirical covariates
vars.empirical <- setdiff(names(dat5_all), names(dat_tr5_all))
length(vars.empirical)
# 84
dat5_all <- dat5_all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical_le200 <- select_top200_empirical_covariates(dat.all=dat5_all,
                                                           time_var = "time_discont",
                                                           event_var = "sens_discont",
                                                           max_vars = 200, nfolds = 5,
                                                           vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat5_all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical_le200),
                        continuous_vars = c("age_st", "ddose_prev", "year", 
                                            "ret_length", "prior_episode"))



saveRDS(ps_fit, paste0(RES_DIR, "PScore/", "ps_weight_tr5_hdPS.rds"))


## Model fitting ########
risk <- list(
  tr4 = list(
    od = list(), dis = list()
  ),
  tr5 = list(
    od = list(), dis = list()
  )
)

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))


### met 3-4 days #########################
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr4_hdPS.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr4$od <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = "year", 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)


risk$tr4$dis <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = "year",
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr4, paste0(RES_DIR, "Risk/", "risk_tr4_hdPS.rds"))

cat("  ✓ Trial 2 risk estimates saved.\n")


### met 5-14 days #########################

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr5_hdPS.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr5$od <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CI, 
  cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr5$dis <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CI, 
  cluster_var = "studyid_rev",
  samples = 499, 
  conf_int = TRUE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr5, paste0(RES_DIR, "Risk/", "risk_tr5_hdPS.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")








# hdPS with ML analysis  ###########################

## PS fitting ########
### met 3-4 days #########################
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr4_hdPS.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc", "stab_numerator", "ps_trunc",
                                                  "ow_wt", "ow_wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length", "prior_episode"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, paste0(RES_DIR, "PScore/", "ps_weight_tr4_hdPS_ML.rds"))

### met 5-14 days #########################
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr5_hdPS.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length", "prior_episode"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, paste0(RES_DIR, "PScore/", "ps_weight_tr5_hdPS_ML.rds"))


## Model fitting ########
risk <- list(
  tr4 = list(
    od = list(), dis = list()
  ),
  tr5 = list(
    od = list(), dis = list()
  )
)

### met 3-4 days #########################
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr4_hdPS_ML.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr4$od <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = "year", 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)


risk$tr4$dis <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = "year",
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr4, paste0(RES_DIR, "Risk/", "risk_tr4_hdPS_ML.rds"))

cat("  ✓ Trial 2 risk estimates saved.\n")


### met 5-14 days #########################

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr5_hdPS_ML.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr5$od <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CI, 
  cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr5$dis <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CI, 
  cluster_var = "studyid_rev",
  samples = 499, 
  conf_int = TRUE, 
  ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr5, paste0(RES_DIR, "Risk/", "risk_tr5_hdPS_ML.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")



