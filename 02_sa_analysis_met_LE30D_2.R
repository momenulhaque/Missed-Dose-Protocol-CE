# 1. INCIDENT USERS ####
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
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/LE30D/"
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

# Shared pre-processing: collapse all reduction categories into one
clean_trt <- function(dat) {
  dat |>
    mutate(
      trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)
    ) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS))
}
raw_first <- read_sas(file.path(DATA_DIR, "TT_analytic_incid_sens_nw_m.sas7bdat"))

dat_tr2_first <- raw_first |>
  filter(trial_gap %in% 3:4) |>
  clean_trt() |>
  filter(ret_length <= 30)

dat_tr3_first <- raw_first |>
  filter(trial_cat == "5-14") |>
  clean_trt() |>
  filter(ret_length <= 30)

rm(raw_first); gc()

ps_tr2_first <- ps_model_wt(
  dat             = dat_tr2_first,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_first, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))

ps_tr3_first <- ps_model_wt(
  dat             = dat_tr3_first,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr3_first, paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))


# Risk curves and confidence intervals for overdose and discontinuation
cat("\n=== Computing Cumulative Incidence ===\n")

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))

# ── INITIALIZE RESULTS STRUCTURE ──

risk <- list(
  tr2 = list(od = list(), dis = list()),
  tr3 = list(od = list(), dis = list())
)

## MET interruption of 3-4 days #############

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev"), 
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

risk$tr2$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs = c("year", "ddose_prev"),  
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_first.rds"))


cat("  ✓ Trial 2 risk estimates saved.\n")

## MET interruption of 5-14 days #############


ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev"), 
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev"), 
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_first.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")





# 2. Overlap-weights ##################
HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))

# ── INITIALIZE RESULTS STRUCTURE ──

risk <- list(
  tr2 = list(od = list(), dis = list()),
  tr3 = list(od = list(), dis = list())
)

## MET interruption of 3-4 days #############
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))



# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "ow_wt"
)

risk$tr2$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "ow_wt"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)

risk$tr2$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"),  
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "ow_wt"
)

risk$tr2$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "ow_wt"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_ow.rds"))


cat("  ✓ Trial 2 risk estimates saved.\n")

## MET interruption of 5-14 days #############


ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

risk$tr3$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

risk$tr3$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_ow.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")



# 3. hdPS analysis  ###########################
clean_trt_hdPS <- function(dat) {
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
      prior_gap_num = row_number(),
      unique_id = paste0(studyid_rev, "_", time0)
    ) %>%
    ungroup()
  
}

raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_hdPS() |>
  filter(ret_length <= 30)

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_hdPS()|>
  filter(ret_length <= 30)

rm(raw_all); gc()

# hdPS analysis  ###########################
## PS fitting ########
### met 3-4 days #########################

dat23_all <- get_dat_all_prev_m(dat1=dat_tr2_all)
dim(dat23_all)
# empirical covariates
vars.empirical <- setdiff(names(dat23_all), names(dat_tr2_all))
length(vars.empirical)
# 58
dat23_all <- dat23_all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical_le200 <- select_top200_empirical_covariates(dat.all=dat23_all,
                                                           time_var = "time_discont",
                                                           event_var = "sens_discont",
                                                           max_vars = 200, nfolds = 5,
                                                           vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat23_all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical_le200),
                        continuous_vars = c("age_st", "ddose_prev", "year", 
                                            "ret_length", "prior_episode"))



saveRDS(ps_fit, paste0(RES_DIR, "PScore/", "ps_weight_34d_hdPS.rds"))


### met 5-14 days #########################

dat514_all <- get_dat_all_prev_m(dat1=dat_tr3_all)
dim(dat514_all)
# empirical covariates
vars.empirical <- setdiff(names(dat514_all), names(dat_tr3_all))
length(vars.empirical)
# 84
dat514_all <- dat514_all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical_le200 <- select_top200_empirical_covariates(dat.all=dat514_all,
                                                           time_var = "time_discont",
                                                           event_var = "sens_discont",
                                                           max_vars = 200, nfolds = 5,
                                                           vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat514_all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical_le200),
                        continuous_vars = c("age_st", "ddose_prev", "year", 
                                            "ret_length", "prior_episode"))



saveRDS(ps_fit, paste0(RES_DIR, "PScore/", "ps_weight_514d_hdPS.rds"))


## Model fitting ########
risk <- list(
  tr2 = list(
    od = list(), dis = list()
  ),
  tr3 = list(
    od = list(), dis = list()
  )
)

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))


### met 3-4 days #########################
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_34d_hdPS.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"), 
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

risk$tr2$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"),  
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_hdPS.rds"))


cat("  ✓ Trial 2 risk estimates saved.\n")

### met 5-14 days #########################

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_514d_hdPS.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_hdPS.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")



# hdPS with ML analysis  ###########################

## PS fitting ########
### met 3-4 days #########################
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_34d_hdPS.rds"))

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

saveRDS(ps_fit_ml, paste0(RES_DIR, "PScore/", "ps_weight_34d_hdPS_ML.rds"))

### met 5-14 days #########################
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_514d_hdPS.rds"))

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

saveRDS(ps_fit_ml, paste0(RES_DIR, "PScore/", "ps_weight_514d_hdPS_ML.rds"))


## Model fitting ########
risk <- list(
  tr2 = list(
    od = list(), dis = list()
  ),
  tr3 = list(
    od = list(), dis = list()
  )
)

### met 3-4 days #########################
names(ps_fit$data_weights)

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_34d_hdPS_ML.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"), 
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

risk$tr2$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"),  
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_hdPS_ML.rds"))

cat("  ✓ Trial 2 risk estimates saved.\n")


### met 5-14 days #########################

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_514d_hdPS_ML.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()


saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_hdPS_ML.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")


# 5. Per-protocol #############

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
clust_vars <- c("studyid_rev", "episode", "time0")

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
raw_first <- read_sas(file.path(DATA_DIR, "tt_analytic_rprev_sens_nw_m.sas7bdat")) |>
  arrange(studyid_rev, episode, time0) |>
  mutate(
    trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)
  ) |>
  rename(trt = trt_s456) |>
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
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_sa2 <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all() |>
  filter(ret_length <= 30)

names(dat_tr3_sa2)

dat_tr3_sa2 <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all() |>
  filter(ret_length <= 30)

# ── 2. APPLY PER-PROTOCOL CENSORING ───────────────────────────────────────────
apply_pp_censoring <- function(dat_index) {
  
  # Find the exact next gap >= 3 days for each index event
  censor_dates <- dat_index %>%
    left_join(all_gaps, by = c("studyid_rev", "episode"), relationship = "many-to-many") %>%
    filter(next_time0 > time0) %>%          # Must happen AFTER the index reinitiation
    group_by(studyid_rev, episode, time0) %>%
    dplyr::slice_min(order_by = next_time0, n = 1, with_ties = FALSE) %>% # <--- FIXED HERE
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
dat_tr2_pp <- apply_pp_censoring(dat_tr2_sa2)
dat_tr3_pp <- apply_pp_censoring(dat_tr3_sa2)

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


# ── 4. EXECUTION FOR TRIAL 2 (3-4 DAYS) ───────────────────────────────────────
cat("\n--- Running Trial 2 (3-4 days) Per-Protocol ---\n")

# A. Calculate Baseline IPTW (Uses pp_out_vars)
ps_tr2_pp <- ps_model_wt(
  dat = dat_tr2_pp, exposure = "trt", out_vars = pp_out_vars, 
  clust_vars = clust_vars, baselinevars = c(baselinevars, "prior_episode"), 
  continuous_vars = c(continuous_vars, "prior_episode")
)
dat_tr2_weighted <- ps_tr2_pp$data_weights

# B. Calculate IPCW and Combined Weights for Overdose
dat_tr2_weighted$ipcw_od <- add_ipcw(
  dat = dat_tr2_weighted, time_var = "time_od_acute_pp", 
  censor_var = "is_censored_od", covariates = c(baselinevars, "prior_episode")
)

dat_tr2_weighted$wt_pp_od_raw <- dat_tr2_weighted$wt_stab * dat_tr2_weighted$ipcw_od
t_low_od <- quantile(dat_tr2_weighted$wt_pp_od_raw, 0.01, na.rm = TRUE)
t_high_od <- quantile(dat_tr2_weighted$wt_pp_od_raw, 0.99, na.rm = TRUE)
dat_tr2_weighted$wt_pp_od <- pmin(pmax(dat_tr2_weighted$wt_pp_od_raw, t_low_od), t_high_od)

# C. Calculate IPCW and Combined Weights for Discontinuation
dat_tr2_weighted$ipcw_dis <- add_ipcw(
  dat = dat_tr2_weighted, time_var = "time_discont_pp", 
  censor_var = "is_censored_dis", covariates = c(baselinevars, "prior_episode")
)

dat_tr2_weighted$wt_pp_dis_raw <- dat_tr2_weighted$wt_stab * dat_tr2_weighted$ipcw_dis
t_low_dis <- quantile(dat_tr2_weighted$wt_pp_dis_raw, 0.01, na.rm = TRUE)
t_high_dis <- quantile(dat_tr2_weighted$wt_pp_dis_raw, 0.99, na.rm = TRUE)
dat_tr2_weighted$wt_pp_dis <- pmin(pmax(dat_tr2_weighted$wt_pp_dis_raw, t_low_dis), t_high_dis)


# D. Outcome Models
HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))
risk_pp <- list(tr2 = list(od = list(), dis = list()), tr3 = list(od = list(), dis = list()))

# Overdose (Using prior_gap instead of prior_episode)
long_dat_tr2_od <- expand_person_period(dat_tr2_weighted, time_var = "time_od_acute_pp", event_var = "sens_od_acute_pp", interval_length = 14)
risk_pp$tr2$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat_tr2_od, trt_var = "trt", covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE, 
  ref_group = "no reversion", wt_type = "wt_pp_od"
)

# Discontinuation (Using prior_gap)
long_dat_tr2_dis <- expand_person_period(dat_tr2_weighted, time_var = "time_discont_pp", event_var = "sens_discont_pp", interval_length = 14)
risk_pp$tr2$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat_tr2_dis, trt_var = "trt", covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE, 
  ref_group = "no reversion", wt_type = "wt_pp_dis"
)

saveRDS(risk_pp$tr2, paste0(RES_DIR, "risk/risk_tr2_per_protocol.rds"))


# ── 5. EXECUTION FOR TRIAL 3 (5-14 DAYS) ──────────────────────────────────────
cat("\n--- Running Trial 3 (5-14 days) Per-Protocol ---\n")

# A. Calculate Baseline IPTW
ps_tr3_pp <- ps_model_wt(
  dat = dat_tr3_pp, exposure = "trt", out_vars = pp_out_vars, 
  clust_vars = clust_vars, baselinevars = c(baselinevars, "prior_episode"), 
  continuous_vars = c(continuous_vars, "prior_episode")
)
dat_tr3_weighted <- ps_tr3_pp$data_weights

# B. Calculate IPCW and Combined Weights for Overdose
dat_tr3_weighted$ipcw_od <- add_ipcw(
  dat = dat_tr3_weighted, time_var = "time_od_acute_pp", 
  censor_var = "is_censored_od", covariates = c(baselinevars, "prior_episode")
)
dat_tr3_weighted$wt_pp_od_raw <- dat_tr3_weighted$wt_stab * dat_tr3_weighted$ipcw_od
t_low_od_3 <- quantile(dat_tr3_weighted$wt_pp_od_raw, 0.01, na.rm = TRUE)
t_high_od_3 <- quantile(dat_tr3_weighted$wt_pp_od_raw, 0.99, na.rm = TRUE)
dat_tr3_weighted$wt_pp_od <- pmin(pmax(dat_tr3_weighted$wt_pp_od_raw, t_low_od_3), t_high_od_3)

# C. Calculate IPCW and Combined Weights for Discontinuation
dat_tr3_weighted$ipcw_dis <- add_ipcw(
  dat = dat_tr3_weighted, time_var = "time_discont_pp", 
  censor_var = "is_censored_dis", covariates = c(baselinevars, "prior_episode")
)
dat_tr3_weighted$wt_pp_dis_raw <- dat_tr3_weighted$wt_stab * dat_tr3_weighted$ipcw_dis
t_low_dis_3 <- quantile(dat_tr3_weighted$wt_pp_dis_raw, 0.01, na.rm = TRUE)
t_high_dis_3 <- quantile(dat_tr3_weighted$wt_pp_dis_raw, 0.99, na.rm = TRUE)
dat_tr3_weighted$wt_pp_dis <- pmin(pmax(dat_tr3_weighted$wt_pp_dis_raw, t_low_dis_3), t_high_dis_3)

# D. Outcome Models
# Overdose (Using prior_gap)
long_dat_tr3_od <- expand_person_period(dat_tr3_weighted, time_var = "time_od_acute_pp", event_var = "sens_od_acute_pp", interval_length = 14)
risk_pp$tr3$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat_tr3_od, trt_var = "trt", covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE, 
  ref_group = "no reversion", wt_type = "wt_pp_od"
)

# Discontinuation (Using prior_gap)
long_dat_tr3_dis <- expand_person_period(dat_tr3_weighted, time_var = "time_discont_pp", event_var = "sens_discont_pp", interval_length = 14)
risk_pp$tr3$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat_tr3_dis, trt_var = "trt", covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE, 
  ref_group = "no reversion", wt_type = "wt_pp_dis"
)

saveRDS(risk_pp$tr3, paste0(RES_DIR, "risk/risk_tr3_per_protocol.rds"))

cat("\n✓ Per-Protocol Analysis Complete. Results saved.\n")


# 6. Cluster ############
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
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/LE30D/"

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
  "udt", "ddose_prev", "virtual_care", "pss", "comsize", "n_pats"
)
continuous_vars <- c("ddose_prev", "year", "ret_length", "comsize", "n_pats")

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
  prior_episode = "Number of prior episode",
  comsize       = "Number of prescribers shared",
  n_pats        = "Number of patients served"
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

# First-interruption-in-any-episode 
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_cl_sens_nw_m.sas7bdat"))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all() |>
  filter(ret_length <= 30)

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all() |>
  filter(ret_length <= 30)

dat_tr2_all %>% count(com)
dat_tr3_all %>% count(com)

# 2. PROPENSITY SCORE WEIGHTS 

ps_tr2_all <- ps_model_wt(
  dat             = dat_tr2_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_episode"),
  continuous_vars = c(continuous_vars, "prior_episode"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_all, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_cl.rds"))



ps_tr3_all <- ps_model_wt(
  dat             = dat_tr3_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_episode"),
  continuous_vars = c(continuous_vars, "prior_episode"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr3_all, paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_cl.rds"))



# 7. RMST & CUMULATIVE INCIDENCE 
# Risk curves and confidence intervals for overdose and discontinuation



cat("\n=== Computing Cumulative Incidence ===\n")

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))

# ── INITIALIZE RESULTS STRUCTURE ──

risk <- list(
  tr2 = list(od = list(), dis = list()),
  tr3 = list(od = list(), dis = list()))


cat("\nTrial 2 (3-4 days)...\n")


cat("  Computing: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_cl.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"),
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

risk$tr2$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_cl.rds"))


cat("  ✓ Trial 2 risk estimates saved.\n")


cat("\nTrial 3 (5-14 days)...\n")


cat("  Computing: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_cl.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_cl.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")
