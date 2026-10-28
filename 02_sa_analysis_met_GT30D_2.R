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
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/GT30D/"
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
  filter(ret_length > 30)

dat_tr3_first <- raw_first |>
  filter(trial_cat == "5-14") |>
  clean_trt() |>
  filter(ret_length > 30)

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

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss", "prior_episode"
)
continuous_vars <- c("ddose_prev", "year", "ret_length", "prior_episode")



raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_hdPS() |>
  filter(ret_length > 30)

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_hdPS()|>
  filter(ret_length >30)

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



# 4. hdPS with ML analysis  ###########################

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
rm(list = ls()); gc()

library(tidyverse)
library(haven)
library(splines)
library(sandwich)
library(mvtnorm)
library(data.table)
library(Matrix)

source("H:/GitHub/Missed-Dose-Protocols-PRI/00_functions.R")

DATA_DIR <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/GT30D/"
dir.create(file.path(RES_DIR, "PScore"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RES_DIR, "Risk"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RES_DIR, "Tables"), recursive = TRUE, showWarnings = FALSE)

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
CLUST_VARS <- c("studyid_rev", "episode", "time0")
MAX_FOLLOWUP <- 365

baselinevars <- c("sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute",
                  "cci_gt0", "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", "udt",
                  "ddose_prev", "clin_size_b", "virtual_care", "pss")

continuous_vars <- c("ddose_prev", "year", "ret_length")

default_out_vars <- c("sens_discont", "time_discont", "sens_od_acute", "time_od_acute",
                      "sens_all_acute", "time_all_acute")

pp_out_vars <- c(default_out_vars, "time_od_acute_pp", "sens_od_acute_pp", "is_censored_od",
                 "interruption_start_date", "interruption_day3_date", "time_to_od_censor")

## 1. HELPER FUNCTIONS ###########################################################

clean_trt_all <- function(dat) {
  dat |> mutate(trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)) |>
    rename(trt = trt_s456) |> filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS), prior_episode = as.numeric(episode - 1)) |>
    arrange(studyid_rev, episode, time0) |> group_by(studyid_rev, episode) |>
    mutate(prior_gap = lag(cumsum(gap_days), default = 0), prior_gap_num = row_number() - 1L) |> ungroup()
}

## 2. LOAD ALL REINITIATIONS WITHIN EPISODES #####################################

raw_first <- read_sas(file.path(DATA_DIR, "tt_analytic_rprev_sens_nw_m.sas7bdat")) |>
  mutate(time0 = as.Date(time0))

interruption_registry <- raw_first |> filter(trial_gap >= 3) |>
  transmute(studyid_rev, episode, subsequent_restart_date = time0,
            interruption_days = as.integer(trial_gap),
            interruption_start_date = subsequent_restart_date - interruption_days,
            interruption_day3_date = interruption_start_date + 2L) |>
  distinct()



stopifnot(all(interruption_registry$interruption_days >= 3),
          all(interruption_registry$interruption_day3_date == interruption_registry$interruption_start_date + 2L))


raw_index <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat")) |>
  mutate(time0 = as.Date(time0))

dat_tr2 <- raw_index |> filter(trial_gap %in% 3:4, ret_length > 30) |> clean_trt_all()
dat_tr3 <- raw_index |> filter(trial_cat == "5-14", ret_length > 30) |> clean_trt_all()

apply_day3_censoring <- function(dat_index, max_followup = 365) {
  
  next_interruption <- dat_index |> select(studyid_rev, episode, time0) |>
    left_join(interruption_registry, by = c("studyid_rev", "episode"), relationship = "many-to-many") |>
    filter(interruption_start_date > time0, interruption_day3_date <= time0 + max_followup) |>
    group_by(studyid_rev, episode, time0) |>
    slice_min(interruption_start_date, n = 1, with_ties = FALSE) |> ungroup() |>
    select(studyid_rev, episode, time0, interruption_start_date, interruption_day3_date)
  
  out <- dat_index |> left_join(next_interruption, by = c("studyid_rev", "episode", "time0")) |>
    mutate(time_to_od_censor = as.numeric(interruption_day3_date - time0),
           time_to_od_censor = replace_na(time_to_od_censor, as.numeric(max_followup)),
           time_to_od_censor = pmin(time_to_od_censor, as.numeric(max_followup)),
           time_od_acute_pp = pmin(time_od_acute, time_to_od_censor, max_followup, na.rm = TRUE),
           sens_od_acute_pp = as.integer(sens_od_acute == 1 & time_od_acute <= time_to_od_censor & time_od_acute <= max_followup),
           is_censored_od = as.integer(!is.na(interruption_day3_date) & time_to_od_censor < time_od_acute &
                                         time_to_od_censor < max_followup))
  
  stopifnot(all(out$time_to_od_censor >= 0, na.rm = TRUE),
            all(out$time_od_acute_pp <= max_followup, na.rm = TRUE))
  
  out
}

dat_tr2_pp <- apply_day3_censoring(dat_tr2, 365)
dat_tr3_pp <- apply_day3_censoring(dat_tr3, 365)


## 3. AUDIT CENSORING ############################################################

safe_min <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_median <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

censoring_audit <- bind_rows(`3-4 days` = dat_tr2_pp, `5-14 days` = dat_tr3_pp, .id = "trial") |>
  summarise(episodes = n(), n_with_next_interruption = sum(!is.na(interruption_day3_date)),
            n_censored = sum(is_censored_od), pct_censored = 100 * mean(is_censored_od),
            min_days_to_censor = safe_min(time_to_od_censor[is_censored_od == 1]),
            median_days_to_censor = safe_median(time_to_od_censor[is_censored_od == 1]),
            original_events_30d = sum(sens_od_acute == 1 & time_od_acute <= 30, na.rm = TRUE),
            retained_events_30d = sum(sens_od_acute_pp == 1 & time_od_acute_pp <= 30, na.rm = TRUE), .by = trial)

print(censoring_audit, n = Inf)

## 4. BASELINE PROPENSITY-SCORE WEIGHTS ##########################################

ps_tr2_pp <- ps_model_wt(dat = dat_tr2_pp, exposure = "trt", out_vars = pp_out_vars, clust_vars = CLUST_VARS,
                         baselinevars = c(baselinevars, "prior_episode"),
                         continuous_vars = c(continuous_vars, "prior_episode"))

ps_tr3_pp <- ps_model_wt(dat = dat_tr3_pp, exposure = "trt", out_vars = pp_out_vars, clust_vars = CLUST_VARS,
                         baselinevars = c(baselinevars, "prior_episode"),
                         continuous_vars = c(continuous_vars, "prior_episode"))

dat_tr2_weighted <- ps_tr2_pp$data_weights
dat_tr3_weighted <- ps_tr3_pp$data_weights






PP_INTERVAL <- 3

long_tr2_od <- expand_person_period(dat_tr2_weighted, time_var = "time_od_acute_pp",
                                    event_var = "sens_od_acute_pp", interval_length = PP_INTERVAL)

long_tr3_od <- expand_person_period(dat_tr3_weighted, time_var = "time_od_acute_pp",
                                    event_var = "sens_od_acute_pp", interval_length = PP_INTERVAL)





add_censor_event <- function(dat) {
  dat |> mutate(pp_id = interaction(studyid_rev, episode, time0, drop = TRUE)) |>
    arrange(pp_id, interval) |> group_by(pp_id) |>
    mutate(censor_event = as.integer(is_censored_od == 1 & row_number() == n())) |> ungroup()
}

long_tr2_od <- add_censor_event(long_tr2_od)
long_tr3_od <- add_censor_event(long_tr3_od)





add_od_ipcw <- function(long_dat, covariates, ps_weight = "wt_stab", trunc = c(0, .99)) {
  
  den_model <- glm(reformulate(c("trt", "ns(interval, df=3)", covariates), response = "censor_event"),
                   data = long_dat, family = binomial())
  
  num_model <- glm(reformulate(c("trt", "ns(interval, df=3)"), response = "censor_event"),
                   data = long_dat, family = binomial())
  
  out <- long_dat |>
    mutate(p_cens_den = pmin(pmax(predict(den_model, type = "response"), 1e-6), 1 - 1e-6),
           p_cens_num = pmin(pmax(predict(num_model, type = "response"), 1e-6), 1 - 1e-6)) |>
    arrange(pp_id, interval) |> group_by(pp_id) |>
    mutate(prob_uncens_den = lag(cumprod(1 - p_cens_den), default = 1),
           prob_uncens_num = lag(cumprod(1 - p_cens_num), default = 1),
           stab_ipcw_od = prob_uncens_num / pmax(prob_uncens_den, .01),
           wt_pp_od_raw = .data[[ps_weight]] * stab_ipcw_od) |> ungroup()
  
  weight_limits <- quantile(out$wt_pp_od_raw, trunc, na.rm = TRUE)
  
  out |> mutate(wt_pp_od = pmin(pmax(wt_pp_od_raw, weight_limits[1]), weight_limits[2]),
                wt_pp_od_lower = weight_limits[1], wt_pp_od_upper = weight_limits[2])
}

ipcw_covariates <- c(baselinevars, "prior_episode")

long_tr2_od <- add_od_ipcw(long_tr2_od, ipcw_covariates)
long_tr3_od <- add_od_ipcw(long_tr3_od, ipcw_covariates)


weight_ess <- function(w) sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)

audit_weights <- function(dat, trial) {
  last_rows <- dat |> group_by(pp_id) |> slice_tail(n = 1) |> ungroup()
  
  last_rows |> summarise(trial = trial, episodes = n(), censored = sum(is_censored_od),
                         min_ipcw = min(stab_ipcw_od, na.rm = TRUE),
                         mean_ipcw = mean(stab_ipcw_od, na.rm = TRUE),
                         median_ipcw = median(stab_ipcw_od, na.rm = TRUE),
                         p99_ipcw = quantile(stab_ipcw_od, .99, na.rm = TRUE),
                         max_ipcw = max(stab_ipcw_od, na.rm = TRUE),
                         mean_final_weight = mean(wt_pp_od, na.rm = TRUE),
                         p99_final_weight = quantile(wt_pp_od, .99, na.rm = TRUE),
                         max_final_weight = max(wt_pp_od, na.rm = TRUE),
                         effective_sample_size = weight_ess(wt_pp_od))
}

weight_audit <- bind_rows(audit_weights(long_tr2_od, "3-4 days"),
                          audit_weights(long_tr3_od, "5-14 days"))

## 5. SMD LOVE PLOTS ############

library(tidyverse)
library(cobalt)
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)

VAR_LABELS <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than treatment* before OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of treatment episode before time zero",
  home_grp_b    = "Unstable housing in past 5 years",
  incar_grp_bl  = "Incarcerated in the past year",
  drug_acute    = "Drug-related acute care visits in the past month",
  cci_gt0       = "Charlson Comorbidity Index (>0)",
  copd_asthma   = "Asthma or COPD (ever)",
  sed_psy_med   = "Dispensations of psychiatric or sedative medication (past month)",
  po            = "Opioid dispensations other than treatment* (past month)",
  sud           = "Other substance use disorder (ever)",
  aud           = "Alcohol use disorder (ever)",
  smd           = "Serious mental disorder (ever)",
  hcv           = "Hepatitis C (ever)",
  tob           = "Tobacco use disorder (past year)",
  udt           = "Urine drug test in the past month",
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Concurrent opioid and stimulant dispensations (past month)",
  prior_gap     = "Total treatment interruptions since initiation",
  prior_episode = "Number of prior treatment episodes"
)


TRT_PLOT_LABELS <- c("no reversion" = "No dose change", "1-<50%" = "1-<50% reduction",
                     ">=50%/<=40mg" = "≥50% reduction/≤40 mg", "increase" = "Dose increase")

plot_var_labels <- VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")]
fmla <- reformulate(c(baselinevars, "prior_episode"), response = "trt")
x_low <- -0.6; x_high <- 0.6; x_breaks <- seq(-0.4, 0.4, 0.2)

prepare_balance_data <- function(dat) {
  dat |> mutate(trt = recode(as.character(trt), !!!TRT_PLOT_LABELS),
                trt = factor(trt, levels = unname(TRT_PLOT_LABELS)),
                incar_grp_bl = if_else(as.character(incar_grp_bl) == "1_pastyear", 1, 0),
                prior_episode = as.numeric(prior_episode))
}

make_love_panel <- function(dat, comparison, title, show_y = FALSE) {
  d <- dat |> filter(trt %in% comparison) |> mutate(trt = droplevels(trt))
  balance <- bal.tab(fmla, data = d, weights = d$wt_trunc, estimand = "ATE", s.d.denom = "pooled",
                     un = TRUE, stats = "mean.diffs", m.threshold = 0.1)
  
  love.plot(balance, stat = "mean.diffs", stars = "none", thresholds = c(m = 0.1), abs = FALSE,
            var.names = plot_var_labels, sample.names = c("Before weighting", "After weighting"),
            position = "none", line = FALSE, colors = c("black", "red"), shapes = c(17, 16), size = 3.5) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_x_continuous(limits = c(x_low, x_high), breaks = x_breaks, expand = c(0, 0)) +
    labs(x = "Standardized mean difference", y = NULL, title = title) +
    theme(title = element_text(color = "black", size = 12), panel.grid = element_blank(),
          axis.line = element_line(color = "black"), axis.ticks.length = unit(0.2, "cm"),
          axis.ticks = element_line(linewidth = 0.8, color = "black"),
          axis.text.x = element_text(color = "black", size = 11),
          axis.text.y = if (show_y) element_text(color = "black", size = 10) else element_blank(),
          axis.title.x = element_text(size = 11), legend.title = element_blank(),
          legend.text = element_text(size = 11, color = "black"), legend.key = element_blank())
}

make_trial_love_plot <- function(dat, trial_label) {
  dat <- prepare_balance_data(dat)
  
  p1 <- make_love_panel(dat, c("No dose change", "1-<50% reduction"),
                        "1-<50% reduction\nvs no change", TRUE)
  
  p2 <- make_love_panel(dat, c("No dose change", "≥50% reduction/≤40 mg"),
                        "≥50% reduction/≤40 mg\nvs no change")
  
  p3 <- make_love_panel(dat, c("No dose change", "Dose increase"),
                        "Dose increase\nvs no change")
  
  panels <- (p1 + plot_spacer() + p2 + plot_spacer() + p3) +
    plot_layout(widths = c(1, 0.05, 1, 0.05, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_blank(), legend.text = element_text(size = 11, color = "black"),
          legend.key = element_blank())
  
  panels <- plot_spacer() + panels + plot_layout(heights = c(0.03, 0.97))
  
  ggdraw(panels) +
    draw_label(paste0(trial_label, "\nNote: Covariates with |SMD| ≥0.10 were included in the outcome model."),
               x = 0.5, y = 0.005, hjust = 0.5, vjust = 0, size = 11)
}



plot_tr2_pp <- make_trial_love_plot(dat_tr2_weighted, "Interruption duration: 3-4 days")
plot_tr3_pp <- make_trial_love_plot(dat_tr3_weighted, "Interruption duration: 5-14 days")

jpeg(file.path(RES_DIR, "Figures", "smd_love_plots_tr2_pp_od.jpg"), width = 4200, height = 2400, res = 300)
print(plot_tr2_pp)
dev.off()

jpeg(file.path(RES_DIR, "Figures", "smd_love_plots_tr3_pp_od.jpg"), width = 4200, height = 2400, res = 300)
print(plot_tr3_pp)
dev.off()


## 6. Outcome models #############


HORIZONS_CI <- c(14, 28, 30)

risk_pp_od <- list(tr2 = list(od = list()), tr3 = list(od = list()))

risk_pp_od$tr2$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_tr2_od, trt_var = "trt", covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE,
  ref_group = "no reversion", wt_type = "wt_pp_od", rmst_tau = NULL
)

risk_pp_od$tr3$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_tr3_od, trt_var = "trt", covs = c("year", "ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, conf_int = TRUE,
  ref_group = "no reversion", wt_type = "wt_pp_od", rmst_tau = NULL
)




## 7. SAVE FULL MODEL RESULTS ###################################################

saveRDS(risk_pp_od$tr2, file.path(RES_DIR, "Risk", "risk_tr2_day3_pp_od.rds"))
saveRDS(risk_pp_od$tr3, file.path(RES_DIR, "Risk", "risk_tr3_day3_pp_od.rds"))

saveRDS(long_tr2_od, file.path(RES_DIR, "PScore", "long_tr2_day3_pp_od.rds"))
saveRDS(long_tr3_od, file.path(RES_DIR, "PScore", "long_tr3_day3_pp_od.rds"))




TRT_LABELS <- c("no reversion" = "No dose change", "1-<50%" = "<50% reduction",
                ">=50%/<=40mg" = "≥50% reduction or ≤40 mg", "increase" = "Dose increase")

TRT_ORDER <- c("No dose change", "<50% reduction", "≥50% reduction or ≤40 mg", "Dose increase")

extract_pp_od_results <- function(result, trial, horizon = 30) {
  
  risks <- result$risks |> filter(horizon == .env$horizon) |>
    transmute(trt = as.character(trt), risk = 100 * risk, risk_lower = 100 * risk_lower, risk_upper = 100 * risk_upper)
  
  contrasts <- result$contrasts |> filter(horizon == .env$horizon) |>
    transmute(trt = as.character(trt), RD = 100 * RD, RD_lower = 100 * RD_lower, RD_upper = 100 * RD_upper,
              RR = RR, RR_lower = RR_lower, RR_upper = RR_upper)
  
  risks |> left_join(contrasts, by = "trt") |>
    mutate(Trial = trial, Horizon = horizon, Strategy = recode(trt, !!!TRT_LABELS),
           Strategy = factor(Strategy, levels = TRT_ORDER),
           `Risk, % (95% CI)` = sprintf("%.2f (%.2f, %.2f)", risk, risk_lower, risk_upper),
           `Risk difference, % (95% CI)` = if_else(trt == "no reversion", "Reference",
                                                   sprintf("%.2f (%.2f, %.2f)", RD, RD_lower, RD_upper)),
           `Risk ratio (95% CI)` = if_else(trt == "no reversion", "Reference",
                                           sprintf("%.2f (%.2f, %.2f)", RR, RR_lower, RR_upper))) |>
    arrange(Strategy) |>
    select(Trial, Horizon, Strategy, `Risk, % (95% CI)`, `Risk difference, % (95% CI)`,
           `Risk ratio (95% CI)`, risk, risk_lower, risk_upper, RD, RD_lower, RD_upper, RR, RR_lower, RR_upper)
}



pp_od_30d_results <- bind_rows(
  extract_pp_od_results(risk_pp_od$tr2$od$ci, trial = "3-4 days", horizon = 30),
  extract_pp_od_results(risk_pp_od$tr3$od$ci, trial = "5-14 days", horizon = 30)
)

print(pp_od_30d_results, n = Inf)



pp_od_30d_table <- pp_od_30d_results |>
  select(`Interruption duration` = Trial, `Dose-adjustment strategy` = Strategy,
         `30-day standardized risk, % (95% CI)` = `Risk, % (95% CI)`,
         `30-day risk difference, % (95% CI)` = `Risk difference, % (95% CI)`)

write_csv(pp_od_30d_table, file.path(RES_DIR, "Tables", "day3_pp_od_30d_table.csv"))






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

## ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/GT30D/"

## ── Shared settings ────────────────────────────────────────────────────────────
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
  filter(ret_length > 30)

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all() |>
  filter(ret_length > 30)


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




# 7. Clusterl; HSDA ########
## > 30 days ############# 

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")   # ps_model_wt(), predict_cuminc(), expand_person_period()

### ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v3/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/GT30D/"

### ── Shared settings ────────────────────────────────────────────────────────────
TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")

out_vars   <- c("sens_discont", "time_discont",
                "sens_od_acute", "time_od_acute",
                "sens_all_acute", "time_all_acute")

clust_vars <- c("clnt_hsda", "clnt_lha", "studyid_rev", "episode", "time0")

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "virtual_care", "pss", "comsize", "n_pats", "clnt_hsda"
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
  n_pats        = "Number of patients served",
  clnt_hsda       = "Health service delivery areas"
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

raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_cl_sens_nw_m.sas7bdat"))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all() |>
  filter(ret_length > 30)

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all() |>
  filter(ret_length > 30)


# 2. PROPENSITY SCORE WEIGHTS 

ps_tr2_all <- ps_model_gow(
  dat             = dat_tr2_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_episode"),
  continuous_vars = c(continuous_vars, "prior_episode"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_all, paste0(RES_DIR, "PScore/", "ps_weight_tr2_cluster_hsda.rds"))


ps_tr3_all <- ps_model_gow(
  dat             = dat_tr3_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_episode"),
  continuous_vars = c(continuous_vars, "prior_episode"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr3_all, paste0(RES_DIR, "PScore/", "ps_weight_tr3_cluster_hsda.rds"))


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

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_cluster_hsda.rds"))

ps_fit$data_weights |> count(trt)


# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$od$curve <- predict_cuminc_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode", "comsize", "n_pats"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$od$ci <- predict_cuminc_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode", "comsize", "n_pats"),
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

risk$tr2$dis$curve <- predict_cuminc_fast(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode", "comsize", "n_pats"),
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$dis$ci <- predict_cuminc_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ddose_prev", "prior_episode", "comsize", "n_pats"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_cluster_hsda.rds"))


cat("  ✓ Trial 2 risk estimates saved.\n")


cat("\nTrial 3 (5-14 days)...\n")


cat("  Computing: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_cluster_hsda.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$od$curve <- predict_cuminc_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode", "comsize", "n_pats"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$od$ci <- predict_cuminc_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode", "comsize", "n_pats"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$dis$curve <- predict_cuminc_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode", "comsize", "n_pats"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$dis$ci <- predict_cuminc_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("year", "ddose_prev", "prior_episode", "comsize", "n_pats"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_cluster_hsda.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")


