# Load Packages ##################
# 01_main_analysis.R
# Two trials - methadone interruption 3-4 and 5-14 days
# Two analytic cohort:
#   first-interruption-in-first-episode - TT_analytic_incid_sens_nw_m data set
#   first-interruption-in-any-episode - TT_analytic_prev_sens_nw_m data set
###

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
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"
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
  prior_episode = "Total prior episode"
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

clean_trt_all <- function(dat) {
  dat |>
    mutate(
      trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)
    ) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS),
           prior_episode = episode - 1) |>
    group_by(studyid_rev, episode) |>
    mutate(
      prior_gap     = cumsum(gap_days),
      prior_gap_num = row_number()
    ) %>%
    ungroup()
  
}

# 1. LOAD AND FILTER RAW SAS DATA #########

# First-interruption-in-first-episode
raw_first <- read_sas(file.path(DATA_DIR, "TT_analytic_incid_sens_nw_m.sas7bdat"))

dat_tr2_first <- raw_first |>
  filter(trial_gap %in% 3:4) |>
  clean_trt()

dat_tr3_first <- raw_first |>
  filter(trial_cat == "5-14") |>
  clean_trt()

rm(raw_first); gc()

# First-interruption-in-any-episode 
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all() 

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all()
  

rm(raw_all); gc()



# 2. PROPENSITY SCORE WEIGHTS ################

## MET 3-4 days: first-interruption-in-first-episode #####

ps_tr2_first <- ps_model_wt(
  dat             = dat_tr2_first,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_first, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))

## MET 3-4 days: first-interruption-in-any-episode #####

ps_tr2_all <- ps_model_wt(
  dat             = dat_tr2_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_episode"),
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_all, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))




## MET 5-14 days: first-interruption-in-first-episode #####
ps_tr3_first <- ps_model_wt(
  dat             = dat_tr3_first,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr3_first, paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))

## MET 5-14 days: first-interruption-in-any-episode #####
ps_tr3_all <- ps_model_wt(
  dat             = dat_tr3_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_episode"),
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr3_all, paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))




# 3. TABLE 1 ##############
# Produces one table per trial × cohort.
# Rows: covariates. Columns: treatment arms (unweighted) + SMD before/after weighting.

tbl_tr2_first <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"),
  bvars = baselinevars, cont_vars = continuous_vars, clustered = FALSE
)

tbl_tr2_all <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"),
  bvars = c(baselinevars, "prior_episode"), 
  cont_vars = c(continuous_vars, "prior_episode"), clustered = TRUE
)

tbl_tr3_first <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"),
  bvars = baselinevars, cont_vars = continuous_vars, clustered = FALSE
)

tbl_tr3_all <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"),
  bvars = c(baselinevars, "prior_episode"), 
  cont_vars = c(continuous_vars, "prior_episode"), clustered = TRUE
)


library(openxlsx)

wb <- createWorkbook()

for (nm in c("tr2_first", "tr2_all", "tr3_first", "tr3_all")) {
  tbl <- get(paste0("tbl_", nm))
  addWorksheet(wb, sheetName = nm)
  writeData(wb, sheet = nm, x = as.data.frame(tbl), rowNames = TRUE)
}

saveWorkbook(wb, paste0(RES_DIR, "Tables/", "table1.xlsx"), overwrite = TRUE)

## EXPORT SMD for TABLE1 - WIDE FORMAT #############

cat("\n=== Extracting SMD Results and Creating Excel Summary (Wide Format) ===\n")

# Function to extract SMDs from bal.tab object
extract_smd_results <- function(bal_tab_obj, comparison_name) {
  
  # Extract balance table
  balance_df <- bal_tab_obj$Balance
  
  # Get before and after weights SMDs
  results <- data.frame(
    Variable = rownames(balance_df),
    SMD_Before = as.numeric(balance_df$Diff.Un),
    SMD_After = as.numeric(balance_df$Diff.Adj),
    Comparison = comparison_name,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

## MET interruption of 3-4 days #############

cat("\nExtracting SMDs: Trial 2 (3-4 days)...\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))

dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

fmla <- as.formula(paste("trt ~", paste(baselinevars, collapse = " + ")))

# Extract SMDs for all pairwise comparisons - FIRST-IN-FIRST-EPISODE
smd_tr2_first_list <- list()

# 1-<50% vs no change
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr2_first_list[["1-<50% vs no change"]] <- extract_smd_results(prep_1, "1-<50% vs no change")


# >=50% vs no change
d_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr2_first_list[[">=50% vs no change"]] <- extract_smd_results(prep_2, ">=50% vs no change")


# Increase vs no change
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr2_first_list[["Increase vs no change"]] <- extract_smd_results(prep_3, "Increase vs no change")


# Combine into long format
smd_tr2_first_long <- bind_rows(smd_tr2_first_list) 

# Convert to wide format with proper column ordering
smd_tr2_first_wide <- smd_tr2_first_long %>%
  pivot_wider(
    names_from = Comparison,
    values_from = c(SMD_Before, SMD_After),
    names_glue = "{Comparison}_{.value}"
  ) %>%
  # Calculate max SMD per variable (actual signed value with largest absolute magnitude)
  rowwise() %>%
  mutate(
    max_SMD_Before = {
      vals <- c(`1-<50% vs no change_SMD_Before`, 
                `>=50% vs no change_SMD_Before`, 
                `Increase vs no change_SMD_Before`)
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) vals[which.max(abs(vals))] else NA
    },
    max_SMD_After = {
      vals <- c(`1-<50% vs no change_SMD_After`, 
                `>=50% vs no change_SMD_After`, 
                `Increase vs no change_SMD_After`)
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) vals[which.max(abs(vals))] else NA
    }
  ) %>%
  ungroup() %>%
  # Reorder columns: Variable, then pairs of Before/After for each comparison, then max columns
  select(Variable,
         `1-<50% vs no change_SMD_Before`, `1-<50% vs no change_SMD_After`,
         `>=50% vs no change_SMD_Before`, `>=50% vs no change_SMD_After`,
         `Increase vs no change_SMD_Before`, `Increase vs no change_SMD_After`,
         max_SMD_Before, max_SMD_After) %>%
  # Rename columns for cleaner display
  rename(
    `1-<50% Before` = `1-<50% vs no change_SMD_Before`,
    `1-<50% After` = `1-<50% vs no change_SMD_After`,
    `>=50% Before` = `>=50% vs no change_SMD_Before`,
    `>=50% After` = `>=50% vs no change_SMD_After`,
    `Increase Before` = `Increase vs no change_SMD_Before`,
    `Increase After` = `Increase vs no change_SMD_After`
  )

cat("  ✓ Trial 2, First-in-first-episode SMDs extracted (wide format)\n")

# ── Extract SMDs: FIRST-IN-ANY-EPISODE ──

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0),
    prior_episode = as.numeric(prior_episode)
  )

fmla_ep <- as.formula(paste("trt ~", paste(c(baselinevars, "prior_episode"), collapse = " + ")))

smd_tr2_ep_list <- list()

# 1-<50% vs no change
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla_ep, data = d_1, weights = d_1$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr2_ep_list[["1-<50% vs no change"]] <- extract_smd_results(prep_1, "1-<50% vs no change")

# >=50% vs no change
d_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla_ep, data = d_2, weights = d_2$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr2_ep_list[[">=50% vs no change"]] <- extract_smd_results(prep_2, ">=50% vs no change")

# Increase vs no change
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla_ep, data = d_3, weights = d_3$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr2_ep_list[["Increase vs no change"]] <- extract_smd_results(prep_3, "Increase vs no change")

smd_tr2_ep_long <- bind_rows(smd_tr2_ep_list)

# Convert to wide format with proper column ordering
smd_tr2_ep_wide <- smd_tr2_ep_long %>%
  pivot_wider(
    names_from = Comparison,
    values_from = c(SMD_Before, SMD_After),
    names_glue = "{Comparison}_{.value}"
  ) %>%
  # Calculate max SMD per variable
  rowwise() %>%
  mutate(
    max_SMD_Before = {
      vals <- c(`1-<50% vs no change_SMD_Before`, 
                `>=50% vs no change_SMD_Before`, 
                `Increase vs no change_SMD_Before`)
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) vals[which.max(abs(vals))] else NA
    },
    max_SMD_After = {
      vals <- c(`1-<50% vs no change_SMD_After`, 
                `>=50% vs no change_SMD_After`, 
                `Increase vs no change_SMD_After`)
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) vals[which.max(abs(vals))] else NA
    }
  ) %>%
  ungroup() %>%
  select(Variable,
         `1-<50% vs no change_SMD_Before`, `1-<50% vs no change_SMD_After`,
         `>=50% vs no change_SMD_Before`, `>=50% vs no change_SMD_After`,
         `Increase vs no change_SMD_Before`, `Increase vs no change_SMD_After`,
         max_SMD_Before, max_SMD_After) %>%
  rename(
    `1-<50% Before` = `1-<50% vs no change_SMD_Before`,
    `1-<50% After` = `1-<50% vs no change_SMD_After`,
    `>=50% Before` = `>=50% vs no change_SMD_Before`,
    `>=50% After` = `>=50% vs no change_SMD_After`,
    `Increase Before` = `Increase vs no change_SMD_Before`,
    `Increase After` = `Increase vs no change_SMD_After`
  )

cat("  ✓ Trial 2, First-in-any-episode SMDs extracted (wide format)\n")

## MET interruption of 5-14 days #############

cat("\nExtracting SMDs: Trial 3 (5-14 days)...\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

fmla <- as.formula(paste("trt ~", paste(baselinevars, collapse = " + ")))

smd_tr3_first_list <- list()

# 1-<50% vs no change
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr3_first_list[["1-<50% vs no change"]] <- extract_smd_results(prep_1, "1-<50% vs no change")

# >=50% vs no change
d_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr3_first_list[[">=50% vs no change"]] <- extract_smd_results(prep_2, ">=50% vs no change")

# Increase vs no change
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr3_first_list[["Increase vs no change"]] <- extract_smd_results(prep_3, "Increase vs no change")

smd_tr3_first_long <- bind_rows(smd_tr3_first_list)

# Convert to wide format with proper column ordering
smd_tr3_first_wide <- smd_tr3_first_long %>%
  pivot_wider(
    names_from = Comparison,
    values_from = c(SMD_Before, SMD_After),
    names_glue = "{Comparison}_{.value}"
  ) %>%
  # Calculate max SMD per variable
  rowwise() %>%
  mutate(
    max_SMD_Before = {
      vals <- c(`1-<50% vs no change_SMD_Before`, 
                `>=50% vs no change_SMD_Before`, 
                `Increase vs no change_SMD_Before`)
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) vals[which.max(abs(vals))] else NA
    },
    max_SMD_After = {
      vals <- c(`1-<50% vs no change_SMD_After`, 
                `>=50% vs no change_SMD_After`, 
                `Increase vs no change_SMD_After`)
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) vals[which.max(abs(vals))] else NA
    }
  ) %>%
  ungroup() %>%
  select(Variable,
         `1-<50% vs no change_SMD_Before`, `1-<50% vs no change_SMD_After`,
         `>=50% vs no change_SMD_Before`, `>=50% vs no change_SMD_After`,
         `Increase vs no change_SMD_Before`, `Increase vs no change_SMD_After`,
         max_SMD_Before, max_SMD_After) %>%
  rename(
    `1-<50% Before` = `1-<50% vs no change_SMD_Before`,
    `1-<50% After` = `1-<50% vs no change_SMD_After`,
    `>=50% Before` = `>=50% vs no change_SMD_Before`,
    `>=50% After` = `>=50% vs no change_SMD_After`,
    `Increase Before` = `Increase vs no change_SMD_Before`,
    `Increase After` = `Increase vs no change_SMD_After`
  )

cat("  ✓ Trial 3, First-in-first-episode SMDs extracted (wide format)\n")

# ── Extract SMDs: FIRST-IN-ANY-EPISODE ──

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0),
    prior_episode = as.numeric(prior_episode)
  )

fmla_ep <- as.formula(paste("trt ~", paste(c(baselinevars, "prior_episode"), collapse = " + ")))

smd_tr3_ep_list <- list()

# 1-<50% vs no change
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla_ep, data = d_1, weights = d_1$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr3_ep_list[["1-<50% vs no change"]] <- extract_smd_results(prep_1, "1-<50% vs no change")

# >=50% vs no change
d_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla_ep, data = d_2, weights = d_2$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr3_ep_list[[">=50% vs no change"]] <- extract_smd_results(prep_2, ">=50% vs no change")

# Increase vs no change
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla_ep, data = d_3, weights = d_3$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
smd_tr3_ep_list[["Increase vs no change"]] <- extract_smd_results(prep_3, "Increase vs no change")

smd_tr3_ep_long <- bind_rows(smd_tr3_ep_list)

# Convert to wide format with proper column ordering
smd_tr3_ep_wide <- smd_tr3_ep_long %>%
  pivot_wider(
    names_from = Comparison,
    values_from = c(SMD_Before, SMD_After),
    names_glue = "{Comparison}_{.value}"
  ) %>%
  # Calculate max SMD per variable
  rowwise() %>%
  mutate(
    max_SMD_Before = {
      vals <- c(`1-<50% vs no change_SMD_Before`, 
                `>=50% vs no change_SMD_Before`, 
                `Increase vs no change_SMD_Before`)
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) vals[which.max(abs(vals))] else NA
    },
    max_SMD_After = {
      vals <- c(`1-<50% vs no change_SMD_After`, 
                `>=50% vs no change_SMD_After`, 
                `Increase vs no change_SMD_After`)
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) vals[which.max(abs(vals))] else NA
    }
  ) %>%
  ungroup() %>%
  select(Variable,
         `1-<50% vs no change_SMD_Before`, `1-<50% vs no change_SMD_After`,
         `>=50% vs no change_SMD_Before`, `>=50% vs no change_SMD_After`,
         `Increase vs no change_SMD_Before`, `Increase vs no change_SMD_After`,
         max_SMD_Before, max_SMD_After) %>%
  rename(
    `1-<50% Before` = `1-<50% vs no change_SMD_Before`,
    `1-<50% After` = `1-<50% vs no change_SMD_After`,
    `>=50% Before` = `>=50% vs no change_SMD_Before`,
    `>=50% After` = `>=50% vs no change_SMD_After`,
    `Increase Before` = `Increase vs no change_SMD_Before`,
    `Increase After` = `Increase vs no change_SMD_After`
  )

cat("  ✓ Trial 3, First-in-any-episode SMDs extracted (wide format)\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SAVE TO EXCEL
# ═════════════════════════════════════════════════════════════════════════════

library(openxlsx)

wb <- createWorkbook()

# Trial 2, First-in-first-episode
addWorksheet(wb, sheetName = "TR2_First-in-first")
writeData(wb, sheet = "TR2_First-in-first", x = smd_tr2_first_wide)

# Trial 2, First-in-any-episode
addWorksheet(wb, sheetName = "TR2_First-in-any")
writeData(wb, sheet = "TR2_First-in-any", x = smd_tr2_ep_wide)

# Trial 3, First-in-first-episode
addWorksheet(wb, sheetName = "TR3_First-in-first")
writeData(wb, sheet = "TR3_First-in-first", x = smd_tr3_first_wide)

# Trial 3, First-in-any-episode
addWorksheet(wb, sheetName = "TR3_First-in-any")
writeData(wb, sheet = "TR3_First-in-any", x = smd_tr3_ep_wide)

saveWorkbook(wb, paste0(RES_DIR, "Tables/", "smd_results_all_wide.xlsx"), overwrite = TRUE)


# 4. PS DISTRIBUTION PLOTS ####
# Histogram + density of propensity scores by treatment arm,
# faceted by interruption length × cohort design.

COMPARISON_COLORS <- c(
  "No dose change"              = "#000000", 
   "Dose increase"               = "#D62728",  
  "1-<50% reduction"            = "#2CA02C",  
  ">=50% reduction/initial dose*"  = "#FF7F0E"  
)


ps_score_plot <- function(dat,
                          x_lab           = "Propensity Score",
                          y_lab           = "Density",
                          legend_position = "bottom",
                          x_points        = c(0, .25, .5, .75, 1)) {
  
  ggplot(dat, aes(x = ps, colour = group, fill = group)) +
    geom_histogram(aes(y = ..density..),
                   position = "identity", alpha = 0.25, bins = 30) +
    geom_density(alpha = 0.05, size = 0.5) +
    scale_colour_manual(values = COMPARISON_COLORS) +
    scale_fill_manual(values   = COMPARISON_COLORS) +
    facet_wrap(~Trial_design, nrow = 2, scales = "free") +
    labs(x = x_lab, y = y_lab) +
    scale_x_continuous(breaks = x_points) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks.y      = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x      = element_line(linewidth = 1.2, color = "black"),
      axis.text.x       = element_text(color = "black", size = 14),
      axis.text.y       = element_text(color = "black", size = 14),
      axis.title.y      = element_text(size = 14),
      legend.position   = legend_position,
      legend.direction  = "horizontal",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 14, color = "black"),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key        = element_blank(),
      strip.background  = element_rect(fill = "lightgray", color = NA),
      strip.text        = element_text(face = "bold", size = 13)
    )
}

prep_ps_dat <- function(ps_rds, trial_label, design_label) {
  readRDS(ps_rds)$data_weights |>
    mutate(
      group = recode(as.character(trt),
                     "no reversion" = "No dose change",
                     "1-<50%"    = "1-<50% reduction",
                     ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                     "increase"     = "Dose increase"),
      group  = factor(group, levels = c("No dose change", "1-<50% reduction", ">=50% reduction/initial dose*", "Dose increase")),
      Trial  = trial_label,
      Design = design_label
    ) |>
    select(group, ps, Trial, Design)
}


dat_ps <- bind_rows(
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"),
              "MET interruption of 3-4 days",  "First-interruption-in-first-episode"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"),
              "MET interruption of 3-4 days",  "First-interruption-in-any-episode"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"),
              "MET interruption of 5-14 days", "First-interruption-in-first-episode"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"),
              "MET interruption of 5-14 days", "First-interruption-in-any-episode")
) |>
  mutate(
    Trial_design = factor(
      paste0(Trial, "\n", Design),
      levels = c(
        "MET interruption of 3-4 days\nFirst-interruption-in-first-episode",
        "MET interruption of 3-4 days\nFirst-interruption-in-any-episode",
        "MET interruption of 5-14 days\nFirst-interruption-in-first-episode",
        "MET interruption of 5-14 days\nFirst-interruption-in-any-episode"
      )
    )
  )

plot_ps <- ps_score_plot(dat_ps) +
  guides(color = guide_legend(nrow = 2))

jpeg(paste0(RES_DIR, "Figures/", "ps_distribution.jpg"),
     width = 5500, height = 4200, res = 300)
plot_ps
dev.off()




# 5. SMD LOVE PLOTS ############
# Each trial (3-4 days, 5-14 days) has 2 cohorts
#   - First-interruption-in-first-episode
#   - First-interruption-in-any-episode


cat("\n=== Creating SMD plots (2 trials × 2 cohorts = 4 figures) ===\n")

COMPARISON_COLORS <- c(
  "No dose change"              = "#000000", 
  "Dose increase"               = "#D62728",  
  "1-<50% reduction"            = "#2CA02C",  
  ">=50% reduction/initial dose*" = "#FF7F0E"  
)

# Global settings for love plots
x_low    <- -0.6
x_high   <- 0.6
x_breaks <- c(-0.4, -0.2, 0, 0.2, 0.4)

# Helper function to create one love plot
make_love_panel_v2 <- function(dat, bvars, var_labels, title_txt,
                               show_y_labels = TRUE) {
  
  fmla <- as.formula(paste("trt ~", paste(bvars, collapse = " + ")))
  
  prep <- bal.tab(
    fmla,
    data        = dat,
    weights     = dat$wt_trunc,
    estimand    = "ATE",
    m.threshold = 0.1,
    s.d.denom   = "pooled",
    un          = TRUE,
    stats       = "mean.diffs"
  )
  
  n_vars  <- nrow(prep$Balance)
  bg_data <- data.frame(
    var_index  = seq_len(n_vars),
    fill_color = rep(c("gray65", "white"), length.out = n_vars)
  )
  
  plt <- love.plot(
    prep,
    stat         = "mean.diffs",
    stars        = "none",
    thresholds   = c(m = 0.1),
    abs          = FALSE,
    var.names    = var_labels,
    sample.names = c("Before", "After"),
    position     = "none",
    line         = FALSE,
    colors       = c("black", "red"),
    shapes       = c(16, 16),
    title        = title_txt,
    size         = 4
  ) +
    geom_rect(
      data        = bg_data,
      aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
          xmin = -Inf, xmax = Inf, fill = fill_color),
      inherit.aes = FALSE,
      alpha       = 0.2
    ) +
    scale_x_continuous(limits = c(x_low, x_high), breaks = x_breaks, expand = c(0, 0)) +
    scale_fill_identity() +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    labs(x = "SMD", y = "", color = "", shape = "") +
    theme_minimal(base_size = 12) +
    theme(
      title              = element_text(color = "black", size = 13, face = "bold"),
      panel.grid         = element_blank(),
      axis.line          = element_line(color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.text.x        = element_text(color = "black", size = 11),
      axis.ticks.y       = element_line(linewidth = 1, color = "black"),
      axis.text.y        = if (show_y_labels)
        element_text(color = "black", size = 10)
      else
        element_blank(),
      axis.title.x       = element_text(size = 11),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(size = 10, color = "black"),
      legend.key.size    = unit(0.3, "cm"),
      plot.margin        = margin(5, 5, 5, 5, "pt")
    )
  
  plt
}

## MET interruption of 3-4 days #############

cat("\nTrial 2 (3-4 days)...\n")

### First-interruption-in-first-episode #########
cat("  Creating SMD plots: First-interruption-in-first-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

fmla <- as.formula(paste("trt ~", paste(baselinevars, collapse = " + ")))

# Panel 1: 1-<50% vs no change (WITH y-axis labels)
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_1 <- love.plot(prep_1, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% baselinevars],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "1-<50% reduction\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2: >=50% vs no change (NO y-axis labels)
d_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_2 <- love.plot(prep_2, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% baselinevars],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "≥50% reduction/initial dose*\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3: Increase vs no change (NO y-axis labels)
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_3 <- love.plot(prep_3, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% baselinevars],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "Dose increase\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

row_tr2_first <- p_1 + plot_spacer() + p_2 + plot_spacer() + p_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

row_tr2_first <- plot_spacer()/row_tr2_first + plot_layout(heights = c(0.05, 1))

# Save figure for Trial 2, First-interruption-in-first-episode
plot_tr2_first_final <- ggdraw(row_tr2_first) +
  draw_label("A. Methadone interruption of 3-4 days: First-interruption-in-first-episode",
             x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  draw_label("Note: SMD ≥ 0.1 included in outcome models. * ≤40mg methadone",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/", "smd_love_plots_tr2_first.jpg"),
     width = 4200, height = 2200, res = 300)
print(plot_tr2_first_final)
dev.off()

cat("  ✓ Trial 2, First-interruption-in-first-episode SMD figure saved.\n")

### First-interruption-in-any-episode #########
cat("  Creating SMD plots: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds")) 

dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0),
    prior_episode = as.numeric(prior_episode)
  )

fmla <- as.formula(paste("trt ~", paste(c(baselinevars, "prior_episode"), collapse = " + ")))

# Panel 1
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_1 <- love.plot(prep_1, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "1-<50% reduction\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2
d_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_2 <- love.plot(prep_2, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "≥50% reduction/initial dose*\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_3 <- love.plot(prep_3, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "Dose increase\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

row_tr2_all <- p_1 + plot_spacer() + p_2 + plot_spacer() + p_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

row_tr2_all <- plot_spacer() + row_tr2_all + plot_layout(heights = c(0.05, 1))

plot_tr2_all_final <- ggdraw(row_tr2_all) +
  draw_label("B. Methadone interruption of 3-4 days: First-interruption-in-any-episode",
             x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  draw_label("Note: SMD ≥ 0.1 included in outcome models. * ≤40mg methadone",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/", "smd_love_plots_tr2_all.jpg"),
     width = 4200, height = 2200, res = 300)
print(plot_tr2_all_final)
dev.off()

cat("  ✓ Trial 2, First-interruption-in-any-episode SMD figure saved.\n")


## MET interruption of 5-14 days #############

cat("\nTrial 3 (5-14 days)...\n")

### First-interruption-in-first-episode #########
cat("  Creating SMD plots: First-interruption-in-first-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

fmla <- as.formula(paste("trt ~", paste(baselinevars, collapse = " + ")))

# Panel 1: 1-<50% vs no change (WITH y-axis labels)
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_1 <- love.plot(prep_1, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% baselinevars],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "1-<50% reduction\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2: >=50% vs no change (NO y-axis labels)
d_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_2 <- love.plot(prep_2, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% baselinevars],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "≥50% reduction/initial dose*\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3: Increase vs no change (NO y-axis labels)
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_3 <- love.plot(prep_3, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% baselinevars],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "Dose increase\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

row_tr3_first <- p_1 + plot_spacer() + p_2 + plot_spacer() + p_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

row_tr3_first <- plot_spacer()/row_tr3_first + plot_layout(heights = c(0.05, 1))

# Save figure for Trial 2, First-interruption-in-first-episode
plot_tr3_first_final <- ggdraw(row_tr3_first) +
  draw_label("A. Methadone interruption of 5-14 days: First-interruption-in-first-episode",
             x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  draw_label("Note: SMD ≥ 0.1 included in outcome models. * ≤40mg methadone",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/", "smd_love_plots_tr3_first.jpg"),
     width = 4200, height = 2200, res = 300)
print(plot_tr3_first_final)
dev.off()

cat("  ✓ Trial 2, First-interruption-in-first-episode SMD figure saved.\n")

### First-interruption-in-any-episode #########
cat("  Creating SMD plots: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds")) 

dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0),
    prior_episode = as.numeric(prior_episode)
  )

fmla <- as.formula(paste("trt ~", paste(c(baselinevars, "prior_episode"), collapse = " + ")))

# Panel 1
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_1 <- love.plot(prep_1, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "1-<50% reduction\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2
d_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_2 <- love.plot(prep_2, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "≥50% reduction/initial dose*\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$wt_trunc,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_3 <- love.plot(prep_3, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "Dose increase\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

row_tr3_all <- p_1 + plot_spacer() + p_2 + plot_spacer() + p_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

row_tr3_all <- plot_spacer() + row_tr3_all + plot_layout(heights = c(0.05, 1))

plot_tr3_all_final <- ggdraw(row_tr3_all) +
  draw_label("B. Methadone interruption of 5-14 days: First-interruption-in-any-episode",
             x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  draw_label("Note: SMD ≥ 0.1 included in outcome models. * ≤40mg methadone",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/", "smd_love_plots_tr3_all.jpg"),
     width = 4200, height = 2200, res = 300)
print(plot_tr3_all_final)
dev.off()

cat("  ✓ Trial 2, First-interruption-in-any-episode SMD figure saved.\n")



# 6. CUMULATIVE INCIDENCE ####
# Risk curves and confidence intervals for overdose and discontinuation
# Two trials × Two cohorts

cat("\n=== Computing Cumulative Incidence ===\n")

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- c(30, 60, 90, 180, 365)

# ── INITIALIZE RESULTS STRUCTURE ──

risk <- list(
  tr2 = list(
    first = list(od = list(), dis = list()),
    all = list(od = list(), dis = list())
  ),
  tr3 = list(
    first = list(od = list(), dis = list()),
    all = list(od = list(), dis = list())
  )
)

## MET interruption of 3-4 days #############

cat("\nTrial 2 (3-4 days)...\n")

### FIRST-INTERRUPTION-IN-FIRST-EPISODE ####
cat("  Computing: First-interruption-in-first-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)


risk$tr2$first$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr2$first$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr2$first$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr2$first$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

### FIRST-INTERRUPTION-IN-ANY-EPISODE ####
cat("  Computing: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$all$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$all$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"), 
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

risk$tr2$all$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs =c("ddose_prev", "prior_episode"),  
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$all$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_34d.rds"))

cat("  ✓ Trial 2 risk estimates saved.\n")

## MET interruption of 5-14 days #############

cat("\nTrial 3 (5-14 days)...\n")

### FIRST-INTERRUPTION-IN-FIRST-EPISODE ####
cat("  Computing: First-interruption-in-first-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$first$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$first$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$first$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)


risk$tr3$first$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

### FIRST-INTERRUPTION-IN-ANY-EPISODE ####
cat("  Computing: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$all$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$all$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$all$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$all$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_514d.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")

