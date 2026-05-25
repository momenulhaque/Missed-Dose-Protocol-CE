# Load Packages ##################
# 01_main_analysis.R
# Trials 2 and 3 - methadone interruption 3-4d and 5-14d
# Two cohort designs:
#   first_interruption  (formerly incident users)    - TT_analytic_incid_sens_nw_m data set
#   all_interruptions   (formerly prevalent new users)     - TT_analytic_rprev_sens_nw_m data set
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
TRT_LEVELS <- c("no reversion", "reduction", "increase")

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
  prior_gap     = "Total OAT interruptions since initiation"
)

# Shared pre-processing: collapse all reduction categories into one
clean_trt <- function(dat) {
  dat |>
    mutate(
      trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%", ">=50%/<=40mg"), "reduction", trt_s456)
    ) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS))
}

clean_trt_all <- function(dat) {
  dat |>
    mutate(
      trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%", ">=50%/<=40mg"), "reduction", trt_s456)
    ) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS)) |>
    group_by(studyid_rev, episode) |>
    mutate(
      prior_gap     = cumsum(gap_days),
      prior_gap_num = row_number()
    ) %>%
    ungroup()
}

# 1. LOAD AND FILTER RAW SAS DATA #########

# First-interruption cohort (incid)
raw_first <- read_sas(file.path(DATA_DIR, "TT_analytic_incid_sens_nw_m.sas7bdat"))

dat_tr2_first <- raw_first |>
  filter(trial_gap %in% 3:4) |>
  clean_trt()

dat_tr3_first <- raw_first |>
  filter(trial_cat == "5-14") |>
  clean_trt()

rm(raw_first); gc()

# All-interruptions cohort (rprev)
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_rprev_sens_nw_m.sas7bdat"))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all()

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all()

rm(raw_all); gc()

# 2. PROPENSITY SCORE WEIGHTS ################

## MET 3-4 days: first-interruption #####

ps_tr2_first <- ps_model_wt(
  dat             = dat_tr2_first,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_first, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"))

## MET 3-4 days: all-interruptions #####
ps_tr2_all <- ps_model_wt(
  dat             = dat_tr2_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_gap"),
  continuous_vars = c(continuous_vars, "prior_gap"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)
saveRDS(ps_tr2_all, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"))

## MET 5-14 days: first-interruption #####
ps_tr3_first <- ps_model_wt(
  dat             = dat_tr3_first,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)
saveRDS(ps_tr3_first, paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"))

## MET 5-14 days: all-interruptions #####
ps_tr3_all <- ps_model_wt(
  dat             = dat_tr3_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_gap"),
  continuous_vars = c(continuous_vars, "prior_gap"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)
saveRDS(ps_tr3_all, paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"))


# 3. TABLE 1 ##############
# Produces one table per trial × cohort.
# Rows: covariates. Columns: treatment arms (unweighted) + SMD before/after weighting.


tbl_tr2_first <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"),
  bvars = baselinevars, cont_vars = continuous_vars, clustered = FALSE
)

tbl_tr2_all <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"),
  bvars = c(baselinevars, "prior_gap"),
  cont_vars = c(continuous_vars, "prior_gap"), clustered = TRUE
)

tbl_tr3_first <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"),
  bvars = baselinevars, cont_vars = continuous_vars, clustered = FALSE
)

tbl_tr3_all <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"),
  bvars = c(baselinevars, "prior_gap"),
  cont_vars = c(continuous_vars, "prior_gap"), clustered = TRUE
)

library(openxlsx)

wb <- createWorkbook()

for (nm in c("tr2_first", "tr2_all", "tr3_first", "tr3_all")) {
  tbl <- get(paste0("tbl_", nm))
  addWorksheet(wb, sheetName = nm)
  writeData(wb, sheet = nm, x = as.data.frame(tbl), rowNames = TRUE)
}

saveWorkbook(wb, paste0(RES_DIR, "Tables/", "table1_all.xlsx"), overwrite = TRUE)

# 4. PS DISTRIBUTION PLOTS ####
# Histogram + density of propensity scores by treatment arm,
# faceted by interruption length × cohort design.

COMPARISON_COLORS <- c(
  "No dose change" = "#440154",
  "Dose reduction" = "black",
  "Dose increase"  = "#D62728"
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
                     "reduction"    = "Dose reduction",
                     "increase"     = "Dose increase"),
      group  = factor(group, levels = c("No dose change", "Dose reduction", "Dose increase")),
      Trial  = trial_label,
      Design = design_label
    ) |>
    select(group, ps, Trial, Design)
}


dat_ps <- bind_rows(
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"),
              "MET interruption of 3-4 days",  "First-interruption cohort"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"),
              "MET interruption of 3-4 days",  "All-interruptions cohort"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"),
              "MET interruption of 5-14 days", "First-interruption cohort"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"),
              "MET interruption of 5-14 days", "All-interruptions cohort")
) |>
  mutate(
    Trial_design = factor(
      paste0(Trial, "\n", Design),
      levels = c(
        "MET interruption of 3-4 days\nFirst-interruption cohort",
        "MET interruption of 3-4 days\nAll-interruptions cohort",
        "MET interruption of 5-14 days\nFirst-interruption cohort",
        "MET interruption of 5-14 days\nAll-interruptions cohort"
      )
    )
  )


plot_ps <- ps_score_plot(dat_ps) +
  guides(color = guide_legend(nrow = 1))

jpeg(paste0(RES_DIR, "Figures/", "ps_distribution.jpg"),
     width = 4400, height = 3000, res = 300)
plot_ps
dev.off()


# 5. SMD LOVE PLOTS ############
# One figure per trial, two rows (first_interruption / all_interruptions),
# two panels per row (reduction vs no change / increase vs no change).

x_low    <- -0.6;  x_high   <- 0.6
x_breaks <- c(-0.4, -0.2, 0, 0.2, 0.4)


## MET 3-4 days love plot (first / all) #####

row_tr2_first <- make_love_row(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"),
  bvars       = baselinevars,
  var_labels  = VAR_LABELS[names(VAR_LABELS) != "prior_gap"]
)

row_tr2_all <- make_love_row(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"),
  bvars      = c(baselinevars, "prior_gap"),
  var_labels = VAR_LABELS
)

plot_tr2 <- ggdraw(
  plot_spacer() / row_tr2_first / plot_spacer() / row_tr2_all +
    plot_layout(heights = c(0.06, 1, 0.04, 1), guides = "collect") &
    theme(legend.position = "bottom")
) +
  draw_label("A. Methadone interruption of 3-4 days: first-interruption cohort",
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, size = 16) +
  draw_label("B. Methadone interruption of 3-4 days: all-interruptions cohort",
             x = 0.02, y = 0.485, hjust = 0, vjust = 0, size = 16) 

jpeg(paste0(RES_DIR, "Figures/", "smd_tr2_34d.jpg"),
    width = 3800, height = 4400, res = 300)
plot_tr2
dev.off()

## MET 5-14 days love plot (first / all) #####

row_tr3_first <- make_love_row(
  paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"),
  bvars      = baselinevars,
  var_labels = VAR_LABELS[names(VAR_LABELS) != "prior_gap"]
)

row_tr3_all <- make_love_row(
  paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"),
  bvars      = c(baselinevars, "prior_gap"),
  var_labels = VAR_LABELS
)

plot_tr3 <- ggdraw(
  plot_spacer() / row_tr3_first / plot_spacer() / row_tr3_all +
    plot_layout(heights = c(0.06, 1, 0.04, 1), guides = "collect") &
    theme(legend.position = "bottom")
) +
  draw_label("A. Methadone 5-14 days — first-interruption cohort",
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, size = 16) +
  draw_label("B. Methadone 5-14 days — all-interruptions cohort",
             x = 0.02, y = 0.485, hjust = 0, vjust = 0, size = 16) +
  draw_label("Note: Confounders with SMD \u2265 0.1 were included in the outcome models as covariates.",
             x = 0.02, y = 0.01, hjust = 0, vjust = 0, size = 12)

jpeg(paste0(RES_DIR, "Figures/", "smd_tr3.jpg"),
     width = 3800, height = 4400, res = 300)
plot_tr3
dev.off()


# 6. CUMULATIVE INCIDENCE ####
HORIZONS_CURVE <- c(seq(7, 357, by = 7), 365)
HORIZONS_CI    <- c(30, 60, 90, 180, 365)

risk <- list(
  tr2 = list(first = list(od = list(), dis = list()),
             all   = list(od = list(), dis = list())),
  tr3 = list(first = list(od = list(), dis = list()),
             all   = list(od = list(), dis = list()))
)

## MET 3-4 days / first-interruption ######
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"))

## Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 7)

risk$tr2$first$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr2$first$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

## Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 7)

risk$tr2$first$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr2$first$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

## MET 3-4 days / all-interruptions #######
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"))

## Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 7)

risk$tr2$all$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "prior_gap"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 1000, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt_trunc"
)

risk$tr2$all$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "prior_gap"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 1000, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt_trunc"
)

## Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 7)

risk$tr2$all$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr2$all$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_34d.rds"))



## MET 5-14 days / first-interruption ######
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"))

## Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 7)

risk$tr3$first$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ret_length", "ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$first$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ret_length", "ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

## Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 7)

risk$tr3$first$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ret_length", "ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$first$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ret_length", "ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

##  MET 5-14 days / all-interruptions ######
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"))

## Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 7)

risk$tr3$all$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("prior_gap"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 99, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$all$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("prior_gap"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 99, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

## Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 7)

risk$tr3$all$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("prior_gap"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 99, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$all$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("prior_gap"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 99, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3.rds"))
