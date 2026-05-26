# 02_stratified_analysis.R
# Stratified analysis by baseline retention length (ret_length)
# ALL-INTERRUPTIONS COHORT ONLY
# Strata: ≤30 days, 31-90 days, 91+ days
# Uses STANDARD MLE (ps_model_wt) - adequate sample sizes in all strata

library(tidyverse)
library(haven)
library(nnet)
library(splines)
library(sandwich)
library(mvtnorm)
library(patchwork)
library(cowplot)

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"
dir.create(paste0(RES_DIR, "PScore/Stratified"),   showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Risk/Stratified"),     showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Figures/Stratified"),  showWarnings = FALSE)

# ── Shared settings ────────────────────────────────────────────────────────────
TRT_LEVELS <- c("no reversion", "reduction", "increase")

out_vars   <- c("sens_discont", "time_discont",
                "sens_od_acute", "time_od_acute",
                "sens_all_acute", "time_all_acute")

clust_vars <- c("studyid_rev", "episode", "time0")

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length",  # ← continuous within strata
  "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss"
)
continuous_vars <- c("ddose_prev", "year", "ret_length")

# Stratification levels
STRATA_LEVELS <- c("le30", "31_90", "91plus")
STRATA_LABELS <- c(
  "le30"   = "≤30 days",
  "31_90"  = "31-90 days",
  "91plus" = "91+ days"
)

# Pre-processing function
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

# Add stratification variable
add_strata <- function(dat) {
  dat |>
    mutate(
      ret_strat = case_when(
        ret_length <= 30 ~ "le30",
        ret_length <= 90 ~ "31_90",
        TRUE ~ "91plus"
      ),
      ret_strat = factor(ret_strat, levels = STRATA_LEVELS)
    )
}


# 1. LOAD AND STRATIFY DATA ##############

cat("Loading and stratifying all-interruptions cohort data...\n")

raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_rprev_sens_nw_m.sas7bdat"))

dat_tr2_all_strat <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all() |>
  add_strata()

dat_tr3_all_strat <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all() |>
  add_strata()

rm(raw_all); gc()

# Print stratum sizes
cat("\n=== Sample sizes by stratum ===\n")
cat("\nTrial 2 (3-4 days), All-interruptions:\n")
print(table(dat_tr2_all_strat$ret_strat, dat_tr2_all_strat$trt))

cat("\nTrial 3 (5-14 days), All-interruptions:\n")
print(table(dat_tr3_all_strat$ret_strat, dat_tr3_all_strat$trt))


# 2. PROPENSITY SCORE MODELS ###################

ps_strat <- list(
  tr2 = list(all = list()),
  tr3 = list(all = list())
)

## Trial 2 (3-4 days), All-interruptions #################
cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("TRIAL 2 (3-4 days): Fitting PS models per stratum (standard MLE)\n")
cat("════════════════════════════════════════════════════════════════\n")

for (stratum in STRATA_LEVELS) {
  cat(sprintf("\n--- Stratum: %s ---\n", STRATA_LABELS[stratum]))
  
  dat_strat <- dat_tr2_all_strat %>% filter(ret_strat == stratum)
  cat(sprintf("Total N = %d\n", nrow(dat_strat)))
  cat("Treatment distribution:\n")
  print(table(dat_strat$trt))
  
  ps_strat$tr2$all[[stratum]] <- ps_model_wt(
    dat             = dat_strat,
    exposure        = "trt",
    baselinevars    = c(baselinevars, "prior_gap"),  # includes ret_length + prior_gap
    continuous_vars = c(continuous_vars, "prior_gap"),
    out_vars        = out_vars,
    clust_vars      = clust_vars
  )
  
  # Save
  saveRDS(ps_strat$tr2$all[[stratum]], 
          paste0(RES_DIR, "PScore/Stratified/ps_tr2_all_", stratum, ".rds"))
  
  # Check weights
  wt_data <- ps_strat$tr2$all[[stratum]]$data_weights
  cat("\nWeight distribution:\n")
  print(summary(wt_data$wt_stab))
  cat(sprintf("Weights > 50: %d (%.1f%%)\n", 
              sum(wt_data$wt_stab > 50, na.rm = TRUE), 
              100 * mean(wt_data$wt_stab > 50, na.rm = TRUE)))
  cat(sprintf("PS outside [0.01, 0.99]: %d (%.1f%%)\n",
              sum(wt_data$ps < 0.01 | wt_data$ps > 0.99, na.rm = TRUE),
              100 * mean(wt_data$ps < 0.01 | wt_data$ps > 0.99, na.rm = TRUE)))
}


## Trial 3 (5-14 days), All-interruptions #######
cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("TRIAL 3 (5-14 days): Fitting PS models per stratum (standard MLE)\n")
cat("════════════════════════════════════════════════════════════════\n")

for (stratum in STRATA_LEVELS) {
  cat(sprintf("\n--- Stratum: %s ---\n", STRATA_LABELS[stratum]))
  
  dat_strat <- dat_tr3_all_strat %>% filter(ret_strat == stratum)
  cat(sprintf("Total N = %d\n", nrow(dat_strat)))
  cat("Treatment distribution:\n")
  print(table(dat_strat$trt))
  
  ps_strat$tr3$all[[stratum]] <- ps_model_wt(
    dat             = dat_strat,
    exposure        = "trt",
    baselinevars    = c(baselinevars, "prior_gap"),
    continuous_vars = c(continuous_vars, "prior_gap"),
    out_vars        = out_vars,
    clust_vars      = clust_vars
  )
  
  # Save
  saveRDS(ps_strat$tr3$all[[stratum]], 
          paste0(RES_DIR, "PScore/Stratified/ps_tr3_all_", stratum, ".rds"))
  
  # Check weights
  wt_data <- ps_strat$tr3$all[[stratum]]$data_weights
  cat("\nWeight distribution:\n")
  print(summary(wt_data$wt_stab))
  cat(sprintf("Weights > 50: %d (%.1f%%)\n", 
              sum(wt_data$wt_stab > 50, na.rm = TRUE), 
              100 * mean(wt_data$wt_stab > 50, na.rm = TRUE)))
  cat(sprintf("PS outside [0.01, 0.99]: %d (%.1f%%)\n",
              sum(wt_data$ps < 0.01 | wt_data$ps > 0.99, na.rm = TRUE),
              100 * mean(wt_data$ps < 0.01 | wt_data$ps > 0.99, na.rm = TRUE)))
}


# 3. PS DISTRIBUTION PLOTS (Stratified) #########

COMPARISON_COLORS <- c(
  "No dose change" = "#440154",
  "Dose reduction" = "black",
  "Dose increase"  = "#D62728"
)

ps_score_plot_strat <- function(dat, title = "") {
  ggplot(dat, aes(x = ps, colour = group, fill = group)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.25, bins = 30) +
    geom_density(alpha = 0.05, size = 0.5) +
    scale_colour_manual(values = COMPARISON_COLORS) +
    scale_fill_manual(values = COMPARISON_COLORS) +
    facet_wrap(~Stratum, nrow = 1, scales = "free") +
    labs(x = "Propensity Score", y = "Density", title = title) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold", size = 11)
    )
}

prep_ps_strat_dat <- function(strata_list) {
  map_dfr(STRATA_LEVELS, function(stratum) {
    strata_list[[stratum]]$data_weights |>
      mutate(
        group = recode(as.character(trt),
                       "no reversion" = "No dose change",
                       "reduction" = "Dose reduction",
                       "increase" = "Dose increase"),
        group = factor(group, levels = c("No dose change", "Dose reduction", "Dose increase")),
        Stratum = STRATA_LABELS[stratum],
        Stratum = factor(Stratum, levels = STRATA_LABELS)
      ) |>
      select(group, ps, Stratum)
  })
}

## Trial 2
dat_ps_tr2_all <- prep_ps_strat_dat(ps_strat$tr2$all)
plot_ps_tr2_all <- ps_score_plot_strat(
  dat_ps_tr2_all, 
  "Trial 2 (3-4 days): PS Distribution by Retention Length"
) + guides(color = guide_legend(nrow = 1))

jpeg(paste0(RES_DIR, "Figures/Stratified/ps_dist_tr2_all_strat.jpg"),
     width = 3800, height = 1600, res = 300)
print(plot_ps_tr2_all)
dev.off()

## Trial 3
dat_ps_tr3_all <- prep_ps_strat_dat(ps_strat$tr3$all)
plot_ps_tr3_all <- ps_score_plot_strat(
  dat_ps_tr3_all, 
  "Trial 3 (5-14 days): PS Distribution by Retention Length"
) + guides(color = guide_legend(nrow = 1))

jpeg(paste0(RES_DIR, "Figures/Stratified/ps_dist_tr3_all_strat.jpg"),
     width = 3800, height = 1600, res = 300)
print(plot_ps_tr3_all)
dev.off()


# 4. CUMULATIVE INCIDENCE (Stratified) ###########

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- c(30, 60, 90, 180, 365)

risk_strat <- list(
  tr2 = list(all = list()),
  tr3 = list(all = list())
)

## ── Trial 2 (3-4 days) ───────────────────────────────────────────────────────
cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("Computing cumulative incidence: Trial 2 (3-4 days)\n")
cat("════════════════════════════════════════════════════════════════\n")

for (stratum in STRATA_LEVELS) {
  cat(sprintf("\n--- Stratum: %s ---\n", STRATA_LABELS[stratum]))
  
  ps_fit <- ps_strat$tr2$all[[stratum]]
  risk_strat$tr2$all[[stratum]] <- list(od = list(), dis = list())
  
  ## Overdose
  cat("  Computing overdose risk...\n")
  long_dat <- expand_person_period(
    ps_fit$data_weights,
    time_var = "time_od_acute",
    event_var = "sens_od_acute",
    interval_length = 14
  )
  
  risk_strat$tr2$all[[stratum]]$od$curve <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ddose_prev", "prior_gap"),
    horizons = HORIZONS_CURVE,
    cluster_var = "studyid_rev",
    samples = 499,
    conf_int = FALSE,
    ref_group = "no reversion",
    wt_type = "ow_wt_trunc"
  )
  
  risk_strat$tr2$all[[stratum]]$od$ci <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ddose_prev", "prior_gap"),
    horizons = HORIZONS_CI,
    cluster_var = "studyid_rev",
    samples = 499,
    conf_int = TRUE,
    ref_group = "no reversion",
    wt_type = "ow_wt_trunc"
  )
  
  ## Discontinuation
  cat("  Computing discontinuation risk...\n")
  long_dat <- expand_person_period(
    ps_fit$data_weights,
    time_var = "time_discont",
    event_var = "sens_discont",
    interval_length = 14
  )
  
  risk_strat$tr2$all[[stratum]]$dis$curve <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = NULL,
    horizons = HORIZONS_CURVE,
    cluster_var = "studyid_rev",
    samples = 499,
    conf_int = FALSE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  risk_strat$tr2$all[[stratum]]$dis$ci <- predict_cuminc(
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
}

saveRDS(risk_strat$tr2$all, paste0(RES_DIR, "Risk/Stratified/risk_tr2_all_strat.rds"))


## ── Trial 3 (5-14 days) ──────────────────────────────────────────────────────
cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("Computing cumulative incidence: Trial 3 (5-14 days)\n")
cat("════════════════════════════════════════════════════════════════\n")

for (stratum in STRATA_LEVELS) {
  cat(sprintf("\n--- Stratum: %s ---\n", STRATA_LABELS[stratum]))
  
  ps_fit <- ps_strat$tr3$all[[stratum]]
  risk_strat$tr3$all[[stratum]] <- list(od = list(), dis = list())
  
  ## Overdose
  cat("  Computing overdose risk...\n")
  long_dat <- expand_person_period(
    ps_fit$data_weights,
    time_var = "time_od_acute",
    event_var = "sens_od_acute",
    interval_length = 14
  )
  
  risk_strat$tr3$all[[stratum]]$od$curve <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("prior_gap"),
    horizons = HORIZONS_CURVE,
    cluster_var = "studyid_rev",
    samples = 499,
    conf_int = FALSE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  risk_strat$tr3$all[[stratum]]$od$ci <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("prior_gap"),
    horizons = HORIZONS_CI,
    cluster_var = "studyid_rev",
    samples = 499,
    conf_int = TRUE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  ## Discontinuation
  cat("  Computing discontinuation risk...\n")
  long_dat <- expand_person_period(
    ps_fit$data_weights,
    time_var = "time_discont",
    event_var = "sens_discont",
    interval_length = 14
  )
  
  risk_strat$tr3$all[[stratum]]$dis$curve <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("prior_gap"),
    horizons = HORIZONS_CURVE,
    cluster_var = "studyid_rev",
    samples = 499,
    conf_int = FALSE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  risk_strat$tr3$all[[stratum]]$dis$ci <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("prior_gap"),
    horizons = HORIZONS_CI,
    cluster_var = "studyid_rev",
    samples = 499,
    conf_int = TRUE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  rm(long_dat); gc()
}
saveRDS(risk_strat$tr3$all, paste0(RES_DIR, "Risk/Stratified/risk_tr3_all_strat.rds"))


