# 02_stratified_analysis.R
# Stratified analysis by baseline retention length (ret_length)
# FIRST-IN-EPISODE (prev) cohort
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
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results_v1/"
dir.create(paste0(RES_DIR, "PScore/Stratified"),   showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Risk/Stratified"),     showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Figures/Stratified"),  showWarnings = FALSE)

# ── Shared settings ────────────────────────────────────────────────────────────
TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")

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
clean_trt_first_ep <- function(dat) {
  dat |>
    mutate(
      trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)
    ) |>    
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS))
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

cat("Loading and stratifying first-in-episode cohort data...\n")

raw_first_ep <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_first_ep_strat <- raw_first_ep |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_first_ep() |>
  add_strata()

dat_tr3_first_ep_strat <- raw_first_ep |>
  filter(trial_cat == "5-14") |>
  clean_trt_first_ep() |>
  add_strata()

rm(raw_first_ep); gc()

# Print stratum sizes
cat("\n=== Sample sizes by stratum ===\n")
cat("\nTrial 2 (3-4 days), First-in-episode:\n")
print(table(dat_tr2_first_ep_strat$ret_strat, dat_tr2_first_ep_strat$trt))

cat("\nTrial 3 (5-14 days), First-in-episode:\n")
print(table(dat_tr3_first_ep_strat$ret_strat, dat_tr3_first_ep_strat$trt))


# 2. PROPENSITY SCORE MODELS ###################

ps_strat <- list(
  tr2 = list(first_ep = list()),
  tr3 = list(first_ep = list())
)

## Trial 2 (3-4 days), First-in-episode #################
cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("TRIAL 2 (3-4 days): Fitting PS models per stratum (standard MLE)\n")
cat("First-in-episode cohort\n")
cat("════════════════════════════════════════════════════════════════\n")

for (stratum in STRATA_LEVELS) {
  cat(sprintf("\n--- Stratum: %s ---\n", STRATA_LABELS[stratum]))
  
  dat_strat <- dat_tr2_first_ep_strat %>% filter(ret_strat == stratum)
  cat(sprintf("Total N = %d\n", nrow(dat_strat)))
  cat("Treatment distribution:\n")
  print(table(dat_strat$trt))
  
  ps_strat$tr2$first_ep[[stratum]] <- ps_model_wt(
    dat             = dat_strat,
    exposure        = "trt",
    baselinevars    = baselinevars,
    continuous_vars = continuous_vars,
    out_vars        = out_vars,
    clust_vars      = clust_vars
  )
  
  # Save
  saveRDS(ps_strat$tr2$first_ep[[stratum]], 
          paste0(RES_DIR, "PScore/Stratified/ps_tr2_34d_first_ep_", stratum, ".rds"))
  
  # Check weights
  wt_data <- ps_strat$tr2$first_ep[[stratum]]$data_weights
  cat("\nWeight distribution:\n")
  print(summary(wt_data$wt_stab))
  cat(sprintf("Weights > 50: %d (%.1f%%)\n", 
              sum(wt_data$wt_stab > 50, na.rm = TRUE), 
              100 * mean(wt_data$wt_stab > 50, na.rm = TRUE)))
  cat(sprintf("PS outside [0.01, 0.99]: %d (%.1f%%)\n",
              sum(wt_data$ps < 0.01 | wt_data$ps > 0.99, na.rm = TRUE),
              100 * mean(wt_data$ps < 0.01 | wt_data$ps > 0.99, na.rm = TRUE)))
}


## Trial 3 (5-14 days), First-in-episode #######
cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("TRIAL 3 (5-14 days): Fitting PS models per stratum (standard MLE)\n")
cat("First-in-episode cohort\n")
cat("════════════════════════════════════════════════════════════════\n")

for (stratum in STRATA_LEVELS) {
  cat(sprintf("\n--- Stratum: %s ---\n", STRATA_LABELS[stratum]))
  
  dat_strat <- dat_tr3_first_ep_strat %>% filter(ret_strat == stratum)
  cat(sprintf("Total N = %d\n", nrow(dat_strat)))
  cat("Treatment distribution:\n")
  print(table(dat_strat$trt))
  
  ps_strat$tr3$first_ep[[stratum]] <- ps_model_wt(
    dat             = dat_strat,
    exposure        = "trt",
    baselinevars    = baselinevars,
    continuous_vars = continuous_vars,
    out_vars        = out_vars,
    clust_vars      = clust_vars
  )
  
  # Save
  saveRDS(ps_strat$tr3$first_ep[[stratum]], 
          paste0(RES_DIR, "PScore/Stratified/ps_tr3_first_ep_", stratum, ".rds"))
  
  # Check weights
  wt_data <- ps_strat$tr3$first_ep[[stratum]]$data_weights
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
  "No dose change"                = "#000000",    # Black
  "1-<50% reduction"              = "#2CA02C",    # Green
  ">=50% reduction/initial dose*" = "#FF7F0E",    # Orange
  "Dose increase"                 = "#D62728"     # Red
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
      legend.text = element_text(size = 10),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold", size = 11)
    )
}

prep_ps_strat_dat <- function(strata_list) {
  map_dfr(STRATA_LEVELS, function(stratum) {
    strata_list[[stratum]]$data_weights |>
      mutate(
        group = recode(as.character(trt),
                       "no reversion"       = "No dose change",
                       "1-<50%"             = "1-<50% reduction",
                       ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                       "increase"           = "Dose increase"),
        group = factor(group, levels = c("No dose change", "1-<50% reduction", 
                                         ">=50% reduction/initial dose*", "Dose increase")),
        Stratum = STRATA_LABELS[stratum],
        Stratum = factor(Stratum, levels = STRATA_LABELS)
      ) |>
      select(group, ps, Stratum)
  })
}

## Trial 2
cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("Creating PS distribution plots: Trial 2 (stratified)\n")
cat("════════════════════════════════════════════════════════════════\n")

dat_ps_tr2_first_ep <- prep_ps_strat_dat(ps_strat$tr2$first_ep)
plot_ps_tr2_first_ep <- ps_score_plot_strat(
  dat_ps_tr2_first_ep, 
  "MET interruption of 3-4 days: PS Distribution"
) + guides(color = guide_legend(nrow = 2))

jpeg(paste0(RES_DIR, "Figures/Stratified/ps_dist_tr2_first_ep_strat.jpg"),
     width = 4200, height = 2000, res = 300)
print(plot_ps_tr2_first_ep)
dev.off()

cat("✓ Trial 2 PS distribution plot saved\n")

## Trial 3
cat("\nCreating PS distribution plots: Trial 3 (stratified)\n")

dat_ps_tr3_first_ep <- prep_ps_strat_dat(ps_strat$tr3$first_ep)
plot_ps_tr3_first_ep <- ps_score_plot_strat(
  dat_ps_tr3_first_ep, 
  "MET interruption of 5-14 days: PS Distribution"
) + guides(color = guide_legend(nrow = 2))

jpeg(paste0(RES_DIR, "Figures/Stratified/ps_dist_tr3_first_ep_strat.jpg"),
     width = 4200, height = 2000, res = 300)
print(plot_ps_tr3_first_ep)
dev.off()

cat("✓ Trial 3 PS distribution plot saved\n")


# 4. CUMULATIVE INCIDENCE (Stratified) ###########

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- c(30, 60, 90, 180, 365)

risk_strat <- list(
  tr2 = list(first_ep = list()),
  tr3 = list(first_ep = list())
)

## ── Trial 2 (3-4 days) ───────────────────────────────────────────────────────
cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("Computing cumulative incidence: Trial 2 (3-4 days)\n")
cat("First-in-episode cohort, stratified by retention length\n")
cat("════════════════════════════════════════════════════════════════\n")

for (stratum in STRATA_LEVELS) {
  cat(sprintf("\n--- Stratum: %s ---\n", STRATA_LABELS[stratum]))
  
  ps_fit <- ps_strat$tr2$first_ep[[stratum]]
  risk_strat$tr2$first_ep[[stratum]] <- list(od = list(), dis = list())
  
  ## Overdose
  cat("  Computing overdose risk...\n")
  long_dat <- expand_person_period(
    ps_fit$data_weights,
    time_var = "time_od_acute",
    event_var = "sens_od_acute",
    interval_length = 14
  )
  
  risk_strat$tr2$first_ep[[stratum]]$od$curve <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ddose_prev", "ret_length"),
    horizons = HORIZONS_CURVE,
    cluster_var = "studyid_rev",
    samples = 199,
    conf_int = FALSE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  risk_strat$tr2$first_ep[[stratum]]$od$ci <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ddose_prev", "ret_length"),
    horizons = HORIZONS_CI,
    cluster_var = "studyid_rev",
    samples = 199,
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
  
  risk_strat$tr2$first_ep[[stratum]]$dis$curve <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ret_length"),
    horizons = HORIZONS_CURVE,
    cluster_var = "studyid_rev",
    samples = 199,
    conf_int = FALSE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  risk_strat$tr2$first_ep[[stratum]]$dis$ci <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ret_length"),
    horizons = HORIZONS_CI,
    cluster_var = "studyid_rev",
    samples = 199,
    conf_int = TRUE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  rm(long_dat); gc()
}

saveRDS(risk_strat$tr2$first_ep, paste0(RES_DIR, "Risk/Stratified/risk_tr2_34d_first_ep_strat.rds"))
cat("\n✓ Trial 2 risk estimates saved\n")


## ── Trial 3 (5-14 days) ──────────────────────────────────────────────────────
cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("Computing cumulative incidence: Trial 3 (5-14 days)\n")
cat("First-in-episode cohort, stratified by retention length\n")
cat("════════════════════════════════════════════════════════════════\n")

for (stratum in STRATA_LEVELS) {
  cat(sprintf("\n--- Stratum: %s ---\n", STRATA_LABELS[stratum]))
  
  ps_fit <- ps_strat$tr3$first_ep[[stratum]]
  risk_strat$tr3$first_ep[[stratum]] <- list(od = list(), dis = list())
  
  ## Overdose
  cat("  Computing overdose risk...\n")
  long_dat <- expand_person_period(
    ps_fit$data_weights,
    time_var = "time_od_acute",
    event_var = "sens_od_acute",
    interval_length = 14
  )
  
  risk_strat$tr3$first_ep[[stratum]]$od$curve <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ddose_prev", "ret_length"),
    horizons = HORIZONS_CURVE,
    cluster_var = "studyid_rev",
    samples = 199,
    conf_int = FALSE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  risk_strat$tr3$first_ep[[stratum]]$od$ci <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ddose_prev", "ret_length"),
    horizons = HORIZONS_CI,
    cluster_var = "studyid_rev",
    samples = 199,
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
  
  risk_strat$tr3$first_ep[[stratum]]$dis$curve <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ret_length"),
    horizons = HORIZONS_CURVE,
    cluster_var = "studyid_rev",
    samples = 199,
    conf_int = FALSE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  risk_strat$tr3$first_ep[[stratum]]$dis$ci <- predict_cuminc(
    Dat_expanded = long_dat,
    trt_var = "trt",
    covs = c("ret_length"),
    horizons = HORIZONS_CI,
    cluster_var = "studyid_rev",
    samples = 199,
    conf_int = TRUE,
    ref_group = "no reversion",
    wt_type = "wt_trunc"
  )
  
  rm(long_dat); gc()
}

saveRDS(risk_strat$tr3$first_ep, paste0(RES_DIR, "Risk/Stratified/risk_tr3_first_ep_strat.rds"))
cat("\n✓ Trial 3 risk estimates saved\n")

# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS COMPLETE
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("STRATIFIED ANALYSIS COMPLETE\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("\nCohort: First-in-episode (prev)\n")
cat("Stratification variable: Retention length at baseline\n")
cat("  • ≤30 days (early treatment initiates)\n")
cat("  • 31-90 days (mid-treatment)\n")
cat("  • 91+ days (established patients)\n")
cat("\nResults saved to:\n")
cat(sprintf("  PS models:     %sPScore/Stratified/\n", RES_DIR))
cat(sprintf("  Risk estimates: %sRisk/Stratified/\n", RES_DIR))
cat(sprintf("  Figures:       %sFigures/Stratified/\n", RES_DIR))
cat("\nFiles created:\n")
cat("  ✓ 6 PS model files (2 trials × 3 strata)\n")
cat("  ✓ 2 PS distribution figures (stratified)\n")
cat("  ✓ 2 risk estimate RDS files (tr2, tr3 × 3 strata each)\n")
cat("════════════════════════════════════════════════════════════════\n")
