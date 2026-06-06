# 02_stratified_analysis_rmst.R
# Stratified analysis by baseline retention length (ret_length)
# ANY EPISODE cohort only
# Strata: ≤30 days, >30 days

library(tidyverse)
library(haven)
library(nnet)
library(splines)
library(sandwich)
library(mvtnorm)
library(forestploter)
library(grid)

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")

# ── Paths & Settings ──────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"
dir.create(paste0(RES_DIR, "PScore/Stratified"),   showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Risk/Stratified"),     showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Figures/Stratified"),  showWarnings = FALSE)

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
out_vars   <- c("sens_discont", "time_discont", "sens_od_acute", "time_od_acute")
clust_vars <- c("studyid_rev", "episode", "time0")

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length", "prior_episode", # ← Key variables for 'any episode'
  "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss"
)

continuous_vars <- c("ddose_prev", "year", "ret_length", "prior_episode")

# Updated Stratification levels
STRATA_LEVELS <- c("le30", "gt30")
STRATA_LABELS <- c("le30" = "≤30 days", "gt30" = ">30 days")

# Pre-processing function

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



add_strata <- function(dat) {
  dat |>
    mutate(
      ret_strat = case_when(
        ret_length <= 30 ~ "le30",
        TRUE ~ "gt30"
      ),
      ret_strat = factor(ret_strat, levels = STRATA_LEVELS)
    )
}

# 1. LOAD AND STRATIFY DATA ##############
cat("Loading and stratifying 'Any episode' cohort data...\n")

raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_strat <- raw_all |> filter(trial_gap %in% 3:4) |> clean_trt_all() |> add_strata()
dat_tr3_strat <- raw_all |> filter(trial_cat == "5-14") |> clean_trt_all() |> add_strata()
rm(raw_all); gc()

# 2. PROPENSITY SCORE MODELS ###################
ps_strat <- list(tr2 = list(), tr3 = list())

fit_ps_models <- function(dat_stratified, trial_name) {
  res <- list()
  for (stratum in STRATA_LEVELS) {
    cat(sprintf("\nFitting %s - Stratum: %s\n", trial_name, STRATA_LABELS[stratum]))
    dat_sub <- dat_stratified %>% filter(ret_strat == stratum)
    
    res[[stratum]] <- ps_model_wt(
      dat = dat_sub, exposure = "trt", baselinevars = baselinevars,
      continuous_vars = continuous_vars, out_vars = out_vars, clust_vars = clust_vars
    )
    saveRDS(res[[stratum]], paste0(RES_DIR, "PScore/Stratified/ps_", trial_name, "_", stratum, "_ow.rds"))
  }
  return(res)
}

ps_strat$tr2 <- fit_ps_models(dat_tr2_strat, "tr2_34d")
ps_strat$tr3 <- fit_ps_models(dat_tr3_strat, "tr3_514d")

# 3. CUMULATIVE INCIDENCE & RMST ###########
HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))
risk_strat <- list(tr2 = list(), tr3 = list())

compute_risks <- function(ps_list, trial_name) {
  res <- list()
  for (stratum in STRATA_LEVELS) {
    cat(sprintf("\nComputing Risk for %s - Stratum: %s\n", trial_name, STRATA_LABELS[stratum]))
    ps_fit <- ps_list[[stratum]]
    res[[stratum]] <- list(od = list(), dis = list())
    
    # Overdose
    long_od <- expand_person_period(ps_fit$data_weights, "time_od_acute", "sens_od_acute", 14)
    res[[stratum]]$od$ci <- predict_cuminc_rmst(
      Dat_expanded = long_od, trt_var = "trt", covs = c("ddose_prev", "prior_episode"),
      horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, 
      conf_int = TRUE, ref_group = "no reversion", wt_type = "wt_trunc", rmst_tau = 365
    )
    rm(long_od); gc()
    
    # Discontinuation
    long_dis <- expand_person_period(ps_fit$data_weights, "time_discont", "sens_discont", 14)
    res[[stratum]]$dis$ci <- predict_cuminc_rmst(
      Dat_expanded = long_dis, trt_var = "trt", covs = c("ddose_prev", "prior_episode"),
      horizons = HORIZONS_CI, cluster_var = "studyid_rev", samples = 499, 
      conf_int = TRUE, ref_group = "no reversion", wt_type = "wt_trunc", rmst_tau = 365
    )
    rm(long_dis); gc()
  }
  return(res)
}

risk_strat$tr2 <- compute_risks(ps_strat$tr2, "tr2")
risk_strat$tr3 <- compute_risks(ps_strat$tr3, "tr3")

saveRDS(risk_strat$tr2, paste0(RES_DIR, "Risk/Stratified/risk_tr2_strat_rmst.rds"))
saveRDS(risk_strat$tr3, paste0(RES_DIR, "Risk/Stratified/risk_tr3_strat_rmst.rds"))

# 4. EXTRACTION FUNCTIONS ###################

extract_rd_data <- function(risk_obj, dat_weights, max_time, outcome = "od", 
                            ref_name = "no reversion", stratum_label, int_label) {
  
  time_v  <- if(outcome == "dis") "time_discont" else "time_od_acute"
  event_v <- if(outcome == "dis") "sens_discont" else "sens_od_acute"
  
  counts <- compute_summary_stats_timewin(
    dat_weights, time_var = time_v, event_var = event_v, 
    weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
  ) %>% dplyr::select(trt_group = trt, Episode)
  
  risks <- risk_obj$risks %>% filter(horizon == max_time) %>%
    mutate(
      risk_pct = risk * 100, risk_lower_pct = risk_lower * 100, risk_upper_pct = risk_upper * 100,
      abs_risk_disp = ifelse(abs(risk_pct) >= 10,
                             sprintf("%.1f (%.1f, %.1f)", risk_pct, risk_lower_pct, risk_upper_pct),
                             sprintf("%.2f (%.2f, %.2f)", risk_pct, risk_lower_pct, risk_upper_pct))
    ) %>% dplyr::select(trt_group = trt, abs_risk_disp)
  
  contrasts <- risk_obj$contrasts %>% filter(horizon == max_time) %>%
    mutate(
      rd_pct = RD * 100, rd_lower_pct = RD_lower * 100, rd_upper_pct = RD_upper * 100,
      rd_disp = ifelse(abs(rd_pct) >= 10,
                       sprintf("%.1f (%.1f, %.1f)", rd_pct, rd_lower_pct, rd_upper_pct),
                       sprintf("%.2f (%.2f, %.2f)", rd_pct, rd_lower_pct, rd_upper_pct))
    ) %>% dplyr::select(trt_group = trt, rd_pct, rd_lower_pct, rd_upper_pct, rd_disp)
  
  counts %>%
    left_join(risks, by = "trt_group") %>%
    left_join(contrasts, by = "trt_group") %>%
    mutate(
      rd_disp = if_else(trt_group == ref_name, "Ref.", rd_disp),
      rd_pct = if_else(trt_group == ref_name, 0, rd_pct),
      rd_lower_pct = if_else(trt_group == ref_name, 0, rd_lower_pct),
      rd_upper_pct = if_else(trt_group == ref_name, 0, rd_upper_pct),
      Interruption = int_label,
      Stratum = stratum_label
    )
}

extract_rmst_data <- function(risk_obj, dat_weights, max_time, outcome = "dis", 
                              ref_name = "no reversion", stratum_label, int_label) {
  
  time_v  <- if(outcome == "dis") "time_discont" else "time_od_acute"
  event_v <- if(outcome == "dis") "sens_discont" else "sens_od_acute"
  
  counts <- compute_summary_stats_timewin(
    dat_weights, time_var = time_v, event_var = event_v, 
    weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
  ) %>% dplyr::select(trt_group = trt, Episode)
  
  rmst_abs <- risk_obj$rmst %>%
    mutate(rmst_disp = sprintf("%.1f (%.1f, %.1f)", RMST, RMST_lower, RMST_upper)) %>% 
    dplyr::select(trt_group = trt, rmst_disp)
  
  rmst_diffs <- risk_obj$rmst_contrasts %>%
    mutate(diff_disp = sprintf("%.1f (%.1f, %.1f)", RMST_diff, RMST_diff_lower, RMST_diff_upper)) %>% 
    dplyr::select(trt_group = trt, RMST_diff, RMST_diff_lower, RMST_diff_upper, diff_disp)
  
  counts %>%
    left_join(rmst_abs, by = "trt_group") %>%
    left_join(rmst_diffs, by = "trt_group") %>%
    mutate(
      diff_disp = if_else(trt_group == ref_name, "Ref.", diff_disp),
      RMST_diff = if_else(trt_group == ref_name, 0, RMST_diff),
      RMST_diff_lower = if_else(trt_group == ref_name, 0, RMST_diff_lower),
      RMST_diff_upper = if_else(trt_group == ref_name, 0, RMST_diff_upper),
      Interruption = int_label,
      Stratum = stratum_label
    )
}


# 5. FIGURE PREPARATION ###################

compile_plot_data <- function(metric_type = "rd", outcome = "od", risk_at = 30) {
  extr_func <- if(metric_type == "rd") extract_rd_data else extract_rmst_data
  
  d1 <- extr_func(risk_strat$tr2$le30[[outcome]]$ci, ps_strat$tr2$le30$data_weights, risk_at, outcome, 
                  stratum_label = "≤30 days", int_label = "3-4 days")
  d2 <- extr_func(risk_strat$tr2$gt30[[outcome]]$ci, ps_strat$tr2$gt30$data_weights, risk_at, outcome, 
                  stratum_label = ">30 days", int_label = "3-4 days")
  d3 <- extr_func(risk_strat$tr3$le30[[outcome]]$ci, ps_strat$tr3$le30$data_weights, risk_at, outcome, 
                  stratum_label = "≤30 days", int_label = "5-14 days")
  d4 <- extr_func(risk_strat$tr3$gt30[[outcome]]$ci, ps_strat$tr3$gt30$data_weights, risk_at, outcome, 
                  stratum_label = ">30 days", int_label = "5-14 days")
  
  bind_rows(d1, d2, d3, d4) %>%
    mutate(
      trt_group = recode(trt_group,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                         "increase" = "Dose increase"),
      Interruption = ifelse(duplicated(Interruption), "", Interruption),
      # Only print the stratum label on the first row of its block
      Stratum = ifelse(trt_group == "No dose change", Stratum, ""),
      ` ` = paste(rep(" ", 16), collapse = " ") 
    ) %>%
    rename(`Dose adjustment` = "trt_group", `Episode at risk\nat reinitiation, n` = "Episode")
}

# --- Figure 2A: RD Overdose (30 days) ---
data_2a <- compile_plot_data("rd", "od", 30) %>%
  rename(`% risk at 1 month\n(95% CI)` = "abs_risk_disp", `Risk Difference (%)\n(95% CI)` = "rd_disp")

p_2a <- forest(
  data_2a[, c("Interruption", "Stratum", "Dose adjustment", "Episode at risk\nat reinitiation, n", "% risk at 1 month\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est = data_2a$rd_pct, lower = data_2a$rd_lower_pct, upper = data_2a$rd_upper_pct,
  ci_column = 7, ref_line = 0, x_trans = "none", xlim = c(-2.5, 2.5), ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk difference (%)",
  title = "A. % difference in risk of overdose-related acute care visits\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7, refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"), xlab_adjust = c("refline", "center"))
)


p_2a_fin <- add_grob(p_2a, row = nrow(data_2a)+2, col = 1:2, order = "background", gb_fn = gridtext::richtext_grob, 
                     text = "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals", gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))
jpeg(paste0(RES_DIR, "Figures/Stratified/Figure2_A_OD_RD.jpg"), width = 3800, height = 1800, res = 300)
p_2a_fin
dev.off()


# --- Figure 2B: RD Discontinuation (365 days) ---
data_2b <- compile_plot_data("rd", "dis", 365) %>%
  rename(`% risk at 12 months\n(95% CI)` = "abs_risk_disp", `Risk Difference (%)\n(95% CI)` = "rd_disp")

p_2b <- forest(
  data_2b[, c("Interruption", "Stratum", "Dose adjustment", "Episode at risk\nat reinitiation, n", "% risk at 12 months\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est = data_2b$rd_pct, lower = data_2b$rd_lower_pct, upper = data_2b$rd_upper_pct,
  ci_column = 7, ref_line = 0, x_trans = "none", xlim = c(-12.5, 12.5), ticks_at = c(-10, -5, 0, 5, 10),
  xlab = "Risk difference (%)",
  title = "B. % difference in risk of treatment discontinuation (12 months)\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7, refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"), xlab_adjust = c("refline", "center"))
)
p_2b_fin <- add_grob(p_2b, row = nrow(data_2b)+2, col = 1:2, order = "background", gb_fn = gridtext::richtext_grob, 
                     text = "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals", gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))
jpeg(paste0(RES_DIR, "Figures/Stratified/Figure2_B_Discont_RD.jpg"), width = 3800, height = 1800, res = 300)
p_2b_fin
dev.off()


# --- Figure 3A: RMST Overdose (365 days) ---
data_3a <- compile_plot_data("rmst", "od", 365) %>%
  rename(`Event-free days\nat 12 months (95% CI)` = "rmst_disp", `Difference in days\n(95% CI)` = "diff_disp")

p_3a <- forest(
  data_3a[, c("Interruption", "Stratum", "Dose adjustment", "Episode at risk\nat reinitiation, n", "Event-free days\nat 12 months (95% CI)", "Difference in days\n(95% CI)", " ")],
  est = data_3a$RMST_diff, lower = data_3a$RMST_diff_lower, upper = data_3a$RMST_diff_upper,
  ci_column = 7, ref_line = 0, x_trans = "none", xlim = c(-5, 5), ticks_at = c(-4, -2, 0, 2, 4), # Narrow bounds for rare OD event
  xlab = "Difference in event-free days",
  title = "A. Restricted mean survival time for overdose-related acute care visits\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7, refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"), xlab_adjust = c("refline", "center"))
)
p_3a_fin <- add_grob(p_3a, row = nrow(data_3a)+2, col = 1:2, order = "background", gb_fn = gridtext::richtext_grob, 
                     text = "Positive difference indicates more event-free days.<br>CI: compatibility ('confidence') intervals", gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))
jpeg(paste0(RES_DIR, "Figures/Stratified/Figure3_A_OD_RMST.jpg"), width = 3800, height = 1800, res = 300)
p_3a_fin
dev.off()


# --- Figure 3B: RMST Discontinuation (365 days) ---
data_3b <- compile_plot_data("rmst", "dis", 365) %>%
  rename(`Event-free days\nat 12 months (95% CI)` = "rmst_disp", `Difference in days\n(95% CI)` = "diff_disp")

p_3b <- forest(
  data_3b[, c("Interruption", "Stratum", "Dose adjustment", "Episode at risk\nat reinitiation, n", "Event-free days\nat 12 months (95% CI)", "Difference in days\n(95% CI)", " ")],
  est = data_3b$RMST_diff, lower = data_3b$RMST_diff_lower, upper = data_3b$RMST_diff_upper,
  ci_column = 7, ref_line = 0, x_trans = "none", xlim = c(-45, 45), ticks_at = c(-40, -20, 0, 20, 40),
  xlab = "Difference in event-free days",
  title = "B. Restricted mean survival time for treatment discontinuation (12-month horizon)\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7, refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"), xlab_adjust = c("refline", "center"))
)
p_3b_fin <- add_grob(p_3b, row = nrow(data_3b)+2, col = 1:2, order = "background", gb_fn = gridtext::richtext_grob, 
                     text = "Positive difference indicates more event-free days.<br>CI: compatibility ('confidence') intervals", gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))
jpeg(paste0(RES_DIR, "Figures/Stratified/Figure3_B_Discont_RMST.jpg"), width = 3800, height = 1800, res = 300)
p_3b_fin
dev.off()