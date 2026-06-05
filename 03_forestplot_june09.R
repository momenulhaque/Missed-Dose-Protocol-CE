library(tidyverse)
library(forestploter)
library(grid)
# FOREST PLOT PREPARATION: Cumulative Incidence & Rates
# Missed-dose trials: MET interruption of 3-4 days and 5-14 days
# Outcomes: Overdose (1 month) and Discontinuation (12 months)

library(tidyverse)

# LOAD RISK RESULTS

RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

# UNIFIED EXTRACTION: Absolute Risk and Risk Difference at time t



risk_tr2 <- readRDS(paste0(RES_DIR, "Risk/risk_tr2_34d_rmst.rds"))

risk_tr3 <- readRDS(paste0(RES_DIR, "Risk/risk_tr3_514d_rmst.rds"))

dat_tr2_first <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))
dat_tr2_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
dat_tr3_first <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))
dat_tr3_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))



library(tidyverse)
library(forestploter)
library(grid)

# ── 1. UPDATED EXTRACTION FUNCTIONS (Cohort Dropped) ──────────────────────────

extract_rd_data <- function(risk_obj, dat_weights, max_time, outcome = "od", 
                            ref_name = "no reversion", int_label) {
  
  time_v  <- if(outcome == "dis") "time_discont" else "time_od_acute"
  event_v <- if(outcome == "dis") "sens_discont" else "sens_od_acute"
  
  counts <- compute_summary_stats_timewin(
    dat_weights, time_var = time_v, event_var = event_v, 
    weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
  ) %>% dplyr::select(trt_group = trt, Episode)
  
  risks <- risk_obj$risks %>%
    filter(horizon == max_time) %>%
    mutate(
      risk_pct = risk * 100,
      risk_lower_pct = risk_lower * 100,
      risk_upper_pct = risk_upper * 100,
      abs_risk_disp = ifelse(
        abs(risk_pct) >= 10,
        sprintf("%.1f (%.1f, %.1f)", risk_pct, risk_lower_pct, risk_upper_pct),
        sprintf("%.2f (%.2f, %.2f)", risk_pct, risk_lower_pct, risk_upper_pct)
      )
    ) %>% dplyr::select(trt_group = trt, abs_risk_disp)
  
  contrasts <- risk_obj$contrasts %>%
    filter(horizon == max_time) %>%
    mutate(
      rd_pct = RD * 100,
      rd_lower_pct = RD_lower * 100,
      rd_upper_pct = RD_upper * 100,
      rd_disp = ifelse(
        abs(rd_pct) >= 10,
        sprintf("%.1f (%.1f, %.1f)", rd_pct, rd_lower_pct, rd_upper_pct),
        sprintf("%.2f (%.2f, %.2f)", rd_pct, rd_lower_pct, rd_upper_pct)
      )
    ) %>% dplyr::select(trt_group = trt, rd_pct, rd_lower_pct, rd_upper_pct, rd_disp)
  
  res <- counts %>%
    left_join(risks, by = "trt_group") %>%
    left_join(contrasts, by = "trt_group") %>%
    mutate(
      rd_disp = if_else(trt_group == ref_name, "Ref.", rd_disp),
      rd_pct = if_else(trt_group == ref_name, 0, rd_pct),
      rd_lower_pct = if_else(trt_group == ref_name, 0, rd_lower_pct),
      rd_upper_pct = if_else(trt_group == ref_name, 0, rd_upper_pct),
      Interruption = int_label
    )
  
  return(res)
}


extract_rmst_data <- function(risk_obj, dat_weights, max_time, outcome = "dis", 
                              ref_name = "no reversion", int_label) {
  
  time_v  <- if(outcome == "dis") "time_discont" else "time_od_acute"
  event_v <- if(outcome == "dis") "sens_discont" else "sens_od_acute"
  
  counts <- compute_summary_stats_timewin(
    dat_weights, time_var = time_v, event_var = event_v, 
    weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
  ) %>% dplyr::select(trt_group = trt, Episode)
  
  rmst_abs <- risk_obj$rmst %>%
    mutate(
      rmst_disp = sprintf("%.1f (%.1f, %.1f)", RMST, RMST_lower, RMST_upper)
    ) %>% dplyr::select(trt_group = trt, rmst_disp)
  
  rmst_diffs <- risk_obj$rmst_contrasts %>%
    mutate(
      diff_disp = sprintf("%.1f (%.1f, %.1f)", RMST_diff, RMST_diff_lower, RMST_diff_upper)
    ) %>% dplyr::select(trt_group = trt, RMST_diff, RMST_diff_lower, RMST_diff_upper, diff_disp)
  
  res <- counts %>%
    left_join(rmst_abs, by = "trt_group") %>%
    left_join(rmst_diffs, by = "trt_group") %>%
    mutate(
      diff_disp = if_else(trt_group == ref_name, "Ref.", diff_disp),
      RMST_diff = if_else(trt_group == ref_name, 0, RMST_diff),
      RMST_diff_lower = if_else(trt_group == ref_name, 0, RMST_diff_lower),
      RMST_diff_upper = if_else(trt_group == ref_name, 0, RMST_diff_upper),
      Interruption = int_label
    )
  return(res)
}


# ── 2. FIGURE 2: RISK DIFFERENCES (Any Episode Cohort) ────────────────────────

# --- FIGURE 2A: Overdose Risk Difference (90 Days) ---
risk_at_od <- 90

tr2_od_rd <- extract_rd_data(risk_tr2$all$od$ci, dat_tr2_all$data_weights, risk_at_od, 
                             outcome = "od", int_label = "3-4 days")
tr3_od_rd <- extract_rd_data(risk_tr3$all$od$ci, dat_tr3_all$data_weights, risk_at_od, 
                             outcome = "od", int_label = "5-14 days")

plot_data_2a <- bind_rows(tr2_od_rd, tr3_od_rd) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct),
    Interruption = ifelse(duplicated(Interruption), "", Interruption),
    ` ` = paste(rep(" ", 16), collapse = " ") 
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `% risk at 3 months\n(95% CI)` = "abs_risk_disp",
    `Risk Difference (%)\n(95% CI)` = "rd_disp"
  )

p_2a <- forest(
  plot_data_2a[, c("Interruption", "Dose adjustment\nafter interruption", "Episode, n", 
                   "% risk at 3 months\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est = plot_data_2a$rd_pct, lower = plot_data_2a$rd_lower_pct, upper = plot_data_2a$rd_upper_pct,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-2.5, 2.5), ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk difference (%)",
  title = "A. Absolute difference in risk of overdose-related acute care visits (3 months)\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

p_2a_final <- add_grob(p_2a, row = nrow(plot_data_2a) + 2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = "Risk: cumulative incidence<br>CI: compatibility intervals", 
                       gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, 
                       x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/June09/", "Figure2_A_OD_RD.jpg"), width = 3500, height = 1400, res = 300)
p_2a_final
dev.off()


# --- FIGURE 2B: Discontinuation Risk Difference (365 Days) ---
risk_at_dis <- 365

tr2_dis_rd <- extract_rd_data(risk_tr2$all$dis$ci, dat_tr2_all$data_weights, risk_at_dis, 
                              outcome = "dis", int_label = "3-4 days")
tr3_dis_rd <- extract_rd_data(risk_tr3$all$dis$ci, dat_tr3_all$data_weights, risk_at_dis, 
                              outcome = "dis", int_label = "5-14 days")

plot_data_2b <- bind_rows(tr2_dis_rd, tr3_dis_rd) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct),
    Interruption = ifelse(duplicated(Interruption), "", Interruption),
    ` ` = paste(rep(" ", 16), collapse = " ") 
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `% risk at 12 months\n(95% CI)` = "abs_risk_disp",
    `Risk Difference (%)\n(95% CI)` = "rd_disp"
  )

p_2b <- forest(
  plot_data_2b[, c("Interruption", "Dose adjustment\nafter interruption", "Episode, n", 
                   "% risk at 12 months\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est = plot_data_2b$rd_pct, lower = plot_data_2b$rd_lower_pct, upper = plot_data_2b$rd_upper_pct,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-12.5, 12.5), ticks_at = c(-10, -5, 0, 5, 10),
  xlab = "Risk difference (%)",
  title = "B. Absolute difference in risk of treatment discontinuation (12 months)\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

p_2b_final <- add_grob(p_2b, row = nrow(plot_data_2b) + 2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = "Risk: cumulative incidence<br>CI: compatibility intervals", 
                       gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, 
                       x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/June09/", "Figure2_B_Discont_RD.jpg"), width = 3500, height = 1400, res = 300)
p_2b_final
dev.off()


# ── 3. FIGURE 3: RMST DIFFERENCES (Any Episode Cohort) ────────────────────────

# --- FIGURE 3A: Overdose RMST ---
# Assumes predict_cuminc_rmst was run with rmst_tau = 365 for OD
risk_at_rmst_od <- 365 

tr2_od_rmst <- extract_rmst_data(risk_tr2$all$od$ci, dat_tr2_all$data_weights, risk_at_rmst_od, 
                                 outcome = "od", int_label = "3-4 days")
tr3_od_rmst <- extract_rmst_data(risk_tr3$all$od$ci, dat_tr3_all$data_weights, risk_at_rmst_od, 
                                 outcome = "od", int_label = "5-14 days")

plot_data_3a <- bind_rows(tr2_od_rmst, tr3_od_rmst) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    RMST_diff = ifelse(trt_group == "No dose change", NA, RMST_diff),
    RMST_diff_lower = ifelse(trt_group == "No dose change", NA, RMST_diff_lower),
    RMST_diff_upper = ifelse(trt_group == "No dose change", NA, RMST_diff_upper),
    Interruption = ifelse(duplicated(Interruption), "", Interruption),
    ` ` = paste(rep(" ", 16), collapse = " ")
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `Event-free days\nat 12 months (95% CI)` = "rmst_disp",
    `Difference in days\n(95% CI)` = "diff_disp"
  )

p_3a <- forest(
  plot_data_3a[, c("Interruption", "Dose adjustment\nafter interruption", "Episode, n", 
                   "Event-free days\nat 12 months (95% CI)", "Difference in days\n(95% CI)", " ")],
  est = plot_data_3a$RMST_diff, lower = plot_data_3a$RMST_diff_lower, upper = plot_data_3a$RMST_diff_upper,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-5, 5), ticks_at = c(-4, -2, 0, 2, 4), # Narrower bounds expected for rare OD
  xlab = "Difference in event-free days",
  title = "A. Restricted mean survival time for overdose-related acute care visits (12-month horizon)\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

txt_rmst <- "Positive difference indicates more event-free days.<br>CI: compatibility intervals"
p_3a_final <- add_grob(p_3a, row = nrow(plot_data_3a)+2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = txt_rmst, gp = gpar(fontsize = 11), 
                       hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/June09/", "Figure3_A_OD_RMST.jpg"), width = 3500, height = 1400, res = 300)
p_3a_final
dev.off()


# --- FIGURE 3B: Discontinuation RMST (365 Days) ---
risk_at_rmst_dis <- 365

tr2_dis_rmst <- extract_rmst_data(risk_tr2$all$dis$ci, dat_tr2_all$data_weights, risk_at_rmst_dis, 
                                  outcome = "dis", int_label = "3-4 days")
tr3_dis_rmst <- extract_rmst_data(risk_tr3$all$dis$ci, dat_tr3_all$data_weights, risk_at_rmst_dis, 
                                  outcome = "dis", int_label = "5-14 days")

plot_data_3b <- bind_rows(tr2_dis_rmst, tr3_dis_rmst) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    RMST_diff = ifelse(trt_group == "No dose change", NA, RMST_diff),
    RMST_diff_lower = ifelse(trt_group == "No dose change", NA, RMST_diff_lower),
    RMST_diff_upper = ifelse(trt_group == "No dose change", NA, RMST_diff_upper),
    Interruption = ifelse(duplicated(Interruption), "", Interruption),
    ` ` = paste(rep(" ", 16), collapse = " ")
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `Event-free days\nat 12 months (95% CI)` = "rmst_disp",
    `Difference in days\n(95% CI)` = "diff_disp"
  )

p_3b <- forest(
  plot_data_3b[, c("Interruption", "Dose adjustment\nafter interruption", "Episode, n", 
                   "Event-free days\nat 12 months (95% CI)", "Difference in days\n(95% CI)", " ")],
  est = plot_data_3b$RMST_diff, lower = plot_data_3b$RMST_diff_lower, upper = plot_data_3b$RMST_diff_upper,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-45, 45), ticks_at = c(-40, -20, 0, 20, 40),
  xlab = "Difference in event-free days",
  title = "B. Restricted mean survival time for treatment discontinuation (12-month horizon)\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

p_3b_final <- add_grob(p_3b, row = nrow(plot_data_3b)+2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = txt_rmst, gp = gpar(fontsize = 11), 
                       hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/June09/", "Figure3_B_Discont_RMST.jpg"), width = 3500, height = 1400, res = 300)
p_3b_final
dev.off()