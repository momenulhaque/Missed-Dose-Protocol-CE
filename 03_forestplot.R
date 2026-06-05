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



compute_summary_stats_timewin <- function(Dat, 
                                          time_var = "time_od_acute", 
                                          event_var = "sens_od_acute",
                                          weight_var = "wt_trunc",
                                          trt_var = "trt",
                                          max_time = 90) {
  
  library(dplyr)
  library(rlang)
  
  Dat %>%
    # Restrict to events within time window
    mutate(
      event_in_window = ifelse(!!sym(time_var) <= max_time, !!sym(event_var), 0)
    ) %>%
    group_by(!!sym(trt_var)) %>%
    summarise(
      Episode = n(),
      Weighted_Episode = sum(!!sym(weight_var), na.rm = TRUE),
      Event = sum(event_in_window, na.rm = TRUE),  # Events within time window
      Weighted_event = sum(!!sym(weight_var) * event_in_window, na.rm = TRUE),
      Total_PY = sum(!!sym(weight_var) * !!sym(time_var), na.rm = TRUE) / 52.17857,
      Rate_PY = Weighted_event / Total_PY,
      Rate_100PY = Rate_PY * 100,
      .groups = 'drop'
    ) %>%
    select(trt_var, Episode, Event, Rate_100PY)
}


extract_rr_data <- function(risk_obj, 
                            dat_weights, 
                            max_time, 
                            outcome = "od",
                            ref_name = "no reversion", 
                            cohort_label, 
                            int_label) {
  # 1. Counts
  if(outcome == "dis"){
    counts <- compute_summary_stats_timewin(
      dat_weights, time_var = "time_discont", event_var = "sens_discont", 
      weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
    ) %>% dplyr::select(trt_group = trt, Episode)
  }
  
  if(outcome == "od"){
    counts <- compute_summary_stats_timewin(
      dat_weights, time_var = "time_od_acute", event_var = "sens_od_acute", 
      weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
    ) %>% dplyr::select(trt_group = trt, Episode)
  }  
  
  # 2. Absolute Risk
  risks <- risk_obj$risks %>%
    filter(horizon == max_time) %>%
    mutate(
      risk_pct = risk * 100,
      abs_risk_disp = sprintf("%.2f (%.2f, %.2f)", risk * 100, risk_lower * 100, risk_upper * 100)
    ) %>% select(trt_group = trt, abs_risk_disp)
  
  # 3. Contrasts (RR)
  contrasts <- risk_obj$contrasts %>%
    filter(horizon == max_time) %>%
    mutate(
      rr_disp = sprintf("%.2f (%.2f, %.2f)", RR, RR_lower, RR_upper)
    ) %>% select(trt_group = trt, RR, RR_lower, RR_upper, rr_disp)
  
  # Merge
  res <- counts %>%
    left_join(risks, by = "trt_group") %>%
    left_join(contrasts, by = "trt_group") %>%
    mutate(
      rr_disp = if_else(trt_group == ref_name, "Ref.", rr_disp),
      RR = if_else(trt_group == ref_name, 1, RR),
      RR_lower = if_else(trt_group == ref_name, 1, RR_lower),
      RR_upper = if_else(trt_group == ref_name, 1, RR_upper),
      Cohort = cohort_label,
      Interruption = int_label
    )
  return(res)
}




extract_rd_data <- function(risk_obj, 
                            dat_weights, 
                            max_time, 
                            outcome = "od",
                            ref_name = "no reversion", 
                            cohort_label, 
                            int_label) {
  
  # 1. Counts
  if(outcome == "dis"){
    counts <- compute_summary_stats_timewin(
      dat_weights, time_var = "time_discont", event_var = "sens_discont", 
      weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
    ) %>% dplyr::select(trt_group = trt, Episode)
  }
  
  if(outcome == "od"){
    counts <- compute_summary_stats_timewin(
      dat_weights, time_var = "time_od_acute", event_var = "sens_od_acute", 
      weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
    ) %>% dplyr::select(trt_group = trt, Episode)
  }
  
  
  # 2. Absolute Risk (converted to %)
  risks <- risk_obj$risks %>%
    filter(horizon == max_time) %>%
    mutate(
      risk_pct = risk * 100,
      risk_lower_pct = risk_lower * 100,
      risk_upper_pct = risk_upper * 100,
      # Conditional formatting: 1 decimal if >= 10%, 2 decimals if < 10%
      abs_risk_disp = ifelse(
        abs(risk_pct) >= 10,
        sprintf("%.1f (%.1f, %.1f)", risk_pct, risk_lower_pct, risk_upper_pct),
        sprintf("%.2f (%.2f, %.2f)", risk_pct, risk_lower_pct, risk_upper_pct)
      )
    ) %>% dplyr::select(trt_group = trt, abs_risk_disp)
  
  # 3. Contrasts (RD converted to %)
  contrasts <- risk_obj$contrasts %>%
    filter(horizon == max_time) %>%
    mutate(
      rd_pct = RD * 100,
      rd_lower_pct = RD_lower * 100,
      rd_upper_pct = RD_upper * 100,
      # Conditional formatting for Risk Difference
      rd_disp = ifelse(
        abs(rd_pct) >= 10,
        sprintf("%.1f (%.1f, %.1f)", rd_pct, rd_lower_pct, rd_upper_pct),
        sprintf("%.2f (%.2f, %.2f)", rd_pct, rd_lower_pct, rd_upper_pct)
      )
    ) %>% dplyr::select(trt_group = trt, rd_pct, rd_lower_pct, rd_upper_pct, rd_disp)
  
  # 4. Merge and finalize reference group
  res <- counts %>%
    left_join(risks, by = "trt_group") %>%
    left_join(contrasts, by = "trt_group") %>%
    mutate(
      rd_disp = if_else(trt_group == ref_name, "Ref.", rd_disp),
      rd_pct = if_else(trt_group == ref_name, 0, rd_pct),
      rd_lower_pct = if_else(trt_group == ref_name, 0, rd_lower_pct),
      rd_upper_pct = if_else(trt_group == ref_name, 0, rd_upper_pct),
      Cohort = cohort_label,
      Interruption = int_label
    )
  
  return(res)
}



# --- HELPER: Extract RMST Differences ---
extract_rmst_data <- function(risk_obj, 
                              dat_weights, 
                              max_time, 
                              outcome = "dis",
                              ref_name = "no reversion", 
                              cohort_label, 
                              int_label) {
  # 1. Counts
  
  # 1. Counts
  if(outcome == "dis"){
    counts <- compute_summary_stats_timewin(
      dat_weights, time_var = "time_discont", event_var = "sens_discont", 
      weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
    ) %>% dplyr::select(trt_group = trt, Episode)
  }
  
  if(outcome == "od"){
    counts <- compute_summary_stats_timewin(
      dat_weights, time_var = "time_od_acute", event_var = "sens_od_acute", 
      weight_var = "wt_trunc", trt_var = "trt", max_time = max_time
    ) %>% dplyr::select(trt_group = trt, Episode)
  }
  
  
  # 2. Absolute RMST
  rmst_abs <- risk_obj$rmst %>%
    mutate(
      rmst_disp = sprintf("%.1f (%.1f, %.1f)", RMST, RMST_lower, RMST_upper)
    ) %>% select(trt_group = trt, rmst_disp)
  
  # 3. RMST Contrasts
  rmst_diffs <- risk_obj$rmst_contrasts %>%
    mutate(
      diff_disp = sprintf("%.1f (%.1f, %.1f)", RMST_diff, RMST_diff_lower, RMST_diff_upper)
    ) %>% select(trt_group = trt, RMST_diff, RMST_diff_lower, RMST_diff_upper, diff_disp)
  
  # Merge
  res <- counts %>%
    left_join(rmst_abs, by = "trt_group") %>%
    left_join(rmst_diffs, by = "trt_group") %>%
    mutate(
      diff_disp = if_else(trt_group == ref_name, "Ref.", diff_disp),
      RMST_diff = if_else(trt_group == ref_name, 0, RMST_diff),
      RMST_diff_lower = if_else(trt_group == ref_name, 0, RMST_diff_lower),
      RMST_diff_upper = if_else(trt_group == ref_name, 0, RMST_diff_upper),
      Cohort = cohort_label,
      Interruption = int_label
    )
  return(res)
}



risk_tr2 <- readRDS(paste0(RES_DIR, "Risk/risk_tr2_34d_rmst.rds"))

risk_tr3 <- readRDS(paste0(RES_DIR, "Risk/risk_tr3_514d_rmst.rds"))

dat_tr2_first <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))
dat_tr2_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
dat_tr3_first <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))
dat_tr3_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))


# od #########
## Figure 2A: met 3-4: od: 90-day-RR ###########

risk_at <- 90

# Build Data for Trial 2
tr2_first_rr <- extract_rr_data(risk_tr2$first$od$ci, outcome="od", 
                                dat_tr2_first$data_weights, risk_at, 
                                cohort_label = "First episode", int_label = "3-4 days")
tr2_all_rr   <- extract_rr_data(risk_tr2$all$od$ci,  outcome="od", 
                                dat_tr2_all$data_weights, risk_at, 
                                cohort_label = "Any episode", int_label = "3-4 days")

plot_data_rr <- bind_rows(tr2_first_rr, tr2_all_rr) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    RR = ifelse(trt_group == "No dose change", NA, RR),
    RR_lower = ifelse(trt_group == "No dose change", NA, RR_lower),
    RR_upper = ifelse(trt_group == "No dose change", NA, RR_upper),
    Cohort = ifelse(duplicated(Cohort), "", Cohort),
    ` ` = paste(rep(" ", 16), collapse = " ") # Empty canvas column
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `% risk at 3 months\n(95% CI)` = "abs_risk_disp",
    `Risk Ratio\n(95% CI)` = "rr_disp"
  )

# Generate Forest Plot
p_rr <- forest(
  plot_data_rr[, c("Cohort", "Dose adjustment\nafter interruption", "Episode, n", 
                   "% risk at 3 months\n(95% CI)", "Risk Ratio\n(95% CI)", " ")],
  est = plot_data_rr$RR, lower = plot_data_rr$RR_lower, upper = plot_data_rr$RR_upper,
  ci_column = 6, 
  ref_line = 1,          # Multiplicative null
  x_trans = "log",       # Log scale for ratios
  xlim = c(0.25, 4),     # Adjust bounds if CI exceeds 4
  ticks_at = c(0.25, 0.5, 1, 2, 4),
  xlab = "Risk ratio",
  title = "A. Relative risk of overdose-related acute care visits among dose adjustment groups\n    after first interruption of 3-4 days\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

# Add Legend
txt_rr <- "Risk: cumulative incidence<br>RR: Ratio of cumulative incidences<br>CI: compatibility intervals"
p_rr_final <- add_grob(p_rr, row = nrow(plot_data_rr)+2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = txt_rr, gp = gpar(fontsize = 11), 
                       hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/", "Figure2_A_3m_RR.jpg"), width = 3500, height = 1400, res = 300)
p_rr_final
dev.off()

## Figure 2A: met 3-4: od: 90-day-RD ###########

risk_at <- 90

# Build Data for Trial 2 (3-4 days interruption)
tr2_first_rd <- extract_rd_data(risk_tr2$first$od$ci, outcome = "od",
                                dat_tr2_first$data_weights, risk_at, 
                                cohort_label = "First episode", int_label = "3-4 days")
tr2_all_rd   <- extract_rd_data(risk_tr2$all$od$ci, outcome = "od",
                                dat_tr2_all$data_weights, risk_at, 
                                cohort_label = "Any episode", int_label = "3-4 days")

# Combine and Format for Plotting
plot_data_rd <- bind_rows(tr2_first_rd, tr2_all_rd) %>%
  mutate(
    # Clean treatment names
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    
    # NA out the numeric values for the reference group so the plot doesn't draw a box at 0
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct),
    
    # Clean repetitive cohort labels
    Cohort = ifelse(duplicated(Cohort), "", Cohort),
    
    # Create the empty string column to act as the canvas for the forest plot visual
    ` ` = paste(rep(" ", 16), collapse = " ") 
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `% risk at 3 months\n(95% CI)` = "abs_risk_disp",
    `Risk Difference (%)\n(95% CI)` = "rd_disp"
  )

# Generate Forest Plot
library(forestploter)
library(grid)

p_rd <- forest(
  plot_data_rd[, c("Cohort", "Dose adjustment\nafter interruption", "Episode, n", 
                   "% risk at 3 months\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est       = plot_data_rd$rd_pct, 
  lower     = plot_data_rd$rd_lower_pct, 
  upper     = plot_data_rd$rd_upper_pct,
  ci_column = 6, 
  ref_line  = 0,          # Additive null for Risk Difference
  x_trans   = "none",     # Linear scale
  xlim      = c(-2.5, 2.5),  # Adjust bounds based on your max expected difference (e.g., +/- 2.5%)
  ticks_at  = c(-2, -1, 0, 1, 2),
  xlab      = "Risk difference (%)",
  title     = "A. % difference in risk of overdose-related acute care visits among dose adjustment groups\n    after first interruption of 3-4 days\n",
  theme     = forest_theme(
    base_size   = 14, 
    ci_pch      = 15, 
    ci_lwd      = 1.7,
    refline_gp  = gpar(lwd = 1, lty = "dashed", col = "red"),
    xlab_adjust = c("refline", "center")
  )
)

# Add Legend/Footnote
txt_rd <- "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals"
p_rd_final <- add_grob(
  p_rd, 
  row = nrow(plot_data_rd) + 2, 
  col = 1:2, 
  order = "background", 
  gb_fn = gridtext::richtext_grob, 
  text = txt_rd, 
  gp = gpar(fontsize = 11), 
  hjust = 0, vjust = 1, halign = 0, valign = 1, 
  x = unit(0, "npc"), y = unit(1, "npc")
)

# Save output
jpeg(paste0(RES_DIR, "Figures/", "Figure2_A_3m_RD.jpg"), width = 3500, height = 1400, res = 300)
p_rd_final
dev.off()



## Figure 2A: met 3-4: od: 30-day-RR ###########

risk_at <- 30

# Build Data for Trial 2
tr2_first_rr <- extract_rr_data(risk_tr2$first$od$ci, outcome="od", 
                                dat_tr2_first$data_weights, risk_at, 
                                cohort_label = "First episode", int_label = "3-4 days")
tr2_all_rr   <- extract_rr_data(risk_tr2$all$od$ci,  outcome="od", 
                                dat_tr2_all$data_weights, risk_at, 
                                cohort_label = "Any episode", int_label = "3-4 days")

plot_data_rr <- bind_rows(tr2_first_rr, tr2_all_rr) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    RR = ifelse(trt_group == "No dose change", NA, RR),
    RR_lower = ifelse(trt_group == "No dose change", NA, RR_lower),
    RR_upper = ifelse(trt_group == "No dose change", NA, RR_upper),
    Cohort = ifelse(duplicated(Cohort), "", Cohort),
    ` ` = paste(rep(" ", 16), collapse = " ") # Empty canvas column
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `% risk at 1 month\n(95% CI)` = "abs_risk_disp",
    `Risk Ratio\n(95% CI)` = "rr_disp"
  )

# Generate Forest Plot
p_rr <- forest(
  plot_data_rr[, c("Cohort", "Dose adjustment\nafter interruption", "Episode, n", 
                   "% risk at 1 month\n(95% CI)", "Risk Ratio\n(95% CI)", " ")],
  est = plot_data_rr$RR, lower = plot_data_rr$RR_lower, upper = plot_data_rr$RR_upper,
  ci_column = 6, 
  ref_line = 1,          # Multiplicative null
  x_trans = "log",       # Log scale for ratios
  xlim = c(0.25, 4),     # Adjust bounds if CI exceeds 4
  ticks_at = c(0.25, 0.5, 1, 2, 4),
  xlab = "Risk ratio",
  title = "A. Relative risk of overdose-related acute care visits among dose adjustment groups\n    after first interruption of 3-4 days\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

# Add Legend
txt_rr <- "Risk: cumulative incidence<br>RR: Ratio of cumulative incidences<br>CI: compatibility intervals"
p_rr_final <- add_grob(p_rr, row = nrow(plot_data_rr)+2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = txt_rr, gp = gpar(fontsize = 11), 
                       hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/", "Figure2_A_1m_RR.jpg"), width = 3500, height = 1400, res = 300)
p_rr_final
dev.off()

## Figure 2A: met 3-4: od: 30-day-RD ###########

risk_at <- 30

# Build Data for Trial 2 (3-4 days interruption)
tr2_first_rd <- extract_rd_data(risk_tr2$first$od$ci, outcome = "od",
                                dat_tr2_first$data_weights, risk_at, 
                                cohort_label = "First episode", int_label = "3-4 days")
tr2_all_rd   <- extract_rd_data(risk_tr2$all$od$ci, outcome = "od",
                                dat_tr2_all$data_weights, risk_at, 
                                cohort_label = "Any episode", int_label = "3-4 days")

# Combine and Format for Plotting
plot_data_rd <- bind_rows(tr2_first_rd, tr2_all_rd) %>%
  mutate(
    # Clean treatment names
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    
    # NA out the numeric values for the reference group so the plot doesn't draw a box at 0
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct),
    
    # Clean repetitive cohort labels
    Cohort = ifelse(duplicated(Cohort), "", Cohort),
    
    # Create the empty string column to act as the canvas for the forest plot visual
    ` ` = paste(rep(" ", 16), collapse = " ") 
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `% risk at 1 month\n(95% CI)` = "abs_risk_disp",
    `Risk Difference (%)\n(95% CI)` = "rd_disp"
  )

# Generate Forest Plot
library(forestploter)
library(grid)

p_rd <- forest(
  plot_data_rd[, c("Cohort", "Dose adjustment\nafter interruption", "Episode, n", 
                   "% risk at 1 month\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est       = plot_data_rd$rd_pct, 
  lower     = plot_data_rd$rd_lower_pct, 
  upper     = plot_data_rd$rd_upper_pct,
  ci_column = 6, 
  ref_line  = 0,          # Additive null for Risk Difference
  x_trans   = "none",     # Linear scale
  xlim      = c(-2.5, 2.5),  # Adjust bounds based on your max expected difference (e.g., +/- 2.5%)
  ticks_at  = c(-2, -1, 0, 1, 2),
  xlab      = "Risk difference (%)",
  title     = "A. % difference in risk of overdose-related acute care visits among dose adjustment groups\n    after first interruption of 3-4 days\n",
  theme     = forest_theme(
    base_size   = 14, 
    ci_pch      = 15, 
    ci_lwd      = 1.7,
    refline_gp  = gpar(lwd = 1, lty = "dashed", col = "red"),
    xlab_adjust = c("refline", "center")
  )
)

# Add Legend/Footnote
txt_rd <- "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals"
p_rd_final <- add_grob(
  p_rd, 
  row = nrow(plot_data_rd) + 2, 
  col = 1:2, 
  order = "background", 
  gb_fn = gridtext::richtext_grob, 
  text = txt_rd, 
  gp = gpar(fontsize = 11), 
  hjust = 0, vjust = 1, halign = 0, valign = 1, 
  x = unit(0, "npc"), y = unit(1, "npc")
)

# Save output
jpeg(paste0(RES_DIR, "Figures/", "Figure2_A_1m_RD.jpg"), width = 3500, height = 1400, res = 300)
p_rd_final
dev.off()



# dis #########
## Figure 3A: met 3-4: dis: 12-months-RD ###########
risk_at <- 365

# Build Data for Trial 2 (3-4 days interruption)
tr2_first_rd <- extract_rd_data(risk_tr2$first$dis$ci, outcome = "dis", 
                                dat_tr2_first$data_weights, risk_at, 
                                cohort_label = "First episode", int_label = "3-4 days")
tr2_all_rd   <- extract_rd_data(risk_tr2$all$dis$ci, outcome = "dis", 
                                dat_tr2_all$data_weights, risk_at, 
                                cohort_label = "Any episode", int_label = "3-4 days")

# Combine and Format for Plotting
plot_data_rd <- bind_rows(tr2_first_rd, tr2_all_rd) %>%
  mutate(
    # Clean treatment names
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    
    # NA out the numeric values for the reference group so the plot doesn't draw a box at 0
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct),
    
    # Clean repetitive cohort labels
    Cohort = ifelse(duplicated(Cohort), "", Cohort),
    
    # Create the empty string column to act as the canvas for the forest plot visual
    ` ` = paste(rep(" ", 16), collapse = " ") 
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `% risk at 12 months\n(95% CI)` = "abs_risk_disp",
    `Risk Difference (%)\n(95% CI)` = "rd_disp"
  )

# Generate Forest Plot
library(forestploter)
library(grid)

p_rd <- forest(
  plot_data_rd[, c("Cohort", "Dose adjustment\nafter interruption", "Episode, n", 
                   "% risk at 12 months\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est       = plot_data_rd$rd_pct, 
  lower     = plot_data_rd$rd_lower_pct, 
  upper     = plot_data_rd$rd_upper_pct,
  ci_column = 6, 
  ref_line  = 0,          # Additive null for Risk Difference
  x_trans   = "none",     # Linear scale
  xlim = c(-12.5, 12.5),
  ticks_at = c(-10, -5, 0, 5, 10),
  xlab      = "Risk difference (%)",
  title     = "A. % difference in risk of treatment discontinuation among dose adjustment groups\n    after first interruption of 3-4 days\n",
  theme     = forest_theme(
    base_size   = 14, 
    ci_pch      = 15, 
    ci_lwd      = 1.7,
    refline_gp  = gpar(lwd = 1, lty = "dashed", col = "red"),
    xlab_adjust = c("refline", "center")
  )
)

# Add Legend/Footnote
txt_rd <- "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals"
p_rd_final <- add_grob(
  p_rd, 
  row = nrow(plot_data_rd) + 2, 
  col = 1:2, 
  order = "background", 
  gb_fn = gridtext::richtext_grob, 
  text = txt_rd, 
  gp = gpar(fontsize = 11), 
  hjust = 0, vjust = 1, halign = 0, valign = 1, 
  x = unit(0, "npc"), y = unit(1, "npc")
)

# Save output
jpeg(paste0(RES_DIR, "Figures/", "Figure3_A_12m_RD.jpg"), width = 3500, height = 1400, res = 300)
p_rd_final
dev.off()



## Figure 3B: met 3-4: dis: 12-months-RD ###########

risk_at <- 365

# Build Data for Trial 2 (5-14 days interruption)
tr3_first_rd <- extract_rd_data(risk_tr3$first$dis$ci, outcome = "dis", 
                                dat_tr3_first$data_weights, risk_at, 
                                cohort_label = "First episode", int_label = "5-14 days")
tr3_all_rd   <- extract_rd_data(risk_tr3$all$dis$ci, outcome = "dis", 
                                dat_tr3_all$data_weights, risk_at, 
                                cohort_label = "Any episode", int_label = "5-14 days")

# Combine and Format for Plotting
plot_data_rd <- bind_rows(tr3_first_rd, tr3_all_rd) %>%
  mutate(
    # Clean treatment names
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    
    # NA out the numeric values for the reference group so the plot doesn't draw a box at 0
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct),
    
    # Clean repetitive cohort labels
    Cohort = ifelse(duplicated(Cohort), "", Cohort),
    
    # Create the empty string column to act as the canvas for the forest plot visual
    ` ` = paste(rep(" ", 16), collapse = " ") 
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `% risk at 12 months\n(95% CI)` = "abs_risk_disp",
    `Risk Difference (%)\n(95% CI)` = "rd_disp"
  )

# Generate Forest Plot
library(forestploter)
library(grid)

p_rd <- forest(
  plot_data_rd[, c("Cohort", "Dose adjustment\nafter interruption", "Episode, n", 
                   "% risk at 12 months\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est       = plot_data_rd$rd_pct, 
  lower     = plot_data_rd$rd_lower_pct, 
  upper     = plot_data_rd$rd_upper_pct,
  ci_column = 6, 
  ref_line  = 0,          # Additive null for Risk Difference
  x_trans   = "none",     # Linear scale
  xlim = c(-12.5, 12.5),
  ticks_at = c(-10, -5, 0, 5, 10),
  xlab      = "Risk difference (%)",
  title     = "B. % difference in risk of treatment discontinuation among dose adjustment groups\n    after first interruption of 5-14 days\n",
  theme     = forest_theme(
    base_size   = 14, 
    ci_pch      = 15, 
    ci_lwd      = 1.7,
    refline_gp  = gpar(lwd = 1, lty = "dashed", col = "red"),
    xlab_adjust = c("refline", "center")
  )
)

# Add Legend/Footnote
txt_rd <- "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals"
p_rd_final <- add_grob(
  p_rd, 
  row = nrow(plot_data_rd) + 2, 
  col = 1:2, 
  order = "background", 
  gb_fn = gridtext::richtext_grob, 
  text = txt_rd, 
  gp = gpar(fontsize = 11), 
  hjust = 0, vjust = 1, halign = 0, valign = 1, 
  x = unit(0, "npc"), y = unit(1, "npc")
)


# Save output
jpeg(paste0(RES_DIR, "Figures/", "Figure3_B_12m_RD.jpg"), width = 3500, height = 1400, res = 300)
p_rd_final
dev.off()

# rmst #################
## eFigure A: met 3-4: dis: rmst at 365 day ###########
risk_at <- 365

# Build Data for Trial 2
tr2_first_rmst <- extract_rmst_data(risk_tr2$first$dis$ci, outcome = "dis",
                                    dat_tr2_first$data_weights, risk_at, 
                                    cohort_label = "First episode", int_label = "3-4 days")
tr2_all_rmst   <- extract_rmst_data(risk_tr2$all$dis$ci, outcome = "dis", 
                                    dat_tr2_all$data_weights, risk_at, 
                                    cohort_label = "Any episode", int_label = "3-4 days")

plot_data_rmst <- bind_rows(tr2_first_rmst, tr2_all_rmst) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    RMST_diff = ifelse(trt_group == "No dose change", NA, RMST_diff),
    RMST_diff_lower = ifelse(trt_group == "No dose change", NA, RMST_diff_lower),
    RMST_diff_upper = ifelse(trt_group == "No dose change", NA, RMST_diff_upper),
    Cohort = ifelse(duplicated(Cohort), "", Cohort),
    ` ` = paste(rep(" ", 16), collapse = " ") # Empty canvas column
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `Event-free days\nat 12 months (95% CI)` = "rmst_disp",
    `Difference in days\n(95% CI)` = "diff_disp"
  )

# Generate Forest Plot
p_rmst <- forest(
  plot_data_rmst[, c("Cohort", "Dose adjustment\nafter interruption", "Episode, n", 
                     "Event-free days\nat 12 months (95% CI)", "Difference in days\n(95% CI)", " ")],
  est = plot_data_rmst$RMST_diff, lower = plot_data_rmst$RMST_diff_lower, upper = plot_data_rmst$RMST_diff_upper,
  ci_column = 6, 
  ref_line = 0,          # Additive null for difference
  x_trans = "none",
  xlim = c(-45, 45),     # IMPORTANT: Adjust this wider if any RMST bounds exceed +/- 45 days
  ticks_at = c(-40, -20, 0, 20, 40),
  xlab = "Difference in event-free days",
  title = "A. Restricted mean survival time for treatment discontinuation\n    after first interruption of 3-4 days\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

# Add Legend
txt_rmst <- "Positive difference indicates more days retained in treatment.<br>CI: compatibility intervals"
p_rmst_final <- add_grob(p_rmst, row = nrow(plot_data_rmst)+2, col = 1:2, order = "background", 
                         gb_fn = gridtext::richtext_grob, text = txt_rmst, gp = gpar(fontsize = 11), 
                         hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/", "eFigure_A_12m_RMST.jpg"), width = 3500, height = 1400, res = 300)
p_rmst_final
dev.off()


## eFigure B: met 5-14: dis: rmst at 365 day ###########
risk_at <- 365

# Build Data for Trial 2
tr3_first_rmst <- extract_rmst_data(risk_tr3$first$dis$ci, outcome = "dis",
                                    dat_tr3_first$data_weights, risk_at, 
                                    cohort_label = "First episode", int_label = "5-14 days")
tr3_all_rmst   <- extract_rmst_data(risk_tr3$all$dis$ci, outcome = "dis", 
                                    dat_tr3_all$data_weights, risk_at, 
                                    cohort_label = "Any episode", int_label = "5-14 days")

plot_data_rmst <- bind_rows(tr3_first_rmst, tr3_all_rmst) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    RMST_diff = ifelse(trt_group == "No dose change", NA, RMST_diff),
    RMST_diff_lower = ifelse(trt_group == "No dose change", NA, RMST_diff_lower),
    RMST_diff_upper = ifelse(trt_group == "No dose change", NA, RMST_diff_upper),
    Cohort = ifelse(duplicated(Cohort), "", Cohort),
    ` ` = paste(rep(" ", 16), collapse = " ") # Empty canvas column
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode, n` = "Episode",
    `Event-free days\nat 12 months (95% CI)` = "rmst_disp",
    `Difference in days\n(95% CI)` = "diff_disp"
  )

# Generate Forest Plot
p_rmst <- forest(
  plot_data_rmst[, c("Cohort", "Dose adjustment\nafter interruption", "Episode, n", 
                     "Event-free days\nat 12 months (95% CI)", "Difference in days\n(95% CI)", " ")],
  est = plot_data_rmst$RMST_diff, lower = plot_data_rmst$RMST_diff_lower, upper = plot_data_rmst$RMST_diff_upper,
  ci_column = 6, 
  ref_line = 0,          # Additive null for difference
  x_trans = "none",
  xlim = c(-45, 45),     # IMPORTANT: Adjust this wider if any RMST bounds exceed +/- 45 days
  ticks_at = c(-40, -20, 0, 20, 40),
  xlab = "Difference in event-free days",
  title = "B. Restricted mean survival time for treatment discontinuation\n    after first interruption of 5-14 days\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

# Add Legend
txt_rmst <- "Positive difference indicates more days retained in treatment.<br>CI: compatibility ('confidence') intervals"
p_rmst_final <- add_grob(p_rmst, row = nrow(plot_data_rmst)+2, col = 1:2, order = "background", 
                         gb_fn = gridtext::richtext_grob, text = txt_rmst, gp = gpar(fontsize = 11), 
                         hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/", "eFigure_B_12m_RMST.jpg"), width = 3500, height = 1400, res = 300)
p_rmst_final
dev.off()




