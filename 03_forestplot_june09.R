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


# 2. FIGURE 2: RISK DIFFERENCES #################

risk_at_od <- 30

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
    `Episode at risk\nat reinitiation, n` = "Episode",
    `% risk at 1 month\n(95% CI)` = "abs_risk_disp",
    `Risk Difference (%)\n(95% CI)` = "rd_disp"
  )

p_2a <- forest(
  plot_data_2a[, c("Interruption", "Dose adjustment\nafter interruption", "Episode at risk\nat reinitiation, n", 
                   "% risk at 1 month\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est = plot_data_2a$rd_pct, lower = plot_data_2a$rd_lower_pct, upper = plot_data_2a$rd_upper_pct,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-2.5, 2.5), ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk difference (%)",
  title = "A. Difference in risk of overdose-related acute care visits\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

p_2a_final <- add_grob(p_2a, row = nrow(plot_data_2a) + 2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals", 
                       gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, 
                       x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/June09/", "Figure2_A_OD_RD.jpg"), width = 3500, height = 1400, res = 300)
p_2a_final
dev.off()


# FIGURE 2B: Discontinuation Risk Difference (365 Days) 
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
    `Episode at risk\nat reinitiation, n` = "Episode",
    `% risk at 12 months\n(95% CI)` = "abs_risk_disp",
    `Risk Difference (%)\n(95% CI)` = "rd_disp"
  )

p_2b <- forest(
  plot_data_2b[, c("Interruption", "Dose adjustment\nafter interruption", "Episode at risk\nat reinitiation, n", 
                   "% risk at 12 months\n(95% CI)", "Risk Difference (%)\n(95% CI)", " ")],
  est = plot_data_2b$rd_pct, lower = plot_data_2b$rd_lower_pct, upper = plot_data_2b$rd_upper_pct,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-12.5, 12.5), ticks_at = c(-10, -5, 0, 5, 10),
  xlab = "Risk difference (%)",
  title = "B. Difference in risk of treatment discontinuation\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

p_2b_final <- add_grob(p_2b, row = nrow(plot_data_2b) + 2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals", 
                       gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, 
                       x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/June09/", "Figure2_B_Discont_RD.jpg"), width = 3500, height = 1400, res = 300)
p_2b_final
dev.off()


# 3. FIGURE 3: RMST DIFFERENCES ###############
# --- FIGURE 3A: Overdose RMST ---

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
    `Episode at risk\nat reinitiation, n` = "Episode",
    `Event-free days\nat 12 months (95% CI)` = "rmst_disp",
    `Difference in days\n(95% CI)` = "diff_disp"
  )

p_3a <- forest(
  plot_data_3a[, c("Interruption", "Dose adjustment\nafter interruption", "Episode at risk\nat reinitiation, n", 
                   "Event-free days\nat 12 months (95% CI)", "Difference in days\n(95% CI)", " ")],
  est = plot_data_3a$RMST_diff, lower = plot_data_3a$RMST_diff_lower, upper = plot_data_3a$RMST_diff_upper,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-5, 5), ticks_at = c(-4, -2, 0, 2, 4), # Narrower bounds expected for rare OD
  xlab = "Difference in event-free days",
  title = "A. Restricted mean survival time for overdose-related acute care visits\n",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

txt_rmst <- "Positive difference indicates more event-free days.<br>CI: compatibility ('confidence') intervals"
p_3a_final <- add_grob(p_3a, row = nrow(plot_data_3a)+2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = txt_rmst, gp = gpar(fontsize = 11), 
                       hjust = 0, vjust = 1, halign = 0, valign = 1, x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/June09/", "Figure3_A_OD_RMST.jpg"), width = 3500, height = 1400, res = 300)
p_3a_final
dev.off()


# FIGURE 3B: Discontinuation RMST (365 Days) 
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
    `Episode at risk\nat reinitiation, n` = "Episode",
    `Event-free days\nat 12 months (95% CI)` = "rmst_disp",
    `Difference in days\n(95% CI)` = "diff_disp"
  )

p_3b <- forest(
  plot_data_3b[, c("Interruption", "Dose adjustment\nafter interruption", "Episode at risk\nat reinitiation, n", 
                   "Event-free days\nat 12 months (95% CI)", "Difference in days\n(95% CI)", " ")],
  est = plot_data_3b$RMST_diff, lower = plot_data_3b$RMST_diff_lower, upper = plot_data_3b$RMST_diff_upper,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-45, 45), ticks_at = c(-40, -20, 0, 20, 40),
  xlab = "Difference in event-free days",
  title = "Restricted mean survival time for treatment discontinuation\n",
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




# % for the text ##########
## overdose # #############

# Load survival package for unadjusted KM estimates
library(survival)
library(dplyr)
library(tidyr)

# Helper function to extract raw KM estimates and RMST
# Load survival package for unadjusted KM estimates
library(survival)
library(dplyr)
library(tidyr)

# Helper function to extract raw KM estimates, event counts, and RMST
compute_raw_km_stats <- function(dat, time_var, event_var, trt_var = "trt", tau = 365) {
  
  # Ensure treatment is a factor
  dat[[trt_var]] <- as.factor(dat[[trt_var]])
  
  # Formula for survival model
  form <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", trt_var))
  
  # Fit unadjusted Kaplan-Meier
  fit <- survfit(form, data = dat)
  
  # 1. Extract Raw Cumulative Incidence (1 - Survival) at 30 and 365 days
  sum_fit <- summary(fit, times = c(30, 365), extend = TRUE)
  
  risks <- data.frame(
    Treatment = gsub(paste0(trt_var, "="), "", sum_fit$strata),
    Time = sum_fit$time,
    Risk_pct = (1 - sum_fit$surv) * 100
  ) %>%
    pivot_wider(names_from = Time, values_from = Risk_pct, names_prefix = "Risk_pct_") %>%
    rename(
      `Raw Risk at 30 days (%)` = Risk_pct_30,
      `Raw Risk at 365 days (%)` = Risk_pct_365
    )
  
  # 2. Extract Raw Event Counts for the specific time windows
  event_counts <- dat %>%
    group_by(!!sym(trt_var)) %>%
    summarise(
      `Events at 30 days` = sum(!!sym(time_var) <= 30 & !!sym(event_var) == 1, na.rm = TRUE),
      `Events at 365 days` = sum(!!sym(time_var) <= 365 & !!sym(event_var) == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(Treatment = !!sym(trt_var)) %>%
    mutate(Treatment = as.character(Treatment))
  
  # 3. Extract Raw RMST (Restricted Mean Survival Time) at tau
  km_table <- as.data.frame(summary(fit, rmean = tau)$table)
  
  rmst_data <- data.frame(
    Treatment = gsub(paste0(trt_var, "="), "", rownames(km_table)),
    `Total Patients` = km_table$records,
    `Total Events (Any Time)` = km_table$events,
    `Raw RMST (Days)` = km_table$rmean,
    `RMST SE` = km_table$`se(rmean)`,
    check.names = FALSE
  )
  
  # 4. Combine and reorder into one neat table
  final_stats <- rmst_data %>%
    left_join(event_counts, by = "Treatment") %>%
    left_join(risks, by = "Treatment") %>%
    mutate(
      Treatment = recode(Treatment,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                         "increase" = "Dose increase")
    ) %>%
    # Reorder columns so events and percentages are side-by-side
    select(Treatment, `Total Patients`, 
           `Events at 30 days`, `Raw Risk at 30 days (%)`, 
           `Events at 365 days`, `Raw Risk at 365 days (%)`, 
           `Total Events (Any Time)`, `Raw RMST (Days)`, `RMST SE`)
  
  return(final_stats)
}
cat("\n=================================================================\n")
cat(" RAW (UNADJUSTED) OVERDOSE STATS: 'ANY EPISODE' COHORT\n")
cat("=================================================================\n")

cat("\n--- Trial 2: Methadone Interruption of 3-4 days ---\n")
raw_od_tr2 <- compute_raw_km_stats(dat_tr2_all$data_weights, 
                                   time_var = "time_od_acute", 
                                   event_var = "sens_od_acute")
print(raw_od_tr2, digits = 3, row.names = FALSE)

cat("\n--- Trial 3: Methadone Interruption of 5-14 days ---\n")
raw_od_tr3 <- compute_raw_km_stats(dat_tr3_all$data_weights, 
                                   time_var = "time_od_acute", 
                                   event_var = "sens_od_acute")
print(raw_od_tr3, digits = 3, row.names = FALSE)



## KM curves for overdose ######

library(tidyverse)
library(survival)
library(broom)
library(patchwork)

# 1. Define color palette
trt_colors_km <- c(
  "No dose change" = "#000000",                   
  "<50% reduction" = "#2CA02C",                   
  "\u226550% reduction/\u226440mg" = "#FF7F0E",  
  "Dose increase" = "#D62728"                     
)

# 2. Helper function: Extract and format Kaplan-Meier data (Main Plot)
build_km_df <- function(dat, time_var = "time_od_acute", event_var = "sens_od_acute") {
  dat <- dat %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "\u226550% reduction/\u226440mg",
                         "increase" = "Dose increase"),
      trt_label = factor(trt_label, levels = names(trt_colors_km))
    )
  
  form <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ trt_label"))
  fit <- survfit(form, data = dat)
  
  fit_df <- tidy(fit) %>%
    mutate(
      cum_inc = 1 - estimate,
      trt_label = gsub("trt_label=", "", strata), 
      trt_label = factor(trt_label, levels = names(trt_colors_km))
    ) %>% select(time, cum_inc, trt_label)
  
  t0 <- expand.grid(time = 0, cum_inc = 0, trt_label = names(trt_colors_km))
  
  fit_df <- bind_rows(t0, fit_df) %>% arrange(trt_label, time)
  return(fit_df)
}

# 3. Helper function: Generate Cumulative Event Table
build_event_table <- function(dat, time_var = "time_od_acute", event_var = "sens_od_acute") {
  dat <- dat %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "\u226550% reduction/\u226440mg",
                         "increase" = "Dose increase"),
      # Reverse factor levels so "No dose change" sits at the top of the table
      trt_label = factor(trt_label, levels = rev(names(trt_colors_km))) 
    )
  
  breaks_x <- c(0, 30, 90, 180, 270, 365)
  
  # Calculate cumulative events at exactly the specified breaks
  table_df <- purrr::map_dfr(breaks_x, function(t) {
    dat %>%
      group_by(trt_label) %>%
      summarise(
        time = t,
        # Explicitly force 0 events at t=0
        raw_events = ifelse(t == 0, 0, sum(!!sym(time_var) <= t & !!sym(event_var) == 1, na.rm = TRUE)),
        .groups = "drop"
      )
  }) %>%
    mutate(
      # Safely apply privacy suppression: Keep 0 at day 0. Mask 1-9. Show >=10.
      events_label = case_when(
        time == 0 ~ "0",
        raw_events > 0 & raw_events < 10 ~ "<10",
        TRUE ~ as.character(raw_events)
      )
    )
  
  # Create text-based ggplot table
  ggplot(table_df, aes(x = time, y = trt_label, label = events_label)) +
    geom_text(size = 4.5, color = "black") +
    scale_x_continuous(limits = c(0, 365), breaks = breaks_x) +
    labs(y = NULL, x = "Cumulative events") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_blank(),
      axis.ticks   = element_blank(),
      axis.text.x  = element_blank(), # Hide x-axis text since main plot has it
      #axis.title.x = element_blank(),
      axis.text.y  = element_text(color = "black", size = 12),
      axis.title.y = element_text(color = "black", size = 12, face = "bold", angle = 0, vjust = 0.5),
      plot.margin  = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")
    )
}

# 4. Helper function: Plot the main KM curve
plot_km_curve <- function(df, title) {
  ggplot(df, aes(x = time, y = cum_inc, color = trt_label)) +
    geom_step(linewidth = 1) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "gray60", linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 365), breaks = c(0, 30, 90, 180, 270, 365)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                       limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
    scale_color_manual(values = trt_colors_km) +
    labs(title = title, x = "Days since time zero\n", y = "Unadjusted Cumulative Incidence (%)") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(color = "black", size = 13),
      legend.position   = "none", # Legend will be handled by patchwork
      plot.margin       = margin(t = 10, r = 5, b = 0, l = 5, unit = "pt")
    )
}

# 5. Generate data, curves, and tables for both trials
df_km_tr2 <- build_km_df(dat_tr2_all$data_weights)
df_km_tr3 <- build_km_df(dat_tr3_all$data_weights)

p_km_tr2  <- plot_km_curve(df_km_tr2, "A. Methadone interruption of 3-4 days")
p_tab_tr2 <- build_event_table(dat_tr2_all$data_weights)

p_km_tr3  <- plot_km_curve(df_km_tr3, "B. Methadone interruption of 5-14 days")
p_tab_tr3 <- build_event_table(dat_tr3_all$data_weights)

# 6. Extract legend from one of the plots to place at the very bottom
dummy_plot <- p_km_tr2 + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size = 13), legend.key.width = unit(2, "cm"))
shared_legend <- cowplot::get_legend(dummy_plot)

# 7. Assemble the vertical layout using Patchwork (FLATTENED)
final_km_composite <- p_km_tr2 / 
  p_tab_tr2 / 
  plot_spacer() / 
  p_km_tr3 / 
  p_tab_tr3 / 
  shared_legend +
  plot_layout(heights = c(5, 1, 0.5, 5, 1, 0.5)) 

# 8. Save output
jpeg(paste0(RES_DIR, "Figures/June09/eFigure4_Unadjusted_OD_KM_Vertical.jpg"), 
     width = 3000, height = 4000, res = 300)
print(final_km_composite)
dev.off()

## discontinuation ##########







cat("\n=================================================================\n")
cat(" RAW (UNADJUSTED) DISCONTINUATION STATS: 'ANY EPISODE' COHORT\n")
cat("=================================================================\n")

cat("\n--- Trial 2: Methadone Interruption of 3-4 days ---\n")
raw_dis_tr2 <- compute_raw_km_stats(dat_tr2_all$data_weights, 
                                    time_var = "time_discont", 
                                    event_var = "sens_discont")
print(raw_dis_tr2, digits = 4, row.names = FALSE)

cat("\n--- Trial 3: Methadone Interruption of 5-14 days ---\n")
raw_dis_tr3 <- compute_raw_km_stats(dat_tr3_all$data_weights, 
                                    time_var = "time_discont", 
                                    event_var = "sens_discont")
print(raw_dis_tr3, digits = 4, row.names = FALSE)


library(purrr)

cat("\n========================================================\n")
cat(" MANUSCRIPT TEXT EXTRACTION (COPY AND PASTE INTO WORD)\n")
cat("========================================================\n")

# Helper function to generate sentences
print_manuscript_text <- function(df, metric_name, estimate_col) {
  df %>%
    # Filter out the reference group so we only get the active comparisons
    filter(`Dose adjustment` != "No dose change") %>%
    mutate(
      # Clean up the stratum labels for text
      clean_stratum = ifelse(Stratum == "", NA, Stratum)
    ) %>%
    fill(clean_stratum, .direction = "down") %>%
    mutate(
      sentence = sprintf(
        "For patients with %s retention, following a %s, the %s for the '%s' group compared to no dose change was %s.",
        clean_stratum, 
        tolower(Interruption), 
        metric_name, 
        `Dose adjustment`, 
        !!sym(estimate_col)
      )
    ) %>%
    pull(sentence) %>%
    walk(~cat("  -", ., "\n"))
}

cat("\n--- 1. OVERDOSE RISK DIFFERENCES (90 DAYS) ---\n")
print_manuscript_text(
  data_2a, 
  metric_name = "absolute risk difference", 
  estimate_col = "Risk Difference (%)\n(95% CI)"
)

cat("\n--- 2. DISCONTINUATION RISK DIFFERENCES (365 DAYS) ---\n")
print_manuscript_text(
  data_2b, 
  metric_name = "absolute risk difference", 
  estimate_col = "Risk Difference (%)\n(95% CI)"
)

cat("\n--- 3. DISCONTINUATION RMST DIFFERENCES (365 DAYS) ---\n")
print_manuscript_text(
  data_3b, 
  metric_name = "difference in event-free days", 
  estimate_col = "Difference in days\n(95% CI)"
)
cat("\n========================================================\n")


