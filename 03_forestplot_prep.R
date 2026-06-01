# FOREST PLOT PREPARATION: Cumulative Incidence & Rates
# Missed-dose trials: MET interruption of 3-4 days and 5-14 days
# Outcomes: Overdose (1 month) and Discontinuation (12 months)

library(tidyverse)

# LOAD RISK RESULTS

RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

# UNIFIED EXTRACTION: Absolute Risk and Risk Difference at time t

extract_risk_and_rd <- function(risk_obj, timepoint = 90, ref_name = "no reversion") {
  
  # 1. Extract absolute risks
  risks_t <- risk_obj$ci$risks %>%
    filter(horizon == timepoint) %>%
    mutate(
      risk_pct = risk * 100,
      risk_lower_pct = risk_lower * 100,
      risk_upper_pct = risk_upper * 100,
      # Format dynamically based on magnitude of the point estimate
      abs_risk_display = ifelse(
        abs(risk_pct) >= 10,
        sprintf("%.1f (%.1f, %.1f)", risk_pct, risk_lower_pct, risk_upper_pct),
        sprintf("%.2f (%.2f, %.2f)", risk_pct, risk_lower_pct, risk_upper_pct)
      )
    ) %>%
    select(trt, abs_risk_display, risk_pct)
  
  # 2. Extract contrasts (Risk Differences)
  contrasts_t <- risk_obj$ci$contrasts %>%
    filter(horizon == timepoint) %>%
    mutate(
      rd_pct = RD * 100,
      rd_lower_pct = RD_lower * 100,
      rd_upper_pct = RD_upper * 100,
      # Use absolute value to correctly format negative risk differences
      rd_display = ifelse(
        abs(rd_pct) >= 10,
        sprintf("%.1f (%.1f, %.1f)", rd_pct, rd_lower_pct, rd_upper_pct),
        sprintf("%.2f (%.2f, %.2f)", rd_pct, rd_lower_pct, rd_upper_pct)
      )
    ) %>%
    select(trt, rd_pct, rd_lower_pct, rd_upper_pct, rd_display)
  
  # 3. Merge and handle the reference group
  merged_data <- risks_t %>%
    left_join(contrasts_t, by = "trt") %>%
    mutate(
      # The reference group will have NA for RD from the left join. Replace with "Ref."
      rd_display = if_else(trt == ref_name, "Ref.", rd_display),
      rd_pct = if_else(trt == ref_name, 0, rd_pct),
      rd_lower_pct = if_else(trt == ref_name, 0, rd_lower_pct),
      rd_upper_pct = if_else(trt == ref_name, 0, rd_upper_pct)
    ) %>%
    rename(trt_group = trt)
  
  return(merged_data)
}


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



risk_tr2 <- readRDS(paste0(RES_DIR, "Risk/risk_tr2_34d.rds"))
risk_tr3 <- readRDS(paste0(RES_DIR, "Risk/risk_tr3_514d.rds"))

dat_tr2_first <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))
dat_tr2_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
dat_tr3_first <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))
dat_tr3_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))


# A. MET interruption of 3-4 days ##############
## at 30 days: od ##################
risk_at = 30

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr2_first <- extract_risk_and_rd(risk_tr2$first$od, timepoint = risk_at, ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr2_first <- compute_summary_stats_timewin(
  dat_tr2_first$data_weights,
  time_var = "time_od_acute",
  event_var = "sens_od_acute",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)

# Final preparatory dataset for this stratum
tr2_plot_first <- counts_tr2_first %>%
  left_join(risk_rd_tr2_first, by = "trt_group") %>%
  mutate(
    Cohort = "First episode", 
    outcome = "Overdose", 
    Interruption = "3-4 days"
  )

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr2_all <- extract_risk_and_rd(risk_tr2$all$od, 
                                       timepoint = risk_at, 
                                       ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr2_all <- compute_summary_stats_timewin(
  dat_tr2_all$data_weights,
  time_var = "time_od_acute",
  event_var = "sens_od_acute",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)



# Final preparatory dataset for this stratum
tr2_plot_all <- counts_tr2_all %>%
  left_join(risk_rd_tr2_all, by = "trt_group") %>%
  mutate(
    Cohort = "Any episode", 
    outcome = "Overdose", 
    Interruption = "3-4 days"
  )


library(stringr)

# 1. Combine and format
tr2_plot <- bind_rows(tr2_plot_first, tr2_plot_all) %>% 
  select(Cohort, Interruption, trt_group, Episode, abs_risk_display, 
         rd_display, rd_pct, rd_lower_pct, rd_upper_pct) %>% 
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct)
  ) %>% 
  rename(
    `Dose adjustment\nafter interruption` = "trt_group", 
    `Episode, n` = "Episode",
    `% risk at 1 month\n(95% CI)` = "abs_risk_display",
    `% difference in 1-month\nrisk (95% CI)` = "rd_display",
    `First interruption` = "Interruption"
  )

# 2. Add an explicit, empty string column for the forest plot space
tr2_plot$` ` <- paste(rep(" ", 16), collapse = " ")

# 3. Clean repetitive labels for the table display
tr2_plot$`Cohort` <- ifelse(duplicated(tr2_plot$`Cohort`), "", tr2_plot$`Cohort`)
tr2_plot$`First interruption` <- ifelse(duplicated(tr2_plot$`First interruption`), "", tr2_plot$`First interruption`)


# 4. Generate the plot
p_dis_tr2 <- forest(
  tr2_plot[, c("Cohort", "Dose adjustment\nafter interruption", 
               "Episode, n", "% risk at 1 month\n(95% CI)", "% difference in 1-month\nrisk (95% CI)", " ")],
  est   = tr2_plot$rd_pct,
  lower = tr2_plot$rd_lower_pct,
  upper = tr2_plot$rd_upper_pct,
  ci_column = 6, # Points to the blank string column
  ref_line = 0,
  sizes = 0.6,
  x_trans = "none",
  xlim = c(-2.5, 2.5),
  ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk difference (%)",
  arrow_lab = NULL, # More precise clinical language
  title = "A. Comparative risk of overdose-related acute care visits among dose adjustment groups\n    after first interruption of 3-4 days\n",
  nudge_y = 0.02,
  theme = forest_theme(
    base_size = 14,
    ci_pch = 15,
    ci_lwd = 1.7,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
    arrow_length = 0.05,
    arrow_gp = gpar(fontsize = 11, lwd = 0.7),
    xlab_adjust = c("refline", "center")
  )
)



txt <- "Risk. cumulative incidence<br>CI. compatibility ('confidence') intervals"

p_dis_tr2_1 <- add_grob(p_dis_tr2, 
                        row = nrow(tr2_plot)+2, 
                        col = 1:2,
                        order = "background",
                        gb_fn = gridtext::richtext_grob,
                        text = txt,
                        gp = gpar(fontsize = 11),
                        hjust = 0, vjust = 1, halign = 0, valign = 1,
                        x = unit(0, "npc"), y = unit(1, "npc"))


jpeg(paste0(RES_DIR, "Figures/", "Figure2_A_1m.jpg"),
     width = 3500,  height = 1400, res = 300)
p_dis_tr2_1
dev.off()




## at 90 days: od ##################
risk_at = 90

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr2_first <- extract_risk_and_rd(risk_tr2$first$od, 
                                         timepoint = risk_at, 
                                         ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr2_first <- compute_summary_stats_timewin(
  dat_tr2_first$data_weights,
  time_var = "time_od_acute",
  event_var = "sens_od_acute",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)

# Final preparatory dataset for this stratum
tr2_plot_first <- counts_tr2_first %>%
  left_join(risk_rd_tr2_first, by = "trt_group") %>%
  mutate(
    Cohort = "First episode", 
    outcome = "Overdose", 
    Interruption = "3-4 days"
  )

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr2_all <- extract_risk_and_rd(risk_tr2$all$od, 
                                       timepoint = risk_at, 
                                       ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr2_all <- compute_summary_stats_timewin(
  dat_tr2_all$data_weights,
  time_var = "time_od_acute",
  event_var = "sens_od_acute",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)



# Final preparatory dataset for this stratum
tr2_plot_all <- counts_tr2_all %>%
  left_join(risk_rd_tr2_all, by = "trt_group") %>%
  mutate(
    Cohort = "Any episode", 
    outcome = "Overdose", 
    Interruption = "3-4 days"
  )


library(stringr)

# 1. Combine and format
tr2_plot <- bind_rows(tr2_plot_first, tr2_plot_all) %>% 
  select(Cohort, Interruption, trt_group, Episode, abs_risk_display, 
         rd_display, rd_pct, rd_lower_pct, rd_upper_pct) %>% 
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct)
  ) %>% 
  rename(
    `Dose adjustment\nafter interruption` = "trt_group", 
    `Episode, n` = "Episode",
    `% risk at 3 months\n(95% CI)` = "abs_risk_display",
    `% difference in 3-month\nrisk (95% CI)` = "rd_display",
    `First interruption` = "Interruption"
  )

# 2. Add an explicit, empty string column for the forest plot space
tr2_plot$` ` <- paste(rep(" ", 16), collapse = " ")

# 3. Clean repetitive labels for the table display
tr2_plot$`Cohort` <- ifelse(duplicated(tr2_plot$`Cohort`), "", tr2_plot$`Cohort`)
tr2_plot$`First interruption` <- ifelse(duplicated(tr2_plot$`First interruption`), "", tr2_plot$`First interruption`)


# 4. Generate the plot
p_dis_tr2 <- forest(
  tr2_plot[, c("Cohort", "Dose adjustment\nafter interruption", 
               "Episode, n", "% risk at 3 months\n(95% CI)", "% difference in 3-month\nrisk (95% CI)", " ")],
  est   = tr2_plot$rd_pct,
  lower = tr2_plot$rd_lower_pct,
  upper = tr2_plot$rd_upper_pct,
  ci_column = 6, # Points to the blank string column
  ref_line = 0,
  sizes = 0.6,
  x_trans = "none",
  xlim = c(-2.5, 2.5),
  ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk difference (%)",
  arrow_lab = NULL, # More precise clinical language
  title = "A. Comparative risk of overdose-related acute care visits among dose adjustment groups\n    after first interruption of 3-4 days\n",
  nudge_y = 0.02,
  theme = forest_theme(
    base_size = 14,
    ci_pch = 15,
    ci_lwd = 1.7,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
    arrow_length = 0.05,
    arrow_gp = gpar(fontsize = 11, lwd = 0.7),
    xlab_adjust = c("refline", "center")
  )
)



txt <- "Risk. cumulative incidence<br>CI. compatibility intervals"

p_dis_tr2_1 <- add_grob(p_dis_tr2, 
                        row = nrow(tr2_plot)+2, 
                        col = 1:2,
                        order = "background",
                        gb_fn = gridtext::richtext_grob,
                        text = txt,
                        gp = gpar(fontsize = 11),
                        hjust = 0, vjust = 1, halign = 0, valign = 1,
                        x = unit(0, "npc"), y = unit(1, "npc"))


jpeg(paste0(RES_DIR, "Figures/", "Figure2_A_3m.jpg"),
     width = 3500,  height = 1400, res = 300)
p_dis_tr2_1
dev.off()



# B. MET interruption of 5-14 days ##############
## at 30 days: od ##################
risk_at = 30

### MET 5-14d: First epiode ###########
risk_rd_tr3_first <- extract_risk_and_rd(risk_tr3$first$od, timepoint = risk_at, ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr3_first <- compute_summary_stats_timewin(
  dat_tr3_first$data_weights,
  time_var = "time_od_acute",
  event_var = "sens_od_acute",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)

# Final preparatory dataset for this stratum
tr3_plot_first <- counts_tr3_first %>%
  left_join(risk_rd_tr3_first, by = "trt_group") %>%
  mutate(
    Cohort = "First episode", 
    outcome = "Overdose", 
    Interruption = "5-14 days"
  )

### MET 5-14d: First interruption in any epiode ###########
risk_rd_tr3_all <- extract_risk_and_rd(risk_tr3$all$od, 
                                       timepoint = risk_at, 
                                       ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr3_all <- compute_summary_stats_timewin(
  dat_tr3_all$data_weights,
  time_var = "time_od_acute",
  event_var = "sens_od_acute",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)



# Final preparatory dataset for this stratum
tr3_plot_all <- counts_tr3_all %>%
  left_join(risk_rd_tr3_all, by = "trt_group") %>%
  mutate(
    Cohort = "Any episode", 
    outcome = "Overdose", 
    Interruption = "5-14 days"
  )


library(stringr)

# 1. Combine and format
tr3_plot <- bind_rows(tr3_plot_first, tr3_plot_all) %>% 
  select(Cohort, Interruption, trt_group, Episode, abs_risk_display, 
         rd_display, rd_pct, rd_lower_pct, rd_upper_pct) %>% 
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct)
  ) %>% 
  rename(
    `Dose adjustment\nafter interruption` = "trt_group", 
    `Episode, n` = "Episode",
    `% risk at 1 month\n(95% CI)` = "abs_risk_display",
    `% difference in 1-month\nrisk (95% CI)` = "rd_display",
    `First interruption` = "Interruption"
  )

# 2. Add an explicit, empty string column for the forest plot space
tr3_plot$` ` <- paste(rep(" ", 16), collapse = " ")

# 3. Clean repetitive labels for the table display
tr3_plot$`Cohort` <- ifelse(duplicated(tr3_plot$`Cohort`), "", tr3_plot$`Cohort`)
tr3_plot$`First interruption` <- ifelse(duplicated(tr3_plot$`First interruption`), "", tr3_plot$`First interruption`)


# 4. Generate the plot
p_dis_tr3 <- forest(
  tr3_plot[, c("Cohort", "Dose adjustment\nafter interruption", 
               "Episode, n", "% risk at 1 month\n(95% CI)", "% difference in 1-month\nrisk (95% CI)", " ")],
  est   = tr3_plot$rd_pct,
  lower = tr3_plot$rd_lower_pct,
  upper = tr3_plot$rd_upper_pct,
  ci_column = 6, # Points to the blank string column
  ref_line = 0,
  sizes = 0.6,
  x_trans = "none",
  xlim = c(-2.5, 2.5),
  ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk difference (%)",
  arrow_lab = NULL, # More precise clinical language
  title = "B. Comparative risk of overdose-related acute care visits among dose adjustment groups\n    after first interruption of 5-14 days\n",
  nudge_y = 0.02,
  theme = forest_theme(
    base_size = 14,
    ci_pch = 15,
    ci_lwd = 1.7,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
    arrow_length = 0.05,
    arrow_gp = gpar(fontsize = 11, lwd = 0.7),
    xlab_adjust = c("refline", "center")
  )
)



txt <- "Risk. cumulative incidence<br>CI. compatibility ('confidence') intervals"

p_dis_tr3_1 <- add_grob(p_dis_tr3, 
                        row = nrow(tr3_plot)+2, 
                        col = 1:2,
                        order = "background",
                        gb_fn = gridtext::richtext_grob,
                        text = txt,
                        gp = gpar(fontsize = 11),
                        hjust = 0, vjust = 1, halign = 0, valign = 1,
                        x = unit(0, "npc"), y = unit(1, "npc"))


jpeg(paste0(RES_DIR, "Figures/", "Figure2_B_1m.jpg"),
     width = 3500,  height = 1400, res = 300)
p_dis_tr3_1
dev.off()




## at 90 days: od ##################
risk_at = 90

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr3_first <- extract_risk_and_rd(risk_tr3$first$od, 
                                         timepoint = risk_at, 
                                         ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr3_first <- compute_summary_stats_timewin(
  dat_tr3_first$data_weights,
  time_var = "time_od_acute",
  event_var = "sens_od_acute",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)

# Final preparatory dataset for this stratum
tr3_plot_first <- counts_tr3_first %>%
  left_join(risk_rd_tr3_first, by = "trt_group") %>%
  mutate(
    Cohort = "First episode", 
    outcome = "Overdose", 
    Interruption = "5-14 days"
  )

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr3_all <- extract_risk_and_rd(risk_tr3$all$od, 
                                       timepoint = risk_at, 
                                       ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr3_all <- compute_summary_stats_timewin(
  dat_tr3_all$data_weights,
  time_var = "time_od_acute",
  event_var = "sens_od_acute",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)



# Final preparatory dataset for this stratum
tr3_plot_all <- counts_tr3_all %>%
  left_join(risk_rd_tr3_all, by = "trt_group") %>%
  mutate(
    Cohort = "Any episode", 
    outcome = "Overdose", 
    Interruption = "5-14 days"
  )


library(stringr)

# 1. Combine and format
tr3_plot <- bind_rows(tr3_plot_first, tr3_plot_all) %>% 
  select(Cohort, Interruption, trt_group, Episode, abs_risk_display, 
         rd_display, rd_pct, rd_lower_pct, rd_upper_pct) %>% 
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct)
  ) %>% 
  rename(
    `Dose adjustment\nafter interruption` = "trt_group", 
    `Episode, n` = "Episode",
    `% risk at 3 months\n(95% CI)` = "abs_risk_display",
    `% difference in 3-month\nrisk (95% CI)` = "rd_display",
    `First interruption` = "Interruption"
  )

# 2. Add an explicit, empty string column for the forest plot space
tr3_plot$` ` <- paste(rep(" ", 16), collapse = " ")

# 3. Clean repetitive labels for the table display
tr3_plot$`Cohort` <- ifelse(duplicated(tr3_plot$`Cohort`), "", tr3_plot$`Cohort`)
tr3_plot$`First interruption` <- ifelse(duplicated(tr3_plot$`First interruption`), "", tr3_plot$`First interruption`)


# 4. Generate the plot
p_dis_tr3 <- forest(
  tr3_plot[, c("Cohort", "Dose adjustment\nafter interruption", 
               "Episode, n", "% risk at 3 months\n(95% CI)", "% difference in 3-month\nrisk (95% CI)", " ")],
  est   = tr3_plot$rd_pct,
  lower = tr3_plot$rd_lower_pct,
  upper = tr3_plot$rd_upper_pct,
  ci_column = 6, # Points to the blank string column
  ref_line = 0,
  sizes = 0.6,
  x_trans = "none",
  xlim = c(-2.5, 2.5),
  ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk difference (%)",
  arrow_lab = NULL, # More precise clinical language
  title = "B. Comparative risk of overdose-related acute care visits among dose adjustment groups\n    after first interruption of 5-14 days\n",
  nudge_y = 0.02,
  theme = forest_theme(
    base_size = 14,
    ci_pch = 15,
    ci_lwd = 1.7,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
    arrow_length = 0.05,
    arrow_gp = gpar(fontsize = 11, lwd = 0.7),
    xlab_adjust = c("refline", "center")
  )
)



txt <- "Risk. cumulative incidence<br>CI. compatibility intervals"

p_dis_tr3_1 <- add_grob(p_dis_tr3, 
                        row = nrow(tr3_plot)+2, 
                        col = 1:2,
                        order = "background",
                        gb_fn = gridtext::richtext_grob,
                        text = txt,
                        gp = gpar(fontsize = 11),
                        hjust = 0, vjust = 1, halign = 0, valign = 1,
                        x = unit(0, "npc"), y = unit(1, "npc"))


jpeg(paste0(RES_DIR, "Figures/", "Figure2_B_3m.jpg"),
     width = 3500,  height = 1400, res = 300)
p_dis_tr3_1
dev.off()



# A. MET interruption of 3-4 days ################
## at 365 days: dis ##################
risk_at = 365

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr2_first <- extract_risk_and_rd(risk_tr2$first$dis, 
                                         timepoint = risk_at, 
                                         ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr2_first <- compute_summary_stats_timewin(
  dat_tr2_first$data_weights,
  time_var = "time_discont",
  event_var = "sens_discont",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)

# Final preparatory dataset for this stratum
tr2_plot_first <- counts_tr2_first %>%
  left_join(risk_rd_tr2_first, by = "trt_group") %>%
  mutate(
    Cohort = "First episode", 
    outcome = "Overdose", 
    Interruption = "3-4 days"
  )

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr2_all <- extract_risk_and_rd(risk_tr2$all$dis, 
                                       timepoint = risk_at, 
                                       ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr2_all <- compute_summary_stats_timewin(
  dat_tr2_all$data_weights,
  time_var = "time_discont",
  event_var = "sens_discont",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)



# Final preparatory dataset for this stratum
tr2_plot_all <- counts_tr2_all %>%
  left_join(risk_rd_tr2_all, by = "trt_group") %>%
  mutate(
    Cohort = "Any episode", 
    outcome = "Overdose", 
    Interruption = "3-4 days"
  )


library(stringr)

# 1. Combine and format
tr2_plot <- bind_rows(tr2_plot_first, tr2_plot_all) %>% 
  select(Cohort, Interruption, trt_group, Episode, abs_risk_display, 
         rd_display, rd_pct, rd_lower_pct, rd_upper_pct) %>% 
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct)
  ) %>% 
  rename(
    `Dose adjustment\nafter interruption` = "trt_group", 
    `Episode, n` = "Episode",
    `% risk at 12 months\n(95% CI)` = "abs_risk_display",
    `% difference in 12-month\nrisk (95% CI)` = "rd_display",
    `First interruption` = "Interruption"
  )

# 2. Add an explicit, empty string column for the forest plot space
tr2_plot$` ` <- paste(rep(" ", 16), collapse = " ")

# 3. Clean repetitive labels for the table display
tr2_plot$`Cohort` <- ifelse(duplicated(tr2_plot$`Cohort`), "", tr2_plot$`Cohort`)
tr2_plot$`First interruption` <- ifelse(duplicated(tr2_plot$`First interruption`), "", tr2_plot$`First interruption`)

library(forestploter)

# 4. Generate the plot
p_dis_tr2 <- forest(
  tr2_plot[, c("Cohort", "Dose adjustment\nafter interruption", 
               "Episode, n", "% risk at 12 months\n(95% CI)", "% difference in 12-month\nrisk (95% CI)", " ")],
  est   = tr2_plot$rd_pct,
  lower = tr2_plot$rd_lower_pct,
  upper = tr2_plot$rd_upper_pct,
  ci_column = 6, # Points to the blank string column
  ref_line = 0,
  sizes = 0.6,
  x_trans = "none",
  xlim = c(-12.5, 12.5),
  ticks_at = c(-10, -5, 0, 5, 10),
  xlab = "Risk difference (%)",
  arrow_lab = NULL, # More precise clinical language
  title = "A. Comparative risk of treatment discontinuation among dose adjustment groups\n    after first interruption of 3-4 days\n",
  nudge_y = 0.02,
  theme = forest_theme(
    base_size = 14,
    ci_pch = 15,
    ci_lwd = 1.7,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
    arrow_length = 0.05,
    arrow_gp = gpar(fontsize = 11, lwd = 0.7),
    xlab_adjust = c("refline", "center")
  )
)



txt <- "Risk. cumulative incidence<br>CI. compatibility intervals"

p_dis_tr2_1 <- add_grob(p_dis_tr2, 
                        row = nrow(tr2_plot)+2, 
                        col = 1:2,
                        order = "background",
                        gb_fn = gridtext::richtext_grob,
                        text = txt,
                        gp = gpar(fontsize = 11),
                        hjust = 0, vjust = 1, halign = 0, valign = 1,
                        x = unit(0, "npc"), y = unit(1, "npc"))


jpeg(paste0(RES_DIR, "Figures/", "Figure3_A_12m.jpg"),
     width = 3500,  height = 1400, res = 300)
p_dis_tr2_1
dev.off()




# B. MET interruption of 5-14 days ################
## at 365 days: dis ##################
risk_at = 365

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr3_first <- extract_risk_and_rd(risk_tr3$first$dis, 
                                         timepoint = risk_at, 
                                         ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr3_first <- compute_summary_stats_timewin(
  dat_tr3_first$data_weights,
  time_var = "time_discont",
  event_var = "sens_discont",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)

# Final preparatory dataset for this stratum
tr3_plot_first <- counts_tr3_first %>%
  left_join(risk_rd_tr3_first, by = "trt_group") %>%
  mutate(
    Cohort = "First episode", 
    outcome = "Overdose", 
    Interruption = "5-14 days"
  )

### MET 3-4d: First interruption in any epiode ###########
risk_rd_tr3_all <- extract_risk_and_rd(risk_tr3$all$dis, 
                                       timepoint = risk_at, 
                                       ref_name = "no reversion")

# Get baseline counts ONLY (drop the confusing Rate_100PY)
counts_tr3_all <- compute_summary_stats_timewin(
  dat_tr3_all$data_weights,
  time_var = "time_discont",
  event_var = "sens_discont",
  weight_var = "wt_trunc",
  trt_var = "trt",
  max_time = risk_at
) %>% 
  select(trt_group = trt, Episode, Event)



# Final preparatory dataset for this stratum
tr3_plot_all <- counts_tr3_all %>%
  left_join(risk_rd_tr3_all, by = "trt_group") %>%
  mutate(
    Cohort = "Any episode", 
    outcome = "Overdose", 
    Interruption = "5-14 days"
  )


library(stringr)

# 1. Combine and format
tr3_plot <- bind_rows(tr3_plot_first, tr3_plot_all) %>% 
  select(Cohort, Interruption, trt_group, Episode, abs_risk_display, 
         rd_display, rd_pct, rd_lower_pct, rd_upper_pct) %>% 
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct)
  ) %>% 
  rename(
    `Dose adjustment\nafter interruption` = "trt_group", 
    `Episode, n` = "Episode",
    `% risk at 12 months\n(95% CI)` = "abs_risk_display",
    `% difference in 12-month\nrisk (95% CI)` = "rd_display",
    `First interruption` = "Interruption"
  )

# 2. Add an explicit, empty string column for the forest plot space
tr3_plot$` ` <- paste(rep(" ", 16), collapse = " ")

# 3. Clean repetitive labels for the table display
tr3_plot$`Cohort` <- ifelse(duplicated(tr3_plot$`Cohort`), "", tr3_plot$`Cohort`)
tr3_plot$`First interruption` <- ifelse(duplicated(tr3_plot$`First interruption`), "", tr3_plot$`First interruption`)

library(forestploter)

# 4. Generate the plot
p_dis_tr3 <- forest(
  tr3_plot[, c("Cohort", "Dose adjustment\nafter interruption", 
               "Episode, n", "% risk at 12 months\n(95% CI)", "% difference in 12-month\nrisk (95% CI)", " ")],
  est   = tr3_plot$rd_pct,
  lower = tr3_plot$rd_lower_pct,
  upper = tr3_plot$rd_upper_pct,
  ci_column = 6, # Points to the blank string column
  ref_line = 0,
  sizes = 0.6,
  x_trans = "none",
  xlim = c(-12.5, 12.5),
  ticks_at = c(-10, -5, 0, 5, 10),
  xlab = "Risk difference (%)",
  arrow_lab = NULL, # More precise clinical language
  title = "B. Comparative risk of treatment discontinuation among dose adjustment groups\n    after first interruption of 5-14 days\n",
  nudge_y = 0.02,
  theme = forest_theme(
    base_size = 14,
    ci_pch = 15,
    ci_lwd = 1.7,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
    arrow_length = 0.05,
    arrow_gp = gpar(fontsize = 11, lwd = 0.7),
    xlab_adjust = c("refline", "center")
  )
)



txt <- "Risk. cumulative incidence<br>CI. compatibility intervals"

p_dis_tr3_1 <- add_grob(p_dis_tr3, 
                        row = nrow(tr3_plot)+2, 
                        col = 1:2,
                        order = "background",
                        gb_fn = gridtext::richtext_grob,
                        text = txt,
                        gp = gpar(fontsize = 11),
                        hjust = 0, vjust = 1, halign = 0, valign = 1,
                        x = unit(0, "npc"), y = unit(1, "npc"))


jpeg(paste0(RES_DIR, "Figures/", "Figure3_B_12m.jpg"),
     width = 3500,  height = 1400, res = 300)
p_dis_tr3_1
dev.off()



