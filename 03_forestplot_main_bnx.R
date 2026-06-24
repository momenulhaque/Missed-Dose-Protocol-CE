library(tidyverse)
library(forestploter)
library(grid)
# FOREST PLOT PREPARATION: Cumulative Incidence & Rates
# Missed-dose trials: BNX interruption of 3-5 days and 6-14 days
# Outcomes: Overdose (1 month) and Discontinuation (12 months)

library(tidyverse)

# LOAD RISK RESULTS



RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/BNX/"

# UNIFIED EXTRACTION: Absolute Risk and Risk Difference at time t



risk_tr4 <- readRDS(paste0(RES_DIR, "Risk/risk_tr4_rmst.rds"))

risk_tr5 <- readRDS(paste0(RES_DIR, "Risk/risk_tr5_rmst.rds"))

dat_tr4_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr4_all.rds"))
dat_tr5_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr5_all.rds"))



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
      rd_disp = if_else(trt_group == ref_name, "Reference", rd_disp),
      rd_pct = if_else(trt_group == ref_name, 0, rd_pct),
      rd_lower_pct = if_else(trt_group == ref_name, 0, rd_lower_pct),
      rd_upper_pct = if_else(trt_group == ref_name, 0, rd_upper_pct),
      Interruption = int_label
    )
  
  return(res)
}


trt_colors <- c(
  "No dose change" = "black", 
  "<50% reduction" = "#2CA02C", 
  "≥50% reduction/≤4mg" = "#FF7F0E", 
  "Dose increase" = "#D62728"
)


# 2. FIGURE 2: OVERDOSE RISK DIFFERENCE (30 Days) #######

risk_at_od <- 30

tr4_od_rd <- extract_rd_data(risk_tr4$all$od$ci, dat_tr4_all$data_weights, risk_at_od, 
                             outcome = "od", int_label = "3-5 days")
tr5_od_rd <- extract_rd_data(risk_tr5$all$od$ci, dat_tr5_all$data_weights, risk_at_od, 
                             outcome = "od", int_label = "6-14 days")


plot_data_4a <- bind_rows(tr4_od_rd, tr5_od_rd) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=4mg" = "≥50% reduction/≤4mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct),
    Interruption = ifelse(duplicated(Interruption), "", Interruption),
    # Assign the new header to the plot canvas column
    `% Difference in 1-month Risk (95% CI)` = paste(rep(" ", 16), collapse = " ") 
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode at risk\nat reinitiation, n` = "Episode",
    `% risk at 1 month\n(95% CI)` = "abs_risk_disp",
    # Clear the header for the numeric column to allow spanning
    ` ` = "rd_disp"
  )




p_4a <- forest(
  plot_data_4a[, c("Interruption", "Dose adjustment\nafter interruption", "Episode at risk\nat reinitiation, n", 
                   "% risk at 1 month\n(95% CI)", " ", "% Difference in 1-month Risk (95% CI)")],
  est = plot_data_4a$rd_pct, lower = plot_data_4a$rd_lower_pct, upper = plot_data_4a$rd_upper_pct,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-2.5, 2.5), ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk difference (%)",
  title = "",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       ci_Theight = 0.2, # <--- ADDS THE CI EDGE BARS
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

# Alignments for Figure 2A
p_4a <- edit_plot(p_4a, col = 1:2, part = "header", hjust = 0, x = unit(0, "npc"))
p_4a <- edit_plot(p_4a, col = 1:2, part = "body", hjust = 0, x = unit(0, "npc"))
p_4a <- edit_plot(p_4a, col = 3:6, part = "header", hjust = 0.5, x = unit(0.5, "npc"))
p_4a <- edit_plot(p_4a, col = 3:5, part = "body", hjust = 0.5, x = unit(0.5, "npc"))
p_4a <- edit_plot(p_4a, col = 6, part = "header", hjust = 0.5, x = unit(0.15, "npc"))

# ---> NEW: Apply Colors and Bigger Point Estimates <---
ci_col_4a <- trt_colors[plot_data_4a$`Dose adjustment\nafter interruption`]

for(i in 1:nrow(plot_data_4a)) {
  if(!is.na(plot_data_4a$rd_pct[i])) { 
    current_color <- ci_col_4a[i]
    p_4a <- edit_plot(p_4a, row = i, col = 6, which = "ci",
                      gp = gpar(col = current_color, fill = current_color, cex = 1.2)) # cex=1.8 makes the box bigger
  }
}

# Add the existing left-aligned footnote
p_4a_final <- add_grob(p_4a, row = nrow(plot_data_4a) + 2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, 
                       text = "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals", 
                       gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, 
                       x = unit(0, "npc"), y = unit(1, "npc"))


jpeg(paste0(RES_DIR, "Figures/", "Figure2_OD_RD.jpg"), width = 5000, height = 1800, res = 400)
p_4a_final
dev.off()


# 3. FIGURE 3: DISCONTINUATION RISK DIFFERENCE (365 Days) ########

risk_at_dis <- 365

tr4_dis_rd <- extract_rd_data(risk_tr4$all$dis$ci, dat_tr4_all$data_weights, risk_at_dis, 
                              outcome = "dis", int_label = "3-5 days")
tr5_dis_rd <- extract_rd_data(risk_tr5$all$dis$ci, dat_tr5_all$data_weights, risk_at_dis, 
                              outcome = "dis", int_label = "6-14 days")

plot_data_4b <- bind_rows(tr4_dis_rd, tr5_dis_rd) %>%
  mutate(
    trt_group = recode(trt_group,
                       "no reversion" = "No dose change",
                       "1-<50%" = "<50% reduction",
                       ">=50%/<=4mg" = "≥50% reduction/≤4mg",
                       "increase" = "Dose increase"),
    rd_pct = ifelse(trt_group == "No dose change", NA, rd_pct),
    rd_lower_pct = ifelse(trt_group == "No dose change", NA, rd_lower_pct),
    rd_upper_pct = ifelse(trt_group == "No dose change", NA, rd_upper_pct),
    Interruption = ifelse(duplicated(Interruption), "", Interruption),
    # Assign the new header to the plot canvas column
    `% Difference in 12-months Risk (95% CI)` = paste(rep(" ", 16), collapse = " ") 
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = "trt_group",
    `Episode at risk\nat reinitiation, n` = "Episode",
    `% risk at 12 months\n(95% CI)` = "abs_risk_disp",
    # Clear the header for the numeric column to allow spanning
    ` ` = "rd_disp"
  )

p_4b <- forest(
  plot_data_4b[, c("Interruption", "Dose adjustment\nafter interruption", "Episode at risk\nat reinitiation, n", 
                   "% risk at 12 months\n(95% CI)", " ", "% Difference in 12-months Risk (95% CI)")],
  est = plot_data_4b$rd_pct, lower = plot_data_4b$rd_lower_pct, upper = plot_data_4b$rd_upper_pct,
  ci_column = 6, ref_line = 0, x_trans = "none", xlim = c(-12.5, 12.5), ticks_at = c(-10, -5, 0, 5, 10),
  xlab = "Risk difference (%)",
  title = "",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       ci_Theight = 0.2, # <--- ADDS THE CI EDGE BARS
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

# Alignments for Figure 2B
p_4b <- edit_plot(p_4b, col = 1:2, part = "header", hjust = 0, x = unit(0, "npc"))
p_4b <- edit_plot(p_4b, col = 1:2, part = "body", hjust = 0, x = unit(0, "npc"))
p_4b <- edit_plot(p_4b, col = 3:6, part = "header", hjust = 0.5, x = unit(0.5, "npc"))
p_4b <- edit_plot(p_4b, col = 3:5, part = "body", hjust = 0.5, x = unit(0.5, "npc"))
p_4b <- edit_plot(p_4b, col = 6, part = "header", hjust = 0.5, x = unit(0.15, "npc"))

# ---> NEW: Apply Colors and Bigger Point Estimates <---
ci_col_4b <- trt_colors[plot_data_4b$`Dose adjustment\nafter interruption`]

for(i in 1:nrow(plot_data_4b)) {
  if(!is.na(plot_data_4b$rd_pct[i])) { 
    current_color <- ci_col_4b[i]
    p_4b <- edit_plot(p_4b, row = i, col = 6, which = "ci",
                      gp = gpar(col = current_color, fill = current_color, cex = 1.2)) # cex=1.8 makes the box bigger
  }
}

# Add the existing left-aligned footnote
p_4b_final <- add_grob(p_4b, row = nrow(plot_data_4b) + 2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, 
                       text = "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals", 
                       gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, 
                       x = unit(0, "npc"), y = unit(1, "npc"))



jpeg(paste0(RES_DIR, "Figures/", "Figure3_Discont_RD.jpg"), width = 5200, height = 1800, res = 400)
p_4b_final
dev.off()