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
  title = "",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

p_2a_final <- add_grob(p_2a, row = nrow(plot_data_2a) + 2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals", 
                       gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, 
                       x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/June09/", "Figure2_A_OD_RD.jpg"), width = 3500, height = 1300, res = 300)
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
  title = "",
  theme = forest_theme(base_size = 14, ci_pch = 15, ci_lwd = 1.7,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       xlab_adjust = c("refline", "center"))
)

p_2b_final <- add_grob(p_2b, row = nrow(plot_data_2b) + 2, col = 1:2, order = "background", 
                       gb_fn = gridtext::richtext_grob, text = "Risk: cumulative incidence<br>CI: compatibility ('confidence') intervals", 
                       gp = gpar(fontsize = 11), hjust = 0, vjust = 1, halign = 0, valign = 1, 
                       x = unit(0, "npc"), y = unit(1, "npc"))

jpeg(paste0(RES_DIR, "Figures/June09/", "Figure2_B_Discont_RD.jpg"), width = 3500, height = 1300, res = 300)
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

# 4. CURVES PLOT ############

library(tidyverse)
library(patchwork)
library(cowplot)

## 1. SETUP & DATA LOADING ------------------------------------------------------
RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

# Load Risk Objects (Adjusted Curves)
risk_tr2 <- readRDS(paste0(RES_DIR, "Risk/risk_tr2_34d_rmst.rds"))
risk_tr3 <- readRDS(paste0(RES_DIR, "Risk/risk_tr3_514d_rmst.rds"))

# Load PS Data (For Raw Risk Tables)
dat_tr2_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
dat_tr3_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

# Define Shared Aesthetics
trt_colors_km <- c(
  "No dose change" = "#000000",                   
  "<50% reduction" = "#2CA02C",                   
  "\u226550% reduction/\u226440mg" = "#FF7F0E",  
  "Dose increase" = "#D62728"                     
)

## 2. DATA PREP FUNCTIONS -------------------------------------------------------

# Helper to extract and format the adjusted risk curve
prep_adj_curve <- function(risk_obj) {
  risk_obj$risks %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "\u226550% reduction/\u226440mg",
                         "increase" = "Dose increase"),
      trt_label = factor(trt_label, levels = names(trt_colors_km))
    ) %>%
    # Add Time 0 with 0% risk for all treatments to anchor the plot
    bind_rows(
      expand.grid(
        horizon = 0, 
        risk = 0, 
        trt_label = factor(names(trt_colors_km), levels = names(trt_colors_km))
      )
    ) %>%
    arrange(trt_label, horizon)
}

# Helper to generate the At-Risk Table (Raw Denominator)
build_at_risk_table <- function(dat, time_var) {
  dat_clean <- dat$data_weights %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "\u226550% reduction/\u226440mg",
                         "increase" = "Dose increase")
    )
  
  breaks_x <- c(0, 30, 90, 180, 270, 365)
  
  at_risk_df <- purrr::map_dfr(breaks_x, function(t) {
    # 1. Calculate counts for each specific treatment group
    counts <- dat_clean %>%
      group_by(trt_label) %>%
      summarise(
        time = t,
        raw_count = sum(!!sym(time_var) >= t, na.rm = TRUE),
        .groups = "drop"
      )
    
    # 2. Calculate the aggregate "Total episodes" row
    total <- tibble(
      trt_label = "Total episodes",
      time = t,
      raw_count = sum(dat_clean[[time_var]] >= t, na.rm = TRUE)
    )
    
    # Bind them together for this time point
    bind_rows(total, counts)
  }) %>%
    mutate(
      # Place "Total episodes" at the top of the table 
      # (which means making it the last factor level because ggplot builds bottom-to-top)
      trt_label = factor(trt_label, levels = rev(c("Total episodes", names(trt_colors_km)))),
      label_text = as.character(raw_count)
    )
  
  # Create a vector to make the "Total episodes" label bold on the y-axis
  y_face <- ifelse(levels(at_risk_df$trt_label) == "Total episodes", "bold", "plain")
  
  ggplot(at_risk_df, aes(x = time, y = trt_label, label = label_text)) +
    # Bold the numbers for the Total row, keep others plain
    geom_text(aes(fontface = ifelse(trt_label == "Total episodes", "bold", "plain")), 
              size = 4.5, color = "black") +
    scale_x_continuous(limits = c(0, 365), breaks = breaks_x) +
    labs(x = NULL, y = NULL, title = "Episodes at risk") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_blank(),
      axis.ticks   = element_blank(),
      axis.text.x  = element_blank(),
      axis.title.y = element_blank(),
      # Apply the dynamic bolding to the y-axis text
      axis.text.y  = element_text(color = "black", size = 12, face = y_face),
      plot.title   = element_text(size = 12, face = "bold", hjust = -0.05),
      plot.margin  = ggplot2::margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")
    )
}

# Helper to plot the main Adjusted Curve (Smooth Line Version)
plot_adj_curve <- function(df, title, y_max, show_30d_line = TRUE) {
  p <- ggplot(df, aes(x = horizon, y = risk, color = trt_label)) +
    
    # Smooth continuous line instead of geom_step
    geom_line(linewidth = 1) + 
    
    scale_x_continuous(limits = c(0, 365), breaks = c(0, 30, 90, 180, 270, 365)) +
    scale_color_manual(values = trt_colors_km) +
    labs(title = title, x = "Days since time zero\n", y = "Adjusted Cumulative Incidence (%)") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(color = "black", size = 13),
      legend.position   = "none",
      plot.margin       = ggplot2::margin(t = 10, r = 5, b = 0, l = 5, unit = "pt") 
    )
  
  if (show_30d_line) {
    p <- p + geom_vline(xintercept = 30, linetype = "dashed", color = "gray60", linewidth = 0.8)
  }
  
  if (y_max <= 0.20) {
    p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                limits = c(0, y_max), breaks = seq(0, y_max, 0.05))
  } else {
    p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                limits = c(0, y_max))
  }
  return(p)
}

## 3. GENERATE FIGURE #####
cat("Generating Adjusted Overdose Curves...\n")

df_od_tr2 <- prep_adj_curve(risk_tr2$all$od$curve)
df_od_tr3 <- prep_adj_curve(risk_tr3$all$od$curve)

# Generate the updated Overdose plots
p_od_tr2  <- plot_adj_curve(df_od_tr2, "A. Methadone interruption of 3-4 days", 
                            y_max = 0.6, show_30d_line = TRUE)
t_od_tr2  <- build_at_risk_table(dat_tr2_all, "time_od_acute")

p_od_tr3  <- plot_adj_curve(df_od_tr3, "B. Methadone interruption of 5-14 days", 
                            y_max = 0.6, show_30d_line = TRUE)
t_od_tr3  <- build_at_risk_table(dat_tr3_all, "time_od_acute")

dummy_leg <- get_legend(p_od_tr2 + theme(legend.position = "bottom", legend.title = element_blank(), 
                                         legend.text = element_text(size = 13), legend.key.width = unit(2, "cm")))

# Assemble OD
od_composite <- (p_od_tr2 / t_od_tr2 / plot_spacer() / p_od_tr3 / t_od_tr3) / dummy_leg +
  plot_layout(heights = c(5, 1.2, 0.5, 5, 1.2, 0.5)) # Increased table height to 1.2

jpeg(paste0(RES_DIR, "Figures/June09/Adjusted_OD_Line_Vertical.jpg"), width = 3200, height = 4000, res = 300)
print(od_composite)
dev.off()


# Generate the updated Discontinuation plots
p_dis_tr2  <- plot_adj_curve(df_dis_tr2, "A. Methadone interruption of 3-4 days", y_max = 0.6, show_30d_line = FALSE)
t_dis_tr2  <- build_at_risk_table(dat_tr2_all, "time_discont")

p_dis_tr3  <- plot_adj_curve(df_dis_tr3, "B. Methadone interruption of 5-14 days", y_max = 0.6, show_30d_line = FALSE)
t_dis_tr3  <- build_at_risk_table(dat_tr3_all, "time_discont")

# Assemble Discont
dis_composite <- (p_dis_tr2 / t_dis_tr2 / plot_spacer() / p_dis_tr3 / t_dis_tr3) / dummy_leg +
  plot_layout(heights = c(5, 1.2, 0.5, 5, 1.2, 0.5))

jpeg(paste0(RES_DIR, "Figures/June09/Adjusted_Discont_Line_Vertical.jpg"), width = 3200, height = 4000, res = 300)
print(dis_composite)
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




## discontinuation ##########

cat("\n--- Trial 2: Discontinuation (Raw) ---\n")
raw_dis_tr2 <- compute_raw_km_stats(dat_tr2_all$data_weights, 
                                    time_var = "time_discont", 
                                    event_var = "sens_discont")
print(raw_dis_tr2, digits = 4, row.names = FALSE)

cat("\n--- Trial 3: Discontinuation (Raw) ---\n")
raw_dis_tr3 <- compute_raw_km_stats(dat_tr3_all$data_weights, 
                                    time_var = "time_discont", 
                                    event_var = "sens_discont")
print(raw_dis_tr3, digits = 4, row.names = FALSE)


