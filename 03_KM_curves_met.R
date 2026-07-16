library(tidyverse)
library(survival)
library(broom)
library(patchwork)

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

# 1. Define color palette
trt_colors_km <- c(
  "No dose change" = "#000000",                   
  "<50% reduction" = "#2CA02C",                   
  "\u226550% reduction/\u226440mg" = "#FF7F0E",  
  "Dose increase" = "#D62728"                     
)

# OVERDOSE #########

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

build_event_table <- function(dat, time_var = "time_od_acute", event_var = "sens_od_acute") {
  dat <- dat %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "\u226550% reduction/\u226440mg",
                         "increase" = "Dose increase")
    )
  
  breaks_x <- c(0, 30, 90, 180, 270, 365)
  
  events_df <- purrr::map_dfr(breaks_x, function(t) {
    dat %>%
      group_by(trt_label) %>%
      summarise(
        time = t,
        raw_count = ifelse(t == 0, 0, sum(!!sym(time_var) <= t & !!sym(event_var) == 1, na.rm = TRUE)),
        .groups = "drop"
      )
  }) %>% mutate(label_text = case_when(time == 0 ~ "0", raw_count > 0 & raw_count < 10 ~ "<10", TRUE ~ as.character(raw_count)))
  
  at_risk_df <- purrr::map_dfr(breaks_x, function(t) {
    dat %>% summarise(time = t, raw_count = sum(!!sym(time_var) >= t, na.rm = TRUE), .groups = "drop")
  }) %>% mutate(trt_label = "Total episodes at risk", label_text = as.character(raw_count))
  
  table_levels <- c(names(trt_colors_km), "Total episodes at risk")
  table_df <- bind_rows(events_df, at_risk_df) %>% mutate(trt_label = factor(trt_label, levels = rev(table_levels)))
  
  ggplot(table_df, aes(x = time, y = trt_label, label = label_text)) +
    geom_text(size = 4.5, color = "black") +
    scale_x_continuous(limits = c(0, 365), breaks = breaks_x) +
    labs(x = NULL, y = NULL, title = "Risk table") + 
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x  = element_blank(), axis.title.y = element_text(color = "black", size = 13, angle=90),
      axis.text.y  = element_text(color = "black", size = 12),
      plot.margin  = ggplot2::margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")
    )
}

plot_km_curve <- function(df) {
  ggplot(df, aes(x = time, y = cum_inc, color = trt_label)) +
    geom_step(linewidth = 1) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "gray60", linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 365), breaks = c(0, 30, 90, 180, 270, 365)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.8)) +
    scale_color_manual(values = trt_colors_km) + 
    labs(title = NULL, x = "Days since time zero\n", y = "Unadjusted Cumulative Incidence (%)") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor  = element_blank(), axis.line = element_line(color = "black"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8), axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(color = "black", size = 13), legend.position = "none", 
      plot.margin       = ggplot2::margin(t = 10, r = 5, b = 0, l = 5, unit = "pt") 
    )
}

dat_tr2_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
dat_tr3_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

df_km_tr2 <- build_km_df(dat_tr2_all$data_weights)
df_km_tr3 <- build_km_df(dat_tr3_all$data_weights)

# Only pass the DF, no titles
p_km_tr2  <- plot_km_curve(df_km_tr2)
p_tab_tr2 <- build_event_table(dat_tr2_all$data_weights)

p_km_tr3  <- plot_km_curve(df_km_tr3)
p_tab_tr3 <- build_event_table(dat_tr3_all$data_weights)

dummy_plot <- p_km_tr2 + theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 13), legend.key.width = unit(2, "cm"))
shared_legend <- cowplot::get_legend(dummy_plot)

# SAVE TRIAL 2 OVERDOSE
tr2_km_od_composite <- p_km_tr2 / p_tab_tr2 / shared_legend + plot_layout(heights = c(5, 1.2, 0.5)) 
jpeg(paste0(RES_DIR, "Figures/June09/eFigure4_Unadjusted_OD_KM_TR2.jpg"), width = 3200, height = 2200, res = 300)
print(tr2_km_od_composite)
dev.off()

# SAVE TRIAL 3 OVERDOSE
tr3_km_od_composite <- p_km_tr3 / p_tab_tr3 / shared_legend + plot_layout(heights = c(5, 1.2, 0.5)) 
jpeg(paste0(RES_DIR, "Figures/June09/eFigure4_Unadjusted_OD_KM_TR3.jpg"), width = 3200, height = 2200, res = 300)
print(tr3_km_od_composite)
dev.off()

# DISCONT ##########
library(tidyverse)
library(survival)
library(broom)
library(patchwork)

# 1. Define color palette (Asterisk included to perfectly match the recode step)
trt_colors_km_discont <- c(
  "No dose change" = "#000000",                   
  "<50% reduction" = "#2CA02C",                   
  "\u226550% reduction/\u226440mg*" = "#FF7F0E",  
  "Dose increase" = "#D62728"                     
)

# 2. Helper function: Extract and format Kaplan-Meier data 
build_km_df_discont <- function(dat, time_var = "time_discont", event_var = "sens_discont") {
  dat <- dat %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "\u226550% reduction/\u226440mg*",
                         "increase" = "Dose increase"),
      trt_label = factor(trt_label, levels = names(trt_colors_km_discont))
    )
  
  form <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ trt_label"))
  fit <- survfit(form, data = dat)
  
  fit_df <- tidy(fit) %>%
    mutate(
      cum_inc = 1 - estimate,
      trt_label = gsub("trt_label=", "", strata), 
      trt_label = factor(trt_label, levels = names(trt_colors_km_discont))
    ) %>% select(time, cum_inc, trt_label)
  
  t0 <- expand.grid(time = 0, cum_inc = 0, trt_label = names(trt_colors_km_discont))
  
  fit_df <- bind_rows(t0, fit_df) %>% arrange(trt_label, time)
  return(fit_df)
}

# 3. Helper function: Generate Cumulative Event & At-Risk Table
build_event_table_discont <- function(dat, time_var = "time_discont", event_var = "sens_discont") {
  dat <- dat %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "\u226550% reduction/\u226440mg*",
                         "increase" = "Dose increase")
    )
  
  breaks_x <- c(0, 30, 90, 180, 270, 365)
  
  events_df <- purrr::map_dfr(breaks_x, function(t) {
    dat %>%
      group_by(trt_label) %>%
      summarise(
        time = t,
        raw_count = ifelse(t == 0, 0, sum(!!sym(time_var) <= t & !!sym(event_var) == 1, na.rm = TRUE)),
        .groups = "drop"
      )
  }) %>%
    mutate(
      label_text = case_when(
        time == 0 ~ "0",
        raw_count > 0 & raw_count < 10 ~ "<10",
        TRUE ~ as.character(raw_count)
      )
    )
  
  at_risk_df <- purrr::map_dfr(breaks_x, function(t) {
    dat %>%
      summarise(
        time = t,
        raw_count = sum(!!sym(time_var) >= t, na.rm = TRUE),
        .groups = "drop"
      )
  }) %>%
    mutate(
      trt_label = "Total episodes at risk",
      label_text = as.character(raw_count) 
    )
  
  table_levels <- c(names(trt_colors_km_discont), "Total episodes at risk")
  
  table_df <- bind_rows(events_df, at_risk_df) %>%
    mutate(
      trt_label = factor(trt_label, levels = rev(table_levels)) 
    )
  
  ggplot(table_df, aes(x = time, y = trt_label, label = label_text)) +
    geom_text(size = 4.5, color = "black") +
    scale_x_continuous(limits = c(0, 365), breaks = breaks_x) +
    labs(x = NULL, y = NULL, title = "Risk of discontinuations table") + 
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_blank(),
      axis.ticks   = element_blank(),
      axis.text.x  = element_blank(), 
      axis.title.y = element_text(color = "black", size = 13, angle = 90),
      axis.text.y  = element_text(color = "black", size = 12),
      plot.title   = element_text(size = 12, face = "bold", hjust = -0.15),
      plot.margin  = ggplot2::margin(t = 0, r = 5, b = 5, l = 5, unit = "pt") 
    )
}

# 4. Helper function: Plot the main KM curve
plot_km_curve_discont <- function(df) {
  ggplot(df, aes(x = time, y = cum_inc, color = trt_label)) +
    geom_step(linewidth = 1) +
    scale_x_continuous(limits = c(0, 365), breaks = c(0, 30, 90, 180, 270, 365)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.8)) +
    scale_color_manual(values = trt_colors_km_discont) + 
    labs(
      title = NULL, 
      x = "Days since time zero\n", 
      y = "Unadjusted Cumulative Incidence (%)"
    ) +
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
}

# 5. Generate data, curves, and tables for both trials
df_km_tr2_dis <- build_km_df_discont(dat_tr2_all$data_weights)
df_km_tr3_dis <- build_km_df_discont(dat_tr3_all$data_weights)

p_km_tr2_dis  <- plot_km_curve_discont(df_km_tr2_dis)
p_tab_tr2_dis <- build_event_table_discont(dat_tr2_all$data_weights)

p_km_tr3_dis  <- plot_km_curve_discont(df_km_tr3_dis)
p_tab_tr3_dis <- build_event_table_discont(dat_tr3_all$data_weights)

# 6. Extract legend
dummy_plot <- p_km_tr2_dis + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size = 13), legend.key.width = unit(2, "cm"))
shared_legend <- cowplot::get_legend(dummy_plot)

# 7A. Assemble and Save Trial 2 (3-4 Days)
tr2_km_composite_discont <- p_km_tr2_dis / 
  p_tab_tr2_dis / 
  shared_legend +
  plot_layout(heights = c(5, 1.2, 0.5)) 

jpeg(paste0(RES_DIR, "Figures/June09/eFigure5_Unadjusted_Discont_KM_TR2.jpg"), 
     width = 3200, height = 2200, res = 300)
print(tr2_km_composite_discont)
dev.off()

# 7B. Assemble and Save Trial 3 (5-14 Days)
tr3_km_composite_discont <- p_km_tr3_dis / 
  p_tab_tr3_dis / 
  shared_legend +
  plot_layout(heights = c(5, 1.2, 0.5)) 

jpeg(paste0(RES_DIR, "Figures/June09/eFigure5_Unadjusted_Discont_KM_TR3.jpg"), 
     width = 3200, height = 2200, res = 300)
print(tr3_km_composite_discont)
dev.off()
















library(tidyverse)
library(survival)
library(broom)
library(patchwork)

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

# 1. Define color palette
trt_colors_km <- c(
  "No dose change" = "#000000",                   
  "<50% reduction" = "#2CA02C",                   
  "\u226550% reduction/\u226440mg" = "#FF7F0E",  
  "Dose increase" = "#D62728"                     
)

# OVERDOSE: LE30D vs GT30D #########

#  UNADJUSTED OVERDOSE KM CURVES
#  4 separate figures: TR2 × {LE30D, GT30D} and TR3 × {LE30D, GT30D}
library(tidyverse)
library(survival)
library(broom)
library(patchwork)
library(cowplot)

# Color mapping
trt_colors_km <- c(
  "No dose change"                 = "black",
  "<50% reduction"                 = "#2CA02C",
  "\u226550% reduction/\u226440mg" = "#FF7F0E",
  "Dose increase"                  = "#D62728"
)

# ---- HELPER: build KM cumulative-incidence data 
build_km_df <- function(dat, time_var = "time_od_acute", event_var = "sens_od_acute") {
  dat <- dat %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion"      = "No dose change",
                         "1-<50%"            = "<50% reduction",
                         ">=50%/<=40mg"      = "\u226550% reduction/\u226440mg",
                         "increase"          = "Dose increase"),
      trt_label = factor(trt_label, levels = names(trt_colors_km))
    )
  
  form <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ trt_label"))
  fit  <- survfit(form, data = dat)
  
  fit_df <- tidy(fit) %>%
    mutate(
      cum_inc   = 1 - estimate,
      trt_label = gsub("trt_label=", "", strata),
      trt_label = factor(trt_label, levels = names(trt_colors_km))
    ) %>%
    select(time, cum_inc, trt_label)
  
  t0 <- expand.grid(time = 0, cum_inc = 0, trt_label = names(trt_colors_km))
  bind_rows(t0, fit_df) %>% arrange(trt_label, time)
}

# ---- HELPER: event / at-risk table 
build_event_table <- function(dat, time_var = "time_od_acute",
                              event_var = "sens_od_acute") {
  dat <- dat %>%
    mutate(trt_label = recode(trt,
                              "no reversion"      = "No dose change",
                              "1-<50%"            = "<50% reduction",
                              ">=50%/<=40mg"      = "\u226550% reduction/\u226440mg",
                              "increase"          = "Dose increase"))
  
  breaks_x <- c(0, 30, 90, 180, 270, 365)
  
  events_df <- purrr::map_dfr(breaks_x, function(t) {
    dat %>%
      group_by(trt_label) %>%
      summarise(
        time = t,
        raw_count = ifelse(t == 0, 0,
                           sum(!!sym(time_var) <= t & !!sym(event_var) == 1, na.rm = TRUE)),
        .groups = "drop"
      )
  }) %>%
    mutate(label_text = case_when(
      time == 0                        ~ "0",
      raw_count > 0 & raw_count < 10   ~ "<10",
      TRUE                             ~ as.character(raw_count)
    ))
  
  at_risk_df <- purrr::map_dfr(breaks_x, function(t) {
    dat %>% summarise(time = t,
                      raw_count = sum(!!sym(time_var) >= t, na.rm = TRUE),
                      .groups = "drop")
  }) %>%
    mutate(trt_label  = "Total episodes at risk",
           label_text = as.character(raw_count))
  
  table_levels <- c(names(trt_colors_km), "Total episodes at risk")
  table_df <- bind_rows(events_df, at_risk_df) %>%
    mutate(trt_label = factor(trt_label, levels = rev(table_levels)))
  
  ggplot(table_df, aes(x = time, y = trt_label, label = label_text)) +
    geom_text(size = 4.5, color = "black") +
    scale_x_continuous(limits = c(0, 365), breaks = breaks_x) +
    labs(x = NULL, y = NULL, title = "Risk table") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(), axis.line = element_blank(),
      axis.ticks   = element_blank(), axis.text.x = element_blank(),
      axis.title.y = element_text(color = "black", size = 13, angle = 90),
      axis.text.y  = element_text(color = "black", size = 12),
      plot.margin  = ggplot2::margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")
    )
}

# ---- HELPER: KM curve 
plot_km_curve <- function(df) {
  ggplot(df, aes(x = time, y = cum_inc, color = trt_label)) +
    geom_step(linewidth = 1) +
    geom_vline(xintercept = 30, linetype = "dashed",
               color = "gray60", linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 365),
                       breaks = c(0, 30, 90, 180, 270, 365)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 0.8)) +
    scale_color_manual(values = trt_colors_km) +
    labs(title = NULL,
         x = "Days since time zero\n",
         y = "Unadjusted Cumulative Incidence (%)") +
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
}

#  1. LOAD PS OBJECTS FOR BOTH STRATA
RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/"
FIG_DIR <- paste0(RES_DIR, "LE30D/Figures/")

ps_tr2_le30d <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr2_34d_all.rds"))
ps_tr2_gt30d <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr2_34d_all.rds"))

ps_tr3_le30d <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr3_all.rds"))
ps_tr3_gt30d <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr3_all.rds"))

#  2. BUILD KM DATA + TABLES + CURVES PER TRIAL × STRATUM
# ---- TR2 × LE30D 
df_km_tr2_le  <- build_km_df(ps_tr2_le30d$data_weights)
p_km_tr2_le   <- plot_km_curve(df_km_tr2_le)
p_tab_tr2_le  <- build_event_table(ps_tr2_le30d$data_weights)

# ---- TR2 × GT30D 
df_km_tr2_gt  <- build_km_df(ps_tr2_gt30d$data_weights)
p_km_tr2_gt   <- plot_km_curve(df_km_tr2_gt)
p_tab_tr2_gt  <- build_event_table(ps_tr2_gt30d$data_weights)

# ---- TR3 × LE30D 
df_km_tr3_le  <- build_km_df(ps_tr3_le30d$data_weights)
p_km_tr3_le   <- plot_km_curve(df_km_tr3_le)
p_tab_tr3_le  <- build_event_table(ps_tr3_le30d$data_weights)

# ---- TR3 × GT30D 
df_km_tr3_gt  <- build_km_df(ps_tr3_gt30d$data_weights)
p_km_tr3_gt   <- plot_km_curve(df_km_tr3_gt)
p_tab_tr3_gt  <- build_event_table(ps_tr3_gt30d$data_weights)

# ---- Shared legend 
dummy_plot <- p_km_tr2_le +
  theme(
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 13),
    legend.key.width  = unit(2, "cm")
  )
shared_legend <- cowplot::get_legend(dummy_plot)

#  3. SAVE 4 SEPARATE JPEGS
save_km_figure <- function(p_km, p_tab, shared_legend, out_path,
                           width = 3200, height = 2200, res = 300) {
  composite <- p_km / p_tab / shared_legend +
    plot_layout(heights = c(5, 1.2, 0.5))
  jpeg(out_path, width = width, height = height, res = res)
  print(composite)
  dev.off()
}

# ---- TR2 LE30D 
save_km_figure(p_km_tr2_le,  p_tab_tr2_le,  shared_legend,
               paste0(FIG_DIR, "eFigure4_Unadjusted_OD_KM_TR2_LE30D.jpg"))

# ---- TR2 GT30D 
save_km_figure(p_km_tr2_gt,  p_tab_tr2_gt,  shared_legend,
               paste0(FIG_DIR, "eFigure4_Unadjusted_OD_KM_TR2_GT30D.jpg"))

# ---- TR3 LE30D 
save_km_figure(p_km_tr3_le,  p_tab_tr3_le,  shared_legend,
               paste0(FIG_DIR, "eFigure4_Unadjusted_OD_KM_TR3_LE30D.jpg"))

# ---- TR3 GT30D 
save_km_figure(p_km_tr3_gt,  p_tab_tr3_gt,  shared_legend,
               paste0(FIG_DIR, "eFigure4_Unadjusted_OD_KM_TR3_GT30D.jpg"))

cat("\n\u2713 All four KM figures saved to:\n  ", FIG_DIR, "\n", sep = "")
