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



















# R2R version: OVERDOSE & DISCONTINUATION ######
## Strata: LE30D vs GT30D ##########

library(tidyverse)
library(survival)
library(broom)
library(patchwork)
library(cowplot)

source("H:/GitHub/Missed-Dose-Protocols-PRI/00_functions.R")

## ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/"
FIG_DIR   <- paste0(RES_DIR, "LE30D/Figures/")


## ── 1. Define Palettes & Line Types ───────────────────────────────────────────
trt_colors_km <- c(
  "No dose change"                 = "black",
  "<50% reduction"                 = "#2CA02C",
  "\u226550% reduction or \u226440mg" = "#FF7F0E",
  "Dose increase"                  = "#D62728"
)

trt_linetypes_km <- c(
  "No dose change"                 = "solid",
  "<50% reduction"                 = "dashed",
  "\u226550% reduction or \u226440mg" = "dotted",
  "Dose increase"                  = "dotdash"
)


## ── 2. Unified Helper Functions ───────────────────────────────────────────────

build_km_df <- function(dat, time_var, event_var) {
  dat <- dat %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion"      = "No dose change",
                         "1-<50%"            = "<50% reduction",
                         ">=50%/<=40mg"      = "\u226550% reduction or \u226440mg",
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

build_event_table <- function(dat, time_var, event_var) {
  dat <- dat %>%
    mutate(trt_label = recode(trt,
                              "no reversion"      = "No dose change",
                              "1-<50%"            = "<50% reduction",
                              ">=50%/<=40mg"      = "\u226550% reduction or \u226440mg",
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
      raw_count >= 0 & raw_count < 10   ~ "<10",
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
    labs(x = NULL, y = NULL, title = "Number of events/episodes during follow-up") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(), axis.line = element_blank(),
      axis.ticks   = element_blank(), axis.text.x = element_blank(),
      axis.title.y = element_text(color = "black", size = 13, angle = 90),
      axis.text.y  = element_text(color = "black", size = 12),
      plot.margin  = ggplot2::margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")
    )
}

plot_km_curve <- function(df, y_max) {
  ggplot(df, aes(x = time, y = cum_inc, color = trt_label, linetype = trt_label)) +
    geom_step(linewidth = 1.5) + # Updated linewidth
    geom_vline(xintercept = 30, linetype = "dashed",
               color = "gray60", linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 365),
                       breaks = c(0, 30, 90, 180, 270, 365)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, y_max)) +
    scale_color_manual(values = trt_colors_km) +
    scale_linetype_manual(values = trt_linetypes_km) +
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


plot_km_curve_dis <- function(df, y_max) {
  ggplot(df, aes(x = time, y = cum_inc, color = trt_label, linetype = trt_label)) +
    geom_step(linewidth = 1.5) + # Updated linewidth
    scale_x_continuous(limits = c(0, 365),
                       breaks = c(0, 30, 90, 180, 270, 365)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, y_max)) +
    scale_color_manual(values = trt_colors_km) +
    scale_linetype_manual(values = trt_linetypes_km) +
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

save_km_figure <- function(p_km, p_tab, shared_legend, out_path,
                           width = 3200, height = 2200, res = 300) {
  composite <- p_km / p_tab / shared_legend +
    plot_layout(heights = c(5, 1.2, 0.5))
  jpeg(out_path, width = width, height = height, res = res)
  print(composite)
  dev.off()
}


## ── 3. Load PS Objects ────────────────────────────────────────────────────────

ps_tr2_le30d <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr2_34d_all.rds"))
ps_tr2_gt30d <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr2_34d_all.rds"))

ps_tr3_le30d <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr3_all.rds"))
ps_tr3_gt30d <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr3_all.rds"))

df_tr2_le <- ps_tr2_le30d$data_weights
df_tr2_gt <- ps_tr2_gt30d$data_weights
df_tr3_le <- ps_tr3_le30d$data_weights
df_tr3_gt <- ps_tr3_gt30d$data_weights


## ── 4. Generate Shared Legend (Fixed Thickness & Bounds) ──────────────────────
dummy_data <- expand.grid(time = 1:2, cum_inc = c(0, 0.1), trt_label = names(trt_colors_km))
dummy_plot <- ggplot(dummy_data, aes(x = time, y = cum_inc, color = trt_label, linetype = trt_label)) +
  geom_step(linewidth = 1.5) +
  scale_color_manual(values = trt_colors_km) +
  scale_linetype_manual(values = trt_linetypes_km) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.5)),
    linetype = guide_legend(override.aes = list(linewidth = 1.5))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 13),
    legend.key.width  = unit(2.2, "cm"), # Adjusted to prevent text/border clipping
    legend.key.height = unit(0.8, "cm")
  )
shared_legend <- cowplot::get_legend(dummy_plot)


## ── 5. Generate & Save OVERDOSE Figures (Y-Max: 0.25) ─────────────────────────

save_km_figure(
  p_km = plot_km_curve(build_km_df(df_tr2_le, "time_od_acute", "sens_od_acute"), 0.25),
  p_tab = build_event_table(df_tr2_le, "time_od_acute", "sens_od_acute"),
  shared_legend = shared_legend,
  out_path = paste0(FIG_DIR, "eFigure5_Unadjusted_OD_KM_TR2_LE30D.jpg")
)

save_km_figure(
  p_km = plot_km_curve(build_km_df(df_tr2_gt, "time_od_acute", "sens_od_acute"), 0.25),
  p_tab = build_event_table(df_tr2_gt, "time_od_acute", "sens_od_acute"),
  shared_legend = shared_legend,
  out_path = paste0(FIG_DIR, "eFigure5_Unadjusted_OD_KM_TR2_GT30D.jpg")
)

save_km_figure(
  p_km = plot_km_curve(build_km_df(df_tr3_le, "time_od_acute", "sens_od_acute"), 0.25),
  p_tab = build_event_table(df_tr3_le, "time_od_acute", "sens_od_acute"),
  shared_legend = shared_legend,
  out_path = paste0(FIG_DIR, "eFigure6_Unadjusted_OD_KM_TR3_LE30D.jpg")
)

save_km_figure(
  p_km = plot_km_curve(build_km_df(df_tr3_gt, "time_od_acute", "sens_od_acute"), 0.25),
  p_tab = build_event_table(df_tr3_gt, "time_od_acute", "sens_od_acute"),
  shared_legend = shared_legend,
  out_path = paste0(FIG_DIR, "eFigure6_Unadjusted_OD_KM_TR3_GT30D.jpg")
)


## ── 6. Generate & Save DISCONTINUATION Figures (Y-Max: 1.0) ───────────────────

save_km_figure(
  p_km = plot_km_curve_dis(build_km_df(df_tr2_le, "time_discont", "sens_discont"), 1.0),
  p_tab = build_event_table(df_tr2_le, "time_discont", "sens_discont"),
  shared_legend = shared_legend,
  out_path = paste0(FIG_DIR, "eFigure5_Unadjusted_DIS_KM_TR2_LE30D.jpg")
)

save_km_figure(
  p_km = plot_km_curve_dis(build_km_df(df_tr2_gt, "time_discont", "sens_discont"), 1.0),
  p_tab = build_event_table(df_tr2_gt, "time_discont", "sens_discont"),
  shared_legend = shared_legend,
  out_path = paste0(FIG_DIR, "eFigure5_Unadjusted_DIS_KM_TR2_GT30D.jpg")
)

save_km_figure(
  p_km = plot_km_curve_dis(build_km_df(df_tr3_le, "time_discont", "sens_discont"), 1.0),
  p_tab = build_event_table(df_tr3_le, "time_discont", "sens_discont"),
  shared_legend = shared_legend,
  out_path = paste0(FIG_DIR, "eFigure6_Unadjusted_DIS_KM_TR3_LE30D.jpg")
)

save_km_figure(
  p_km = plot_km_curve_dis(build_km_df(df_tr3_gt, "time_discont", "sens_discont"), 1.0),
  p_tab = build_event_table(df_tr3_gt, "time_discont", "sens_discont"),
  shared_legend = shared_legend,
  out_path = paste0(FIG_DIR, "eFigure6_Unadjusted_DIS_KM_TR3_GT30D.jpg")
)






## ── 7. Extract Unadjusted Risk at 30 and 365 Days to CSV ─────────────────────

get_unadjusted_risk <- function(dat, time_var, event_var, outcome_name, trt_strata, timing_strata) {
  dat <- dat %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion"      = "No dose change",
                         "1-<50%"            = "<50% reduction",
                         ">=50%/<=40mg"      = "\u226550% reduction or \u226440mg",
                         "increase"          = "Dose increase"),
      trt_label = factor(trt_label, levels = names(trt_colors_km))
    )
  
  form <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ trt_label"))
  fit  <- survfit(form, data = dat)
  
  # Get KM estimates at exactly 30 and 365 days
  sum_fit <- summary(fit, times = c(30, 365))
  
  tibble(
    outcome          = outcome_name,
    interruption     = trt_strata,
    timing           = timing_strata,
    treatment_arm    = gsub("trt_label=", "", sum_fit$strata),
    time_days        = sum_fit$time,
    n_at_risk        = sum_fit$n.risk,
    cum_events       = sum_fit$n.event,
    cum_risk_pct     = round((1 - sum_fit$surv) * 100, 2),
    lower_95_ci_pct  = round((1 - sum_fit$upper) * 100, 2),
    upper_95_ci_pct  = round((1 - sum_fit$lower) * 100, 2),
    risk_with_95ci   = sprintf("%.2f%% (%.2f%% - %.2f%%)", 
                               (1 - sum_fit$surv) * 100, 
                               (1 - sum_fit$upper) * 100, 
                               (1 - sum_fit$lower) * 100)
  )
}

# Define evaluation scenarios
scenarios <- list(
  # Overdose
  list(dat = df_tr2_le, t_var = "time_od_acute", e_var = "sens_od_acute", outcome = "Overdose", trt = "3-4 days (TR2)", timing = "\u226430 days (LE30D)"),
  list(dat = df_tr2_gt, t_var = "time_od_acute", e_var = "sens_od_acute", outcome = "Overdose", trt = "3-4 days (TR2)", timing = ">30 days (GT30D)"),
  list(dat = df_tr3_le, t_var = "time_od_acute", e_var = "sens_od_acute", outcome = "Overdose", trt = "5-14 days (TR3)", timing = "\u226430 days (LE30D)"),
  list(dat = df_tr3_gt, t_var = "time_od_acute", e_var = "sens_od_acute", outcome = "Overdose", trt = "5-14 days (TR3)", timing = ">30 days (GT30D)"),
  
  # Discontinuation
  list(dat = df_tr2_le, t_var = "time_discont",  e_var = "sens_discont",  outcome = "Discontinuation", trt = "3-4 days (TR2)", timing = "\u226430 days (LE30D)"),
  list(dat = df_tr2_gt, t_var = "time_discont",  e_var = "sens_discont",  outcome = "Discontinuation", trt = "3-4 days (TR2)", timing = ">30 days (GT30D)"),
  list(dat = df_tr3_le, t_var = "time_discont",  e_var = "sens_discont",  outcome = "Discontinuation", trt = "5-14 days (TR3)", timing = "\u226430 days (LE30D)"),
  list(dat = df_tr3_gt, t_var = "time_discont",  e_var = "sens_discont",  outcome = "Discontinuation", trt = "5-14 days (TR3)", timing = ">30 days (GT30D)")
)

# Extract and combine into one table
unadjusted_risk_table <- purrr::map_dfr(scenarios, function(s) {
  get_unadjusted_risk(s$dat, s$t_var, s$e_var, s$outcome, s$trt, s$timing)
})

# Save to CSV
write_csv(unadjusted_risk_table, paste0(FIG_DIR, "Unadjusted_Risk_30_365_Days.csv"))

head(unadjusted_risk_table)




unadjusted_risk_table |> filter(outcome == "Discontinuation", 
                                time_days == 365,
                                interruption == "3-4 days (TR2)") |> select(c(1:4, "cum_risk_pct"))



unadjusted_risk_table |> filter(outcome == "Discontinuation", 
                                time_days == 365,
                                interruption == "5-14 days (TR3)") |> select(c(1:4, "cum_risk_pct"))
