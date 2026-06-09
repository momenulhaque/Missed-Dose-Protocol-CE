# overdose #########
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
                         "increase" = "Dose increase")
    )
  
  breaks_x <- c(0, 30, 90, 180, 270, 365)
  
  # A. Calculate cumulative events per treatment (with <10 suppression)
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
  
  # B. Calculate total episodes at risk at time t
  at_risk_df <- purrr::map_dfr(breaks_x, function(t) {
    dat %>%
      summarise(
        time = t,
        # At risk means they have not experienced the event or been censored prior to time t
        raw_count = sum(!!sym(time_var) >= t, na.rm = TRUE),
        .groups = "drop"
      )
  }) %>%
    mutate(
      trt_label = "Total episodes at risk",
      label_text = as.character(raw_count) # We don't apply <10 to the total at risk
    )
  
  # C. Combine and order the factors so "Total at risk" goes to the bottom
  table_levels <- c(names(trt_colors_km), "Total episodes at risk")
  
  table_df <- bind_rows(events_df, at_risk_df) %>%
    mutate(
      trt_label = factor(trt_label, levels = rev(table_levels)) 
    )
  
  # Create text-based ggplot table
  ggplot(table_df, aes(x = time, y = trt_label, label = label_text)) +
    geom_text(size = 4.5, color = "black") +
    scale_x_continuous(limits = c(0, 365), breaks = breaks_x) +
    labs(x = NULL, y = NULL, title = "Risk table") + # Removed y-axis title since row labels are self-explanatory
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_blank(),
      axis.ticks   = element_blank(),
      axis.text.x  = element_blank(), 
      axis.title.y = element_text(color = "black", size = 13, angle=90),
      axis.text.y  = element_text(color = "black", size = 12),
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
# Heights adjusted to 1.2 for the tables to accommodate the new At-Risk row
final_km_composite <- p_km_tr2 / 
  p_tab_tr2 / 
  plot_spacer() / 
  p_km_tr3 / 
  p_tab_tr3 / 
  shared_legend +
  plot_layout(heights = c(5, 1.2, 0.5, 5, 1.2, 0.5)) 

# 8. Save output
jpeg(paste0(RES_DIR, "Figures/June09/eFigure4_Unadjusted_OD_KM_Vertical.jpg"), 
     width = 3200, height = 4000, res = 300)
print(final_km_composite)
dev.off()



# discontinuation ####
library(tidyverse)
library(survival)
library(broom)
library(patchwork)

# 1. Define color palette (Consistent across all figures)
trt_colors_km <- c(
  "No dose change" = "#000000",                   
  "<50% reduction" = "#2CA02C",                   
  "\u226550% reduction/\u226440mg*" = "#FF7F0E",  
  "Dose increase" = "#D62728"                     
)

# 2. Helper function: Extract and format Kaplan-Meier data (Main Plot)
build_km_df_discont <- function(dat, time_var = "time_discont", event_var = "sens_discont") {
  dat <- dat %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion" = "No dose change",
                         "1-<50%" = "<50% reduction",
                         ">=50%/<=40mg" = "\u226550% reduction/\u226440mg*",
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
  
  # A. Calculate cumulative events per treatment (with <10 suppression)
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
  
  # B. Calculate total episodes at risk at time t
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
      label_text = as.character(raw_count) # We don't apply <10 to the total at risk
    )
  
  # C. Combine and order the factors so "Total at risk" goes to the bottom
  table_levels <- c(names(trt_colors_km), "Total episodes at risk")
  
  table_df <- bind_rows(events_df, at_risk_df) %>%
    mutate(
      trt_label = factor(trt_label, levels = rev(table_levels)) 
    )
  
  # Create text-based ggplot table
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
      plot.margin  = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")
    )
}

# 4. Helper function: Plot the main KM curve
plot_km_curve_discont <- function(df, title) {
  ggplot(df, aes(x = time, y = cum_inc, color = trt_label)) +
    geom_step(linewidth = 1) +
    
    # We do NOT use a 30-day dashed line here since the clinical focus is 12-months (365 days)
    
    scale_x_continuous(limits = c(0, 365), breaks = c(0, 30, 90, 180, 270, 365)) +
    
    # Discontinuation risks are much higher than overdose, so we don't cap the y-axis at 15%
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.8)) +
    scale_color_manual(values = trt_colors_km) +
    
    labs(
      title = title, 
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
      legend.position   = "none", # Legend handled by patchwork
      plot.margin       = margin(t = 10, r = 5, b = 0, l = 5, unit = "pt")
    )
}

# 5. Generate data, curves, and tables for both trials
df_km_tr2_dis <- build_km_df_discont(dat_tr2_all$data_weights)
df_km_tr3_dis <- build_km_df_discont(dat_tr3_all$data_weights)

p_km_tr2_dis  <- plot_km_curve_discont(df_km_tr2_dis, "A. Methadone interruption of 3-4 days")
p_tab_tr2_dis <- build_event_table_discont(dat_tr2_all$data_weights)

p_km_tr3_dis  <- plot_km_curve_discont(df_km_tr3_dis, "B. Methadone interruption of 5-14 days")
p_tab_tr3_dis <- build_event_table_discont(dat_tr3_all$data_weights)

# 6. Extract legend from one of the plots to place at the very bottom
dummy_plot <- p_km_tr2_dis + 
  theme(legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size = 13), legend.key.width = unit(2, "cm"))
shared_legend <- cowplot::get_legend(dummy_plot)

# 7. Assemble the vertical layout using Patchwork
final_km_composite_discont <- p_km_tr2_dis / 
  p_tab_tr2_dis / 
  plot_spacer() / 
  p_km_tr3_dis / 
  p_tab_tr3_dis / 
  shared_legend +
  plot_layout(heights = c(5, 1.2, 0.5, 5, 1.2, 0.5)) 

# 8. Save output
jpeg(paste0(RES_DIR, "Figures/June09/eFigure5_Unadjusted_Discont_KM_Vertical.jpg"), 
     width = 3200, height = 4000, res = 300)
print(final_km_composite_discont)
dev.off()