# CURVES PLOT ############

library(tidyverse)
library(patchwork)
library(cowplot)

## 1. SETUP & DATA LOADING ------------------------------------------------------
RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/"

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
    labs(x = NULL, y = NULL, title = "Number of episodes at risk during follow-up") +
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
plot_adj_curve <- function(df, y_max, show_30d_line = TRUE, y_lab="Adjusted Cumulative Incidence (%)") {
  p <- ggplot(df, aes(x = horizon, y = risk, color = trt_label)) +
    
    # Smooth continuous line instead of geom_step
    geom_line(linewidth = 1) + 
    
    scale_x_continuous(limits = c(0, 365), breaks = c(0, 30, 90, 180, 270, 365)) +
    scale_color_manual(values = trt_colors_km) +
    labs(title = NULL, x = "Days since time zero\n", y = y_lab) +
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
p_od_tr2  <- plot_adj_curve(df_od_tr2, y_max = 0.8, show_30d_line = TRUE,
                            y_lab = "Adjusted 1-month risk (%)")
t_od_tr2  <- build_at_risk_table(dat_tr2_all, "time_od_acute")

p_od_tr3  <- plot_adj_curve(df_od_tr3, y_max = 0.8, show_30d_line = TRUE,
                            y_lab = "Adjusted 1-month risk (%)")
t_od_tr3  <- build_at_risk_table(dat_tr3_all, "time_od_acute")

dummy_leg <- get_legend(p_od_tr2 + theme(legend.position = "bottom", legend.title = element_blank(), 
                                         legend.text = element_text(size = 13), legend.key.width = unit(2, "cm")))

# Assemble OD
od_composite_tr2 <- (p_od_tr2 / t_od_tr2) / dummy_leg +
  plot_layout(heights = c(5, 1.2, 0.5)) # Increased table height to 1.2

jpeg(paste0(RES_DIR, "Figures/June09/Adjusted_Risk_OD_Line_TR2.jpg"), width = 3200, height = 2000, res = 300)
print(od_composite_tr2)
dev.off()


od_composite_tr3 <- (p_od_tr3 / t_od_tr3) / dummy_leg +
  plot_layout(heights = c(5, 1.2, 0.5)) # Increased table height to 1.2


jpeg(paste0(RES_DIR, "Figures/June09/Adjusted_Risk_OD_Line_TR3.jpg"), width = 3200, height = 2000, res = 300)
print(od_composite_tr3)
dev.off()

# Generate the updated Discontinuation plots
p_dis_tr2  <- plot_adj_curve(df_dis_tr2, y_max = 0.8, show_30d_line = FALSE,
                             y_lab = "Adjusted 12-month risk (%)")
t_dis_tr2  <- build_at_risk_table(dat_tr2_all, "time_discont")

p_dis_tr3  <- plot_adj_curve(df_dis_tr3, y_max = 0.8, show_30d_line = FALSE,
                             y_lab = "Adjusted 12-month risk (%)")
t_dis_tr3  <- build_at_risk_table(dat_tr3_all, "time_discont")

# Assemble Discont
dis_composite_tr2 <- (p_dis_tr2 / t_dis_tr2) / dummy_leg +
  plot_layout(heights = c(5, 1.2, 0.5))

jpeg(paste0(RES_DIR, "Figures/June09/Adjusted_Risk_Discont_Line_TR2.jpg"), width = 3200, height = 2000, res = 300)
print(dis_composite_tr2)
dev.off()


# Assemble Discont
dis_composite_tr3 <- (p_dis_tr3 / t_dis_tr3) / dummy_leg +
  plot_layout(heights = c(5, 1.2, 0.5))

jpeg(paste0(RES_DIR, "Figures/June09/Adjusted_Risk_Discont_Line_TR3.jpg"), width = 3200, height = 2000, res = 300)
print(dis_composite_tr3)
dev.off()



# CURVES PLOT: LE30D vs GT30D #########

#  ADJUSTED RISK CURVES: LE30D vs GT30D
#  8 separate figures: {OD, Discont} × {TR2, TR3} × {LE30D, GT30D}
library(tidyverse)
library(patchwork)
library(cowplot)

## 1. SETUP & DATA LOADING ----------------------------------------------------
RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/"
FIG_DIR <- paste0(RES_DIR, "LE30D/Figures/")

# Risk objects (adjusted curves) for both strata
risk_tr2_le30d <- readRDS(paste0(RES_DIR, "LE30D/Risk/risk_tr2.rds"))
risk_tr2_gt30d <- readRDS(paste0(RES_DIR, "GT30D/Risk/risk_tr2.rds"))
risk_tr3_le30d <- readRDS(paste0(RES_DIR, "LE30D/Risk/risk_tr3.rds"))
risk_tr3_gt30d <- readRDS(paste0(RES_DIR, "GT30D/Risk/risk_tr3.rds"))

# PS-weighted data objects (for at-risk tables)
dat_tr2_le30d <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr2_34d_all.rds"))
dat_tr2_gt30d <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr2_34d_all.rds"))
dat_tr3_le30d <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr3_all.rds"))
dat_tr3_gt30d <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr3_all.rds"))

trt_colors_km <- c(
  "No dose change"                 = "#000000",
  "<50% reduction"                 = "#2CA02C",
  "\u226550% reduction/\u226440mg" = "#FF7F0E",
  "Dose increase"                  = "#D62728"
)

## 2. DATA PREP FUNCTIONS ------------------------------------------------------
prep_adj_curve <- function(risk_obj) {
  risk_obj$risks %>%
    mutate(
      trt_label = recode(trt,
                         "no reversion"     = "No dose change",
                         "1-<50%"           = "<50% reduction",
                         ">=50%/<=40mg"     = "\u226550% reduction/\u226440mg",
                         "increase"         = "Dose increase"),
      trt_label = factor(trt_label, levels = names(trt_colors_km))
    ) %>%
    bind_rows(
      expand.grid(
        horizon = 0, risk = 0,
        trt_label = factor(names(trt_colors_km), levels = names(trt_colors_km))
      )
    ) %>%
    arrange(trt_label, horizon)
}

build_at_risk_table <- function(dat, time_var) {
  dat_clean <- dat$data_weights %>%
    mutate(trt_label = recode(trt,
                              "no reversion"     = "No dose change",
                              "1-<50%"           = "<50% reduction",
                              ">=50%/<=40mg"     = "\u226550% reduction/\u226440mg",
                              "increase"         = "Dose increase"))
  
  breaks_x <- c(0, 30, 90, 180, 270, 365)
  
  at_risk_df <- purrr::map_dfr(breaks_x, function(t) {
    counts <- dat_clean %>%
      group_by(trt_label) %>%
      summarise(time = t,
                raw_count = sum(!!sym(time_var) >= t, na.rm = TRUE),
                .groups = "drop")
    total <- tibble(
      trt_label = "Total episodes",
      time = t,
      raw_count = sum(dat_clean[[time_var]] >= t, na.rm = TRUE)
    )
    bind_rows(total, counts)
  }) %>%
    mutate(
      trt_label = factor(trt_label,
                         levels = rev(c("Total episodes", names(trt_colors_km)))),
      label_text = as.character(raw_count)
    )
  
  y_face <- ifelse(levels(at_risk_df$trt_label) == "Total episodes", "bold", "plain")
  
  ggplot(at_risk_df, aes(x = time, y = trt_label, label = label_text)) +
    geom_text(aes(fontface = ifelse(trt_label == "Total episodes", "bold", "plain")),
              size = 4.5, color = "black") +
    scale_x_continuous(limits = c(0, 365), breaks = breaks_x) +
    labs(x = NULL, y = NULL, title = "Number of episodes at risk during follow-up") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(), axis.line = element_blank(),
      axis.ticks   = element_blank(), axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(color = "black", size = 12, face = y_face),
      plot.title   = element_text(size = 12, face = "bold", hjust = -0.05),
      plot.margin  = ggplot2::margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")
    )
}

plot_adj_curve <- function(df, y_max, show_30d_line = TRUE,
                           y_lab = "Adjusted Cumulative Incidence (%)") {
  p <- ggplot(df, aes(x = horizon, y = risk, color = trt_label)) +
    geom_line(linewidth = 1) +
    scale_x_continuous(limits = c(0, 365),
                       breaks = c(0, 30, 90, 180, 270, 365)) +
    scale_color_manual(values = trt_colors_km) +
    labs(title = NULL, x = "Days since time zero\n", y = y_lab) +
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
    p <- p + geom_vline(xintercept = 30, linetype = "dashed",
                        color = "gray60", linewidth = 0.8)
  }
  
  p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                limits = c(0, 0.8))
  return(p)
}

## 3. COMPOSITE + SAVE HELPER --------------------------------------------------
save_curve_figure <- function(p_curve, p_table, shared_legend, out_path,
                              width = 3200, height = 2000, res = 300) {
  composite <- (p_curve / p_table) / shared_legend +
    plot_layout(heights = c(5, 1.2, 0.5))
  jpeg(out_path, width = width, height = height, res = res)
  print(composite)
  dev.off()
}

## 4. BUILD DATA + CURVES + TABLES FOR ALL 8 FIGURES --------------------------
# ---- OVERDOSE ----
df_od_tr2_le  <- prep_adj_curve(risk_tr2_le30d$od$curve)
df_od_tr2_gt  <- prep_adj_curve(risk_tr2_gt30d$od$curve)
df_od_tr3_le  <- prep_adj_curve(risk_tr3_le30d$od$curve)
df_od_tr3_gt  <- prep_adj_curve(risk_tr3_gt30d$od$curve)

p_od_tr2_le   <- plot_adj_curve(df_od_tr2_le,  y_max = 0.05, show_30d_line = TRUE, y_lab = "Adjusted 1-month risk (%)")
p_od_tr2_gt   <- plot_adj_curve(df_od_tr2_gt,  y_max = 0.05, show_30d_line = TRUE, y_lab = "Adjusted 1-month risk (%)")
p_od_tr3_le   <- plot_adj_curve(df_od_tr3_le,  y_max = 0.05, show_30d_line = TRUE, y_lab = "Adjusted 1-month risk (%)")
p_od_tr3_gt   <- plot_adj_curve(df_od_tr3_gt,  y_max = 0.05, show_30d_line = TRUE, y_lab = "Adjusted 1-month risk (%)")

t_od_tr2_le   <- build_at_risk_table(dat_tr2_le30d, "time_od_acute")
t_od_tr2_gt   <- build_at_risk_table(dat_tr2_gt30d, "time_od_acute")
t_od_tr3_le   <- build_at_risk_table(dat_tr3_le30d, "time_od_acute")
t_od_tr3_gt   <- build_at_risk_table(dat_tr3_gt30d, "time_od_acute")

# ---- DISCONTINUATION ----
df_dis_tr2_le <- prep_adj_curve(risk_tr2_le30d$dis$curve)
df_dis_tr2_gt <- prep_adj_curve(risk_tr2_gt30d$dis$curve)
df_dis_tr3_le <- prep_adj_curve(risk_tr3_le30d$dis$curve)
df_dis_tr3_gt <- prep_adj_curve(risk_tr3_gt30d$dis$curve)

p_dis_tr2_le  <- plot_adj_curve(df_dis_tr2_le, y_max = 0.9, show_30d_line = FALSE, y_lab = "Adjusted 12-month risk (%)")
p_dis_tr2_gt  <- plot_adj_curve(df_dis_tr2_gt, y_max = 0.9, show_30d_line = FALSE, y_lab = "Adjusted 12-month risk (%)")
p_dis_tr3_le  <- plot_adj_curve(df_dis_tr3_le, y_max = 0.9, show_30d_line = FALSE, y_lab = "Adjusted 12-month risk (%)")
p_dis_tr3_gt  <- plot_adj_curve(df_dis_tr3_gt, y_max = 0.9, show_30d_line = FALSE, y_lab = "Adjusted 12-month risk (%)")

t_dis_tr2_le  <- build_at_risk_table(dat_tr2_le30d, "time_discont")
t_dis_tr2_gt  <- build_at_risk_table(dat_tr2_gt30d, "time_discont")
t_dis_tr3_le  <- build_at_risk_table(dat_tr3_le30d, "time_discont")
t_dis_tr3_gt  <- build_at_risk_table(dat_tr3_gt30d, "time_discont")

# ---- Shared legend (from any curve) ----
shared_legend <- get_legend(
  p_od_tr2_le +
    theme(legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = 13),
          legend.key.width  = unit(2, "cm"))
)

## 5. SAVE 8 SEPARATE JPEGS ----------------------------------------------------
# --- OD ---
save_curve_figure(p_od_tr2_le, t_od_tr2_le, shared_legend,
                  paste0(FIG_DIR, "Adjusted_Risk_OD_Line_TR2_LE30D.jpg"))
save_curve_figure(p_od_tr2_gt, t_od_tr2_gt, shared_legend,
                  paste0(FIG_DIR, "Adjusted_Risk_OD_Line_TR2_GT30D.jpg"))
save_curve_figure(p_od_tr3_le, t_od_tr3_le, shared_legend,
                  paste0(FIG_DIR, "Adjusted_Risk_OD_Line_TR3_LE30D.jpg"))
save_curve_figure(p_od_tr3_gt, t_od_tr3_gt, shared_legend,
                  paste0(FIG_DIR, "Adjusted_Risk_OD_Line_TR3_GT30D.jpg"))

# --- Discontinuation ---
save_curve_figure(p_dis_tr2_le, t_dis_tr2_le, shared_legend,
                  paste0(FIG_DIR, "Adjusted_Risk_Discont_Line_TR2_LE30D.jpg"))
save_curve_figure(p_dis_tr2_gt, t_dis_tr2_gt, shared_legend,
                  paste0(FIG_DIR, "Adjusted_Risk_Discont_Line_TR2_GT30D.jpg"))
save_curve_figure(p_dis_tr3_le, t_dis_tr3_le, shared_legend,
                  paste0(FIG_DIR, "Adjusted_Risk_Discont_Line_TR3_LE30D.jpg"))
save_curve_figure(p_dis_tr3_gt, t_dis_tr3_gt, shared_legend,
                  paste0(FIG_DIR, "Adjusted_Risk_Discont_Line_TR3_GT30D.jpg"))

cat("\n\u2713 All eight adjusted risk-curve figures saved to:\n  ", FIG_DIR, "\n", sep = "")

