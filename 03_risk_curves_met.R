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
### ---- OVERDOSE ----
df_od_tr2_le  <- prep_adj_curve(risk_tr2_le30d$od$curve)
df_od_tr2_gt  <- prep_adj_curve(risk_tr2_gt30d$od$curve)
df_od_tr3_le  <- prep_adj_curve(risk_tr3_le30d$od$curve)
df_od_tr3_gt  <- prep_adj_curve(risk_tr3_gt30d$od$curve)

p_od_tr2_le   <- plot_adj_curve(df_od_tr2_le,  y_max = 0.05, show_30d_line = TRUE, y_lab = "Adjusted risk (%)")
p_od_tr2_gt   <- plot_adj_curve(df_od_tr2_gt,  y_max = 0.05, show_30d_line = TRUE, y_lab = "Adjusted risk (%)")
p_od_tr3_le   <- plot_adj_curve(df_od_tr3_le,  y_max = 0.05, show_30d_line = TRUE, y_lab = "Adjusted risk (%)")
p_od_tr3_gt   <- plot_adj_curve(df_od_tr3_gt,  y_max = 0.05, show_30d_line = TRUE, y_lab = "Adjusted risk (%)")

t_od_tr2_le   <- build_at_risk_table(dat_tr2_le30d, "time_od_acute")
t_od_tr2_gt   <- build_at_risk_table(dat_tr2_gt30d, "time_od_acute")
t_od_tr3_le   <- build_at_risk_table(dat_tr3_le30d, "time_od_acute")
t_od_tr3_gt   <- build_at_risk_table(dat_tr3_gt30d, "time_od_acute")

### ---- DISCONTINUATION ----
df_dis_tr2_le <- prep_adj_curve(risk_tr2_le30d$dis$curve)
df_dis_tr2_gt <- prep_adj_curve(risk_tr2_gt30d$dis$curve)
df_dis_tr3_le <- prep_adj_curve(risk_tr3_le30d$dis$curve)
df_dis_tr3_gt <- prep_adj_curve(risk_tr3_gt30d$dis$curve)


p_dis_tr2_le  <- plot_adj_curve(df_dis_tr2_le, y_max = 0.9, show_30d_line = FALSE, y_lab = "Adjusted risk (%)")
p_dis_tr2_gt  <- plot_adj_curve(df_dis_tr2_gt, y_max = 0.9, show_30d_line = FALSE, y_lab = "Adjusted risk (%)")
p_dis_tr3_le  <- plot_adj_curve(df_dis_tr3_le, y_max = 0.9, show_30d_line = FALSE, y_lab = "Adjusted risk (%)")
p_dis_tr3_gt  <- plot_adj_curve(df_dis_tr3_gt, y_max = 0.9, show_30d_line = FALSE, y_lab = "Adjusted risk (%)")

t_dis_tr2_le  <- build_at_risk_table(dat_tr2_le30d, "time_discont")
t_dis_tr2_gt  <- build_at_risk_table(dat_tr2_gt30d, "time_discont")
t_dis_tr3_le  <- build_at_risk_table(dat_tr3_le30d, "time_discont")
t_dis_tr3_gt  <- build_at_risk_table(dat_tr3_gt30d, "time_discont")

### ---- Shared legend (from any curve) ----
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


# COMBINED ADJUSTED RISK CURVES: LE30D vs GT30D ########
#
# Creates two high-resolution JPEG figures:
#   1. Adjusted_Risk_OD_Combined.jpg
#   2. Adjusted_Risk_Discont_Combined.jpg
#
# Each figure contains four panels:
#   Top-left:     3-4-day interruption, within <=30 days of OAT initiation
#   Top-right:    3-4-day interruption, after >30 days of OAT initiation
#   Bottom-left:  5-14-day interruption, within <=30 days of OAT initiation
#   Bottom-right: 5-14-day interruption, after >30 days of OAT initiation
#
# Each panel includes:
#   - Adjusted cumulative-incidence curves
#   - Number-of-episodes-at-risk table
#
# Important correction:
#   coord_cartesian() is used for y-axis limits so estimates outside the visible
#   range are not removed before geom_line() draws the curve.

options(error = NULL)

library(tidyverse)
library(patchwork)
library(cowplot)
library(grid)


# 1. PATHS AND OUTPUT FILES ########

RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/"

FIG_DIR <- paste0(
  RES_DIR,
  "LE30D/Figures/"
)

dir.create(
  FIG_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

OD_FILE <- paste0(
  FIG_DIR,
  "Adjusted_Risk_OD_Combined.jpg"
)

DIS_FILE <- paste0(
  FIG_DIR,
  "Adjusted_Risk_Discont_Combined.jpg"
)


# 2. LOAD ADJUSTED-RISK OBJECTS ###########


risk_tr2_le30d <- readRDS(
  paste0(
    RES_DIR,
    "LE30D/Risk/risk_tr2.rds"
  )
)

risk_tr2_gt30d <- readRDS(
  paste0(
    RES_DIR,
    "GT30D/Risk/risk_tr2.rds"
  )
)

risk_tr3_le30d <- readRDS(
  paste0(
    RES_DIR,
    "LE30D/Risk/risk_tr3.rds"
  )
)

risk_tr3_gt30d <- readRDS(
  paste0(
    RES_DIR,
    "GT30D/Risk/risk_tr3.rds"
  )
)


# 3. LOAD PS-WEIGHTED DATA OBJECTS FOR RISK TABLES

dat_tr2_le30d <- readRDS(
  paste0(
    RES_DIR,
    "LE30D/PScore/ps_weight_tr2_34d_all.rds"
  )
)

dat_tr2_gt30d <- readRDS(
  paste0(
    RES_DIR,
    "GT30D/PScore/ps_weight_tr2_34d_all.rds"
  )
)

dat_tr3_le30d <- readRDS(
  paste0(
    RES_DIR,
    "LE30D/PScore/ps_weight_tr3_all.rds"
  )
)

dat_tr3_gt30d <- readRDS(
  paste0(
    RES_DIR,
    "GT30D/PScore/ps_weight_tr3_all.rds"
  )
)


# 4. TREATMENT SETTINGS

TRT_LEVELS <- c(
  "no reversion",
  "1-<50%",
  ">=50%/<=40mg",
  "increase"
)

TRT_LABELS <- c(
  "No dose change",
  "<50% reduction",
  "\u226550% reduction/\u226440mg",
  "Dose increase"
)

trt_colors_km <- c(
  "No dose change"                 = "#000000",
  "<50% reduction"                 = "#2CA02C",
  "\u226550% reduction/\u226440mg" = "#FF7F0E",
  "Dose increase"                  = "#D62728"
)

breaks_x <- c(
  0,
  30,
  90,
  180,
  270,
  365
)


# 5. PUBLICATION-QUALITY FIGURE SETTINGS

FIG_FONT_FAMILY <- "sans"

PANEL_TITLE_SIZE <- 20
AXIS_TITLE_SIZE  <- 17
AXIS_TEXT_SIZE   <- 15

TABLE_TITLE_SIZE <- 15
TABLE_LABEL_SIZE <- 14
TABLE_VALUE_SIZE <- 6

LEGEND_TEXT_SIZE <- 15

CURVE_LINE_WIDTH <- 1.40
AXIS_LINE_WIDTH  <- 0.85
GRID_LINE_WIDTH  <- 0.40
REF_LINE_WIDTH   <- 0.90

CURVE_TABLE_SPACER <- 0.22

# Extended horizontal display range creates room for values at days 0 and 365.
X_DISPLAY_MIN <- -18
X_DISPLAY_MAX <- 383

# Outcome-specific visible scales.
OD_Y_MAX      <- 0.25
OD_Y_INTERVAL <- 0.05

DIS_Y_MAX      <- 0.60
DIS_Y_INTERVAL <- 0.10


# 6. PREPARE ADJUSTED CURVE DATA #########

prep_adj_curve <- function(risk_obj) {
  
  if (
    is.null(risk_obj) ||
    is.null(risk_obj$risks)
  ) {
    stop(
      "The supplied risk object does not contain a $risks element."
    )
  }
  
  required_vars <- c(
    "horizon",
    "risk",
    "trt"
  )
  
  missing_vars <- setdiff(
    required_vars,
    names(risk_obj$risks)
  )
  
  if (length(missing_vars) > 0) {
    stop(
      "Missing variable(s) in risk_obj$risks: ",
      paste(
        missing_vars,
        collapse = ", "
      )
    )
  }
  
  curve_df <- risk_obj$risks %>%
    transmute(
      horizon = as.numeric(horizon),
      
      risk = as.numeric(risk),
      
      trt_label = recode(
        as.character(trt),
        
        "no reversion" =
          "No dose change",
        
        "1-<50%" =
          "<50% reduction",
        
        ">=50%/<=40mg" =
          "\u226550% reduction/\u226440mg",
        
        "increase" =
          "Dose increase",
        
        .default = NA_character_
      ),
      
      trt_label = factor(
        trt_label,
        levels = TRT_LABELS
      )
    ) %>%
    filter(
      !is.na(horizon),
      !is.na(trt_label)
    )
  
  # Add a day-zero anchor only for strategies without a day-zero estimate.
  groups_with_day0 <- curve_df %>%
    filter(horizon == 0) %>%
    pull(trt_label) %>%
    as.character() %>%
    unique()
  
  groups_without_day0 <- setdiff(
    TRT_LABELS,
    groups_with_day0
  )
  
  day0_rows <- tibble(
    horizon = 0,
    
    risk = 0,
    
    trt_label = factor(
      groups_without_day0,
      levels = TRT_LABELS
    )
  )
  
  bind_rows(
    curve_df,
    day0_rows
  ) %>%
    distinct(
      trt_label,
      horizon,
      .keep_all = TRUE
    ) %>%
    arrange(
      trt_label,
      horizon
    )
}


# 7. EXTRACT DATA FROM PS OBJECT #########

extract_data_weights <- function(dat_object) {
  
  if (
    is.list(dat_object) &&
    "data_weights" %in% names(dat_object)
  ) {
    return(
      dat_object$data_weights
    )
  }
  
  if (is.data.frame(dat_object)) {
    return(dat_object)
  }
  
  stop(
    "The supplied PS object is neither a data frame nor an object ",
    "containing a data_weights element."
  )
}


# 8. NUMBER-AT-RISK TABLE #######

build_at_risk_table <- function(
    dat,
    time_var,
    show_strategy_names = TRUE,
    show_table_title = TRUE
) {
  
  dat_clean <- extract_data_weights(
    dat
  )
  
  required_vars <- c(
    "trt",
    time_var
  )
  
  missing_vars <- setdiff(
    required_vars,
    names(dat_clean)
  )
  
  if (length(missing_vars) > 0) {
    stop(
      "Missing variable(s) for risk table: ",
      paste(
        missing_vars,
        collapse = ", "
      )
    )
  }
  
  dat_clean <- dat_clean %>%
    filter(
      !is.na(trt),
      !is.na(.data[[time_var]])
    ) %>%
    mutate(
      trt_label = recode(
        as.character(trt),
        
        "no reversion" =
          "No dose change",
        
        "1-<50%" =
          "<50% reduction",
        
        ">=50%/<=40mg" =
          "\u226550% reduction/\u226440mg",
        
        "increase" =
          "Dose increase",
        
        .default = NA_character_
      ),
      
      trt_label = factor(
        trt_label,
        levels = TRT_LABELS
      )
    ) %>%
    filter(
      !is.na(trt_label)
    )
  
  arm_counts <- purrr::map_dfr(
    breaks_x,
    function(t) {
      
      dat_clean %>%
        group_by(trt_label) %>%
        summarise(
          time = t,
          
          raw_count = if (t == 0) {
            dplyr::n()
          } else {
            sum(
              .data[[time_var]] >= t,
              na.rm = TRUE
            )
          },
          
          .groups = "drop"
        )
    }
  )
  
  total_counts <- purrr::map_dfr(
    breaks_x,
    function(t) {
      
      tibble(
        trt_label = "Total episodes",
        
        time = t,
        
        raw_count = if (t == 0) {
          nrow(dat_clean)
        } else {
          sum(
            dat_clean[[time_var]] >= t,
            na.rm = TRUE
          )
        }
      )
    }
  )
  
  table_levels <- c(
    TRT_LABELS,
    "Total episodes"
  )
  
  at_risk_df <- bind_rows(
    arm_counts,
    total_counts
  ) %>%
    mutate(
      trt_label = factor(
        as.character(trt_label),
        levels = rev(table_levels)
      ),
      
      label_text = case_when(
        is.na(raw_count) ~ "",
        
        raw_count == 0 ~ "0",
        
        raw_count > 0 &
          raw_count < 10 ~ "<10",
        
        TRUE ~ scales::comma(
          raw_count,
          accuracy = 1
        )
      ),
      
      total_row =
        as.character(trt_label) ==
        "Total episodes"
    )
  
  ggplot(
    at_risk_df,
    aes(
      x = time,
      y = trt_label,
      label = label_text
    )
  ) +
    geom_text(
      aes(
        fontface = ifelse(
          total_row,
          "bold",
          "plain"
        )
      ),
      
      size = TABLE_VALUE_SIZE,
      
      family = FIG_FONT_FAMILY,
      
      color = "black",
      
      show.legend = FALSE
    ) +
    
    scale_x_continuous(
      breaks = breaks_x,
      
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    
    scale_y_discrete(
      drop = FALSE,
      
      expand = expansion(
        mult = c(
          0.008,
          0.008
        )
      )
    ) +
    
    coord_cartesian(
      xlim = c(
        X_DISPLAY_MIN,
        X_DISPLAY_MAX
      ),
      
      clip = "off"
    ) +
    
    labs(
      x = NULL,
      
      y = NULL,
      
      title = if (show_table_title) {
        "Number of episodes at risk during follow-up"
      } else {
        NULL
      }
    ) +
    
    theme_minimal(
      base_size = 12,
      
      base_family = FIG_FONT_FAMILY
    ) +
    
    theme(
      panel.grid =
        element_blank(),
      
      axis.line =
        element_blank(),
      
      axis.ticks =
        element_blank(),
      
      axis.text.x =
        element_blank(),
      
      axis.title.x =
        element_blank(),
      
      axis.title.y =
        element_blank(),
      
      axis.text.y = if (show_strategy_names) {
        element_text(
          color = "black",
          
          size = TABLE_LABEL_SIZE,
          
          family = FIG_FONT_FAMILY,
          
          margin = ggplot2::margin(
            r = 5
          )
        )
      } else {
        element_blank()
      },
      
      plot.title =
        element_text(
          color = "black",
          
          size = TABLE_TITLE_SIZE,
          
          face = "bold",
          
          family = FIG_FONT_FAMILY,
          
          hjust = 0,
          
          margin = ggplot2::margin(
            t = 0,
            r = 0,
            b = 10,
            l = 0,
            unit = "pt"
          )
        ),
      
      plot.margin =
        ggplot2::margin(
          t = 0,
          r = 22,
          b = 0,
          l = if (show_strategy_names) {
            22
          } else {
            10
          },
          unit = "pt"
        )
    )
}


# 9. ADJUSTED CURVE PLOT #########
#
# The visible plotting limits are controlled with coord_cartesian().
# This prevents ggplot2 from deleting observations above y_max before
# geom_line() joins the curve points.

plot_adj_curve <- function(
    df,
    y_max,
    y_break_interval,
    show_30d_line = TRUE,
    y_lab = "Adjusted risk (%)",
    title_nm = NULL
) {
  
  required_columns <- c(
    "horizon",
    "risk",
    "trt_label"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(df)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Curve data are missing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  maximum_observed_risk <- suppressWarnings(
    max(
      df$risk,
      na.rm = TRUE
    )
  )
  
  if (
    is.finite(maximum_observed_risk) &&
    maximum_observed_risk > y_max
  ) {
    warning(
      sprintf(
        paste0(
          "Observed maximum risk is %.3f, but y_max is %.3f. ",
          "The curve will be visually truncated. Increase y_max if the ",
          "complete curve should be shown."
        ),
        maximum_observed_risk,
        y_max
      )
    )
  }
  
  y_axis_breaks <- seq(
    from = 0,
    
    to = y_max,
    
    by = y_break_interval
  )
  
  y_axis_breaks <- sort(
    unique(
      c(
        y_axis_breaks,
        y_max
      )
    )
  )
  
  y_axis_breaks <- y_axis_breaks[
    y_axis_breaks >= 0 &
      y_axis_breaks <= y_max
  ]
  
  p <- ggplot(
    df,
    aes(
      x = horizon,
      y = risk,
      color = trt_label,
      group = trt_label
    )
  ) +
    
    geom_line(
      linewidth = CURVE_LINE_WIDTH,
      lineend = "round",
      linejoin = "round",
      
      # Do not silently remove missing values.
      na.rm = FALSE
    ) +
    
    scale_x_continuous(
      breaks = breaks_x,
      
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    
    # No scale limits here, because scale limits remove data.
    scale_y_continuous(
      labels = scales::percent_format(
        accuracy = 1
      ),
      
      breaks = y_axis_breaks,
      
      expand = expansion(
        mult = c(0, 0.02)
      )
    ) +
    
    scale_color_manual(
      values = trt_colors_km,
      
      breaks = TRT_LABELS,
      
      limits = TRT_LABELS,
      
      drop = FALSE
    ) +
    
    coord_cartesian(
      xlim = c(
        X_DISPLAY_MIN,
        X_DISPLAY_MAX
      ),
      
      ylim = c(
        0,
        y_max
      ),
      
      clip = "off"
    ) +
    
    labs(
      title = title_nm,
      
      x = "Days since time zero",
      
      y = y_lab,
      
      color = NULL
    ) +
    
    theme_minimal(
      base_size = 15,
      
      base_family = FIG_FONT_FAMILY
    ) +
    
    theme(
      plot.background =
        element_rect(
          fill = "white",
          color = NA
        ),
      
      panel.background =
        element_rect(
          fill = "white",
          color = NA
        ),
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major.x =
        element_line(
          color = "grey90",
          linewidth = GRID_LINE_WIDTH
        ),
      
      panel.grid.major.y =
        element_line(
          color = "grey90",
          linewidth = GRID_LINE_WIDTH
        ),
      
      axis.line =
        element_line(
          color = "black",
          linewidth = AXIS_LINE_WIDTH
        ),
      
      axis.ticks =
        element_line(
          color = "black",
          linewidth = AXIS_LINE_WIDTH
        ),
      
      axis.ticks.length =
        grid::unit(
          0.20,
          "cm"
        ),
      
      axis.text =
        element_text(
          color = "black",
          size = AXIS_TEXT_SIZE,
          family = FIG_FONT_FAMILY
        ),
      
      axis.title.x =
        element_text(
          color = "black",
          size = AXIS_TITLE_SIZE,
          family = FIG_FONT_FAMILY,
          margin = ggplot2::margin(
            t = 8
          )
        ),
      
      axis.title.y =
        element_text(
          color = "black",
          size = AXIS_TITLE_SIZE,
          family = FIG_FONT_FAMILY,
          margin = ggplot2::margin(
            r = 8
          )
        ),
      
      plot.title =
        element_text(
          color = "black",
          size = PANEL_TITLE_SIZE,
          face = "bold",
          family = FIG_FONT_FAMILY,
          hjust = 0.5,
          lineheight = 1.05,
          margin = ggplot2::margin(
            t = 2,
            r = 0,
            b = 9,
            l = 0,
            unit = "pt"
          )
        ),
      
      legend.position =
        "none",
      
      plot.margin =
        ggplot2::margin(
          t = 12,
          r = 22,
          b = 0,
          l = 22,
          unit = "pt"
        )
    )
  
  if (show_30d_line) {
    p <- p +
      geom_vline(
        xintercept = 30,
        linetype = "dashed",
        color = "grey50",
        linewidth = REF_LINE_WIDTH
      )
  }
  
  p
}


# 10. PREPARE CURVE DATA #######

# Overdose
df_od_tr2_le <- prep_adj_curve(
  risk_tr2_le30d$od$curve
)

df_od_tr2_gt <- prep_adj_curve(
  risk_tr2_gt30d$od$curve
)

df_od_tr3_le <- prep_adj_curve(
  risk_tr3_le30d$od$curve
)

df_od_tr3_gt <- prep_adj_curve(
  risk_tr3_gt30d$od$curve
)


# Discontinuation
df_dis_tr2_le <- prep_adj_curve(
  risk_tr2_le30d$dis$curve
)

df_dis_tr2_gt <- prep_adj_curve(
  risk_tr2_gt30d$dis$curve
)

df_dis_tr3_le <- prep_adj_curve(
  risk_tr3_le30d$dis$curve
)

df_dis_tr3_gt <- prep_adj_curve(
  risk_tr3_gt30d$dis$curve
)


# 11. CURVE COMPLETENESS DIAGNOSTICS #########

check_curve_completeness <- function(
    df,
    panel_name,
    selected_y_max
) {
  
  df %>%
    group_by(trt_label) %>%
    summarise(
      panel = panel_name,
      
      n_rows = n(),
      
      n_horizons =
        n_distinct(horizon),
      
      minimum_horizon =
        suppressWarnings(
          min(
            horizon,
            na.rm = TRUE
          )
        ),
      
      maximum_horizon =
        suppressWarnings(
          max(
            horizon,
            na.rm = TRUE
          )
        ),
      
      missing_risk =
        sum(
          is.na(risk)
        ),
      
      maximum_risk =
        suppressWarnings(
          max(
            risk,
            na.rm = TRUE
          )
        ),
      
      selected_y_max =
        selected_y_max,
      
      exceeds_axis =
        maximum_risk >
        selected_y_max,
      
      .groups = "drop"
    )
}


curve_completeness <- bind_rows(
  
  check_curve_completeness(
    df_od_tr2_le,
    "OD: 3-4 days, within 30 days",
    OD_Y_MAX
  ),
  
  check_curve_completeness(
    df_od_tr2_gt,
    "OD: 3-4 days, after 30 days",
    OD_Y_MAX
  ),
  
  check_curve_completeness(
    df_od_tr3_le,
    "OD: 5-14 days, within 30 days",
    OD_Y_MAX
  ),
  
  check_curve_completeness(
    df_od_tr3_gt,
    "OD: 5-14 days, after 30 days",
    OD_Y_MAX
  ),
  
  check_curve_completeness(
    df_dis_tr2_le,
    "DIS: 3-4 days, within 30 days",
    DIS_Y_MAX
  ),
  
  check_curve_completeness(
    df_dis_tr2_gt,
    "DIS: 3-4 days, after 30 days",
    DIS_Y_MAX
  ),
  
  check_curve_completeness(
    df_dis_tr3_le,
    "DIS: 5-14 days, within 30 days",
    DIS_Y_MAX
  ),
  
  check_curve_completeness(
    df_dis_tr3_gt,
    "DIS: 5-14 days, after 30 days",
    DIS_Y_MAX
  )
)

cat(
  "\n=== Curve completeness and visible-range diagnostic ===\n"
)

print(
  curve_completeness,
  n = Inf
)

write_csv(
  curve_completeness,
  paste0(
    FIG_DIR,
    "Adjusted_Curve_Completeness_Diagnostic.csv"
  )
)


# 12. BUILD OVERDOSE CURVES #########

p_od_tr2_le <- plot_adj_curve(
  df = df_od_tr2_le,
  
  y_max = OD_Y_MAX,
  
  y_break_interval = OD_Y_INTERVAL,
  
  show_30d_line = TRUE,
  
  y_lab = "Adjusted risk (%)",
  
  title_nm =
    "Interruption of 3-4 days within 30 days of OAT"
)


p_od_tr2_gt <- plot_adj_curve(
  df = df_od_tr2_gt,
  
  y_max = OD_Y_MAX,
  
  y_break_interval = OD_Y_INTERVAL,
  
  show_30d_line = TRUE,
  
  y_lab = NULL,
  
  title_nm =
    "Interruption of 3-4 days after 30 days of OAT"
)


p_od_tr3_le <- plot_adj_curve(
  df = df_od_tr3_le,
  
  y_max = OD_Y_MAX,
  
  y_break_interval = OD_Y_INTERVAL,
  
  show_30d_line = TRUE,
  
  y_lab = "Adjusted risk (%)",
  
  title_nm =
    "Interruption of 5-14 days within 30 days of OAT"
)


p_od_tr3_gt <- plot_adj_curve(
  df = df_od_tr3_gt,
  
  y_max = OD_Y_MAX,
  
  y_break_interval = OD_Y_INTERVAL,
  
  show_30d_line = TRUE,
  
  y_lab = NULL,
  
  title_nm =
    "Interruption of 5-14 days after 30 days of OAT"
)


# 13. BUILD OVERDOSE RISK TABLES ##########

t_od_tr2_le <- build_at_risk_table(
  dat = dat_tr2_le30d,
  
  time_var = "time_od_acute",
  
  show_strategy_names = TRUE,
  
  show_table_title = TRUE
)


t_od_tr2_gt <- build_at_risk_table(
  dat = dat_tr2_gt30d,
  
  time_var = "time_od_acute",
  
  show_strategy_names = FALSE,
  
  show_table_title = FALSE
)


t_od_tr3_le <- build_at_risk_table(
  dat = dat_tr3_le30d,
  
  time_var = "time_od_acute",
  
  show_strategy_names = TRUE,
  
  show_table_title = TRUE
)


t_od_tr3_gt <- build_at_risk_table(
  dat = dat_tr3_gt30d,
  
  time_var = "time_od_acute",
  
  show_strategy_names = FALSE,
  
  show_table_title = FALSE
)


# 14. BUILD DISCONTINUATION CURVES ##########

p_dis_tr2_le <- plot_adj_curve(
  df = df_dis_tr2_le,
  
  y_max = DIS_Y_MAX,
  
  y_break_interval = DIS_Y_INTERVAL,
  
  show_30d_line = FALSE,
  
  y_lab = "Adjusted risk (%)",
  
  title_nm =
    "Interruption of 3-4 days within 30 days of OAT"
)


p_dis_tr2_gt <- plot_adj_curve(
  df = df_dis_tr2_gt,
  
  y_max = DIS_Y_MAX,
  
  y_break_interval = DIS_Y_INTERVAL,
  
  show_30d_line = FALSE,
  
  y_lab = NULL,
  
  title_nm =
    "Interruption of 3-4 days after 30 days of OAT"
)


p_dis_tr3_le <- plot_adj_curve(
  df = df_dis_tr3_le,
  
  y_max = DIS_Y_MAX,
  
  y_break_interval = DIS_Y_INTERVAL,
  
  show_30d_line = FALSE,
  
  y_lab = "Adjusted risk (%)",
  
  title_nm =
    "Interruption of 5-14 days within 30 days of OAT"
)


p_dis_tr3_gt <- plot_adj_curve(
  df = df_dis_tr3_gt,
  
  y_max = DIS_Y_MAX,
  
  y_break_interval = DIS_Y_INTERVAL,
  
  show_30d_line = FALSE,
  
  y_lab = NULL,
  
  title_nm =
    "Interruption of 5-14 days after 30 days of OAT"
)


# 15. BUILD DISCONTINUATION RISK TABLES ########

t_dis_tr2_le <- build_at_risk_table(
  dat = dat_tr2_le30d,
  
  time_var = "time_discont",
  
  show_strategy_names = TRUE,
  
  show_table_title = TRUE
)


t_dis_tr2_gt <- build_at_risk_table(
  dat = dat_tr2_gt30d,
  
  time_var = "time_discont",
  
  show_strategy_names = FALSE,
  
  show_table_title = FALSE
)


t_dis_tr3_le <- build_at_risk_table(
  dat = dat_tr3_le30d,
  
  time_var = "time_discont",
  
  show_strategy_names = TRUE,
  
  show_table_title = TRUE
)


t_dis_tr3_gt <- build_at_risk_table(
  dat = dat_tr3_gt30d,
  
  time_var = "time_discont",
  
  show_strategy_names = FALSE,
  
  show_table_title = FALSE
)


# 16. CREATE SHARED LEGEND ##########

legend_plot <- ggplot(
  df_od_tr2_le,
  aes(
    x = horizon,
    y = risk,
    color = trt_label,
    group = trt_label
  )
) +
  
  geom_line(
    linewidth = CURVE_LINE_WIDTH,
    lineend = "round"
  ) +
  
  scale_color_manual(
    values = trt_colors_km,
    breaks = TRT_LABELS,
    limits = TRT_LABELS,
    drop = FALSE
  ) +
  
  labs(
    color = NULL
  ) +
  
  theme_void(
    base_family = FIG_FONT_FAMILY
  ) +
  
  theme(
    legend.position =
      "bottom",
    
    legend.justification =
      "center",
    
    legend.box =
      "horizontal",
    
    legend.title =
      element_blank(),
    
    legend.text =
      element_text(
        color = "black",
        size = LEGEND_TEXT_SIZE,
        family = FIG_FONT_FAMILY
      ),
    
    legend.key.width =
      grid::unit(
        2.4,
        "cm"
      ),
    
    legend.key.height =
      grid::unit(
        0.45,
        "cm"
      ),
    
    legend.spacing.x =
      grid::unit(
        0.35,
        "cm"
      ),
    
    legend.margin =
      ggplot2::margin(
        t = 5,
        r = 0,
        b = 0,
        l = 0,
        unit = "pt"
      )
  )

shared_legend <- cowplot::get_legend(
  legend_plot
)


# 17. BUILD ONE CURVE-TABLE PANEL ########

build_curve_table_panel <- function(
    p_curve,
    p_table,
    spacer_height = CURVE_TABLE_SPACER
) {
  
  p_spacer <- ggplot() +
    theme_void() +
    theme(
      plot.margin = ggplot2::margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0,
        unit = "pt"
      )
    )
  
  patchwork::wrap_plots(
    list(
      p_curve,
      p_spacer,
      p_table
    ),
    
    ncol = 1,
    
    heights = c(
      5.15,
      spacer_height,
      1.25
    )
  )
}


# 18. BUILD OVERDOSE CURVE-TABLE PANELS #########

panel_od_tr2_le <- build_curve_table_panel(
  p_curve = p_od_tr2_le,
  p_table = t_od_tr2_le
)

panel_od_tr2_gt <- build_curve_table_panel(
  p_curve = p_od_tr2_gt,
  p_table = t_od_tr2_gt
)

panel_od_tr3_le <- build_curve_table_panel(
  p_curve = p_od_tr3_le,
  p_table = t_od_tr3_le
)

panel_od_tr3_gt <- build_curve_table_panel(
  p_curve = p_od_tr3_gt,
  p_table = t_od_tr3_gt
)


# 19. BUILD DISCONTINUATION CURVE-TABLE PANELS ########

panel_dis_tr2_le <- build_curve_table_panel(
  p_curve = p_dis_tr2_le,
  p_table = t_dis_tr2_le
)

panel_dis_tr2_gt <- build_curve_table_panel(
  p_curve = p_dis_tr2_gt,
  p_table = t_dis_tr2_gt
)

panel_dis_tr3_le <- build_curve_table_panel(
  p_curve = p_dis_tr3_le,
  p_table = t_dis_tr3_le
)

panel_dis_tr3_gt <- build_curve_table_panel(
  p_curve = p_dis_tr3_gt,
  p_table = t_dis_tr3_gt
)


# 20. COMBINE FOUR PANELS SAFELY USING COWPLOT ########

build_combined_figure <- function(
    panel1,
    panel2,
    panel3,
    panel4,
    shared_legend
) {
  
  top_row <- cowplot::plot_grid(
    panel1,
    panel2,
    
    ncol = 2,
    
    align = "h",
    
    axis = "tb",
    
    rel_widths = c(
      1,
      1
    ),
    
    labels = NULL
  )
  
  bottom_row <- cowplot::plot_grid(
    panel3,
    panel4,
    
    ncol = 2,
    
    align = "h",
    
    axis = "tb",
    
    rel_widths = c(
      1,
      1
    ),
    
    labels = NULL
  )
  
  main_grid <- cowplot::plot_grid(
    top_row,
    bottom_row,
    
    ncol = 1,
    
    align = "v",
    
    axis = "lr",
    
    rel_heights = c(
      1,
      1
    ),
    
    labels = NULL
  )
  
  cowplot::plot_grid(
    main_grid,
    shared_legend,
    
    ncol = 1,
    
    rel_heights = c(
      1,
      0.07
    ),
    
    labels = NULL
  )
}


# 21. BUILD COMBINED FIGURES ########

od_combined <- build_combined_figure(
  panel1 = panel_od_tr2_le,
  panel2 = panel_od_tr2_gt,
  panel3 = panel_od_tr3_le,
  panel4 = panel_od_tr3_gt,
  shared_legend = shared_legend
)


dis_combined <- build_combined_figure(
  panel1 = panel_dis_tr2_le,
  panel2 = panel_dis_tr2_gt,
  panel3 = panel_dis_tr3_le,
  panel4 = panel_dis_tr3_gt,
  shared_legend = shared_legend
)


### SAVE HIGH-RESOLUTION JPEG FILES ##########

jpeg(
  filename = OD_FILE,
  width = 5600,
  height = 4200,
  res = 300,
  quality = 75,
  bg = "white"
)

print(od_combined)

dev.off()


jpeg(
  filename = DIS_FILE,
  width = 5600,
  height = 4200,
  res = 300,
  quality = 75,
  bg = "white"
)

print(dis_combined)

dev.off()
