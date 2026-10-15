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




options(error = NULL)
library(tidyverse)
library(patchwork)
library(cowplot)
library(grid)

# 1. PATHS & OUTPUT FILES
RES_DIR <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/"
FIG_DIR <- paste0(RES_DIR, "LE30D/Figures/")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

OD_FILE  <- paste0(FIG_DIR, "Adjusted_Risk_OD_Combined.jpg")
DIS_FILE <- paste0(FIG_DIR, "Adjusted_Risk_Discont_Combined.jpg")

# 2. LOAD DATA OBJECTS
risk_tr2_le30d <- readRDS(paste0(RES_DIR, "LE30D/Risk/risk_tr2.rds"))
risk_tr2_gt30d <- readRDS(paste0(RES_DIR, "GT30D/Risk/risk_tr2.rds"))
risk_tr3_le30d <- readRDS(paste0(RES_DIR, "LE30D/Risk/risk_tr3.rds"))
risk_tr3_gt30d <- readRDS(paste0(RES_DIR, "GT30D/Risk/risk_tr3.rds"))

dat_tr2_le30d <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr2_34d_all.rds"))
dat_tr2_gt30d <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr2_34d_all.rds"))
dat_tr3_le30d <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr3_all.rds"))
dat_tr3_gt30d <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr3_all.rds"))

# 3. SETTINGS & CONSTANTS
TRT_LABELS <- c("No dose change", "<50% reduction", "\u226550% reduction or \u226440mg", "Dose increase")
trt_colors_km    <- c("No dose change" = "#000000", "<50% reduction" = "#2CA02C", "\u226550% reduction or \u226440mg" = "#FF7F0E", "Dose increase" = "#D62728")
trt_linetypes_km <- c("No dose change" = "solid",   "<50% reduction" = "dashed",  "\u226550% reduction or \u226440mg" = "dotted",  "Dose increase" = "dotdash")

breaks_x <- c(0, 30, 90, 180, 270, 365)
FIG_FONT <- "sans"
CURVE_TABLE_SPACER <- 0.22
X_MIN <- -28; X_MAX <- 383
OD_Y_MAX <- 0.25; OD_Y_INT <- 0.05
DIS_Y_MAX <- 0.60; DIS_Y_INT <- 0.10

# 4. HELPER FUNCTIONS
prep_adj_curve <- function(risk_obj) {
  if (is.null(risk_obj$risks)) stop("Supplied object lacks $risks.")
  curve_df <- risk_obj$risks %>%
    transmute(
      horizon = as.numeric(horizon), risk = as.numeric(risk),
      trt_label = factor(recode(as.character(trt),
                                "no reversion" = "No dose change", "1-<50%" = "<50% reduction",
                                ">=50%/<=40mg" = "\u226550% reduction or \u226440mg", "increase" = "Dose increase"
      ), levels = TRT_LABELS)
    ) %>% filter(!is.na(horizon), !is.na(trt_label))
  
  missing_day0 <- setdiff(TRT_LABELS, curve_df %>% filter(horizon == 0) %>% pull(trt_label) %>% as.character())
  day0_rows <- tibble(horizon = 0, risk = 0, trt_label = factor(missing_day0, levels = TRT_LABELS))
  
  bind_rows(curve_df, day0_rows) %>% distinct(trt_label, horizon, .keep_all = TRUE) %>% arrange(trt_label, horizon)
}

extract_data_weights <- function(dat_obj) {
  if (is.list(dat_obj) && "data_weights" %in% names(dat_obj)) return(dat_obj$data_weights)
  if (is.data.frame(dat_obj)) return(dat_obj)
  stop("Invalid PS object.")
}

build_at_risk_table <- function(dat, time_var, show_strategy_names = TRUE, show_table_title = TRUE) {
  dat_clean <- extract_data_weights(dat) %>%
    filter(!is.na(trt), !is.na(.data[[time_var]])) %>%
    mutate(trt_label = factor(recode(as.character(trt),
                                     "no reversion" = "No dose change", "1-<50%" = "<50% reduction",
                                     ">=50%/<=40mg" = "\u226550% reduction or \u226440mg", "increase" = "Dose increase"
    ), levels = TRT_LABELS)) %>% filter(!is.na(trt_label))
  
  arm_counts <- map_dfr(breaks_x, ~ dat_clean %>% group_by(trt_label) %>%
                          summarise(time = .x, raw_count = if (.x == 0) n() else sum(.data[[time_var]] >= .x, na.rm = TRUE), .groups = "drop"))
  
  total_counts <- map_dfr(breaks_x, ~ tibble(trt_label = "Total episodes", time = .x,
                                             raw_count = if (.x == 0) nrow(dat_clean) else sum(dat_clean[[time_var]] >= .x, na.rm = TRUE)))
  
  at_risk_df <- bind_rows(arm_counts, total_counts) %>%
    mutate(
      trt_label = factor(as.character(trt_label), levels = rev(c(TRT_LABELS, "Total episodes"))),
      label_text = case_when(is.na(raw_count) ~ "", raw_count == 0 ~ "0", raw_count < 10 ~ "<10", TRUE ~ scales::comma(raw_count, accuracy = 1)),
      total_row = (as.character(trt_label) == "Total episodes")
    )
  
  ggplot(at_risk_df, aes(x = time, y = trt_label, label = label_text)) +
    geom_text(aes(fontface = ifelse(total_row, "bold", "plain"),
                  hjust = ifelse(time == 0, 0.8, ifelse(time == 30, 0.4, 0.5))),
              size = 5, family = FIG_FONT, color = "black") +
    scale_x_continuous(breaks = breaks_x, expand = c(0, 0)) +
    scale_y_discrete(
      drop = FALSE, 
      expand = c(0.008, 0.008),
      labels = function(labs) {
        lapply(labs, function(x) if (x == "Total episodes") expression(bold("Total episodes"))[[1]] else x)
      }
    ) +
    coord_cartesian(xlim = c(X_MIN, X_MAX), clip = "off") +
    labs(
      x = NULL, y = NULL, 
      title = if (show_table_title) "Number of episodes at risk during follow-up" else " "
    ) +
    theme_minimal(base_size = 12, base_family = FIG_FONT) +
    theme(
      panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_blank(), axis.title = element_blank(),
      axis.text.y = if (show_strategy_names) element_text(color = "black", size = 15, family = FIG_FONT, margin = ggplot2::margin(r = 5, unit = "pt")) else element_blank(),
      plot.title = element_text(
        color = if (show_table_title) "black" else "transparent", 
        size = 15, face = "bold", family = FIG_FONT, hjust = 0, 
        margin = ggplot2::margin(b = 10, unit = "pt")
      ),
      plot.margin = ggplot2::margin(t = 0, r = 22, b = 0, l = if (show_strategy_names) 22 else 10, unit = "pt")
    )
}

plot_adj_curve <- function(df, y_max, y_interval, show_30d = TRUE, y_lab = "Standardized risk (%)", title_nm = NULL) {
  y_breaks <- sort(unique(c(seq(0, y_max, by = y_interval), y_max)))
  
  p <- ggplot(df, aes(x = horizon, y = risk, color = trt_label, linetype = trt_label, group = trt_label)) +
    geom_line(linewidth = 1.40, lineend = "round", linejoin = "round") +
    scale_x_continuous(breaks = breaks_x, expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = y_breaks, expand = expansion(mult = c(0, 0.02))) +
    scale_color_manual(values = trt_colors_km, breaks = TRT_LABELS, limits = TRT_LABELS, drop = FALSE) +
    scale_linetype_manual(values = trt_linetypes_km, breaks = TRT_LABELS, limits = TRT_LABELS, drop = FALSE) +
    coord_cartesian(xlim = c(X_MIN, X_MAX), ylim = c(0, y_max), clip = "off") +
    labs(title = title_nm, x = "Days since time zero", y = y_lab, color = NULL, linetype = NULL) +
    theme_minimal(base_size = 15, base_family = FIG_FONT) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.40),
      axis.line = element_line(color = "black", linewidth = 0.85),
      axis.ticks = element_line(color = "black", linewidth = 0.85),
      axis.ticks.length = grid::unit(0.20, "cm"),
      axis.text = element_text(color = "black", size = 15, family = FIG_FONT),
      axis.title.x = element_text(color = "black", size = 17, family = FIG_FONT, margin = ggplot2::margin(t = 8, unit = "pt")),
      axis.title.y = element_text(color = "black", size = 17, family = FIG_FONT, margin = ggplot2::margin(r = 8, unit = "pt")),
      plot.title = element_text(color = "black", size = 18, face = "bold", family = FIG_FONT, hjust = 0.5, lineheight = 1.05, margin = ggplot2::margin(t = 2, b = 9, unit = "pt")),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 12, r = 22, b = 0, l = 22, unit = "pt")
    )
  if (show_30d) p <- p + geom_vline(xintercept = 30, linetype = "dashed", color = "grey50", linewidth = 0.90)
  p
}

# 5. PREPARE CURVE DATA
df_od_tr2_le  <- prep_adj_curve(risk_tr2_le30d$od$curve)
df_od_tr2_gt  <- prep_adj_curve(risk_tr2_gt30d$od$curve)
df_od_tr3_le  <- prep_adj_curve(risk_tr3_le30d$od$curve)
df_od_tr3_gt  <- prep_adj_curve(risk_tr3_gt30d$od$curve)

df_dis_tr2_le <- prep_adj_curve(risk_tr2_le30d$dis$curve)
df_dis_tr2_gt <- prep_adj_curve(risk_tr2_gt30d$dis$curve)
df_dis_tr3_le <- prep_adj_curve(risk_tr3_le30d$dis$curve)
df_dis_tr3_gt <- prep_adj_curve(risk_tr3_gt30d$dis$curve)

# 6. DIAGNOSTICS
check_curve_completeness <- function(df, panel_nm, selected_y_max) {
  df %>% group_by(trt_label) %>% summarise(
    panel = panel_nm, n_rows = n(), n_horizons = n_distinct(horizon),
    minimum_horizon = suppressWarnings(min(horizon, na.rm = TRUE)),
    maximum_horizon = suppressWarnings(max(horizon, na.rm = TRUE)),
    missing_risk = sum(is.na(risk)), maximum_risk = suppressWarnings(max(risk, na.rm = TRUE)),
    selected_y_max = selected_y_max, exceeds_axis = maximum_risk > selected_y_max, .groups = "drop"
  )
}
curve_completeness <- bind_rows(
  check_curve_completeness(df_od_tr2_le,  "OD: 3-4 days, <=30d",  OD_Y_MAX),
  check_curve_completeness(df_od_tr2_gt,  "OD: 3-4 days, >30d",   OD_Y_MAX),
  check_curve_completeness(df_od_tr3_le,  "OD: 5-14 days, <=30d", OD_Y_MAX),
  check_curve_completeness(df_od_tr3_gt,  "OD: 5-14 days, >30d",  OD_Y_MAX),
  check_curve_completeness(df_dis_tr2_le, "DIS: 3-4 days, <=30d", DIS_Y_MAX),
  check_curve_completeness(df_dis_tr2_gt, "DIS: 3-4 days, >30d",  DIS_Y_MAX),
  check_curve_completeness(df_dis_tr3_le, "DIS: 5-14 days, <=30d",DIS_Y_MAX),
  check_curve_completeness(df_dis_tr3_gt, "DIS: 5-14 days, >30d", DIS_Y_MAX)
)
print(curve_completeness, n = Inf)
write_csv(curve_completeness, paste0(FIG_DIR, "Adjusted_Curve_Completeness_Diagnostic.csv"))

# 7. GENERATE CURVES & RISK TABLES
p_od_tr2_le <- plot_adj_curve(df_od_tr2_le, OD_Y_MAX, OD_Y_INT, TRUE, "Standardized risk (%)", "Interruptions of 3-4 days, within 30 days of initiation")
p_od_tr2_gt <- plot_adj_curve(df_od_tr2_gt, OD_Y_MAX, OD_Y_INT, TRUE, NULL,                 "Interruptions of 3-4 days, after 30 days of initiation")
p_od_tr3_le <- plot_adj_curve(df_od_tr3_le, OD_Y_MAX, OD_Y_INT, TRUE, "Standardized risk (%)", "Interruptions of 5-14 days, within 30 days of initiation")
p_od_tr3_gt <- plot_adj_curve(df_od_tr3_gt, OD_Y_MAX, OD_Y_INT, TRUE, NULL,                 "Interruptions of 5-14 days, after 30 days of initiation")

t_od_tr2_le <- build_at_risk_table(dat_tr2_le30d, "time_od_acute", TRUE,  TRUE)
t_od_tr2_gt <- build_at_risk_table(dat_tr2_gt30d, "time_od_acute", FALSE, FALSE)
t_od_tr3_le <- build_at_risk_table(dat_tr3_le30d, "time_od_acute", TRUE,  TRUE)
t_od_tr3_gt <- build_at_risk_table(dat_tr3_gt30d, "time_od_acute", FALSE, FALSE)

p_dis_tr2_le <- plot_adj_curve(df_dis_tr2_le, DIS_Y_MAX, DIS_Y_INT, FALSE, "Standardized risk (%)", "Interruptions of 3-4 days, within 30 days of initiation")
p_dis_tr2_gt <- plot_adj_curve(df_dis_tr2_gt, DIS_Y_MAX, DIS_Y_INT, FALSE, NULL,                 "Interruptions of 3-4 days, after 30 days of initiation")
p_dis_tr3_le <- plot_adj_curve(df_dis_tr3_le, DIS_Y_MAX, DIS_Y_INT, FALSE, "Standardized risk (%)", "Interruptions of 5-14 days, within 30 days of initiation")
p_dis_tr3_gt <- plot_adj_curve(df_dis_tr3_gt, DIS_Y_MAX, DIS_Y_INT, FALSE, NULL,                 "Interruptions of 5-14 days, after 30 days of initiation")

t_dis_tr2_le <- build_at_risk_table(dat_tr2_le30d, "time_discont", TRUE,  TRUE)
t_dis_tr2_gt <- build_at_risk_table(dat_tr2_gt30d, "time_discont", FALSE, FALSE)
t_dis_tr3_le <- build_at_risk_table(dat_tr3_le30d, "time_discont", TRUE,  TRUE)
t_dis_tr3_gt <- build_at_risk_table(dat_tr3_gt30d, "time_discont", FALSE, FALSE)

# 8. SHARED LEGEND
legend_plot <- ggplot(df_od_tr2_le, aes(x = horizon, y = risk, color = trt_label, linetype = trt_label, group = trt_label)) +
  geom_line(linewidth = 1.40, lineend = "round") +
  scale_color_manual(values = trt_colors_km, breaks = TRT_LABELS, limits = TRT_LABELS, drop = FALSE) +
  scale_linetype_manual(values = trt_linetypes_km, breaks = TRT_LABELS, limits = TRT_LABELS, drop = FALSE) +
  labs(color = NULL, linetype = NULL) + theme_void(base_family = FIG_FONT) +
  theme(
    legend.position = "bottom", legend.justification = "center", legend.box = "horizontal",
    legend.title = element_blank(), legend.text = element_text(color = "black", size = 15, family = FIG_FONT),
    legend.key.width = grid::unit(2.4, "cm"), legend.key.height = grid::unit(0.45, "cm"),
    legend.spacing.x = grid::unit(0.35, "cm"), legend.margin = ggplot2::margin(t = 5, unit = "pt")
  )
shared_legend <- cowplot::get_legend(legend_plot)

# 9. BUILD SUBPANELS & COMBINED FIGURES
build_curve_table_panel <- function(p_curve, p_table, spacer_height = CURVE_TABLE_SPACER) {
  p_spacer <- ggplot() + theme_void() + theme(plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  patchwork::wrap_plots(list(p_curve, p_spacer, p_table), ncol = 1, heights = c(5.15, spacer_height, 1.25))
}

panel_od_tr2_le <- build_curve_table_panel(p_od_tr2_le, t_od_tr2_le)
panel_od_tr2_gt <- build_curve_table_panel(p_od_tr2_gt, t_od_tr2_gt)
panel_od_tr3_le <- build_curve_table_panel(p_od_tr3_le, t_od_tr3_le)
panel_od_tr3_gt <- build_curve_table_panel(p_od_tr3_gt, t_od_tr3_gt)

panel_dis_tr2_le <- build_curve_table_panel(p_dis_tr2_le, t_dis_tr2_le)
panel_dis_tr2_gt <- build_curve_table_panel(p_dis_tr2_gt, t_dis_tr2_gt)
panel_dis_tr3_le <- build_curve_table_panel(p_dis_tr3_le, t_dis_tr3_le)
panel_dis_tr3_gt <- build_curve_table_panel(p_dis_tr3_gt, t_dis_tr3_gt)

build_combined_figure <- function(panel1, panel2, panel3, panel4, legend, rel_widths = c(1.28, 1)) {
  top_row    <- cowplot::plot_grid(panel1, panel2, ncol = 2, align = "h", axis = "tb", rel_widths = rel_widths)
  bottom_row <- cowplot::plot_grid(panel3, panel4, ncol = 2, align = "h", axis = "tb", rel_widths = rel_widths)
  main_grid  <- cowplot::plot_grid(top_row, bottom_row, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1))
  
  cowplot::plot_grid(main_grid, legend, ncol = 1, rel_heights = c(1, 0.07))
}

od_combined  <- build_combined_figure(panel_od_tr2_le, panel_od_tr2_gt, panel_od_tr3_le, panel_od_tr3_gt, shared_legend)
dis_combined <- build_combined_figure(panel_dis_tr2_le, panel_dis_tr2_gt, panel_dis_tr3_le, panel_dis_tr3_gt, shared_legend)

save_figure <- function(plot_obj, filepath) {
  jpeg(filename = filepath, width = 5600, height = 4300, res = 300, quality = 75, bg = "white")
  print(plot_obj)
  dev.off()
}

save_figure(od_combined, OD_FILE)
save_figure(dis_combined, DIS_FILE)