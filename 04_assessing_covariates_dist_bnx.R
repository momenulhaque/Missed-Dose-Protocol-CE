library(tidyverse)
library(patchwork)
library(ggh4x) 

# 1. Define colors AND explicitly set the requested column order
COMPARISON_COLORS <- c(
  "No dose change"       = "#000000", 
  "1-<50% reduction"     = "#2CA02C",  
  "≥50% reduction/≤4mg" = "#FF7F0E",
  "Dose increase"        = "#D62728"  
)

# Define shapes for the weighting schemes
WEIGHT_SHAPES <- c(
  "Unweighted" = 16, # Circle
  "IPTW"       = 17, # Triangle
  "OW"         = 15  # Square
)

# 2. The Formatted Grid Plotting Function
build_master_grid_plot <- function(ps_obj, trial_title, bvars, cont_vars, var_labels) {
  
  # A. Prep Dataset
  dat <- ps_obj$data_weights %>%
    mutate(
      trt = recode(trt,
                   "no reversion" = "No dose change",
                   "1-<50%" = "1-<50% reduction",
                   ">=50%/<=4mg" = "≥50% reduction/≤4mg",
                   "increase" = "Dose increase"),
      # This locks the column order to match the COMPARISON_COLORS list exactly
      trt = factor(trt, levels = names(COMPARISON_COLORS)),
      wt_unw = 1 
    )
  
  if("incar_grp_bl" %in% colnames(dat)) dat$incar_grp_bl <- ifelse(dat$incar_grp_bl == "1_pastyear", 1, 0)
  dat <- dat %>% mutate(across(all_of(bvars), ~ as.numeric(as.character(.))))
  
  # Calculate Overall Unweighted Means (The "Clinical Middle")
  overall_means <- dat %>%
    summarise(across(all_of(bvars), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Covariate_Raw", values_to = "Overall_Mean")
  
  # B. Calculate Weighted Summaries by Treatment
  get_weighted_summaries <- function(weight_col, scheme_name) {
    dat %>%
      group_by(trt) %>%
      summarise(across(all_of(bvars), ~ weighted.mean(.x, w = !!sym(weight_col), na.rm = TRUE)), .groups = "drop") %>%
      mutate(Scheme = scheme_name)
  }
  
  df_plot <- bind_rows(
    get_weighted_summaries("wt_unw", "Unweighted"),
    get_weighted_summaries("wt_trunc", "IPTW"),
    get_weighted_summaries("ow_wt", "OW")
  ) %>%
    pivot_longer(cols = all_of(bvars), names_to = "Covariate_Raw", values_to = "Value") %>%
    left_join(overall_means, by = "Covariate_Raw")
  
  # C1. Format Values and Clean Names
  df_plot <- df_plot %>%
    mutate(
      Value = ifelse(Covariate_Raw %in% cont_vars, Value, Value * 100),
      Overall_Mean = ifelse(Covariate_Raw %in% cont_vars, Overall_Mean, Overall_Mean * 100),
      Clean_Name = var_labels[Covariate_Raw],
      Clean_Name = ifelse(Covariate_Raw %in% cont_vars, Clean_Name, paste0(Clean_Name, " (%)")),
      Clean_Name = stringr::str_wrap(Clean_Name, width = 30),
      Scheme = factor(Scheme, levels = rev(c("Unweighted", "IPTW", "OW")))
    )
  
  # C2. Enforce the exact row order from var_labels
  ordered_raw_names <- names(var_labels)[names(var_labels) %in% bvars]
  ordered_clean_names <- sapply(ordered_raw_names, function(x) {
    c_name <- var_labels[x]
    c_name <- ifelse(x %in% cont_vars, c_name, paste0(c_name, " (%)"))
    stringr::str_wrap(c_name, width = 30)
  })
  df_plot$Clean_Name <- factor(df_plot$Clean_Name, levels = unique(ordered_clean_names))
  
  # D. CREATE GEOM_BLANK LIMITS
  df_limits <- df_plot %>%
    group_by(Clean_Name) %>%
    summarise(
      min_val = min(c(Value, Overall_Mean), na.rm = TRUE),
      max_val = max(c(Value, Overall_Mean), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      range = max_val - min_val,
      range = ifelse(range == 0, 1, range),
      min_val = min_val - 0.1 * range,
      max_val = max_val + 0.1 * range
    ) %>%
    crossing(trt = factor(levels(df_plot$trt), levels = levels(df_plot$trt))) %>%
    pivot_longer(cols = c(min_val, max_val), values_to = "Value") %>%
    mutate(Scheme = factor("Unweighted", levels = levels(df_plot$Scheme))) 
  
  # E. Build the Ultimate Grid Plot
  p_grid <- ggplot(df_plot, aes(x = Value, y = Scheme, color = trt, shape = Scheme)) +
    geom_blank(data = df_limits, aes(x = Value, y = Scheme)) +
    geom_vline(aes(xintercept = Overall_Mean), linetype = "dashed", color = "gray50", linewidth = 0.8) +
    
    geom_point(size = 4.5, alpha = 0.9) +
    
    ggh4x::facet_grid2(Clean_Name ~ trt, scales = "free_x", independent = "x", switch = "y") +
    
    scale_color_manual(values = COMPARISON_COLORS) +
    scale_shape_manual(values = WEIGHT_SHAPES) + 
    
    labs(
      title = NULL,
      x = "Absolute Mean or Prevalence (%)",
      y = NULL,
      caption = "Note: Dashed vertical line represents the overall unweighted population mean. \nShapes represent weighting schemes: Circle = Unweighted, Triangle = IPTW, Square = OW."
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none", 
      strip.background = element_rect(fill = "gray95", color = "black"),
      strip.text.x = element_text(face = "bold", size = 12), 
      strip.text.y.left = element_text(face = "bold", size = 10, angle = 0, hjust = 1),
      strip.placement = "outside", 
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "gray80"),
      axis.text.y = element_text(face = "bold", color = "black", size = 10),
      axis.text.x = element_text(color = "black", size = 10),
      panel.spacing.y = unit(0.5, "lines"), 
      plot.caption = element_text(hjust = 0, size = 12, color = "black", margin = ggplot2::margin(t = 15))
    )
  
  return(p_grid)
}

# ==============================================================================
# EXECUTION: ALL COVARIATES
# ==============================================================================
full_bvars <- c(names(VAR_LABELS)[names(VAR_LABELS) %in% baselinevars], "prior_episode")
full_cont_vars <- c(continuous_vars, "prior_episode")

# Load Models
dat_tr4_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr4_all.rds"))
dat_tr5_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr5_all.rds"))



# Trial 4
p_grid_tr4_all <- build_master_grid_plot(dat_tr4_all, "3-5 Days Interruption", 
                                         full_bvars, full_cont_vars, VAR_LABELS)


jpeg(paste0(RES_DIR, "Figures/Diagnostic_TR4_MasterGrid.jpg"), 
     width = 5500, height = 7500, res = 300)
print(p_grid_tr4_all)
dev.off()


# Trial 5
p_grid_tr5_all <- build_master_grid_plot(ps_tr5_all, "6-14 Days Interruption", 
                                         full_bvars, full_cont_vars, VAR_LABELS)



jpeg(paste0(RES_DIR, "Figures/Diagnostic_TR5_MasterGrid.jpg"), 
     width = 5500, height = 7500, res = 300)
print(p_grid_tr5_all)
dev.off()

