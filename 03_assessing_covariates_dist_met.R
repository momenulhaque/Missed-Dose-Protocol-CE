# FIGURE: ALL COVARIATES #########

library(tidyverse)
library(patchwork)
library(ggh4x) 

# 1. Define colors AND explicitly set the requested column order
COMPARISON_COLORS <- c(
  "No dose change"       = "#000000", 
  "1-<50% reduction"     = "#2CA02C",  
  "≥50% reduction/≤40mg" = "#FF7F0E",
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
                   ">=50%/<=40mg" = "≥50% reduction/≤40mg",
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


full_bvars <- c(names(VAR_LABELS)[names(VAR_LABELS) %in% baselinevars], "prior_episode")
full_cont_vars <- c(continuous_vars, "prior_episode")

# Load Models
ps_tr2_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
ps_tr3_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

# Trial 2
p_grid_tr2_all <- build_master_grid_plot(ps_tr2_all, "3-4 Days Interruption", 
                                         full_bvars, full_cont_vars, VAR_LABELS)

jpeg(paste0(RES_DIR, "Figures/June09/Diagnostic_TR2_MasterGrid.jpg"), 
     width = 5500, height = 7400, res = 300)
print(p_grid_tr2_all)
dev.off()

# Trial 3
p_grid_tr3_all <- build_master_grid_plot(ps_tr3_all, "5-14 Days Interruption", 
                                         full_bvars, full_cont_vars, VAR_LABELS)

jpeg(paste0(RES_DIR, "Figures/June09/Diagnostic_TR3_MasterGrid.jpg"), 
     width = 5500, height = 7400, res = 300)
print(p_grid_tr3_all)
dev.off()



# FIGURE: STRATIFIED SOME COVS ########

# Purpose:
# Compare the overall overlap-weighted population with the original
# unweighted population and the IPTW-weighted population using only:
#
#   1. Duration of OAT episode before interruption
#   2. Last daily methadone dose before interruption
#
# Four cohorts:
#   - 3-4-day interruption, <=30 days
#   - 3-4-day interruption, >30 days
#   - 5-14-day interruption, <=30 days
#   - 5-14-day interruption, >30 days
#
# Outputs:
#   1. One compact figure
#   2. One CSV with means and exact differences
#
# Importantly, summaries are calculated over the complete weighted population,
# not separately by treatment strategy. This directly characterizes the
# population targeted by overlap weighting.

options(error = NULL)

library(tidyverse)
library(ggh4x)
library(grid)


## 1. PATHS ######


RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/"

FIG_DIR <- paste0(
  RES_DIR,
  "LE30D/Figures/"
)

TABLE_DIR <- paste0(
  RES_DIR,
  "LE30D/Tables/"
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


## 2. LOAD FOUR STRATIFIED PS OBJECTS ########

ps_tr2_le30d <- readRDS(
  paste0(
    RES_DIR,
    "LE30D/PScore/ps_weight_tr2_34d_all.rds"
  )
)

ps_tr2_gt30d <- readRDS(
  paste0(
    RES_DIR,
    "GT30D/PScore/ps_weight_tr2_34d_all.rds"
  )
)

ps_tr3_le30d <- readRDS(
  paste0(
    RES_DIR,
    "LE30D/PScore/ps_weight_tr3_all.rds"
  )
)

ps_tr3_gt30d <- readRDS(
  paste0(
    RES_DIR,
    "GT30D/PScore/ps_weight_tr3_all.rds"
  )
)


## 3. FOCUSED VARIABLES ########

FOCUSED_VARIABLES <- c(
  "ret_length",
  "ddose_prev"
)

FOCUSED_LABELS <- c(
  ret_length =
    "Duration of OAT episode before interruption, days",
  
  ddose_prev =
    "Last daily methadone dose before interruption, mg"
)

SCHEME_LEVELS <- c(
  "Unweighted",
  "IPTW",
  "Overlap weighting"
)

SCHEME_COLORS <- c(
  "Unweighted"       = "#000000",
  "IPTW"             = "#4472C4",
  "Overlap weighting" = "#D62728"
)

SCHEME_SHAPES <- c(
  "Unweighted"        = 16,
  "IPTW"              = 17,
  "Overlap weighting" = 15
)


## 4. DATA AND WEIGHT HELPERS #######

extract_data_weights <- function(ps_obj) {
  
  if (
    is.list(ps_obj) &&
    "data_weights" %in% names(ps_obj)
  ) {
    return(ps_obj$data_weights)
  }
  
  if (is.data.frame(ps_obj)) {
    return(ps_obj)
  }
  
  stop(
    "The supplied object is neither a data frame nor a PS object ",
    "containing data_weights."
  )
}


find_weight_column <- function(
    dat,
    candidates,
    weight_description
) {
  
  available <- candidates[
    candidates %in% names(dat)
  ]
  
  if (length(available) == 0) {
    
    stop(
      "Could not find the ",
      weight_description,
      " column. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  available[1]
}


safe_weighted_mean <- function(
    x,
    w
) {
  
  x <- suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
  
  w <- suppressWarnings(
    as.numeric(
      as.character(w)
    )
  )
  
  valid <- (
    !is.na(x) &
      is.finite(x) &
      !is.na(w) &
      is.finite(w) &
      w >= 0
  )
  
  if (
    sum(valid) == 0 ||
    sum(w[valid]) <= 0
  ) {
    return(NA_real_)
  }
  
  weighted.mean(
    x = x[valid],
    w = w[valid],
    na.rm = TRUE
  )
}


effective_sample_size <- function(w) {
  
  w <- suppressWarnings(
    as.numeric(
      as.character(w)
    )
  )
  
  valid <- (
    !is.na(w) &
      is.finite(w) &
      w >= 0
  )
  
  w <- w[valid]
  
  if (
    length(w) == 0 ||
    sum(w^2) <= 0
  ) {
    return(NA_real_)
  }
  
  sum(w)^2 / sum(w^2)
}


## 5. SUMMARIZE ONE STRATIFIED COHORT #######

summarize_overlap_population <- function(
    ps_obj,
    cohort
) {
  
  dat <- extract_data_weights(
    ps_obj
  )
  
  missing_variables <- setdiff(
    FOCUSED_VARIABLES,
    names(dat)
  )
  
  if (length(missing_variables) > 0) {
    
    stop(
      cohort,
      " is missing: ",
      paste(
        missing_variables,
        collapse = ", "
      )
    )
  }
  
  # The code detects likely weight names rather than assuming a single name.
  iptw_column <- find_weight_column(
    dat = dat,
    candidates = c(
      "wt_trunc",
      "wt_stab",
      "iptw"
    ),
    weight_description = "IPTW"
  )
  
  overlap_column <- find_weight_column(
    dat = dat,
    candidates = c(
      "ow_wt",
      "wt_ow",
      "overlap_weight",
      "overlap_weights"
    ),
    weight_description = "overlap-weight"
  )
  
  dat <- dat %>%
    mutate(
      across(
        all_of(FOCUSED_VARIABLES),
        ~ suppressWarnings(
          as.numeric(
            as.character(.x)
          )
        )
      ),
      
      wt_unweighted = 1
    )
  
  means_long <- bind_rows(
    
    tibble(
      Cohort = cohort,
      Characteristic_Raw = FOCUSED_VARIABLES,
      Scheme = "Unweighted",
      Mean = map_dbl(
        FOCUSED_VARIABLES,
        function(variable) {
          safe_weighted_mean(
            x = dat[[variable]],
            w = dat$wt_unweighted
          )
        }
      )
    ),
    
    tibble(
      Cohort = cohort,
      Characteristic_Raw = FOCUSED_VARIABLES,
      Scheme = "IPTW",
      Mean = map_dbl(
        FOCUSED_VARIABLES,
        function(variable) {
          safe_weighted_mean(
            x = dat[[variable]],
            w = dat[[iptw_column]]
          )
        }
      )
    ),
    
    tibble(
      Cohort = cohort,
      Characteristic_Raw = FOCUSED_VARIABLES,
      Scheme = "Overlap weighting",
      Mean = map_dbl(
        FOCUSED_VARIABLES,
        function(variable) {
          safe_weighted_mean(
            x = dat[[variable]],
            w = dat[[overlap_column]]
          )
        }
      )
    )
  ) %>%
    mutate(
      Characteristic = unname(
        FOCUSED_LABELS[Characteristic_Raw]
      ),
      
      Scheme = factor(
        Scheme,
        levels = SCHEME_LEVELS
      )
    )
  
  comparison_wide <- means_long %>%
    select(
      Cohort,
      Characteristic_Raw,
      Characteristic,
      Scheme,
      Mean
    ) %>%
    pivot_wider(
      names_from = Scheme,
      values_from = Mean
    ) %>%
    mutate(
      `IPTW minus unweighted` =
        IPTW - Unweighted,
      
      `Overlap minus unweighted` =
        `Overlap weighting` - Unweighted,
      
      `Overlap minus IPTW` =
        `Overlap weighting` - IPTW,
      
      `Overlap vs unweighted, % difference` =
        100 *
        (
          `Overlap weighting` - Unweighted
        ) /
        Unweighted
    )
  
  sample_size <- tibble(
    Cohort = cohort,
    
    Episodes = nrow(dat),
    
    Unweighted_ESS = nrow(dat),
    
    IPTW_ESS = effective_sample_size(
      dat[[iptw_column]]
    ),
    
    Overlap_ESS = effective_sample_size(
      dat[[overlap_column]]
    ),
    
    IPTW_weight_column = iptw_column,
    
    Overlap_weight_column = overlap_column
  )
  
  list(
    means_long = means_long,
    comparison_wide = comparison_wide,
    sample_size = sample_size
  )
}


## 6. SUMMARIZE THE FOUR COHORTS #####

summary_tr2_le30d <- summarize_overlap_population(
  ps_obj = ps_tr2_le30d,
  cohort = "3-4 days, within 30 days"
)

summary_tr2_gt30d <- summarize_overlap_population(
  ps_obj = ps_tr2_gt30d,
  cohort = "3-4 days, after 30 days"
)

summary_tr3_le30d <- summarize_overlap_population(
  ps_obj = ps_tr3_le30d,
  cohort = "5-14 days, within 30 days"
)

summary_tr3_gt30d <- summarize_overlap_population(
  ps_obj = ps_tr3_gt30d,
  cohort = "5-14 days, after 30 days"
)


## 7. COMBINE NUMERIC RESULTS #######

means_all <- bind_rows(
  summary_tr2_le30d$means_long,
  summary_tr2_gt30d$means_long,
  summary_tr3_le30d$means_long,
  summary_tr3_gt30d$means_long
)

comparison_all <- bind_rows(
  summary_tr2_le30d$comparison_wide,
  summary_tr2_gt30d$comparison_wide,
  summary_tr3_le30d$comparison_wide,
  summary_tr3_gt30d$comparison_wide
)

sample_size_all <- bind_rows(
  summary_tr2_le30d$sample_size,
  summary_tr2_gt30d$sample_size,
  summary_tr3_le30d$sample_size,
  summary_tr3_gt30d$sample_size
)


# Order cohorts to match the manuscript presentation.
COHORT_LEVELS <- c(
  "3-4 days, within 30 days",
  "3-4 days, after 30 days",
  "5-14 days, within 30 days",
  "5-14 days, after 30 days"
)

means_all <- means_all %>%
  mutate(
    Cohort = factor(
      Cohort,
      levels = COHORT_LEVELS
    ),
    
    Characteristic = factor(
      Characteristic,
      levels = unname(
        FOCUSED_LABELS[FOCUSED_VARIABLES]
      )
    )
  )

comparison_all <- comparison_all %>%
  mutate(
    Cohort = factor(
      Cohort,
      levels = COHORT_LEVELS
    )
  )


## 8. CREATE ONE FOCUSED FIGURE ########

p_overlap_characterization <- ggplot(
  means_all,
  aes(
    x = Mean,
    y = Scheme,
    color = Scheme,
    shape = Scheme
  )
) +
  
  # Connect the three population summaries within each cohort/characteristic.
  geom_line(
    aes(
      group = 1
    ),
    color = "grey65",
    linewidth = 0.65
  ) +
  
  geom_point(
    size = 4.5,
    stroke = 0.9
  ) +
  
  ggh4x::facet_grid2(
    Characteristic ~ Cohort,
    scales = "free_x",
    independent = "x",
    switch = "y"
  ) +
  
  scale_color_manual(
    values = SCHEME_COLORS,
    breaks = SCHEME_LEVELS,
    drop = FALSE
  ) +
  
  scale_shape_manual(
    values = SCHEME_SHAPES,
    breaks = SCHEME_LEVELS,
    drop = FALSE
  ) +
  
  labs(
    title =
      "Baseline characteristics of the original, IPTW, and overlap-weighted populations",
    
    x =
      "Mean",
    
    y = NULL,
    
    color = NULL,
    
    shape = NULL,
    
    caption = paste0(
      "The overlap-weighted population represents episodes receiving greater ",
      "weight when multiple reinitiation strategies were clinically plausible."
    )
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  theme(
    legend.position = "bottom",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12,
      color = "black"
    ),
    
    legend.key.width = grid::unit(
      1.2,
      "cm"
    ),
    
    plot.title = element_text(
      size = 17,
      face = "bold",
      color = "black",
      hjust = 0
    ),
    
    strip.background = element_rect(
      fill = "grey94",
      color = "grey50"
    ),
    
    strip.text.x = element_text(
      size = 11.5,
      face = "bold",
      color = "black"
    ),
    
    strip.text.y.left = element_text(
      size = 11.5,
      face = "bold",
      angle = 0,
      hjust = 1,
      color = "black"
    ),
    
    strip.placement = "outside",
    
    axis.text.x = element_text(
      size = 10.5,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 12.5,
      color = "black",
      margin = ggplot2::margin(
        t = 8
      )
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.y = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    
    panel.spacing.x = grid::unit(
      0.5,
      "lines"
    ),
    
    panel.spacing.y = grid::unit(
      0.5,
      "lines"
    ),
    
    plot.caption = element_text(
      hjust = 0,
      size = 11,
      color = "black",
      margin = ggplot2::margin(
        t = 12
      )
    ),
    
    plot.margin = ggplot2::margin(
      t = 15,
      r = 15,
      b = 12,
      l = 15,
      unit = "pt"
    )
  )


## 9. SAVE THE SINGLE FOCUSED FIGURE #####

FIG_FILE <- paste0(
  FIG_DIR,
  "OverlapPopulation_Focused_Characterization.jpg"
)

jpeg(
  filename = FIG_FILE,
  width = 5200,
  height = 2200,
  res = 300,
  quality = 100,
  bg = "white"
)

print(
  p_overlap_characterization
)

dev.off()


## 10. SAVE THE SMALL NUMERIC TABLE #########

comparison_for_export <- comparison_all %>%
  mutate(
    across(
      c(
        Unweighted,
        IPTW,
        `Overlap weighting`,
        `IPTW minus unweighted`,
        `Overlap minus unweighted`,
        `Overlap minus IPTW`,
        `Overlap vs unweighted, % difference`
      ),
      ~ round(.x, 1)
    )
  ) %>%
  arrange(
    Cohort,
    Characteristic_Raw
  )

CSV_FILE <- paste0(
  TABLE_DIR,
  "OverlapPopulation_Focused_Characterization.csv"
)

write_csv(
  comparison_for_export,
  CSV_FILE
)

write_csv(
  sample_size_all,
  paste0(
    TABLE_DIR,
    "OverlapPopulation_EffectiveSampleSize.csv"
  )
)


## 11. PRINT THE EXACT VALUES NEEDED FOR THE MANUSCRIPT ######

cat(
  "\n=== Focused overlap-population comparison ===\n"
)

print(
  comparison_for_export %>%
    select(
      Cohort,
      Characteristic,
      Unweighted,
      IPTW,
      `Overlap weighting`,
      `Overlap minus unweighted`,
      `Overlap vs unweighted, % difference`
    ),
  n = Inf
)


cat(
  "\n=== Effective sample sizes ===\n"
)

print(
  sample_size_all,
  n = Inf
)


cat(
  "\nFigure saved:\n",
  FIG_FILE,
  "\n\n",
  sep = ""
)

cat(
  "Numeric table saved:\n",
  CSV_FILE,
  "\n",
  sep = ""
)