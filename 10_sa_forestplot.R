library(tidyverse)
library(forestploter)
library(grid)

# 1. SETUP & MASTER MANIFEST ---------------------------------------------------
RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

# Master Manifest
# Added "Overlap weights based analysis" with obj_type "ow_all"
# Master Manifest (FIXED COLUMN NAMES)
sa_manifest <- tibble::tribble(
  ~category, ~sa_label, ~obj_type, ~risk_file, ~ps_file, # <--- Fixed names here
  
  "0. Primary Analysis", "Primary analysis", "primary_all", "Risk/risk_tr2_34d_rmst.rds", "PScore/ps_weight_tr2_34d_all.rds",
  
  "1. Sample restriction", "First OAT episode", "primary_first", "Risk/risk_tr2_34d_rmst.rds", "PScore/ps_weight_tr2_34d_first.rds",
  "1. Sample restriction", "Completed induction", "standard", "Sensitivities/Risk/risk_tr2_34d_SA02_Induction.rds", "Sensitivities/ps_wt_tr2_34d_SA02_Induction.rds",
  "1. Sample restriction", "Missed at <80mg of methadone", "standard", "Sensitivities/Risk/risk_tr2_34d_SA03_Dose_LT80.rds", "Sensitivities/ps_wt_tr2_34d_SA03_Dose_LT80.rds",
  "1. Sample restriction", "Missed at >=80mg of methadone", "standard", "Sensitivities/Risk/risk_tr2_34d_SA04_Dose_GE80.rds", "Sensitivities/ps_wt_tr2_34d_SA04_Dose_GE80.rds",
  "1. Sample restriction", "Days since OAT initiation: \u226430 days", "standard", "Sensitivities/Risk/risk_tr2_34d_SA05_Ret_LE30.rds", "Sensitivities/ps_wt_tr2_34d_SA05_Ret_LE30.rds",
  "1. Sample restriction", "Days since OAT initiation: >30 days", "standard", "Sensitivities/Risk/risk_tr2_34d_SA06_Ret_GT30.rds", "Sensitivities/ps_wt_tr2_34d_SA06_Ret_GT30.rds",
  "1. Sample restriction", "Restricted to high fentanyl area", "standard", "Sensitivities/Risk/risk_tr2_34d_SA09_HighFentArea.rds", "Sensitivities/ps_wt_tr2_34d_SA09_HighFentArea.rds",
  "1. Sample restriction", "1-year or longer cumulative OAT experience", "standard", "Sensitivities/Risk/risk_tr2_34d_SA11_CumOAT_GE1yr.rds", "Sensitivities/ps_wt_tr2_34d_SA11_CumOAT_GE1yr.rds",
  
  "2. Timeline restriction", "Restricted to post fentanyl era", "standard", "Sensitivities/Risk/risk_tr2_34d_SA10_PostFentEra.rds", "Sensitivities/ps_wt_tr2_34d_SA10_PostFentEra.rds",
  
  "3. Exposure classification", "Alternative time zero: all interruptions", "rprev", "Risk/risk_tr2_34d_rprev_ow.rds", "PScore/ps_weight_tr2_34d_rprev_ow.rds",
  
  "4. Outcome definition", "30 days washout for prevalent users", "standard", "Sensitivities/Risk/risk_tr2_34d_SA13_30dWashout.rds", "Sensitivities/ps_wt_tr2_34d_SA13_30dWashout.rds",
  "4. Outcome definition", "No truncation at 12 months", "standard", "Sensitivities/Risk/risk_tr2_34d_SA12_NoTruncation.rds", "Sensitivities/ps_wt_tr2_34d_SA12_NoTruncation.rds",
  "4. Outcome definition", "Composite outcome: Overdose or death", "standard", "Sensitivities/Risk/risk_tr2_34d_SA16_Composite_OD.rds", "Sensitivities/ps_wt_tr2_34d_SA16_Composite_OD.rds",
  
  "5. Model specification", "High-dimensional propensity score (hdPS)", "hdps", "Risk/risk_hdPS_34d.rds", "PScore/ps_weight_34d_hdPS.rds",
  "5. Model specification", "hdPS with machine learning algorithm", "hdps", "Risk/risk_hdPS_ML_34d.rds", "PScore/ps_weight_34d_hdPS_ML.rds",
  "5. Model specification", "Overlap weights based analysis", "ow_all", "Risk/OW/risk_tr2_34d.rds", "PScore/ps_weight_tr2_34d_all.rds"
)

# Added 'no reversion' to track the reference group's N
trt_levels <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_clean <- c("no reversion" = "No dose change", "1-<50%" = "<50% reduction", 
               ">=50%/<=40mg" = "\u226550% reduction/\u226440mg", "increase" = "Dose increase")

# Define exact colors for the graph
trt_colors <- c("No dose change" = "black", 
                "<50% reduction" = "#2CA02C", 
                "\u226550% reduction/\u226440mg" = "#FF7F0E", 
                "Dose increase" = "#D62728")

# 2. SMART EXTRACTION HELPER ---------------------------------------------------
# 2. SMART EXTRACTION HELPER ---------------------------------------------------
extract_nested_data <- function(risk_file, ps_file, obj_type, outcome, max_time) {
  
  full_risk <- paste0(RES_DIR, risk_file)
  full_ps   <- paste0(RES_DIR, ps_file)
  
  na_block <- tibble(trt = trt_levels, Episode = NA_character_, rd_pct = NA_real_, lci_pct = NA_real_, uci_pct = NA_real_, rd_disp = "")
  
  if (!file.exists(full_risk) || !file.exists(full_ps)) {
    cat("  [MISSING] Data not found for:", risk_file, "\n")
    return(na_block)
  }
  
  risk_obj <- readRDS(full_risk)
  ps_obj   <- readRDS(full_ps)
  
  # --- FIX: Dynamic Variable Mapping for SA12 and SA16 ---
  if (grepl("SA12_NoTruncation", risk_file)) {
    time_v  <- if(outcome == "dis") "t_discont" else "t_od_acute"
    event_v <- if(outcome == "dis") "s_discont" else "s_od_acute"
  } else if (grepl("SA16_Composite_OD", risk_file) && outcome == "od") {
    time_v  <- "time_od_comp_outcome"
    event_v <- "sens_od_comp_outcome"
  } else {
    time_v  <- if(outcome == "dis") "time_discont" else "time_od_acute"
    event_v <- if(outcome == "dis") "sens_discont" else "sens_od_acute"
  }
  # -------------------------------------------------------
  
  # 1. Episode Counts (Dynamically use ow_wt for OW analysis)
  dat_w <- if("data_weights" %in% names(ps_obj)) ps_obj$data_weights else ps_obj
  weight_col <- if(obj_type == "ow_all") "ow_wt" else "wt_trunc"
  
  counts <- tryCatch({
    compute_summary_stats_timewin(
      dat_w, time_var = time_v, event_var = event_v, 
      weight_var = weight_col, trt_var = "trt", max_time = max_time
    ) %>% select(trt, Episode) %>% mutate(Episode = as.character(Episode))
  }, error = function(e) tibble(trt = trt_levels, Episode = NA_character_))
  
  # 2. Risk Differences
  contrasts <- tryCatch({
    if (obj_type %in% c("primary_all", "ow_all")) {
      risk_obj$all[[outcome]]$ci$contrasts
    } else if (obj_type == "primary_first") {
      risk_obj$first[[outcome]]$ci$contrasts
    } else if (obj_type %in% c("hdps", "rprev")) {
      risk_obj[[outcome]]$contrasts
    } else if (obj_type == "standard") {
      risk_obj[[outcome]]$ci$contrasts
    } else { NULL }
  }, error = function(e) NULL)
  
  if (is.null(contrasts) || nrow(contrasts) == 0) {
    cat("  [STRUCTURE ERROR/PENDING] No contrasts in:", risk_file, "\n")
    return(na_block)
  }
  
  rd_data <- contrasts %>%
    filter(horizon == max_time) %>%
    mutate(
      rd_pct = RD * 100, lci_pct = RD_lower * 100, uci_pct = RD_upper * 100,
      rd_disp = sprintf("%.2f (%.2f, %.2f)", rd_pct, lci_pct, uci_pct)
    ) %>% select(trt, rd_pct, lci_pct, uci_pct, rd_disp)
  
  cat("  [SUCCESS] Extracted:", risk_file, "\n")
  
  # Merge, ensuring "no reversion" gets "Ref." assigned
  tibble(trt = trt_levels) %>% 
    left_join(counts, by = "trt") %>% 
    left_join(rd_data, by = "trt") %>%
    mutate(rd_disp = ifelse(trt == "no reversion", "Ref.", rd_disp))
}

# 3. BUILD WIDE DATASET --------------------------------------------------------
# UPDATED HELPER: Now accepts the specific manifest as an argument
build_wide_df <- function(manifest, outcome = "od", max_time = 30) {
  final_rows <- list()
  
  for (i in 1:nrow(manifest)) {
    row_meta <- manifest[i, ]
    trt_data <- extract_nested_data(row_meta$risk_file, row_meta$ps_file, row_meta$obj_type, outcome, max_time)
    
    trt_formatted <- trt_data %>%
      mutate(
        Category = row_meta$category,
        SA_Label = row_meta$sa_label,
        Dose_Adj = recode(trt, !!!trt_clean),
        Episode = replace_na(Episode, ""),
        rd_disp = replace_na(rd_disp, "")
      )
    final_rows[[length(final_rows) + 1]] <- trt_formatted
  }
  
  df_final <- bind_rows(final_rows) %>%
    group_by(Category) %>%
    mutate(Category_Col = ifelse(row_number() == 1, Category, "")) %>%
    group_by(Category, SA_Label) %>%
    mutate(SA_Col = ifelse(row_number() == 1, SA_Label, "")) %>%
    ungroup() %>%
    mutate(` ` = paste(rep(" ", 18), collapse = " ")) %>%
    select(Category_Col, SA_Col, Dose_Adj, Episode, rd_disp, rd_pct, lci_pct, uci_pct, ` `) %>%
    rename(
      `Category` = Category_Col,
      `Sensitivity Analysis` = SA_Col,
      `Dose adjustment` = Dose_Adj,
      `Episodes\nat risk` = Episode,
      `Risk Difference (%)\n(95% CI)` = rd_disp
    )
  return(df_final)
}

# 4. PLOT GENERATION -----------------------------------------------------------
cat("\n--- Extracting Data ---\n")

# FIX: Explicitly passing the Trial 2 manifest to the function
plot_data <- build_wide_df(manifest = sa_manifest, outcome = "od", max_time = 30)

# Create an alternating background color pattern
num_sas <- nrow(sa_manifest)
bg_pattern <- rep(c("#FFFFFF", "#F5F5F5"), length.out = num_sas)
bg_fill_colors <- rep(bg_pattern, each = 4) 

# Dynamically map the colors for the CI plot
ci_col_mapped <- trt_colors[plot_data$`Dose adjustment`]
ci_col_mapped[is.na(ci_col_mapped)] <- "black" # Fallback safety

cat("\n--- Generating Base Plot ---\n")
p <- forest(
  plot_data[, c("Category", "Sensitivity Analysis", "Dose adjustment", "Episodes\nat risk", "Risk Difference (%)\n(95% CI)", " ")],
  est = plot_data$rd_pct,
  lower = plot_data$lci_pct,
  upper = plot_data$uci_pct,
  ci_column = 6, 
  ref_line = 0,
  xlim = c(-2.5, 2.5),
  ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk Difference (%)",
  title = "Methadone Interruption 3-4 Days: Adjusted Overdose Risk Differences at 30 Days\n",
  theme = forest_theme(
    base_size = 12,
    ci_pch = 15,
    ci_col = "black",  
    ci_fill = "black", 
    ci_lwd = 1.8,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "black"),
    core = list(
      bg_params = list(fill = bg_fill_colors)
    )
  )
)

cat("--- Applying Dynamic Row Colors ---\n")
for(i in 1:nrow(plot_data)) {
  row_color <- ci_col_mapped[i]
  
  if(!is.na(row_color) && row_color != "black") {
    p <- edit_plot(p, 
                   row = i, 
                   col = 6, 
                   which = "ci", 
                   gp = gpar(col = row_color, fill = row_color))
  }
}


# 5. SAVE ----------------------------------------------------------------------
jpeg(paste0(RES_DIR, "Figures/June09/SA_Forest_Wide_TR2_OD.jpg"), width = 3800, height = 5500, res = 300)
p
dev.off()

cat("\n--- Plot Saved! ---\n")