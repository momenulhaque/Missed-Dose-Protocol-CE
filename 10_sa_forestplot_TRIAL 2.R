# Trial 2 ##############################

library(tidyverse)
library(forestploter)
library(grid)

## 1. SETUP & MASTER MANIFEST ---------------------------------------------------
RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

# Master Manifest (Updated Labels)
sa_manifest <- tibble::tribble(
  ~category, ~sa_label, ~obj_type, ~risk_file, ~ps_file,
  
  "0. Primary Analysis", "No prior interruption within the episode", "primary_all", "Risk/risk_tr2_34d_rmst.rds", "PScore/ps_weight_tr2_34d_all.rds",
  
  "1. Sample restriction", "No prior interruption", "primary_first", "Risk/risk_tr2_34d_rmst.rds", "PScore/ps_weight_tr2_34d_first.rds",
  "1. Sample restriction", "No interruption in past 30 days within the episode", "standard", "Sensitivities/Risk/risk_tr2_34d_sa2.rds", "Sensitivities/ps_weight_tr2_34d_sa2.rds",
    "1. Sample restriction", "Completed induction", "standard", "Sensitivities/Risk/risk_tr2_34d_SA02_Induction.rds", "Sensitivities/ps_wt_tr2_34d_SA02_Induction.rds",
  "1. Sample restriction", "Missed at <80mg of methadone", "standard", "Sensitivities/Risk/risk_tr2_34d_SA03_Dose_LT80.rds", "Sensitivities/ps_wt_tr2_34d_SA03_Dose_LT80.rds",
  "1. Sample restriction", "Missed at >=80mg of methadone", "standard", "Sensitivities/Risk/risk_tr2_34d_SA04_Dose_GE80.rds", "Sensitivities/ps_wt_tr2_34d_SA04_Dose_GE80.rds",
  "1. Sample restriction", "Days since OAT initiation: \u226430 days", "standard", "Sensitivities/Risk/risk_tr2_34d_SA05_Ret_LE30.rds", "Sensitivities/ps_wt_tr2_34d_SA05_Ret_LE30.rds",
  "1. Sample restriction", "Days since OAT initiation: >30 days", "standard", "Sensitivities/Risk/risk_tr2_34d_SA06_Ret_GT30.rds", "Sensitivities/ps_wt_tr2_34d_SA06_Ret_GT30.rds",
  "1. Sample restriction", "Restricted to high fentanyl area", "standard", "Sensitivities/Risk/risk_tr2_34d_SA09_HighFentArea.rds", "Sensitivities/ps_wt_tr2_34d_SA09_HighFentArea.rds",
  "1. Sample restriction", "1-year or longer cumulative OAT experience", "standard", "Sensitivities/Risk/risk_tr2_34d_SA11_CumOAT_GE1yr.rds", "Sensitivities/ps_wt_tr2_34d_SA11_CumOAT_GE1yr.rds",
  
  "2. Timeline restriction", "Restricted to post-fentanyl era (\u2265 April 1, 2012)", "standard", "Sensitivities/Risk/risk_tr2_34d_SA10_PostFentEra.rds", "Sensitivities/ps_wt_tr2_34d_SA10_PostFentEra.rds",
  
  "3. Exposure classification", "5 category dose adjustment", "standard", "Sensitivities/risk_tr2_34d_sa2.rds", "Sensitivities/ps_weight_tr2_34d_sa2.rds",
  
  "4. Outcome definition", "30 days washout for prevalent users", "standard", "Sensitivities/Risk/risk_tr2_34d_SA13_30dWashout.rds", "Sensitivities/ps_wt_tr2_34d_SA13_30dWashout.rds",
  "4. Outcome definition", "No truncation at 12 months", "standard", "Sensitivities/Risk/risk_tr2_34d_SA12_NoTruncation.rds", "Sensitivities/ps_wt_tr2_34d_SA12_NoTruncation.rds",
  "4. Outcome definition", "Composite outcome: Overdose or death", "standard", "Sensitivities/Risk/risk_tr2_34d_SA16_Composite_OD.rds", "Sensitivities/ps_wt_tr2_34d_SA16_Composite_OD.rds",
  
  "5. Model specification", "High-dimensional propensity score (hdPS)", "hdps", "Risk/risk_hdPS_34d.rds", "PScore/ps_weight_34d_hdPS.rds",
  "5. Model specification", "hdPS with machine learning algorithm", "hdps", "Risk/risk_hdPS_ML_34d.rds", "PScore/ps_weight_34d_hdPS_ML.rds",
  "5. Model specification", "Overlap weights based analysis", "ow_all", "Risk/OW/risk_tr2_34d.rds", "PScore/ps_weight_tr2_34d_all.rds"
)

trt_levels <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_clean <- c("no reversion" = "No dose change", "1-<50%" = "<50% reduction", 
               ">=50%/<=40mg" = "\u226550% reduction/\u226440mg", "increase" = "Dose increase")

trt_colors <- c("No dose change" = "black", 
                "<50% reduction" = "#2CA02C", 
                "\u226550% reduction/\u226440mg" = "#FF7F0E", 
                "Dose increase" = "#D62728")

## 2. SMART EXTRACTION HELPER ---------------------------------------------------
extract_nested_data <- function(risk_file, ps_file, obj_type, outcome, max_time) {
  
  full_risk <- paste0(RES_DIR, risk_file)
  full_ps   <- paste0(RES_DIR, ps_file)
  
  na_block <- tibble(trt = trt_levels, Episode = NA_character_, abs_risk_disp = "", 
                     rd_pct = NA_real_, lci_pct = NA_real_, uci_pct = NA_real_, rd_disp = "")
  
  if (!file.exists(full_risk) || !file.exists(full_ps)) return(na_block)
  
  risk_obj <- readRDS(full_risk)
  ps_obj   <- readRDS(full_ps)
  
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
  
  dat_w <- if("data_weights" %in% names(ps_obj)) ps_obj$data_weights else ps_obj
  
  # 1. Episode Counts (FIX: Explicitly convert trt factor to character to prevent join failures)
  counts <- tryCatch({
    dat_w %>% 
      count(trt, name = "Episode") %>% 
      mutate(trt = as.character(trt), Episode = as.character(Episode))
  }, error = function(e) tibble(trt = trt_levels, Episode = NA_character_))
  
  # 2. Risk Differences & Absolute Risks
  contrasts <- tryCatch({
    if (obj_type %in% c("primary_all", "ow_all")) risk_obj$all[[outcome]]$ci$contrasts
    else if (obj_type == "primary_first") risk_obj$first[[outcome]]$ci$contrasts
    else if (obj_type == "hdps") risk_obj[[outcome]]$contrasts
    else if (obj_type == "standard") risk_obj[[outcome]]$ci$contrasts
    else NULL 
  }, error = function(e) NULL)
  
  risks_df <- tryCatch({
    if (obj_type %in% c("primary_all", "ow_all")) risk_obj$all[[outcome]]$ci$risks
    else if (obj_type == "primary_first") risk_obj$first[[outcome]]$ci$risks
    else if (obj_type == "hdps") risk_obj[[outcome]]$risks
    else if (obj_type == "standard") risk_obj[[outcome]]$ci$risks
    else NULL 
  }, error = function(e) NULL)
  
  if (is.null(contrasts) || nrow(contrasts) == 0 || is.null(risks_df)) return(na_block)
  
  # Absolute Risk Extraction
  abs_data <- risks_df %>%
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
    ) %>% select(trt, abs_risk_disp)
  
  # Risk Difference Extraction
  rd_data <- contrasts %>%
    filter(horizon == max_time) %>%
    mutate(
      rd_pct = RD * 100, lci_pct = RD_lower * 100, uci_pct = RD_upper * 100,
      rd_disp = ifelse(
        abs(rd_pct) >= 10,
        sprintf("%.1f (%.1f, %.1f)", rd_pct, lci_pct, uci_pct),
        sprintf("%.2f (%.2f, %.2f)", rd_pct, lci_pct, uci_pct)
      )
    ) %>% select(trt, rd_pct, lci_pct, uci_pct, rd_disp)
  
  tibble(trt = trt_levels) %>% 
    left_join(counts, by = "trt") %>% 
    left_join(abs_data, by = "trt") %>%
    left_join(rd_data, by = "trt") %>%
    mutate(rd_disp = ifelse(trt == "no reversion", "Ref.", rd_disp))
}

## 3. BUILD WIDE DATASET --------------------------------------------------------
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
        abs_risk_disp = replace_na(abs_risk_disp, ""),
        rd_disp = replace_na(rd_disp, "")
      )
    final_rows[[length(final_rows) + 1]] <- trt_formatted
  }
  
  time_str <- if(max_time == 30) "1 month" else "12 months"
  risk_col_name <- paste0("% risk at ", time_str, "\n(95% CI)")
  
  df_final <- bind_rows(final_rows) %>%
    group_by(Category) %>%
    mutate(Category_Col = ifelse(row_number() == 1, Category, "")) %>%
    group_by(Category, SA_Label) %>%
    mutate(SA_Col = {
      # FIX: Detect "\n" and split the text elegantly across the 4 available rows
      parts <- unlist(strsplit(as.character(SA_Label[1]), "\n"))
      out_col <- rep("", n())
      max_fill <- min(length(parts), n()) # Fill up to 4 rows
      out_col[1:max_fill] <- parts[1:max_fill]
      out_col
    }) %>%
    ungroup() %>%
    mutate(` ` = paste(rep(" ", 18), collapse = " ")) %>%
    select(Category_Col, SA_Col, Dose_Adj, Episode, abs_risk_disp, rd_disp, rd_pct, lci_pct, uci_pct, ` `)
  
  names(df_final)[names(df_final) == "Category_Col"] <- "Category"
  names(df_final)[names(df_final) == "SA_Col"] <- "Sensitivity Analysis"
  names(df_final)[names(df_final) == "Dose_Adj"] <- "Dose adjustment"
  names(df_final)[names(df_final) == "Episode"] <- "Episodes at risk\nat reinitiation"
  names(df_final)[names(df_final) == "abs_risk_disp"] <- risk_col_name
  names(df_final)[names(df_final) == "rd_disp"] <- "Risk Difference (%)\n(95% CI)"
  
  return(df_final)
}

## OVERDOSE #########
### 4. PLOT GENERATION -----------------------------------------------------------
plot_data <- build_wide_df(manifest = sa_manifest, outcome = "od", max_time = 30)

num_sas <- nrow(sa_manifest)
bg_pattern <- rep(c("#FFFFFF", "#F5F5F5"), length.out = num_sas)
bg_fill_colors <- rep(bg_pattern, each = 4) 

ci_col_mapped <- trt_colors[plot_data$`Dose adjustment`]
ci_col_mapped[is.na(ci_col_mapped)] <- "black"

risk_col <- "% risk at 1 month\n(95% CI)"

p <- forest(
  plot_data[, c("Category", "Sensitivity Analysis", "Dose adjustment", "Episodes at risk\nat reinitiation", risk_col, "Risk Difference (%)\n(95% CI)", " ")],
  est = plot_data$rd_pct,
  lower = plot_data$lci_pct,
  upper = plot_data$uci_pct,
  ci_column = 7,
  ref_line = 0,
  xlim = c(-2.5, 2.5),
  ticks_at = c(-2, -1, 0, 1, 2),
  xlab = "Risk Difference (%)",
  title = "Methadone Interruption 3-4 Days: Adjusted Risk Differences at 1 month\n",
  theme = forest_theme(
    base_size = 12,
    ci_pch = 15,
    ci_col = "black",  
    ci_fill = "black", 
    ci_lwd = 1.8,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "black"),
    core = list(bg_params = list(fill = bg_fill_colors)) 
    # Removed the padding override, letting the natural grid handle it seamlessly
  )
)

for(i in 1:nrow(plot_data)) {
  row_color <- ci_col_mapped[i]
  if(!is.na(row_color) && row_color != "black") {
    p <- edit_plot(p, row = i, col = 7, which = "ci", gp = gpar(col = row_color, fill = row_color))
  }
}


### 5. SAVE ----------------------------------------------------------------------
jpeg(paste0(RES_DIR, "Figures/June09/SA_Forest_Wide_TR2_OD.jpg"), width = 4400, height = 5500, res = 300)
p
dev.off()



## DISCONTINUATION #########
### 6. PLOT GENERATION -----------------------------------------------------------
plot_data <- build_wide_df(manifest = sa_manifest, outcome = "dis", max_time = 365)

num_sas <- nrow(sa_manifest)
bg_pattern <- rep(c("#FFFFFF", "#F5F5F5"), length.out = num_sas)
bg_fill_colors <- rep(bg_pattern, each = 4) 

ci_col_mapped <- trt_colors[plot_data$`Dose adjustment`]
ci_col_mapped[is.na(ci_col_mapped)] <- "black"

risk_col <- "% risk at 12 months\n(95% CI)"

p <- forest(
  plot_data[, c("Category", "Sensitivity Analysis", "Dose adjustment", "Episodes at risk\nat reinitiation", risk_col, "Risk Difference (%)\n(95% CI)", " ")],
  est = plot_data$rd_pct,
  lower = plot_data$lci_pct,
  upper = plot_data$uci_pct,
  ci_column = 7,
  ref_line = 0,
  xlim = c(-15, 15),
  ticks_at = c(-10, -5, 0, 5, 10),
  xlab = "Risk Difference (%)",
  title = "Methadone Interruption 3-4 Days: Adjusted Risk Differences at 12 months\n",
  theme = forest_theme(
    base_size = 12,
    ci_pch = 15,
    ci_col = "black",  
    ci_fill = "black", 
    ci_lwd = 1.8,
    refline_gp = gpar(lwd = 1, lty = "dashed", col = "black"),
    core = list(bg_params = list(fill = bg_fill_colors)) 
    # Removed the padding override, letting the natural grid handle it seamlessly
  )
)

for(i in 1:nrow(plot_data)) {
  row_color <- ci_col_mapped[i]
  if(!is.na(row_color) && row_color != "black") {
    p <- edit_plot(p, row = i, col = 7, which = "ci", gp = gpar(col = row_color, fill = row_color))
  }
}


### 7. SAVE ----------------------------------------------------------------------
jpeg(paste0(RES_DIR, "Figures/June09/SA_Forest_Wide_TR2_DIS.jpg"), width = 4400, height = 5500, res = 300)
p
dev.off()