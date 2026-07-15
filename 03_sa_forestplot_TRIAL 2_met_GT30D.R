library(tidyverse)
library(forestploter)
library(grid)

# 1. SETUP & MASTER MANIFEST ---------------------------------------------------
RES_DIR <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/GT30D/"

sa_manifest_od <- tibble::tribble(
  ~category, ~sa_label, ~obj_type, ~risk_file, ~ps_file,
  "0. Primary Analysis", "First interruption in\nany episode", "standard", "Risk/risk_tr2.rds", "PScore/ps_weight_tr2_34d_all.rds",
  "1. Sample restriction", "First interruption in\nfirst episode", "standard", "Risk/risk_tr2_first.rds", "PScore/ps_weight_tr2_34d_first.rds",
  "1. Sample restriction", "Missed at <80mg\nof methadone", "standard", "Sensitivities/Risk/risk_tr2_34d_SA03_Dose_LT80.rds", "Sensitivities/ps_wt_tr2_34d_SA03_Dose_LT80.rds",
  "1. Sample restriction", "Missed at \u226580mg\nof methadone", "standard", "Sensitivities/Risk/risk_tr2_34d_SA04_Dose_GE80.rds", "Sensitivities/ps_wt_tr2_34d_SA04_Dose_GE80.rds",
  "1. Sample restriction", "Restricted to high\nfentanyl area", "standard", "Sensitivities/Risk/risk_tr2_34d_SA09_HighFentArea.rds", "Sensitivities/ps_wt_tr2_34d_SA09_HighFentArea.rds",
  "1. Sample restriction", "1-year or longer\ncumulative OAT", "standard", "Sensitivities/Risk/risk_tr2_34d_SA11_CumOAT_GE1yr.rds", "Sensitivities/ps_wt_tr2_34d_SA11_CumOAT_GE1yr.rds",
  "1. Sample restriction", "30 days of no OAT\ndispensation before initiation", "standard", "Sensitivities/Risk/risk_tr2_34d_SA13_30dWashout.rds", "Sensitivities/ps_wt_tr2_34d_SA13_30dWashout.rds",
  "1. Sample restriction", "Alternative trial:\nInterruption of 4 days", "standard", "Sensitivities/Risk/risk_tr2_4d_rmst.rds", "Sensitivities/ps_weight_tr2_4d_all.rds",
  
  "2. Timeline restriction", "Restricted to post-fentanyl\nera (\u2265 04/01/2012)", "standard", "Sensitivities/Risk/risk_tr2_34d_SA10_PostFentEra.rds", "Sensitivities/ps_wt_tr2_34d_SA10_PostFentEra.rds",
  
  "3. Exposure classification", "5 category dose\nadjustment", "standard", "Sensitivities/Risk/risk_tr2_34d_5cat_od.rds", "Sensitivities/ps_weight_tr2_34d_5cat.rds",
  
  "4. Outcome definition", "Overdose-related acute\ncare visits or death", "standard", "Sensitivities/Risk/risk_tr2_34d_SA16_Composite_OD.rds", "Sensitivities/ps_wt_tr2_34d_SA16_Composite_OD.rds",
  "4. Outcome definition", "All-cause acute\ncare visits", "standard", "Sensitivities/Risk/risk_tr2_34d_all_acute.rds", "PScore/ps_weight_tr2_34d_all.rds",
  
  "5. Causal contrast", "Censoring at subsequent\ninterruption of >2 days", "standard", "Risk/risk_tr2_per_protocol.rds", "PScore/ps_weight_34d_hdPS.rds",
  
  "6. Model specification", "High-dimensional propensity\nscore (hdPS)", "standard", "Risk/risk_tr2_hdPS.rds", "PScore/ps_weight_34d_hdPS.rds",
  "6. Model specification", "hdPS with machine\nlearning algorithm", "standard", "Risk/risk_tr2_hdPS_ML.rds", "PScore/ps_weight_34d_hdPS_ML.rds",
  "6. Model specification", "Overlap weights\nbased analysis", "standard", "Risk/risk_tr2_ow.rds", "PScore/ps_weight_tr2_34d_all.rds",
  "6. Model specification", "PS model with additional\nprovider-level covariates", "standard", "Risk/risk_tr2_cl.rds", "PScore/ps_weight_tr2_34d_all_cl.rds"
)

sa_manifest_dis <- tibble::tribble(
  ~category, ~sa_label, ~obj_type, ~risk_file, ~ps_file,
  "0. Primary Analysis", "First interruption in\nany episode", "standard", "Risk/risk_tr2.rds", "PScore/ps_weight_tr2_34d_all.rds",
  "1. Sample restriction", "First interruption in\nfirst episode", "standard", "Risk/risk_tr2_first.rds", "PScore/ps_weight_tr2_34d_first.rds",
  "1. Sample restriction", "Missed at <80mg\nof methadone", "standard", "Sensitivities/Risk/risk_tr2_34d_SA03_Dose_LT80.rds", "Sensitivities/ps_wt_tr2_34d_SA03_Dose_LT80.rds",
  "1. Sample restriction", "Missed at \u226580mg\nof methadone", "standard", "Sensitivities/Risk/risk_tr2_34d_SA04_Dose_GE80.rds", "Sensitivities/ps_wt_tr2_34d_SA04_Dose_GE80.rds",
  "1. Sample restriction", "Restricted to high\nfentanyl area", "standard", "Sensitivities/Risk/risk_tr2_34d_SA09_HighFentArea.rds", "Sensitivities/ps_wt_tr2_34d_SA09_HighFentArea.rds",
  "1. Sample restriction", "1-year or longer\ncumulative OAT", "standard", "Sensitivities/Risk/risk_tr2_34d_SA11_CumOAT_GE1yr.rds", "Sensitivities/ps_wt_tr2_34d_SA11_CumOAT_GE1yr.rds",
  "1. Sample restriction", "30 days of no OAT\ndispensation before initiation", "standard", "Sensitivities/Risk/risk_tr2_34d_SA13_30dWashout.rds", "Sensitivities/ps_wt_tr2_34d_SA13_30dWashout.rds",
  "1. Sample restriction", "Alternative trial:\nInterruption of 4 days", "standard", "Sensitivities/Risk/risk_tr2_4d_rmst.rds", "Sensitivities/ps_weight_tr2_4d_all.rds",
  
  "2. Timeline restriction", "Restricted to post-fentanyl\nera (\u2265 04/01/2012)", "standard", "Sensitivities/Risk/risk_tr2_34d_SA10_PostFentEra.rds", "Sensitivities/ps_wt_tr2_34d_SA10_PostFentEra.rds",
  
  "3. Exposure classification", "5 category dose\nadjustment", "standard", "Sensitivities/Risk/risk_tr2_34d_5cat_dis.rds", "Sensitivities/ps_weight_tr2_34d_5cat.rds",
  
  "4. Outcome definition", "No dispensation of OAT\nfor 30 consecutive days", "standard", "Sensitivities/Risk/risk_tr2_34d_discont30.rds", "Sensitivities/ps_weight_tr2_34d_discont30.rds",
  
  "5. Causal contrast", "Censoring at subsequent\ninterruption of >2 days", "standard", "Risk/risk_tr2_per_protocol.rds", "PScore/ps_weight_34d_hdPS.rds",
  
  "6. Model specification", "High-dimensional propensity\nscore (hdPS)", "standard", "Risk/risk_tr2_hdPS.rds", "PScore/ps_weight_34d_hdPS.rds",
  "6. Model specification", "hdPS with machine\nlearning algorithm", "standard", "Risk/risk_tr2_hdPS_ML.rds", "PScore/ps_weight_34d_hdPS_ML.rds",
  "6. Model specification", "Overlap weights\nbased analysis", "standard", "Risk/risk_tr2_ow.rds", "PScore/ps_weight_tr2_34d_all.rds",
  "6. Model specification", "PS model with additional\nprovider-level covariates", "standard", "Risk/risk_tr2_cl.rds", "PScore/ps_weight_tr2_34d_all_cl.rds"
)

trt_clean <- c(
  "no reversion"      = "No dose change",
  "1-<50%"            = "<50% reduction",
  ">=50%/<=40mg"      = "\u226550% reduction/\u226440mg",
  "increase"          = "Dose increase",
  "increase <50%"     = "<50% increase",
  "increase >=50%"    = "\u226550% increase"
)

trt_colors <- c(
  "No dose change"                 = "black",
  "<50% reduction"                 = "#2CA02C",
  "\u226550% reduction/\u226440mg" = "#FF7F0E",
  "Dose increase"                  = "#D62728",
  "<50% increase"                  = "#9467BD",
  "\u226550% increase"              = "#A00000"
)

# 2A. EVENT-VARIABLE RESOLVER --------------------------------------------------
resolve_event_var <- function(risk_file, outcome) {
  if (grepl("SA12_NoTruncation", risk_file, fixed = TRUE)) {
    return(if (outcome == "od") "s_od_acute" else "s_discont")
  }
  if (grepl("SA16_Composite_OD", risk_file, fixed = TRUE)) {
    return(if (outcome == "od") "sens_od_comp_outcome" else "sens_discont")
  }
  if (grepl("discont30", risk_file, fixed = TRUE)) {
    return(if (outcome == "dis") "sens_discont30" else "sens_od_acute")
  }
  if (grepl("all_acute", risk_file, fixed = TRUE)) {
    return(if (outcome == "od") "sens_all_acute" else "sens_discont")
  }
  if (outcome == "od")  return("sens_od_acute")
  if (outcome == "dis") return("sens_discont")
  return(NA_character_)
}

# 2B. SMART EXTRACTION HELPER (WHOLE-SA MIN-EVENT SUPPRESSION) -----------------
extract_nested_data <- function(risk_file, ps_file, obj_type, outcome,
                                max_time, sa_label, min_events = 1) {
  
  full_risk <- paste0(RES_DIR, risk_file)
  full_ps   <- paste0(RES_DIR, ps_file)
  
  is_5cat <- grepl("5 category", sa_label)
  local_levels <- if (is_5cat) {
    c("no reversion", "1-<50%", ">=50%/<=40mg", "increase <50%", "increase >=50%")
  } else {
    c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
  }
  
  # --- STEP 1: EPISODE + EVENT COUNTS ---
  counts <- tibble(trt = local_levels,
                   Episode = NA_character_,
                   n_events = NA_integer_)
  
  event_var <- resolve_event_var(risk_file, outcome)
  
  if (!file.exists(full_ps)) {
    message(paste("WARNING: PS weights file not found:", full_ps))
  } else {
    ps_obj <- readRDS(full_ps)
    dat_w  <- if ("data_weights" %in% names(ps_obj)) ps_obj$data_weights else ps_obj
    
    if ("trt" %in% colnames(dat_w)) {
      real_counts <- dat_w %>%
        dplyr::group_by(trt) %>%
        dplyr::summarise(Episode = as.character(dplyr::n()), .groups = "drop") %>%
        dplyr::mutate(trt = as.character(trt))
      
      if (!is.na(event_var) && event_var %in% colnames(dat_w)) {
        real_events <- dat_w %>%
          dplyr::group_by(trt) %>%
          dplyr::summarise(
            n_events = sum(.data[[event_var]] == 1, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          dplyr::mutate(trt = as.character(trt))
      } else {
        message(paste("WARNING: event var", event_var,
                      "missing in", basename(full_ps)))
        real_events <- tibble(trt = character(), n_events = integer())
      }
      
      counts <- tibble(trt = local_levels) %>%
        dplyr::left_join(real_counts, by = "trt") %>%
        dplyr::left_join(real_events, by = "trt") %>%
        dplyr::mutate(
          Episode  = dplyr::coalesce(Episode, "0"),
          n_events = dplyr::coalesce(n_events, 0L)
        )
    } else {
      message(paste("WARNING: 'trt' column missing in", full_ps))
    }
  }
  
  na_block <- counts %>%
    mutate(
      abs_risk_disp = "",
      rd_pct = NA_real_, lci_pct = NA_real_, uci_pct = NA_real_,
      rd_disp = ""
    )
  
  # --- STEP 2: EXTRACT RISKS ---
  if (!file.exists(full_risk)) {
    message(paste("WARNING: Risk file not found:", full_risk))
    attr(na_block, "suppressed")   <- FALSE
    attr(na_block, "min_n_events") <- suppressWarnings(min(counts$n_events, na.rm = TRUE))
    return(na_block)
  }
  
  risk_obj  <- readRDS(full_risk)
  risks_df  <- NULL
  contrasts <- NULL
  
  tryCatch({
    if (obj_type %in% c("primary_all", "ow_all")) {
      out_node  <- if (!is.null(risk_obj$all[[outcome]])) outcome else names(risk_obj$all)[1]
      base_node <- risk_obj$all[[out_node]]
    } else if (obj_type == "primary_first") {
      out_node  <- if (!is.null(risk_obj$first[[outcome]])) outcome else names(risk_obj$first)[1]
      base_node <- risk_obj$first[[out_node]]
    } else {
      out_node  <- if (!is.null(risk_obj[[outcome]])) outcome else names(risk_obj)[1]
      base_node <- risk_obj[[out_node]]
    }
    node      <- if (!is.null(base_node$ci)) base_node$ci else base_node
    risks_df  <- node$risks
    contrasts <- node$contrasts
  }, error = function(e) { NULL })
  
  if (is.null(contrasts) || nrow(contrasts) == 0 || is.null(risks_df)) {
    attr(na_block, "suppressed")   <- FALSE
    attr(na_block, "min_n_events") <- suppressWarnings(min(counts$n_events, na.rm = TRUE))
    return(na_block)
  }
  
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
    ) %>% dplyr::select(trt, abs_risk_disp)
  
  rd_data <- contrasts %>%
    filter(horizon == max_time) %>%
    mutate(
      rd_pct  = RD * 100,
      lci_pct = RD_lower * 100,
      uci_pct = RD_upper * 100,
      rd_disp = ifelse(
        abs(rd_pct) >= 10,
        sprintf("%.1f (%.1f, %.1f)", rd_pct, lci_pct, uci_pct),
        sprintf("%.2f (%.2f, %.2f)", rd_pct, lci_pct, uci_pct)
      )
    ) %>% dplyr::select(trt, rd_pct, lci_pct, uci_pct, rd_disp)
  
  res <- tibble(trt = local_levels) %>%
    dplyr::left_join(counts,   by = "trt") %>%
    dplyr::left_join(abs_data, by = "trt") %>%
    dplyr::left_join(rd_data,  by = "trt") %>%
    dplyr::mutate(rd_disp = ifelse(trt == "no reversion", "Reference", rd_disp))
  
  # --- STEP 3: WHOLE-SA LOW-EVENT SUPPRESSION --------------------------------
  min_n_arm <- suppressWarnings(min(res$n_events, na.rm = TRUE))
  suppress  <- !is.finite(min_n_arm) || min_n_arm < min_events
  
  if (suppress) {
    res <- res %>%
      mutate(
        abs_risk_disp = "",
        rd_disp       = "",
        rd_pct        = NA_real_,
        lci_pct       = NA_real_,
        uci_pct       = NA_real_
      )
  }
  
  attr(res, "suppressed")   <- suppress
  attr(res, "min_n_events") <- min_n_arm
  return(res)
}

# 3. BUILD WIDE DATASET --------------------------------------------------------
build_wide_df <- function(manifest, outcome = "od", max_time = 30,
                          min_events = 1) {
  final_rows <- list()
  
  for (i in 1:nrow(manifest)) {
    row_meta <- manifest[i, ]
    trt_data <- extract_nested_data(
      row_meta$risk_file, row_meta$ps_file, row_meta$obj_type,
      outcome, max_time, row_meta$sa_label,
      min_events = min_events
    )
    suppressed_flag <- isTRUE(attr(trt_data, "suppressed"))
    min_n_events    <- attr(trt_data, "min_n_events")
    
    sa_label_display <- if (suppressed_flag) {
      paste0(row_meta$sa_label, "\u2020")
    } else {
      row_meta$sa_label
    }
    
    trt_formatted <- trt_data %>%
      mutate(
        Category = row_meta$category,
        SA_Label = sa_label_display,
        Dose_Adj = recode(trt, !!!trt_clean),
        Episode = replace_na(Episode, ""),
        abs_risk_disp = replace_na(abs_risk_disp, ""),
        rd_disp = replace_na(rd_disp, ""),
        bg_color = ifelse(i %% 2 == 1, "#FFFFFF", "#F5F5F5"),
        min_arm_events = min_n_events
      )
    final_rows[[length(final_rows) + 1]] <- trt_formatted
  }
  
  time_str <- if (max_time == 30) "1 month" else "12 months"
  risk_col_name <- paste0("% risk at ", time_str, "\n(95% CI)")
  
  df_final <- bind_rows(final_rows) %>%
    group_by(Category) %>%
    mutate(Category_Col = ifelse(row_number() == 1, Category, "")) %>%
    group_by(Category, SA_Label) %>%
    mutate(SA_Col = {
      parts <- unlist(strsplit(as.character(SA_Label[1]), "\n"))
      parts <- paste0("    ", parts)
      out_col <- rep("", n())
      max_fill <- min(length(parts), n())
      out_col[1:max_fill] <- parts[1:max_fill]
      out_col
    }) %>%
    ungroup() %>%
    mutate(`Plot_Canvas` = paste(rep(" ", 45), collapse = ""))
  
  categories <- unique(df_final$Category)
  df_with_headers <- bind_rows(lapply(categories, function(cat) {
    sub_df <- df_final %>% filter(Category == cat)
    
    header_row <- sub_df[1, ]
    header_row[] <- NA
    header_row$Category <- cat
    header_row$SA_Col <- cat
    header_row$Dose_Adj <- ""
    header_row$Episode <- ""
    header_row$abs_risk_disp <- ""
    header_row$rd_disp <- ""
    header_row$`Plot_Canvas` <- ""
    header_row$bg_color <- "#E0E0E0"
    header_row$is_header <- TRUE
    sub_df$is_header <- FALSE
    
    bind_rows(header_row, sub_df)
  }))
  
  df_with_headers <- df_with_headers %>%
    select(is_header, SA_Col, Dose_Adj, Episode,
           abs_risk_disp, rd_disp, rd_pct, lci_pct, uci_pct,
           `Plot_Canvas`, bg_color, min_arm_events)
  
  names(df_with_headers)[names(df_with_headers) == "SA_Col"]        <- "Sensitivity Analysis"
  names(df_with_headers)[names(df_with_headers) == "Dose_Adj"]      <- "Dose adjustment"
  names(df_with_headers)[names(df_with_headers) == "Episode"]       <- "Episodes at risk\nat reinitiation"
  names(df_with_headers)[names(df_with_headers) == "abs_risk_disp"] <- risk_col_name
  if(outcome == "od") names(df_with_headers)[names(df_with_headers) == "Plot_Canvas"] <- "% difference in 1-month risk (95% CI)"
  if(outcome == "dis") names(df_with_headers)[names(df_with_headers) == "Plot_Canvas"] <- "% difference in 12-month risk (95% CI)"
  names(df_with_headers)[names(df_with_headers) == "rd_disp"]   <- " "
  
  return(df_with_headers)
}

# 4A. PLOT GENERATION: OVERDOSE ------------------------------------------------
cat("\n--- Generating Overdose Plot (Trial 3) ---\n")
plot_data_od <- build_wide_df(
  manifest    = sa_manifest_od,
  outcome     = "od",
  max_time    = 30,
  min_events  = 1        # <-- tune: suppress SA if any arm has < this many OD events
)

bg_fill_od <- plot_data_od$bg_color
ci_col_od  <- trt_colors[plot_data_od$`Dose adjustment`]
ci_col_od[is.na(ci_col_od)] <- "black"

risk_col_od   <- "% risk at 1 month\n(95% CI)"
header_idx_od <- which(plot_data_od$is_header)

p_od <- forest(
  plot_data_od[, c("Sensitivity Analysis", "Dose adjustment",
                   "Episodes at risk\nat reinitiation", risk_col_od,
                   " ", "% difference in 1-month risk (95% CI)")],
  est = plot_data_od$rd_pct, lower = plot_data_od$lci_pct, upper = plot_data_od$uci_pct,
  ci_column = 6,
  vert_line = c(-3, -2, -1, 1, 2, 3),
  ref_line = 0, xlim = c(-3.5, 3.5), ticks_at = c(-3, -2, -1, 0, 1, 2, 3),
  xlab = "Risk difference (%)",
  theme = forest_theme(base_size = 12, ci_pch = 15, ci_col = "black", ci_fill = "black", ci_lwd = 1.8,
                       vertline_lty = "dotted",
                       vertline_col = "grey50",
                       vertline_lwd = 1,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
                       core = list(bg_params = list(fill = bg_fill_od),
                                   padding = unit(c(2, 2), "mm")))
)

p_od <- edit_plot(p_od, row = header_idx_od, gp = gpar(fontface = "bold"))
p_od <- edit_plot(p_od, col = 3, part = "body",   hjust = 0.5, x = unit(0.5, "npc"))
p_od <- edit_plot(p_od, col = 3, part = "header", hjust = 0.5, x = unit(0.5, "npc"))
p_od <- edit_plot(p_od, col = 4, part = "header", hjust = 0.5, x = unit(0.5, "npc"))
p_od <- edit_plot(p_od, col = 6, part = "body", hjust = 0, x = unit(.2, "npc"))
p_od <- edit_plot(p_od, col = 6, part = "header", hjust = 0.5, x = unit(0.15, "npc"))



for (i in 1:nrow(plot_data_od)) {
  if (!plot_data_od$is_header[i] && !is.na(ci_col_od[i]) && ci_col_od[i] != "black") {
    p_od <- edit_plot(p_od, row = i, col = 6, which = "ci",
                      gp = gpar(col = ci_col_od[i], fill = ci_col_od[i]))
  }
}


jpeg(paste0(RES_DIR, "Figures/SA_Forest_TR2_OD.jpg"),
     width = 3600, height = 5000, res = 300)
print(p_od)
dev.off()

# 4B. PLOT GENERATION: DISCONTINUATION -----------------------------------------
cat("\n--- Generating Discontinuation Plot (Trial 3) ---\n")
plot_data_dis <- build_wide_df(
  manifest    = sa_manifest_dis,
  outcome     = "dis",
  max_time    = 365,
  min_events  = 5       # <-- tune: suppress SA if any arm has < this many discont events
) %>% filter(!grepl("Composite outcome", `Sensitivity Analysis`))

bg_fill_dis <- plot_data_dis$bg_color
ci_col_dis  <- trt_colors[plot_data_dis$`Dose adjustment`]
ci_col_dis[is.na(ci_col_dis)] <- "black"

risk_col_dis   <- "% risk at 12 months\n(95% CI)"
header_idx_dis <- which(plot_data_dis$is_header)

p_dis <- forest(
  plot_data_dis[, c("Sensitivity Analysis", "Dose adjustment",
                    "Episodes at risk\nat reinitiation", risk_col_dis,
                    " ", "% difference in 12-month risk (95% CI)")],
  est = plot_data_dis$rd_pct, lower = plot_data_dis$lci_pct, upper = plot_data_dis$uci_pct,
  ci_column = 6,
  vert_line = c(-20, -10, 10, 20),
  ref_line = 0, xlim = c(-22, 22), ticks_at = c(-20, -10, 0, 10, 20),
  xlab = "Risk difference (%)",
  theme = forest_theme(base_size = 12, ci_pch = 15, ci_col = "black", ci_fill = "black", ci_lwd = 1.8,
                       vertline_lty = "dotted",
                       vertline_col = "grey50",
                       vertline_lwd = 1,
                       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
                       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
                       core = list(bg_params = list(fill = bg_fill_dis),
                                   padding = unit(c(2, 2), "mm")))
)

p_dis <- edit_plot(p_dis, row = header_idx_dis, gp = gpar(fontface = "bold"))
p_dis <- edit_plot(p_dis, col = 3, part = "body",   hjust = 0.5, x = unit(0.5, "npc"))
p_dis <- edit_plot(p_dis, col = 3, part = "header", hjust = 0.5, x = unit(0.5, "npc"))
p_dis <- edit_plot(p_dis, col = 4, part = "header", hjust = 0.5, x = unit(0.5, "npc"))
p_dis <- edit_plot(p_dis, col = 6, part = "body", hjust = 0, x = unit(.2, "npc"))
p_dis <- edit_plot(p_dis, col = 6, part = "header", hjust = 0.5, x = unit(0.15, "npc"))


for (i in 1:nrow(plot_data_dis)) {
  if (!plot_data_dis$is_header[i] && !is.na(ci_col_dis[i]) && ci_col_dis[i] != "black") {
    p_dis <- edit_plot(p_dis, row = i, col = 6, which = "ci",
                       gp = gpar(col = ci_col_dis[i], fill = ci_col_dis[i]))
  }
}

jpeg(paste0(RES_DIR, "Figures/SA_Forest_TR2_DIS.jpg"),
     width = 3700, height = 4800, res = 300)
print(p_dis)
dev.off()
