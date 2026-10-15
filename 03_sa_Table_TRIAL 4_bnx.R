library(tidyverse)
library(forestploter)
library(grid)

# 1. SETUP & MASTER MANIFEST ---------------------------------------------------
RES_DIR <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/BNX/"

sa_manifest_od <- tibble::tribble(
  ~category, ~sa_label, ~obj_type, ~risk_file, ~ps_file,
  "0. Primary Analysis", "First interruption in\nany episode", "primary_all", "Risk/risk_tr4_rmst.rds", "PScore/ps_weight_tr4_all.rds",
  "1. Sample restriction", "First interruption in\nfirst episode", "primary_first", "Risk/risk_tr4_rmst.rds", "PScore/ps_weight_tr4_first.rds",
  "1. Sample restriction", "Interruption following\ndaily witnessed ingestion", "standard", "Risk/risk_tr4_dwi.rds", "PScore/ps_weight_tr4_dwi.rds",
  "1. Sample restriction", "Interruption following\ntake-home dispensing", "standard", "Risk/risk_tr4_dwi_0.rds", "PScore/ps_weight_tr4_dwi_0.rds",
  "1. Sample restriction", "People who completed\ninduction", "standard", "Sensitivities/Risk/risk_tr4_SA02_Induction.rds", "Sensitivities/ps_wt_tr4_SA02_Induction.rds",
  "1. Sample restriction", "Interruption at <16/1mg of\nbuprenorphine/nalxone", "standard", "Sensitivities/Risk/risk_tr4_SA03_Dose_LT80.rds", "Sensitivities/ps_wt_tr4_SA03_Dose_LT80.rds",
  "1. Sample restriction", "Interruption at \u226516of\nbuprenorphine/nalxone", "standard", "Sensitivities/Risk/risk_tr4_SA04_Dose_GE80.rds", "Sensitivities/ps_wt_tr4_SA04_Dose_GE80.rds",
  "1. Sample restriction", "Days since OAT initiation:\n\u226430 days", "standard", "Sensitivities/Risk/risk_tr4_SA05_Ret_LE30.rds", "Sensitivities/ps_wt_tr4_SA05_Ret_LE30.rds",
  "1. Sample restriction", "Days since OAT initiation:\n>30 days", "standard", "Sensitivities/Risk/risk_tr4_SA06_Ret_GT30.rds", "Sensitivities/ps_wt_tr4_SA06_Ret_GT30.rds",
  "1. Sample restriction", "Restricted to high\nfentanyl area", "standard", "Sensitivities/Risk/risk_tr4_SA09_HighFentArea.rds", "Sensitivities/ps_wt_tr4_SA09_HighFentArea.rds",
  "1. Sample restriction", "1-year or longer\ncumulative OAT", "standard", "Sensitivities/Risk/risk_tr4_SA11_CumOAT_GE1yr.rds", "Sensitivities/ps_wt_tr4_SA11_CumOAT_GE1yr.rds",
  "1. Sample restriction", "30 days of no OAT\ndispensation before initiation", "standard", "Sensitivities/Risk/risk_tr4_SA13_30dWashout.rds", "Sensitivities/ps_wt_tr4_SA13_30dWashout.rds",
  "1. Sample restriction", "Alternative interruption\nthresholds: 3 to 4 days", "standard", "Risk/risk_tr4_g23.rds", "PScore/ps_weight_tr4_g23.rds",
  
    
  "2. Timeline restriction", "Restricted to post-fentanyl\nera (\u2265 04/01/2012)", "standard", "Sensitivities/Risk/risk_tr4_SA10_PostFentEra.rds", "Sensitivities/ps_wt_tr4_SA10_PostFentEra.rds",
  
  "3. Exposure classification", "5 category dose\nadjustment", "standard", "Sensitivities/Risk/risk_tr4_5cat_od.rds", "Sensitivities/ps_weight_tr4_5cat.rds",
  
  "4. Outcome definition", "Overdose-related acute\ncare visits or death", "standard", "Sensitivities/Risk/risk_tr4_SA16_Composite_OD.rds", "Sensitivities/ps_wt_tr4_SA16_Composite_OD.rds",
  "4. Outcome definition", "All-cause acute\ncare visits", "standard", "Sensitivities/Risk/risk_tr4_all_acute.rds", "PScore/ps_weight_tr4_all.rds",
  
  "5. Causal contrast", "Censoring at subsequent\ninterruption of >2 days", "standard", "Sensitivities/Risk/risk_tr4_per_protocol.rds", "PScore/ps_weight_tr4_hdPS.rds",
  
  "5. Model specification", "High-dimensional propensity\nscore (hdPS)", "hdps", "Risk/risk_tr4_hdPS.rds", "PScore/ps_weight_tr4_hdPS.rds",  
  "5. Model specification", "hdPS with machine\nlearning algorithm", "hdps", "Risk/risk_tr4_hdPS_ML.rds", "PScore/ps_weight_tr4_hdPS_ML.rds",
  "5. Model specification", "Overlap weights\nbased analysis", "standard", "Risk/risk_tr4_rmst_ow.rds", "PScore/ps_weight_tr4_all.rds",
  "5. Model specification", "PS model with additional\nprovider-level covariates", "standard", "Risk/risk_tr4_rmst_cl.rds", "PScore/ps_weight_tr4_all_cl.rds",
  "5. Model specification", "PS model with additional\ntake-home dose related\ncovariate ", "standard", "Risk/risk_tr4_carry.rds", "PScore/ps_weight_tr4_carry.rds"
  
)

sa_manifest_dis <- tibble::tribble(
  ~category, ~sa_label, ~obj_type, ~risk_file, ~ps_file,
  "0. Primary Analysis", "First interruption in\nany episode", "primary_all", "Risk/risk_tr4_rmst.rds", "PScore/ps_weight_tr4_all.rds",
  "1. Sample restriction", "First interruption in\nfirst episode", "primary_first", "Risk/risk_tr4_rmst.rds", "PScore/ps_weight_tr4_first.rds",
  "1. Sample restriction", "Interruption following\ndaily witnessed ingestion", "standard", "Risk/risk_tr4_dwi.rds", "PScore/ps_weight_tr4_dwi.rds",
  "1. Sample restriction", "Interruption following\ntake-home dispensing", "standard", "Risk/risk_tr4_dwi_0.rds", "PScore/ps_weight_tr4_dwi_0.rds",
  "1. Sample restriction", "People who completed\ninduction", "standard", "Sensitivities/Risk/risk_tr4_SA02_Induction.rds", "Sensitivities/ps_wt_tr4_SA02_Induction.rds",
  "1. Sample restriction", "Interruption at <16/1mg of\nbuprenorphine/nalxone", "standard", "Sensitivities/Risk/risk_tr4_SA03_Dose_LT80.rds", "Sensitivities/ps_wt_tr4_SA03_Dose_LT80.rds",
  "1. Sample restriction", "Interruption at \u226516of\nbuprenorphine/nalxone", "standard", "Sensitivities/Risk/risk_tr4_SA04_Dose_GE80.rds", "Sensitivities/ps_wt_tr4_SA04_Dose_GE80.rds",
  "1. Sample restriction", "Days since OAT initiation:\n\u226430 days", "standard", "Sensitivities/Risk/risk_tr4_SA05_Ret_LE30.rds", "Sensitivities/ps_wt_tr4_SA05_Ret_LE30.rds",
  "1. Sample restriction", "Days since OAT initiation:\n>30 days", "standard", "Sensitivities/Risk/risk_tr4_SA06_Ret_GT30.rds", "Sensitivities/ps_wt_tr4_SA06_Ret_GT30.rds",
  "1. Sample restriction", "Restricted to high\nfentanyl area", "standard", "Sensitivities/Risk/risk_tr4_SA09_HighFentArea.rds", "Sensitivities/ps_wt_tr4_SA09_HighFentArea.rds",
  "1. Sample restriction", "1-year or longer\ncumulative OAT", "standard", "Sensitivities/Risk/risk_tr4_SA11_CumOAT_GE1yr.rds", "Sensitivities/ps_wt_tr4_SA11_CumOAT_GE1yr.rds",
  "1. Sample restriction", "30 days of no OAT\ndispensation before initiation", "standard", "Sensitivities/Risk/risk_tr4_SA13_30dWashout.rds", "Sensitivities/ps_wt_tr4_SA13_30dWashout.rds",
  "1. Sample restriction", "Alternative interruption\nthresholds: 3 to 4 days", "standard", "Risk/risk_tr4_g23.rds", "PScore/ps_weight_tr4_g23.rds",
  
  
  "2. Timeline restriction", "Restricted to post-fentanyl\nera (\u2265 04/01/2012)", "standard", "Sensitivities/Risk/risk_tr4_SA10_PostFentEra.rds", "Sensitivities/ps_wt_tr4_SA10_PostFentEra.rds",
  
  "3. Exposure classification", "5 category dose\nadjustment", "standard", "Sensitivities/Risk/risk_tr4_5cat_dis.rds", "Sensitivities/ps_weight_tr4_5cat.rds",
  
  #  "4. Follow-up", "No truncation of follow-up\nat 12 months", "standard", "Sensitivities/Risk/risk_tr4_SA12_NoTruncation.rds", "Sensitivities/ps_wt_tr4_SA12_NoTruncation.rds",
  
  "4. Outcome definition", "No dispensation of OAT\nfor 30 consecutive days", "standard", "Sensitivities/Risk/risk_tr4_discont30.rds", "Sensitivities/ps_weight_tr4_discont30.rds",
  
  "5. Causal contrast", "Censoring at subsequent\ninterruption of >2 days", "standard", "Sensitivities/Risk/risk_tr4_per_protocol.rds", "PScore/ps_weight_tr4_hdPS.rds",
  
  "5. Model specification", "High-dimensional propensity\nscore (hdPS)", "hdps", "Risk/risk_tr4_hdPS.rds", "PScore/ps_weight_tr4_hdPS.rds",  
  "5. Model specification", "hdPS with machine\nlearning algorithm", "hdps", "Risk/risk_tr4_hdPS_ML.rds", "PScore/ps_weight_tr4_hdPS_ML.rds",
  "5. Model specification", "Overlap weights\nbased analysis", "standard", "Risk/risk_tr4_rmst_ow.rds", "PScore/ps_weight_tr4_all.rds",
  "5. Model specification", "PS model with additional\nprovider-level covariates", "standard", "Risk/risk_tr4_rmst_cl.rds", "PScore/ps_weight_tr4_all_cl.rds",
  "5. Model specification", "PS model with additional\ntake-home dose related\ncovariate ", "standard", "Risk/risk_tr4_carry.rds", "PScore/ps_weight_tr4_carry.rds"
  
)

trt_clean <- c(
  "no reversion" = "No dose change", 
  "1-<50%" = "<50% reduction", 
  ">=50%/<=4mg" = "\u226550% reduction/\u22644mg", 
  "increase" = "Dose increase",
  "increase <50%" = "<50% increase",
  "increase >=50%"  = "\u226550% increase"
)

trt_colors <- c(
  "No dose change" = "black", 
  "<50% reduction" = "#2CA02C", 
  "\u226550% reduction/\u22644mg" = "#FF7F0E", 
  "Dose increase" = "#D62728",
  "<50% increase" = "#9467BD", 
  "\u226550% increase"  = "#A00000"
)

# 2. SMART EXTRACTION HELPER ---------------------------------------------------
extract_nested_data <- function(risk_file, ps_file, obj_type, outcome, max_time, sa_label) {
  
  full_risk <- paste0(RES_DIR, risk_file)
  full_ps   <- paste0(RES_DIR, ps_file)
  
  is_5cat <- grepl("5 category", sa_label)
  local_levels <- if(is_5cat) {
    c("no reversion", "1-<50%", ">=50%/<=4mg", "increase <50%", "increase >=50%")
  } else {
    c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
  }
  
  na_block <- tibble(trt = local_levels, Episode = NA_character_, abs_risk_disp = "", 
                     rd_pct = NA_real_, lci_pct = NA_real_, uci_pct = NA_real_, rd_disp = "")
  
  if (!file.exists(full_risk) || !file.exists(full_ps)) return(na_block)
  
  risk_obj <- readRDS(full_risk)
  ps_obj   <- readRDS(full_ps)
  
  dat_w <- if("data_weights" %in% names(ps_obj)) ps_obj$data_weights else ps_obj
  
  counts <- tryCatch({
    dat_w %>% 
      count(trt, name = "Episode") %>% 
      mutate(trt = as.character(trt), Episode = as.character(Episode))
  }, error = function(e) tibble(trt = local_levels, Episode = NA_character_))
  
  # --- FIX FOR 5-CATEGORY DISCONTINUATION PLOT ---
  risks_df <- NULL
  contrasts <- NULL
  
  tryCatch({
    # Smart Node Locator: If the exact 'outcome' name isn't found, 
    # it safely defaults to the first available node (which fixes the 5cat 'dis' bug)
    if (obj_type %in% c("primary_all", "ow_all")) {
      out_node <- if (!is.null(risk_obj$all[[outcome]])) outcome else names(risk_obj$all)[1]
      base_node <- risk_obj$all[[out_node]]
    } else if (obj_type == "primary_first") {
      out_node <- if (!is.null(risk_obj$first[[outcome]])) outcome else names(risk_obj$first)[1]
      base_node <- risk_obj$first[[out_node]]
    } else {
      out_node <- if (!is.null(risk_obj[[outcome]])) outcome else names(risk_obj)[1]
      base_node <- risk_obj[[out_node]]
    }
    
    node <- if (!is.null(base_node$ci)) base_node$ci else base_node
    risks_df <- node$risks
    contrasts <- node$contrasts
  }, error = function(e) { NULL })
  
  if (is.null(contrasts) || nrow(contrasts) == 0 || is.null(risks_df)) return(na_block)
  
  abs_data <- risks_df %>%
    filter(horizon == max_time) %>%
    mutate(
      risk_pct = risk * 100, risk_lower_pct = risk_lower * 100, risk_upper_pct = risk_upper * 100,
      abs_risk_disp = ifelse(abs(risk_pct) >= 10,
                             sprintf("%.1f (%.1f, %.1f)", risk_pct, risk_lower_pct, risk_upper_pct),
                             sprintf("%.2f (%.2f, %.2f)", risk_pct, risk_lower_pct, risk_upper_pct))
    ) %>% select(trt, abs_risk_disp)
  
  rd_data <- contrasts %>%
    filter(horizon == max_time) %>%
    mutate(
      rd_pct = RD * 100, lci_pct = RD_lower * 100, uci_pct = RD_upper * 100,
      rd_disp = ifelse(abs(rd_pct) >= 10,
                       sprintf("%.1f (%.1f, %.1f)", rd_pct, lci_pct, uci_pct),
                       sprintf("%.2f (%.2f, %.2f)", rd_pct, lci_pct, uci_pct))
    ) %>% select(trt, rd_pct, lci_pct, uci_pct, rd_disp)
  
  tibble(trt = local_levels) %>% 
    left_join(counts, by = "trt") %>% 
    left_join(abs_data, by = "trt") %>%
    left_join(rd_data, by = "trt") %>%
    mutate(rd_disp = ifelse(trt == "no reversion", "Reference", rd_disp))
}



# 3. EXPORT ESTIMATES TO WORD TABLES -------------------------------------------

library(tidyverse)
library(officer)
library(flextable)

build_word_table <- function(manifest, outcome, max_time) {
  
  purrr::pmap_dfr(manifest, function(category, sa_label, obj_type, risk_file, ps_file) {
    
    extract_nested_data(risk_file, ps_file, obj_type, outcome, max_time, sa_label) %>%
      transmute(
        Category = sub("^\\d+\\.\\s*", "", category),
        `Sensitivity analysis` = str_squish(str_replace_all(sa_label, "\n", " ")),
        `Dose adjustment` = recode(trt, !!!trt_clean),
        `Episodes at risk` = replace_na(Episode, ""),
        `Adjusted risk, % (95% CI)` = replace_na(abs_risk_disp, ""),
        `Risk difference, percentage points (95% CI)` = replace_na(rd_disp, "")
      )
  })
}


format_word_table <- function(dat) {
  
  flextable(dat) %>%
    merge_v(j = c("Category", "Sensitivity analysis")) %>%
    valign(j = c("Category", "Sensitivity analysis"), valign = "top") %>%
    align(j = c("Category", "Sensitivity analysis", "Dose adjustment"), align = "left", part = "all") %>%
    align(j = c("Episodes at risk", "Adjusted risk, % (95% CI)", "Risk difference, percentage points (95% CI)"), align = "center", part = "all") %>%
    bold(part = "header") %>%
    bg(part = "header", bg = "#D9E2F3") %>%
    color(part = "header", color = "#000000") %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 9, part = "body") %>%
    fontsize(size = 9, part = "header") %>%
    padding(padding = 3, part = "all") %>%
    width(j = 1, width = 1.15) %>%
    width(j = 2, width = 2.55) %>%
    width(j = 3, width = 1.45) %>%
    width(j = 4, width = 0.85) %>%
    width(j = 5, width = 1.55) %>%
    width(j = 6, width = 1.85) %>%
    set_table_properties(layout = "fixed", width = 1) %>%
    paginate(init = TRUE, hdr_ftr = TRUE)
}


# 4A. PREPARE OVERDOSE AND DISCONTINUATION TABLES ------------------------------

word_data_od  <- build_word_table(sa_manifest_od,  outcome = "od",  max_time = 30)
word_data_dis <- build_word_table(sa_manifest_dis, outcome = "dis", max_time = 365)

ft_od  <- format_word_table(word_data_od)
ft_dis <- format_word_table(word_data_dis)


# 4B. CREATE LANDSCAPE WORD DOCUMENT -------------------------------------------

landscape_section <- prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.6, bottom = 0.6, left = 0.55, right = 0.55)
)

doc <- read_docx() %>%
  body_set_default_section(landscape_section) %>%
  
  body_add_par(
    "Supplement Table. Sensitivity analyses of buprenorphine/naloxone dose-adjustment strategies",
    style = "heading 1"
  ) %>%
  
  body_add_par(
    "A. Overdose-related acute care visits at 1 month",
    style = "heading 2"
  ) %>%
  
  body_add_flextable(ft_od) %>%
  
  body_add_par(
    "Risk differences are expressed in percentage points. No dose change is the reference strategy. CI indicates confidence interval.",
    style = "Normal"
  ) %>%
  
  body_add_break() %>%
  
  body_add_par(
    "B. OAT discontinuation at 12 months",
    style = "heading 2"
  ) %>%
  
  body_add_flextable(ft_dis) %>%
  
  body_add_par(
    "Risk differences are expressed in percentage points. No dose change is the reference strategy. CI indicates confidence interval; OAT, opioid agonist treatment.",
    style = "Normal"
  )

print(doc, target = paste0(RES_DIR, "Tables/SA_Risk_Estimates_TR4.docx"))