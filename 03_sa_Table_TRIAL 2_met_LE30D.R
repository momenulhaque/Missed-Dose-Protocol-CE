library(tidyverse)
library(gt)
library(webshot2) # Required for high-res gtsave image export

# 1. SETUP & MASTER MANIFEST ---------------------------------------------------
RES_DIR <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/LE30D/"



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
  "6. Model specification", "PS model with additional\nprovider-level covariates", "standard", "Risk/risk_tr2_cl.rds", "PScore/ps_weight_tr2_34d_all_cl.rds",
  "6. Model specification", "PS model with additional\nprovider-level covariates\nplus HSDA", "standard", "Risk/risk_tr2_cluster_hsda.rds", "PScore/ps_weight_tr2_cluster_hsda.rds",
  "6. Model specification", "Non-parametric bootstric\nfor 95% CI", "standard", "Risk/risk_tr2_boot.rds", "PScore/ps_weight_tr2_34d_all.rds"
  
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
  "6. Model specification", "PS model with additional\nprovider-level covariates", "standard", "Risk/risk_tr2_cl.rds", "PScore/ps_weight_tr2_34d_all_cl.rds",
  "6. Model specification", "PS model with additional\nprovider-level covariates\nplus HSDA", "standard", "Risk/risk_tr2_cluster_hsda.rds", "PScore/ps_weight_tr2_cluster_hsda.rds",
  "6. Model specification", "Non-parametric bootstric\nfor 95% CI", "standard", "Risk/risk_tr2_boot.rds", "PScore/ps_weight_tr2_34d_all.rds"
  
)


library(dplyr)
library(stringr)
library(tidyr)
library(openxlsx)
library(officer)
library(flextable)

OUT_DIR <- file.path(RES_DIR, "Tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

clean_sa_table <- function(x, outcome) {
  risk_col <- if (outcome == "od") "% risk at 1 month\n(95% CI)" else "% risk at 12 months\n(95% CI)"
  rd_col <- if (outcome == "od") "% difference in 1-month risk (95% CI)" else "% difference in 12-month risk (95% CI)"

  x %>%
    filter(!is_header) %>%
    transmute(
      category = str_remove(`Sensitivity Analysis`, "^\\s+"),
      sensitivity_analysis = str_squish(str_replace_all(`Sensitivity Analysis`, "\\n", " ")),
      dose_adjustment = `Dose adjustment`,
      episodes_at_risk = suppressWarnings(as.integer(`Episodes at risk\nat reinitiation`)),
      risk_95ci = .data[[risk_col]],
      risk_difference_95ci = .data[[rd_col]],
      min_arm_events = min_arm_events
    ) %>%
    # build_wide_df() places the category only in a separate header row; restore it
    mutate(category = NA_character_) 
}

# More reliable cleaner: use manifest order to reconstruct category and SA labels.
make_export_data <- function(manifest, outcome, max_time, min_events) {
  
  rows <- vector("list", nrow(manifest))
  
  for (i in seq_len(nrow(manifest))) {
    
    m <- manifest[i, ]
    
    z <- extract_nested_data(
      risk_file = m$risk_file,
      ps_file   = m$ps_file,
      obj_type  = m$obj_type,
      outcome   = outcome,
      max_time  = max_time,
      sa_label  = m$sa_label,
      min_events = min_events
    )
    
    # Analysis-level attributes
    suppressed_flag <- isTRUE(attr(z, "suppressed"))
    
    min_events_value <- attr(z, "min_n_events")
    
    if (is.null(min_events_value) ||
        length(min_events_value) == 0L ||
        !is.finite(min_events_value)) {
      min_events_value <- NA_integer_
    } else {
      min_events_value <- as.integer(min_events_value)
    }
    
    # Construct scalar labels before mutate()
    outcome_label <- if (outcome == "od") {
      "Overdose"
    } else {
      "Treatment discontinuation"
    }
    
    category_label <- stringr::str_remove(
      as.character(m$category),
      "^[0-9]+\\.\\s*"
    )
    
    category_order_value <- suppressWarnings(
      as.integer(
        stringr::str_extract(
          as.character(m$category),
          "^[0-9]+"
        )
      )
    )
    
    analysis_label <- stringr::str_squish(
      stringr::str_replace_all(
        as.character(m$sa_label),
        "\n",
        " "
      )
    )
    
    if (suppressed_flag) {
      analysis_label <- paste0(analysis_label, "†")
    }
    
    rows[[i]] <- z %>%
      dplyr::mutate(
        outcome = outcome_label,
        category = category_label,
        category_order = category_order_value,
        analysis_order = i,
        sensitivity_analysis = analysis_label,
        
        dose_adjustment = dplyr::recode(
          as.character(trt),
          !!!trt_clean
        ),
        
        episodes_at_risk = suppressWarnings(
          as.integer(Episode)
        ),
        
        risk_95ci = dplyr::na_if(
          abs_risk_disp,
          ""
        ),
        
        risk_difference_95ci = dplyr::na_if(
          rd_disp,
          ""
        ),
        
        suppressed = suppressed_flag,
        min_arm_events = min_events_value
      ) %>%
      dplyr::select(
        outcome,
        category_order,
        category,
        analysis_order,
        sensitivity_analysis,
        dose_adjustment,
        episodes_at_risk,
        risk_95ci,
        risk_difference_95ci,
        suppressed,
        min_arm_events
      )
  }
  
  dplyr::bind_rows(rows)
}
od_export <- make_export_data(
  manifest = sa_manifest_od,
  outcome = "od",
  max_time = 30,
  min_events = 1
)

dis_export <- make_export_data(
  manifest = sa_manifest_dis,
  outcome = "dis",
  max_time = 365,
  min_events = 5
)

# ---------- Excel: analysis-ready values plus publication display ----------
wb <- createWorkbook()
addWorksheet(wb, "Overdose")
addWorksheet(wb, "Discontinuation")
addWorksheet(wb, "Data dictionary")

writeDataTable(wb, "Overdose", od_export, tableStyle = "TableStyleMedium2")
writeDataTable(wb, "Discontinuation", dis_export, tableStyle = "TableStyleMedium2")
writeData(wb, "Data dictionary", data.frame(
  variable = c("category", "sensitivity_analysis", "dose_adjustment", "episodes_at_risk",
               "risk_95ci", "risk_difference_95ci", "suppressed", "min_arm_events"),
  definition = c("Sensitivity-analysis domain", "Analysis specification", "Dose-adjustment strategy",
                 "Number of episodes at risk at reinitiation", "Absolute risk, percentage (95% CI)",
                 "Risk difference versus no dose change, percentage points (95% CI)",
                 "TRUE when estimates were suppressed by the minimum-event rule",
                 "Minimum event count across included treatment arms")
))

for (s in c("Overdose", "Discontinuation", "Data dictionary")) {
  freezePane(wb, s, firstRow = TRUE)
  setColWidths(wb, s, cols = 1:11, widths = "auto")
}
setColWidths(wb, "Overdose", cols = c(3,5,8,9), widths = c(22,42,24,32))
setColWidths(wb, "Discontinuation", cols = c(3,5,8,9), widths = c(22,42,24,32))

xlsx_path <- file.path(OUT_DIR, "sa_Table_TRIAL_2_met.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# ---------- Word: compact, publication-ready tables ----------
make_word_df <- function(dat, risk_heading, rd_heading) {
  out <- list()
  cats <- unique(dat$category)
  r <- 0L
  for (cat in cats) {
    r <- r + 1L
    out[[r]] <- tibble(
      `Sensitivity analysis` = cat,
      `Dose adjustment` = "",
      `Episodes at risk` = "",
      !!risk_heading := "",
      !!rd_heading := "",
      row_type = "category"
    )
    dcat <- dat %>% filter(category == cat)
    for (a in unique(dcat$analysis_order)) {
      da <- dcat %>% filter(analysis_order == a)
      for (j in seq_len(nrow(da))) {
        r <- r + 1L
        out[[r]] <- tibble(
          `Sensitivity analysis` = if (j == 1L) da$sensitivity_analysis[j] else "",
          `Dose adjustment` = da$dose_adjustment[j],
          `Episodes at risk` = ifelse(is.na(da$episodes_at_risk[j]), "", format(da$episodes_at_risk[j], big.mark = ",")),
          !!risk_heading := coalesce(da$risk_95ci[j], ""),
          !!rd_heading := coalesce(da$risk_difference_95ci[j], ""),
          row_type = "data"
        )
      }
    }
  }
  bind_rows(out)
}

make_ft <- function(dat, risk_heading, rd_heading) {
  
  d <- make_word_df(
    dat = dat,
    risk_heading = risk_heading,
    rd_heading = rd_heading
  )
  
  category_rows <- which(d$row_type == "category")
  
  d <- d %>%
    dplyr::select(-row_type)
  
  ft <- flextable::flextable(d)
  
  ft <- flextable::set_header_labels(
    ft,
    `Sensitivity analysis` = "Sensitivity analysis",
    `Dose adjustment` = "Dose adjustment",
    `Episodes at risk` = "Episodes at risk\nat reinitiation"
  )
  
  ft <- flextable::theme_booktabs(ft)
  
  # General font formatting
  ft <- flextable::font(
    ft,
    fontname = "Arial",
    part = "all"
  )
  
  ft <- flextable::fontsize(
    ft,
    size = 8.5,
    part = "body"
  )
  
  ft <- flextable::fontsize(
    ft,
    size = 8.5,
    part = "header"
  )
  
  # Main column header
  ft <- flextable::bold(
    ft,
    part = "header"
  )
  
  ft <- flextable::color(
    ft,
    color = "white",
    part = "header"
  )
  
  ft <- flextable::bg(
    ft,
    bg = "#1F4E78",
    part = "header"
  )
  
  # Alignment
  ft <- flextable::align(
    ft,
    j = 1:2,
    align = "left",
    part = "body"
  )
  
  ft <- flextable::align(
    ft,
    j = 3:5,
    align = "center",
    part = "body"
  )
  
  ft <- flextable::align(
    ft,
    j = 1:5,
    align = "center",
    part = "header"
  )
  
  ft <- flextable::valign(
    ft,
    valign = "center",
    part = "all"
  )
  
  # Cell padding
  ft <- flextable::padding(
    ft,
    padding.top = 3,
    padding.bottom = 3,
    padding.left = 4,
    padding.right = 4,
    part = "all"
  )
  
  # Column widths
  ft <- flextable::width(
    ft,
    j = 1,
    width = 2.40
  )
  
  ft <- flextable::width(
    ft,
    j = 2,
    width = 1.55
  )
  
  ft <- flextable::width(
    ft,
    j = 3,
    width = 1.05
  )
  
  ft <- flextable::width(
    ft,
    j = 4,
    width = 1.55
  )
  
  ft <- flextable::width(
    ft,
    j = 5,
    width = 2.05
  )
  
  # Format category rows
  if (length(category_rows) > 0L) {
    
    ft <- flextable::bg(
      ft,
      i = category_rows,
      j = 1:5,
      bg = "#D9E2F3",
      part = "body"
    )
    
    ft <- flextable::bold(
      ft,
      i = category_rows,
      j = 1:5,
      bold = TRUE,
      part = "body"
    )
    
    ft <- flextable::color(
      ft,
      i = category_rows,
      j = 1:5,
      color = "#1F4E78",
      part = "body"
    )
    
    # Merge each category row independently
    for (row_i in category_rows) {
      
      ft <- flextable::merge_at(
        ft,
        i = row_i,
        j = 1:5,
        part = "body"
      )
    }
    
    # Align the merged category labels to the left
    ft <- flextable::align(
      ft,
      i = category_rows,
      j = 1,
      align = "left",
      part = "body"
    )
    
    # Add a horizontal border to each category row
    ft <- flextable::hline(
      ft,
      i = category_rows,
      border = officer::fp_border(
        color = "#9EADBA",
        width = 0.7
      ),
      part = "body"
    )
  }
  
  ft <- flextable::set_table_properties(
    ft,
    layout = "fixed",
    width = 1
  )
  
  ft <- flextable::paginate(
    ft,
    init = TRUE,
    hdr_ftr = TRUE
  )
  
  return(ft)
}

ft_od <- make_ft(od_export, "1-month risk, %\n(95% CI)", "1-month risk difference, pp\n(95% CI)")
ft_dis <- make_ft(dis_export, "12-month risk, %\n(95% CI)", "12-month risk difference, pp\n(95% CI)")

landscape <- prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.55, bottom = 0.55, left = 0.55, right = 0.55)
)

doc <- read_docx()
doc <- body_add_par(doc, "Sensitivity analyses: overdose", style = "heading 1")
doc <- body_add_flextable(doc, ft_od)
doc <- body_add_par(doc, "Risk differences are percentage-point contrasts versus no dose change. † Estimates suppressed when the minimum event-count rule was not met.", style = "Normal")
doc <- body_end_section_landscape(doc)
doc <- body_add_par(doc, "Sensitivity analyses: treatment discontinuation", style = "heading 1")
doc <- body_add_flextable(doc, ft_dis)
doc <- body_add_par(doc, "Risk differences are percentage-point contrasts versus no dose change. † Estimates suppressed when fewer than 5 events occurred in at least one treatment arm.", style = "Normal")
doc <- body_set_default_section(doc, landscape)

docx_path <- file.path(OUT_DIR, "sa_Table_TRIAL_2_met.docx")
print(doc, target = docx_path)

message("Created:\n", xlsx_path, "\n", docx_path)

