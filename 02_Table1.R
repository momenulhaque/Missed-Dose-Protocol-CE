#  TABLE 1: LE30D vs GT30D ##############
#  Stratified by retention length (LE30D vs GT30D) for TR2 and TR3
#  NO SMDs in this file.


options(error = NULL)

library(tidyverse)
library(haven)
library(tableone)
library(openxlsx)

# ── Setup ────────────────────────────────────────────────────────────────────
DATA_DIR <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR  <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/"
OUT_XLSX <- paste0(RES_DIR, "Tables/Table1_Baseline_Stratified.xlsx")
dir.create(dirname(OUT_XLSX), showWarnings = FALSE, recursive = TRUE)

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
TRT_LABELS <- c("No dose change, n (%)",
                "<50% reduction, n (%)",
                "\u226550% reduction/\u226440 mg, n (%)",
                "Dose increase, n (%)")

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss", "prior_episode"
)

VAR_LABELS <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "  30-39",
  age_40_49     = "  40-49",
  age_ge50      = "  50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date [Median (Q1-Q3)]",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode prior to the index interruption (days) [Median (Q1-Q3)]",
  home_grp_b    = "Unstable housing in past 5 years",
  incar_grp_bl  = "Incarcerated in the past year",
  drug_acute    = "Drug-related acute care visits in the past month",
  cci_gt0       = "Charlson Comorbidity Index (>0)",
  copd_asthma   = "Asthma or COPD in the past year",
  sed_psy_med   = "Dispensations of psychiatric or sedative medication in the past month",
  po            = "Dispensations opioid other than OAT in the past month",
  sud           = "Other substance use disorder (ever)",
  aud           = "Alcohol use disorder (ever)",
  smd           = "Serious mental disorder (ever)",
  hcv           = "Hepatitis C (ever)",
  tob           = "Tobacco use disorder (past year)",
  udt           = "Urine drug test in the past month",
  ddose_prev    = "Last daily dose [Median (Q1-Q3)]",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_episode = "Number of prior episodes [Median (Q1-Q3)]"
)

# ── Cleaning helper ──────────────────────────────────────────────────────────
clean_trt_all <- function(dat) {
  dat |>
    mutate(trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"),
                              "1-<50%", trt_s456)) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS),
           prior_episode = as.numeric(episode - 1)) |>
    group_by(studyid_rev, episode) |>
    mutate(prior_gap     = cumsum(gap_days),
           prior_gap_num = row_number()) |>
    ungroup()
}

# ── Load cohorts ─────────────────────────────────────────────────────────────
cat("Loading raw data ...\n")
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))
dat_tr2_le <- raw_all |> filter(trial_gap %in% 3:4)  |> clean_trt_all() |> filter(ret_length <= 30)
dat_tr2_gt <- raw_all |> filter(trial_gap %in% 3:4)  |> clean_trt_all() |> filter(ret_length >  30)
dat_tr3_le <- raw_all |> filter(trial_cat == "5-14") |> clean_trt_all() |> filter(ret_length <= 30)
dat_tr3_gt <- raw_all |> filter(trial_cat == "5-14") |> clean_trt_all() |> filter(ret_length >  30)
rm(raw_all); gc()

# ── Post-processing helpers ──────────────────────────────────────────────────
# 1. Replace "median [Q1, Q3]" comma with hyphen for ALL median rows.
comma_to_hyphen <- function(s) {
  # Matches "[num, num]" (with optional decimals) inside median brackets
  gsub("\\[([-0-9.]+),\\s*([-0-9.]+)\\]", "[\\1-\\2]", s)
}

# 2. Round year values to integers: "2017.00 [2013.00-2020.00]" -> "2017 [2013-2020]"
round_year_cell <- function(s) {
  # Empty or non-median cells: return as-is
  if (is.na(s) || s == "") return(s)
  # Format A: "M [Q1-Q3]"
  m <- regmatches(s, regexec("^\\s*([-0-9.]+)\\s*\\[\\s*([-0-9.]+)\\s*-\\s*([-0-9.]+)\\s*\\]\\s*$", s))[[1]]
  if (length(m) == 4) {
    return(sprintf("%d [%d-%d]",
                   round(as.numeric(m[2])),
                   round(as.numeric(m[3])),
                   round(as.numeric(m[4]))))
  }
  # Format B: plain single number (should not occur here but safe)
  m2 <- regmatches(s, regexec("^\\s*([-0-9.]+)\\s*$", s))[[1]]
  if (length(m2) == 2) return(sprintf("%d", round(as.numeric(m2[2]))))
  s
}

# ── Core function: build Table 1 (no SMD) ────────────────────────────────────
build_table1 <- function(dat) {
  cont_vars <- c("year", "ret_length", "ddose_prev", "prior_episode")
  cat_vars  <- setdiff(baselinevars, cont_vars)
  
  tab <- CreateTableOne(vars = baselinevars, strata = "trt",
                        data = dat, factorVars = cat_vars,
                        test = FALSE, smd = FALSE)
  
  mat <- print(tab, showAllLevels = FALSE,
               printToggle = FALSE, noSpaces = TRUE,
               contDigits = 2, catDigits = 1,
               nonnormal = cont_vars)   # forces median [IQR] instead of mean (SD)
  
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  
  # ---- Fix 1: convert "[Q1, Q3]" to "[Q1-Q3]" across ALL cells 
  out[] <- lapply(out, function(col) sapply(col, comma_to_hyphen))
  
  # ---- Fix 2: round year row to integers 
  rn <- rownames(out)
  rn_key <- gsub(" \\(mean \\(SD\\)\\)| = 1 \\(%\\)| \\(median \\[IQR\\]\\)", "", rn)
  year_row <- which(rn_key == "year")
  if (length(year_row) > 0) {
    out[year_row, ] <- lapply(out[year_row, , drop = FALSE], round_year_cell)
  }
  
  # ---- Rename arm columns 
  arm_cols <- match(TRT_LEVELS, colnames(out))
  arm_cols <- arm_cols[!is.na(arm_cols)]
  colnames(out)[arm_cols] <- TRT_LABELS[seq_along(arm_cols)]
  
  # ---- Apply VAR_LABELS to row names
  rn_clean <- rn_key
  for (v in names(VAR_LABELS)) rn_clean[rn_clean == v] <- VAR_LABELS[v]
  rownames(out) <- rn_clean
  
  tibble::rownames_to_column(out, "Characteristic")
}

# ── Build ────────────────────────────────────────────────────────────────────
cat("\nBuilding Table 1 for TR2 LE30D ...\n"); tab_tr2_le <- build_table1(dat_tr2_le)
cat("Building Table 1 for TR2 GT30D ...\n");   tab_tr2_gt <- build_table1(dat_tr2_gt)
cat("Building Table 1 for TR3 LE30D ...\n");   tab_tr3_le <- build_table1(dat_tr3_le)
cat("Building Table 1 for TR3 GT30D ...\n");   tab_tr3_gt <- build_table1(dat_tr3_gt)

# ── Save to Excel ────────────────────────────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "TR2_LE30D_3-4d_within_30d");  writeData(wb, "TR2_LE30D_3-4d_within_30d",  tab_tr2_le)
addWorksheet(wb, "TR2_GT30D_3-4d_after_30d");   writeData(wb, "TR2_GT30D_3-4d_after_30d",   tab_tr2_gt)
addWorksheet(wb, "TR3_LE30D_5-14d_within_30d"); writeData(wb, "TR3_LE30D_5-14d_within_30d", tab_tr3_le)
addWorksheet(wb, "TR3_GT30D_5-14d_after_30d");  writeData(wb, "TR3_GT30D_5-14d_after_30d",  tab_tr3_gt)

header_style <- createStyle(textDecoration = "bold", halign = "center",
                            fgFill = "#D9D9D9", border = "TopBottom",
                            wrapText = TRUE)
for (sh in names(wb)) {
  addStyle(wb, sh, header_style, rows = 1, cols = 1:ncol(tab_tr2_le),
           gridExpand = TRUE)
  setColWidths(wb, sh, cols = 1:ncol(tab_tr2_le),
               widths = c(55, 22, 22, 24, 22))
  freezePane(wb, sh, firstActiveRow = 2, firstActiveCol = 2)
}

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("\n\u2713 Table 1 saved: %s\n", OUT_XLSX))





#  PAIRWISE SMDs: LE30D vs GT30D  ###########
#  Mirrors the love-plot code: 2-arm subsets vs "No dose change" reference,
#  wt_trunc weights, s.d.denom = "pooled".
#  All SMDs shown as signed values to 3 decimals.
#  "Max" column = the SMD with largest |value|, sign preserved.
#  Output: one Excel workbook, 4 sheets (TR2/TR3 × LE30D/GT30D).
options(error = NULL)

library(tidyverse)
library(haven)
library(cobalt)
library(openxlsx)

# ── Setup ────────────────────────────────────────────────────────────────────
DATA_DIR <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR  <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/"
OUT_XLSX <- paste0(RES_DIR, "Tables/Table1_SMDs_Stratified.xlsx")
dir.create(dirname(OUT_XLSX), showWarnings = FALSE, recursive = TRUE)

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")

# Same trt recoding used in love-plot code
TRT_RECODE <- c(
  "no reversion"   = "No dose change",
  "1-<50%"         = "1-<50% reduction",
  ">=50%/<=40mg"   = "\u226550% reduction/\u226440mg",
  "increase"       = "Dose increase"
)

REFERENCE_ARM   <- "No dose change"
COMPARISON_ARMS <- c("1-<50% reduction",
                     "\u226550% reduction/\u226440mg",
                     "Dose increase")

# Short comparison labels for column headers
COMP_SHORT <- c(
  "1-<50% reduction"               = "<50%red vs Ref",
  "\u226550% reduction/\u226440mg" = "\u226550%red/\u226440mg vs Ref",
  "Dose increase"                  = "Increase vs Ref"
)

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss", "prior_episode"
)

VAR_LABELS <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode prior to the index interruption (days)",
  home_grp_b    = "Unstable housing in past 5 years",
  incar_grp_bl  = "Incarcerated in the past year",
  drug_acute    = "Drug-related acute care visits in the past month",
  cci_gt0       = "Charlson Comorbidity Index (>0)",
  copd_asthma   = "Asthma or COPD in the past year",
  sed_psy_med   = "Dispensations of psychiatric or sedative medication in the past month",
  po            = "Dispensations opioid other than OAT in the past month",
  sud           = "Other substance use disorder (ever)",
  aud           = "Alcohol use disorder (ever)",
  smd           = "Serious mental disorder (ever)",
  hcv           = "Hepatitis C (ever)",
  tob           = "Tobacco use disorder (past year)",
  udt           = "Urine drug test in the past month",
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_episode = "Number of prior episodes"
)

# ── Load PS-weighted objects ─────────────────────────────────────────────────
cat("Loading PS-weighted objects ...\n")
ps_tr2_le <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr2_34d_all.rds"))
ps_tr2_gt <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr2_34d_all.rds"))
ps_tr3_le <- readRDS(paste0(RES_DIR, "LE30D/PScore/ps_weight_tr3_all.rds"))
ps_tr3_gt <- readRDS(paste0(RES_DIR, "GT30D/PScore/ps_weight_tr3_all.rds"))

# ── Prep helper — mirrors the love-plot recode ───────────────────────────────
prep_dat <- function(ps_obj) {
  ps_obj$data_weights |>
    mutate(
      trt = recode(trt, !!!TRT_RECODE),
      incar_grp_bl = ifelse(as.character(incar_grp_bl) == "1_pastyear",
                            1, as.numeric(as.character(incar_grp_bl))),
      prior_episode = as.numeric(prior_episode)
    )
}

# ── Helper: pick SMD with largest |value|, preserve original sign ────────────
pick_max_signed <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x[which.max(abs(x))]
}

# ── Extract pairwise SMDs — 2-arm bal.tab per comparison ─────────────────────
compute_pairwise_smds <- function(dat) {
  
  fmla <- as.formula(paste("trt ~", paste(baselinevars, collapse = " + ")))
  
  # Matrices: rows = baseline var, cols = comparison arm
  smd_un  <- matrix(NA_real_, nrow = length(baselinevars),
                    ncol = length(COMPARISON_ARMS),
                    dimnames = list(baselinevars, COMPARISON_ARMS))
  smd_adj <- smd_un
  
  for (arm in COMPARISON_ARMS) {
    cat(sprintf("  Comparing '%s' vs Ref ...\n", arm))
    
    sub <- subset(dat, trt %in% c(REFERENCE_ARM, arm))
    sub$trt <- droplevels(factor(sub$trt, levels = c(REFERENCE_ARM, arm)))
    
    if (any(table(sub$trt) == 0)) {
      cat("    Skipped: one arm empty.\n"); next
    }
    
    bt <- tryCatch(
      bal.tab(
        fmla,
        data       = sub,
        weights    = sub$wt_trunc,
        estimand   = "ATE",
        s.d.denom  = "pooled",
        un         = TRUE,
        stats      = "mean.diffs",
        abs        = FALSE,      # keep signed
        binary     = "std",
        continuous = "std"
      ),
      error = function(e) { cat("    bal.tab error: ", e$message, "\n"); NULL }
    )
    if (is.null(bt)) next
    
    bal_tbl <- bt$Balance
    un_col  <- intersect(colnames(bal_tbl),
                         c("Diff.Un", "M.Threshold.Un", "M.Diff.Un"))[1]
    adj_col <- intersect(colnames(bal_tbl),
                         c("Diff.Adj", "M.Threshold.Adj", "M.Diff.Adj"))[1]
    
    vec_un  <- bal_tbl[, un_col]
    vec_adj <- bal_tbl[, adj_col]
    names(vec_un)  <- rownames(bal_tbl)
    names(vec_adj) <- rownames(bal_tbl)
    
    for (v in baselinevars) {
      if (v %in% names(vec_un)) {
        if (is.finite(vec_un[v]))  smd_un[v,  arm] <- vec_un[v]
        if (is.finite(vec_adj[v])) smd_adj[v, arm] <- vec_adj[v]
      } else {
        # Expanded categorical: keep the level whose |SMD| is largest, preserve sign
        exp_names <- grep(paste0("^", v, "_"), names(vec_un), value = TRUE)
        if (length(exp_names) > 0) {
          exp_un  <- vec_un[exp_names];  exp_un  <- exp_un[is.finite(exp_un)]
          exp_adj <- vec_adj[exp_names]; exp_adj <- exp_adj[is.finite(exp_adj)]
          if (length(exp_un)  > 0) smd_un[v,  arm] <- exp_un[which.max(abs(exp_un))]
          if (length(exp_adj) > 0) smd_adj[v, arm] <- exp_adj[which.max(abs(exp_adj))]
        }
      }
    }
  }
  
  # Max: signed SMD with largest absolute magnitude
  max_signed_un  <- apply(smd_un,  1, pick_max_signed)
  max_signed_adj <- apply(smd_adj, 1, pick_max_signed)
  
  # Format: 3 decimals, blank for NA / non-finite
  fmt <- function(x) ifelse(is.na(x) | !is.finite(x), "", sprintf("%.3f", x))
  
  out <- data.frame(
    Characteristic = VAR_LABELS[baselinevars],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # Before-weighting SMDs (3 comparisons + max)
  for (arm in COMPARISON_ARMS) {
    col_name <- paste0("Before wtg: ", COMP_SHORT[arm])
    out[[col_name]] <- fmt(smd_un[, arm])
  }
  out[["Before wtg: SMD (largest |value|)"]] <- fmt(max_signed_un)
  
  # After-weighting SMDs (3 comparisons + max)
  for (arm in COMPARISON_ARMS) {
    col_name <- paste0("After wtg: ", COMP_SHORT[arm])
    out[[col_name]] <- fmt(smd_adj[, arm])
  }
  out[["After wtg: SMD (largest |value|)"]] <- fmt(max_signed_adj)
  
  out
}

# ── Run all 4 cohorts ────────────────────────────────────────────────────────
cat("\nTR2 LE30D ...\n")
dat_tr2_le <- prep_dat(ps_tr2_le)
smd_tr2_le <- compute_pairwise_smds(dat_tr2_le)

cat("\nTR2 GT30D ...\n")
dat_tr2_gt <- prep_dat(ps_tr2_gt)
smd_tr2_gt <- compute_pairwise_smds(dat_tr2_gt)

cat("\nTR3 LE30D ...\n")
dat_tr3_le <- prep_dat(ps_tr3_le)
smd_tr3_le <- compute_pairwise_smds(dat_tr3_le)

cat("\nTR3 GT30D ...\n")
dat_tr3_gt <- prep_dat(ps_tr3_gt)
smd_tr3_gt <- compute_pairwise_smds(dat_tr3_gt)

# ── Save to Excel ────────────────────────────────────────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "TR2_LE30D_3-4d_within_30d");  writeData(wb, "TR2_LE30D_3-4d_within_30d",  smd_tr2_le)
addWorksheet(wb, "TR2_GT30D_3-4d_after_30d");   writeData(wb, "TR2_GT30D_3-4d_after_30d",   smd_tr2_gt)
addWorksheet(wb, "TR3_LE30D_5-14d_within_30d"); writeData(wb, "TR3_LE30D_5-14d_within_30d", smd_tr3_le)
addWorksheet(wb, "TR3_GT30D_5-14d_after_30d");  writeData(wb, "TR3_GT30D_5-14d_after_30d",  smd_tr3_gt)

header_style <- createStyle(textDecoration = "bold", halign = "center",
                            fgFill = "#D9D9D9", border = "TopBottom",
                            wrapText = TRUE)

# Highlight the "max" columns (bold + subtle yellow fill)
max_style <- createStyle(textDecoration = "bold", fgFill = "#FFF2CC")

for (sh in names(wb)) {
  n_col <- ncol(smd_tr2_le)
  addStyle(wb, sh, header_style, rows = 1, cols = 1:n_col, gridExpand = TRUE)
  # Max columns are 5 (Before) and 9 (After)
  addStyle(wb, sh, max_style, rows = 2:(nrow(smd_tr2_le) + 1),
           cols = c(5, 9), gridExpand = TRUE)
  setColWidths(wb, sh, cols = 1:n_col,
               widths = c(55, 18, 22, 18, 22, 18, 22, 18, 22))
  freezePane(wb, sh, firstActiveRow = 2, firstActiveCol = 2)
}

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("\n\u2713 SMDs saved: %s\n", OUT_XLSX))

