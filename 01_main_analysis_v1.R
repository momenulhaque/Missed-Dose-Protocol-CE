# Load Packages ##################
# 01_main_analysis.R
# Trials 2 and 3 - methadone interruption 3-4d and 5-14d
# Two cohort designs:
#   first_interruption of episode  (formerly removed)    - TT_analytic_prev_sens_nw_m data set
#   all_interruptions of episode  (formerly prevalent new users)     - TT_analytic_rprev_sens_nw_m data set
###

library(tidyverse)
library(haven)
library(nnet)
library(splines)
library(sandwich)
library(mvtnorm)
library(tableone)
library(survey)
library(cobalt)
library(patchwork)
library(cowplot)

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")   # ps_model_wt(), predict_cuminc(), expand_person_period()

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results_v1/"
dir.create(paste0(RES_DIR, "PScore"),   showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Risk"),     showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Tables"),   showWarnings = FALSE)
dir.create(paste0(RES_DIR, "Figures"),  showWarnings = FALSE)


# ── Shared settings ────────────────────────────────────────────────────────────
TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")

out_vars   <- c("sens_discont", "time_discont",
                "sens_od_acute", "time_od_acute",
                "sens_all_acute", "time_all_acute")

clust_vars <- c("studyid_rev", "episode", "time0")

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss"
)
continuous_vars <- c("ddose_prev", "year", "ret_length")

# Human-readable variable labels for Table 1 and love plots
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
  ret_length    = "Duration of OAT episode before time zero",
  home_grp_b    = "Unstable housing in past 5 years",
  incar_grp_bl  = "Incarcerated in the past year",
  drug_acute    = "Drug-related acute care visits in the past month",
  cci_gt0       = "Charlson Comorbidity Index (>0)",
  copd_asthma   = "Asthma or COPD in the past year",
  sed_psy_med   = "Dispensations of psychiatric or sedative medication in the past month",
  po            = "Dispensations opioid other than OAT in the past month",
  sud           = "Other substance use disorder",
  aud           = "Alcohol use disorder",
  smd           = "Serious mental disorder",
  hcv           = "Hepatitis C",
  tob           = "Tobacco use disorder (past year)",
  udt           = "Urine drug test in the past month",
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_gap     = "Total OAT interruptions since initiation"
)

# Shared pre-processing: collapse all reduction categories into one
clean_trt <- function(dat) {
  dat |>
    mutate(
      trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)
    ) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS))
}

clean_trt_all <- function(dat) {
  dat |>
    mutate(
      trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)
    ) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt = factor(trt, levels = TRT_LEVELS)) |>
    group_by(studyid_rev, episode) |>
    mutate(
      prior_gap     = cumsum(gap_days),
      prior_gap_num = row_number()
    ) %>%
    ungroup()
}

# 1. LOAD AND FILTER RAW SAS DATA #########

# First-interruption in life
raw_first <- read_sas(file.path(DATA_DIR, "TT_analytic_incid_sens_nw_m.sas7bdat"))

dat_tr2_first <- raw_first |>
  filter(trial_gap %in% 3:4) |>
  clean_trt()

dat_tr3_first <- raw_first |>
  filter(trial_cat == "5-14") |>
  clean_trt()

rm(raw_first); gc()

# First-interruption in episode
raw_first_ep <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_first_ep <- raw_first_ep |>
  filter(trial_gap %in% 3:4) |>
  clean_trt()

dat_tr3_first_ep <- raw_first_ep |>
  filter(trial_cat == "5-14") |>
  clean_trt()

rm(raw_first_ep); gc()

# All-interruptions cohort 
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_rprev_sens_nw_m.sas7bdat"))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all()

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all()

rm(raw_all); gc()

# 2. PROPENSITY SCORE WEIGHTS ################

## MET 3-4 days: first-interruption in life #####

ps_tr2_first <- ps_model_wt(
  dat             = dat_tr2_first,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_first, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"))



## MET 3-4 days: first-interruption in episode #####

ps_tr2_first_ep <- ps_model_wt(
  dat             = dat_tr2_first_ep,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_first_ep, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ep_ow.rds"))


## MET 3-4 days: all-interruptions #####
ps_tr2_all <- ps_model_wt(
  dat             = dat_tr2_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_gap"),
  continuous_vars = c(continuous_vars, "prior_gap"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)
saveRDS(ps_tr2_all, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"))




## MET 5-14 days: first-interruption in life #####
ps_tr3_first <- ps_model_wt(
  dat             = dat_tr3_first,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)
saveRDS(ps_tr3_first, paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"))

## MET 5-14 days: first-interruption in episode #####
ps_tr3_first_ep <- ps_model_wt(
  dat             = dat_tr3_first_ep,
  exposure        = "trt",
  baselinevars    = baselinevars,
  continuous_vars = continuous_vars,
  out_vars        = out_vars,
  clust_vars      = clust_vars
)
saveRDS(ps_tr3_first_ep, paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ep_ow.rds"))


## MET 5-14 days: all-interruptions #####
ps_tr3_all <- ps_model_wt(
  dat             = dat_tr3_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_gap"),
  continuous_vars = c(continuous_vars, "prior_gap"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)
saveRDS(ps_tr3_all, paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"))




# 3. TABLE 1 ##############
# Produces one table per trial × cohort.
# Rows: covariates. Columns: treatment arms (unweighted) + SMD before/after weighting.


tbl_tr2_first <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"),
  bvars = baselinevars, cont_vars = continuous_vars, clustered = FALSE
)

tbl_tr2_first_ep <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ep_ow.rds"),
  bvars = baselinevars, cont_vars = continuous_vars, clustered = FALSE
)

tbl_tr2_all <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"),
  bvars = c(baselinevars, "prior_gap"),
  cont_vars = c(continuous_vars, "prior_gap"), clustered = TRUE
)



tbl_tr3_first <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"),
  bvars = baselinevars, cont_vars = continuous_vars, clustered = FALSE
)

tbl_tr3_first_ep <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ep_ow.rds"),
  bvars = baselinevars, cont_vars = continuous_vars, clustered = FALSE
)

tbl_tr3_all <- make_table1(
  paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"),
  bvars = c(baselinevars, "prior_gap"),
  cont_vars = c(continuous_vars, "prior_gap"), clustered = TRUE
)


library(openxlsx)

wb <- createWorkbook()

for (nm in c("tr2_first", "tr2_first_ep",  "tr2_all", "tr3_first_ep", "tr3_all")) {
  tbl <- get(paste0("tbl_", nm))
  addWorksheet(wb, sheetName = nm)
  writeData(wb, sheet = nm, x = as.data.frame(tbl), rowNames = TRUE)
}

saveWorkbook(wb, paste0(RES_DIR, "Tables/", "table1_all.xlsx"), overwrite = TRUE)

# 4. PS DISTRIBUTION PLOTS ####
# Histogram + density of propensity scores by treatment arm,
# faceted by interruption length × cohort design.

COMPARISON_COLORS <- c(
  "No dose change"              = "#000000", 
   "Dose increase"               = "#D62728",  
  "1-<50% reduction"            = "#2CA02C",  
  ">=50% reduction/initial dose*"  = "#FF7F0E"  
)


ps_score_plot <- function(dat,
                          x_lab           = "Propensity Score",
                          y_lab           = "Density",
                          legend_position = "bottom",
                          x_points        = c(0, .25, .5, .75, 1)) {
  
  ggplot(dat, aes(x = ps, colour = group, fill = group)) +
    geom_histogram(aes(y = ..density..),
                   position = "identity", alpha = 0.25, bins = 30) +
    geom_density(alpha = 0.05, size = 0.5) +
    scale_colour_manual(values = COMPARISON_COLORS) +
    scale_fill_manual(values   = COMPARISON_COLORS) +
    facet_wrap(~Trial_design, nrow = 2, scales = "free") +
    labs(x = x_lab, y = y_lab) +
    scale_x_continuous(breaks = x_points) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks.y      = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x      = element_line(linewidth = 1.2, color = "black"),
      axis.text.x       = element_text(color = "black", size = 14),
      axis.text.y       = element_text(color = "black", size = 14),
      axis.title.y      = element_text(size = 14),
      legend.position   = legend_position,
      legend.direction  = "horizontal",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 14, color = "black"),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key        = element_blank(),
      strip.background  = element_rect(fill = "lightgray", color = NA),
      strip.text        = element_text(face = "bold", size = 13)
    )
}

prep_ps_dat <- function(ps_rds, trial_label, design_label) {
  readRDS(ps_rds)$data_weights |>
    mutate(
      group = recode(as.character(trt),
                     "no reversion" = "No dose change",
                     "1-<50%"    = "1-<50% reduction",
                     ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                     "increase"     = "Dose increase"),
      group  = factor(group, levels = c("No dose change", "1-<50% reduction", ">=50% reduction/initial dose*", "Dose increase")),
      Trial  = trial_label,
      Design = design_label
    ) |>
    select(group, ps, Trial, Design)
}


dat_ps <- bind_rows(
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"),
              "MET interruption of 3-4 days",  "First-interruption"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ep_ow.rds"),
              "MET interruption of 3-4 days",  "First-interruption in episode"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"),
              "MET interruption of 3-4 days",  "All-interruptions"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"),
              "MET interruption of 5-14 days", "First-interruption"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ep_ow.rds"),
              "MET interruption of 5-14 days", "First-interruption in episode"),
  prep_ps_dat(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"),
              "MET interruption of 5-14 days", "All-interruptions cohort")
) |>
  mutate(
    Trial_design = factor(
      paste0(Trial, "\n", Design),
      levels = c(
        "MET interruption of 3-4 days\nFirst-interruption",
        "MET interruption of 3-4 days\nFirst-interruption in episode",
        "MET interruption of 3-4 days\nAll-interruptions",
        "MET interruption of 5-14 days\nFirst-interruption",
        "MET interruption of 5-14 days\nFirst-interruption in episode",
        "MET interruption of 5-14 days\nAll-interruptions"
      )
    )
  )


plot_ps <- ps_score_plot(dat_ps) +
  guides(color = guide_legend(nrow = 1))

jpeg(paste0(RES_DIR, "Figures/", "ps_distribution.jpg"),
     width = 5500, height = 3000, res = 300)
plot_ps
dev.off()


# 5. SMD LOVE PLOTS ############
cat("\n=== Creating SMD plots (6 separate figures, one per cohort) ===\n")

# Global settings for love plots
x_low    <- -0.6
x_high   <- 0.6
x_breaks <- c(-0.4, -0.2, 0, 0.2, 0.4)

# Helper function to create one love plot with flexible y-axis labels
make_love_panel_v2 <- function(dat, bvars, var_labels, title_txt,
                               show_y_labels = TRUE) {
  
  fmla <- as.formula(paste("trt ~", paste(bvars, collapse = " + ")))
  
  prep <- bal.tab(
    fmla,
    data        = dat,
    weights     = dat$wt_trunc,
    estimand    = "ATE",
    m.threshold = 0.1,
    s.d.denom   = "pooled",
    un          = TRUE,
    stats       = "mean.diffs"
  )
  
  n_vars  <- nrow(prep$Balance)
  bg_data <- data.frame(
    var_index  = seq_len(n_vars),
    fill_color = rep(c("gray65", "white"), length.out = n_vars)
  )
  
  plt <- love.plot(
    prep,
    stat         = "mean.diffs",
    stars        = "none",
    thresholds   = c(m = 0.1),
    abs          = FALSE,
    var.names    = var_labels,
    sample.names = c("Before", "After"),
    position     = "none",
    line         = FALSE,
    colors       = c("black", "red"),
    shapes       = c(16, 16),
    title        = title_txt,
    size         = 4
  ) +
    geom_rect(
      data        = bg_data,
      aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
          xmin = -Inf, xmax = Inf, fill = fill_color),
      inherit.aes = FALSE,
      alpha       = 0.2
    ) +
    scale_x_continuous(limits = c(x_low, x_high), breaks = x_breaks, expand = c(0, 0)) +
    scale_fill_identity() +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    labs(x = "SMD", y = "", color = "", shape = "") +
    theme_minimal(base_size = 12) +
    theme(
      title              = element_text(color = "black", size = 13, face = "bold"),
      panel.grid         = element_blank(),
      axis.line          = element_line(color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.text.x        = element_text(color = "black", size = 11),
      axis.ticks.y       = element_line(linewidth = 1, color = "black"),
      axis.text.y        = if (show_y_labels)
        element_text(color = "black", size = 10)
      else
        element_blank(),
      axis.title.x       = element_text(size = 11),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(size = 10, color = "black"),
      legend.key.size    = unit(0.3, "cm"),
      plot.margin        = margin(5, 5, 5, 5, "pt")
    )
  
  plt
}

# Helper: Create one 3-panel row (1-<50% | >=50% | increase) for one trial × cohort
make_love_row_v2 <- function(ps_rds, bvars, var_labels, show_y_labels = TRUE,
                             trt_recode = c(
                               "no reversion"       = "No dose change",
                               "1-<50%"             = "1-<50%",
                               ">=50%/<=40mg"       = ">=50%*",
                               "increase"           = "Increase")) {
  
  dat <- readRDS(ps_rds)$data_weights |>
    mutate(
      trt          = recode(as.character(trt), !!!trt_recode),
      incar_grp_bl = if_else(incar_grp_bl == "1_pastyear", 1L, 0L)
    )
  
  # 1-<50% vs no change
  d_1 <- subset(dat, trt %in% c("No dose change", "1-<50%"))
  d_1$trt <- droplevels(factor(d_1$trt))
  p_1 <- make_love_panel_v2(d_1, bvars, var_labels,
                            title_txt = "1-<50% vs no change",
                            show_y_labels = show_y_labels)
  
  # >=50% vs no change
  d_2 <- subset(dat, trt %in% c("No dose change", ">=50%*"))
  d_2$trt <- droplevels(factor(d_2$trt))
  p_2 <- make_love_panel_v2(d_2, bvars, var_labels,
                            title_txt = ">=50%* vs no change",
                            show_y_labels = FALSE)  # Never show y-labels for middle/right
  
  # Increase vs no change
  d_3 <- subset(dat, trt %in% c("No dose change", "Increase"))
  d_3$trt <- droplevels(factor(d_3$trt))
  p_3 <- make_love_panel_v2(d_3, bvars, var_labels,
                            title_txt = "Increase vs no change",
                            show_y_labels = FALSE)  # Never show y-labels for right panels
  
  # Stack the three panels horizontally
  p_1 + p_2 + p_3 + plot_layout(widths = c(1.2, 1, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.text = element_text(size = 9))
}

# ── BUILD 2×3 GRID: 2 Columns (Trials) × 3 Rows (Cohorts) ──

cat("\n  Building Trial 2 rows...\n")

# TRIAL 2 - Row 1: First-in-life
row_tr2_first <- make_love_row_v2(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"),
  bvars       = baselinevars,
  var_labels  = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
  show_y_labels = TRUE  # Left column: show labels
)

# TRIAL 2 - Row 2: First-in-episode (safety check)
row_tr2_first_ep <- make_love_row_v2(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ep_ow.rds"),
  bvars       = baselinevars,
  var_labels  = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
  show_y_labels = TRUE
)

# TRIAL 2 - Row 3: All-interruptions
row_tr2_all <- make_love_row_v2(
  paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"),
  bvars      = c(baselinevars, "prior_gap"),
  var_labels = VAR_LABELS,
  show_y_labels = TRUE
)

# ── CREATE INDIVIDUAL COMPARISON PANELS FOR EACH COHORT ──
# Pattern: 3 horizontal panels (1-<50% | >=50% | increase) with spacers between
# First panel shows y-axis labels, others blank

cat("\n  Creating comparison panels...\n")

# ═══════════════════════════════════════════════════════════════════════════════
# TRIAL 2 (3-4 DAYS)
# ═════════════════════════════════════════════════════════════════════════════

# ── FIRST-IN-LIFE ──
cat("    Trial 2, First-in-life...\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

fmla <- as.formula(paste("trt ~", paste(baselinevars, collapse = " + ")))

# Panel 1: 1-<50% vs no change (WITH y-axis labels)
d_tr2_first_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_tr2_first_1$trt <- droplevels(d_tr2_first_1$trt)
prep_tr2_first_1 <- bal.tab(fmla, data = d_tr2_first_1, weights = d_tr2_first_1$wt_trunc,
                            estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                            un = TRUE, stats = "mean.diffs")
p_tr2_first_1 <- love.plot(prep_tr2_first_1, stat = "mean.diffs", stars = "none",
                           thresholds = c(m = 0.1), abs = FALSE,
                           var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                           sample.names = c("Before", "After"), position = "none", line = FALSE,
                           colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2: >=50% vs no change (NO y-axis labels)
d_tr2_first_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_tr2_first_2$trt <- droplevels(d_tr2_first_2$trt)
prep_tr2_first_2 <- bal.tab(fmla, data = d_tr2_first_2, weights = d_tr2_first_2$wt_trunc,
                            estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                            un = TRUE, stats = "mean.diffs")
p_tr2_first_2 <- love.plot(prep_tr2_first_2, stat = "mean.diffs", stars = "none",
                           thresholds = c(m = 0.1), abs = FALSE,
                           var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                           sample.names = c("Before", "After"), position = "none", line = FALSE,
                           colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3: Increase vs no change (NO y-axis labels)
d_tr2_first_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_tr2_first_3$trt <- droplevels(d_tr2_first_3$trt)
prep_tr2_first_3 <- bal.tab(fmla, data = d_tr2_first_3, weights = d_tr2_first_3$wt_trunc,
                            estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                            un = TRUE, stats = "mean.diffs")
p_tr2_first_3 <- love.plot(prep_tr2_first_3, stat = "mean.diffs", stars = "none",
                           thresholds = c(m = 0.1), abs = FALSE,
                           var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                           sample.names = c("Before", "After"), position = "none", line = FALSE,
                           colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Stack panels horizontally for first-in-life
row_tr2_first <- p_tr2_first_1 + plot_spacer() + p_tr2_first_2 + plot_spacer() + p_tr2_first_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

# ── FIRST-IN-EPISODE (safety check) ──
cat("    Trial 2, First-in-episode (safety)...\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ep_ow.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

# Panel 1: (WITH y-axis labels)
d_tr2_ep_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_tr2_ep_1$trt <- droplevels(d_tr2_ep_1$trt)
prep_tr2_ep_1 <- bal.tab(fmla, data = d_tr2_ep_1, weights = d_tr2_ep_1$wt_trunc,
                         estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                         un = TRUE, stats = "mean.diffs")
p_tr2_ep_1 <- love.plot(prep_tr2_ep_1, stat = "mean.diffs", stars = "none",
                        thresholds = c(m = 0.1), abs = FALSE,
                        var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                        sample.names = c("Before", "After"), position = "none", line = FALSE,
                        colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2: (NO y-axis labels)
d_tr2_ep_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_tr2_ep_2$trt <- droplevels(d_tr2_ep_2$trt)
prep_tr2_ep_2 <- bal.tab(fmla, data = d_tr2_ep_2, weights = d_tr2_ep_2$wt_trunc,
                         estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                         un = TRUE, stats = "mean.diffs")
p_tr2_ep_2 <- love.plot(prep_tr2_ep_2, stat = "mean.diffs", stars = "none",
                        thresholds = c(m = 0.1), abs = FALSE,
                        var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                        sample.names = c("Before", "After"), position = "none", line = FALSE,
                        colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3: (NO y-axis labels)
d_tr2_ep_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_tr2_ep_3$trt <- droplevels(d_tr2_ep_3$trt)
prep_tr2_ep_3 <- bal.tab(fmla, data = d_tr2_ep_3, weights = d_tr2_ep_3$wt_trunc,
                         estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                         un = TRUE, stats = "mean.diffs")
p_tr2_ep_3 <- love.plot(prep_tr2_ep_3, stat = "mean.diffs", stars = "none",
                        thresholds = c(m = 0.1), abs = FALSE,
                        var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                        sample.names = c("Before", "After"), position = "none", line = FALSE,
                        colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Stack panels horizontally for first-in-episode
row_tr2_ep <- p_tr2_ep_1 + plot_spacer() + p_tr2_ep_2 + plot_spacer() + p_tr2_ep_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

# ── ALL-INTERRUPTIONS ──
cat("    Trial 2, All-interruptions...\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

fmla_all <- as.formula(paste("trt ~", paste(c(baselinevars, "prior_gap"), collapse = " + ")))

# Panel 1: (WITH y-axis labels)
d_tr2_all_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_tr2_all_1$trt <- droplevels(d_tr2_all_1$trt)
prep_tr2_all_1 <- bal.tab(fmla_all, data = d_tr2_all_1, weights = d_tr2_all_1$wt_trunc,
                          estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                          un = TRUE, stats = "mean.diffs")
p_tr2_all_1 <- love.plot(prep_tr2_all_1, stat = "mean.diffs", stars = "none",
                         thresholds = c(m = 0.1), abs = FALSE,
                         var.names = VAR_LABELS,
                         sample.names = c("Before", "After"), position = "none", line = FALSE,
                         colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "bottom", legend.box = "horizontal")

# Panel 2: (NO y-axis labels)
d_tr2_all_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_tr2_all_2$trt <- droplevels(d_tr2_all_2$trt)
prep_tr2_all_2 <- bal.tab(fmla_all, data = d_tr2_all_2, weights = d_tr2_all_2$wt_trunc,
                          estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                          un = TRUE, stats = "mean.diffs")
p_tr2_all_2 <- love.plot(prep_tr2_all_2, stat = "mean.diffs", stars = "none",
                         thresholds = c(m = 0.1), abs = FALSE,
                         var.names = VAR_LABELS,
                         sample.names = c("Before", "After"), position = "none", line = FALSE,
                         colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3: (NO y-axis labels)
d_tr2_all_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_tr2_all_3$trt <- droplevels(d_tr2_all_3$trt)
prep_tr2_all_3 <- bal.tab(fmla_all, data = d_tr2_all_3, weights = d_tr2_all_3$wt_trunc,
                          estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                          un = TRUE, stats = "mean.diffs")
p_tr2_all_3 <- love.plot(prep_tr2_all_3, stat = "mean.diffs", stars = "none",
                         thresholds = c(m = 0.1), abs = FALSE,
                         var.names = VAR_LABELS,
                         sample.names = c("Before", "After"), position = "bottom", line = FALSE,
                         colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "bottom", legend.box = "horizontal")

# Stack panels horizontally for all-interruptions
row_tr2_all <- p_tr2_all_1 + plot_spacer() + p_tr2_all_2 + plot_spacer() + p_tr2_all_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1), guides = "collect")

# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE TRIAL 2 FIGURE
# ═════════════════════════════════════════════════════════════════════════════

plot_tr2_combined <- plot_spacer() /
  row_tr2_first /
  plot_spacer() /
  row_tr2_ep /
  plot_spacer() /
  row_tr2_all +
  plot_layout(heights = c(0.04, 1, 0.02, 1, 0.02, 1)) &
  theme(legend.position = "bottom", legend.box = "horizontal")

plot_tr2_final <- ggdraw(plot_tr2_combined) +
  draw_label("A. Methadone interruption of 3-4 days: First-in-life",
             x = 0.02, y = 0.97, hjust = 0, vjust = 0, size = 14) +
  draw_label("B. Methadone interruption of 3-4 days: First-in-episode (safety check)",
             x = 0.02, y = 0.63, hjust = 0, vjust = 0, size = 14) +
  draw_label("C. Methadone interruption of 3-4 days: All-interruptions",
             x = 0.02, y = 0.31, hjust = 0, vjust = 0, size = 14) +
  draw_label("Note: SMD ≥ 0.1 included in outcome models. * ≤40mg methadone",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/", "smd_love_plots_tr2.jpg"),
     width = 5200, height = 6800, res = 300)
print(plot_tr2_final)
dev.off()

cat("  ✓ Trial 2 SMD figure saved.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# TRIAL 3 (5-14 DAYS)
# ═════════════════════════════════════════════════════════════════════════════

# ── FIRST-IN-LIFE ──
cat("    Trial 3, First-in-life...\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

fmla <- as.formula(paste("trt ~", paste(baselinevars, collapse = " + ")))

# Panel 1: (WITH y-axis labels)
d_tr3_first_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_tr3_first_1$trt <- droplevels(d_tr3_first_1$trt)
prep_tr3_first_1 <- bal.tab(fmla, data = d_tr3_first_1, weights = d_tr3_first_1$wt_trunc,
                            estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                            un = TRUE, stats = "mean.diffs")
p_tr3_first_1 <- love.plot(prep_tr3_first_1, stat = "mean.diffs", stars = "none",
                           thresholds = c(m = 0.1), abs = FALSE,
                           var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                           sample.names = c("Before", "After"), position = "none", line = FALSE,
                           colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2: (NO y-axis labels)
d_tr3_first_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_tr3_first_2$trt <- droplevels(d_tr3_first_2$trt)
prep_tr3_first_2 <- bal.tab(fmla, data = d_tr3_first_2, weights = d_tr3_first_2$wt_trunc,
                            estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                            un = TRUE, stats = "mean.diffs")
p_tr3_first_2 <- love.plot(prep_tr3_first_2, stat = "mean.diffs", stars = "none",
                           thresholds = c(m = 0.1), abs = FALSE,
                           var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                           sample.names = c("Before", "After"), position = "none", line = FALSE,
                           colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3: (NO y-axis labels)
d_tr3_first_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_tr3_first_3$trt <- droplevels(d_tr3_first_3$trt)
prep_tr3_first_3 <- bal.tab(fmla, data = d_tr3_first_3, weights = d_tr3_first_3$wt_trunc,
                            estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                            un = TRUE, stats = "mean.diffs")
p_tr3_first_3 <- love.plot(prep_tr3_first_3, stat = "mean.diffs", stars = "none",
                           thresholds = c(m = 0.1), abs = FALSE,
                           var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                           sample.names = c("Before", "After"), position = "none", line = FALSE,
                           colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

row_tr3_first <- p_tr3_first_1 + plot_spacer() + p_tr3_first_2 + plot_spacer() + p_tr3_first_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

# ── FIRST-IN-EPISODE (safety check) ──
cat("    Trial 3, First-in-episode (safety)...\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ep_ow.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

# Panel 1: (WITH y-axis labels)
d_tr3_ep_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_tr3_ep_1$trt <- droplevels(d_tr3_ep_1$trt)
prep_tr3_ep_1 <- bal.tab(fmla, data = d_tr3_ep_1, weights = d_tr3_ep_1$wt_trunc,
                         estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                         un = TRUE, stats = "mean.diffs")
p_tr3_ep_1 <- love.plot(prep_tr3_ep_1, stat = "mean.diffs", stars = "none",
                        thresholds = c(m = 0.1), abs = FALSE,
                        var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                        sample.names = c("Before", "After"), position = "none", line = FALSE,
                        colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2: (NO y-axis labels)
d_tr3_ep_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_tr3_ep_2$trt <- droplevels(d_tr3_ep_2$trt)
prep_tr3_ep_2 <- bal.tab(fmla, data = d_tr3_ep_2, weights = d_tr3_ep_2$wt_trunc,
                         estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                         un = TRUE, stats = "mean.diffs")
p_tr3_ep_2 <- love.plot(prep_tr3_ep_2, stat = "mean.diffs", stars = "none",
                        thresholds = c(m = 0.1), abs = FALSE,
                        var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                        sample.names = c("Before", "After"), position = "none", line = FALSE,
                        colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3: (NO y-axis labels)
d_tr3_ep_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_tr3_ep_3$trt <- droplevels(d_tr3_ep_3$trt)
prep_tr3_ep_3 <- bal.tab(fmla, data = d_tr3_ep_3, weights = d_tr3_ep_3$wt_trunc,
                         estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                         un = TRUE, stats = "mean.diffs")
p_tr3_ep_3 <- love.plot(prep_tr3_ep_3, stat = "mean.diffs", stars = "none",
                        thresholds = c(m = 0.1), abs = FALSE,
                        var.names = VAR_LABELS[names(VAR_LABELS) != "prior_gap"],
                        sample.names = c("Before", "After"), position = "none", line = FALSE,
                        colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

row_tr3_ep <- p_tr3_ep_1 + plot_spacer() + p_tr3_ep_2 + plot_spacer() + p_tr3_ep_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

# ── ALL-INTERRUPTIONS ──
cat("    Trial 3, All-interruptions...\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"))
dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = ">=50% reduction/initial dose*",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0)
  )

fmla_all <- as.formula(paste("trt ~", paste(c(baselinevars, "prior_gap"), collapse = " + ")))

# Panel 1: (WITH y-axis labels)
d_tr3_all_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_tr3_all_1$trt <- droplevels(d_tr3_all_1$trt)
prep_tr3_all_1 <- bal.tab(fmla_all, data = d_tr3_all_1, weights = d_tr3_all_1$wt_trunc,
                          estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                          un = TRUE, stats = "mean.diffs")
p_tr3_all_1 <- love.plot(prep_tr3_all_1, stat = "mean.diffs", stars = "none",
                         thresholds = c(m = 0.1), abs = FALSE,
                         var.names = VAR_LABELS,
                         sample.names = c("Before", "After"), position = "none", line = FALSE,
                         colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "bottom", legend.box = "horizontal")

# Panel 2: (NO y-axis labels)
d_tr3_all_2 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
d_tr3_all_2$trt <- droplevels(d_tr3_all_2$trt)
prep_tr3_all_2 <- bal.tab(fmla_all, data = d_tr3_all_2, weights = d_tr3_all_2$wt_trunc,
                          estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                          un = TRUE, stats = "mean.diffs")
p_tr3_all_2 <- love.plot(prep_tr3_all_2, stat = "mean.diffs", stars = "none",
                         thresholds = c(m = 0.1), abs = FALSE,
                         var.names = VAR_LABELS,
                         sample.names = c("Before", "After"), position = "none", line = FALSE,
                         colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3: (NO y-axis labels)
d_tr3_all_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_tr3_all_3$trt <- droplevels(d_tr3_all_3$trt)
prep_tr3_all_3 <- bal.tab(fmla_all, data = d_tr3_all_3, weights = d_tr3_all_3$wt_trunc,
                          estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                          un = TRUE, stats = "mean.diffs")
p_tr3_all_3 <- love.plot(prep_tr3_all_3, stat = "mean.diffs", stars = "none",
                         thresholds = c(m = 0.1), abs = FALSE,
                         var.names = VAR_LABELS,
                         sample.names = c("Before", "After"), position = "bottom", line = FALSE,
                         colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "") +
  theme(title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "bottom", legend.box = "horizontal")

row_tr3_all <- p_tr3_all_1 + plot_spacer() + p_tr3_all_2 + plot_spacer() + p_tr3_all_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1), guides = "collect")

# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE TRIAL 3 FIGURE
# ═════════════════════════════════════════════════════════════════════════════

plot_tr3_combined <- plot_spacer() /
  row_tr3_first /
  plot_spacer() /
  row_tr3_ep /
  plot_spacer() /
  row_tr3_all +
  plot_layout(heights = c(0.04, 1, 0.02, 1, 0.02, 1)) &
  theme(legend.position = "bottom", legend.box = "horizontal")

plot_tr3_final <- ggdraw(plot_tr3_combined) +
  draw_label("A. Methadone interruption of 5-14 days: First-in-life",
             x = 0.02, y = 0.97, hjust = 0, vjust = 0, size = 14) +
  draw_label("B. Methadone interruption of 5-14 days: First-in-episode (safety check)",
             x = 0.02, y = 0.63, hjust = 0, vjust = 0, size = 14) +
  draw_label("C. Methadone interruption of 5-14 days: All-interruptions",
             x = 0.02, y = 0.31, hjust = 0, vjust = 0, size = 14) +
  draw_label("Note: SMD ≥ 0.1 included in outcome models. * ≤40mg methadone",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/", "smd_love_plots_tr3.jpg"),
     width = 5200, height = 6800, res = 300)
print(plot_tr3_final)
dev.off()

cat("  ✓ Trial 3 SMD figure saved.\n")


# 6. CUMULATIVE INCIDENCE ####
HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- c(30, 60, 90, 180, 365)

risk <- list(
  tr2 = list(first = list(od = list(), dis = list()),
             first_ep = list(od = list(), dis = list()),
             all   = list(od = list(), dis = list())),
  tr3 = list(first = list(od = list(), dis = list()),
             first_ep = list(od = list(), dis = list()),
             all   = list(od = list(), dis = list()))
)


## MET 3-4 days / first-interruption ######
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ow.rds"))

## Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr2$first$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr2$first$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

## Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr2$first$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr2$first$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = NULL,
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()


## MET 3-4 days / first-in-episode ######

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first_ep_ow.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$first_ep$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = "ddose_prev", 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 199, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)
risk$tr2$first_ep$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = NULL, 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 199, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)

risk$tr2$first_ep$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs = NULL, 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 199, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)
risk$tr2$first_ep$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = NULL, 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 199, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)
rm(long_dat); gc()



## MET 3-4 days / all-interruptions #######
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all_ow.rds"))

## Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr2$all$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt_trunc"
)

risk$tr2$all$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt_trunc"
)

## Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr2$all$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = "ddose_prev",
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr2$all$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = "ddose_prev",
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_34d.rds"))



## MET 5-14 days / first-interruption ######
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ow.rds"))

## Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$first$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$first$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

## Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$first$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$first$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

## MET 5-14 days / first-in-episode ######
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first_ep_ow.rds"))

## Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$first_ep$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$first_ep$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

## Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 7)

risk$tr3$first_ep$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$first_ep$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)


##  MET 5-14 days / all-interruptions ######
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all_ow.rds"))

## Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$all$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "prior_gap"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$all$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "prior_gap"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

## Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 7)

risk$tr3$all$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "prior_gap"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 199, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$all$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "prior_gap"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 199, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3.rds"))




