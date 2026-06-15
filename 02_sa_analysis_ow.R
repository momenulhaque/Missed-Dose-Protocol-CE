# Load Packages ##################
# Two trials - methadone interruption 3-4 and 5-14 days
#   first-interruption-in-any-episode - TT_analytic_prev_sens_nw_m data set
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
  prior_gap     = "Total OAT interruptions since initiation",
  prior_episode = "Total prior episode"
)

# 5. SMD LOVE PLOTS ############
# Each trial (3-4 days, 5-14 days) one cohort
#   - First-interruption-in-any-episode


cat("\n=== Creating SMD plots (2 trials × 1 cohorts = 2 figures) ===\n")

COMPARISON_COLORS <- c(
  "No dose change"              = "#000000", 
  "Dose increase"               = "#D62728",  
  "1-<50% reduction"            = "#2CA02C",  
  "≥50% reduction/≤40mg" = "#FF7F0E"  
)

# Global settings for love plots
x_low    <- -0.6
x_high   <- 0.6
x_breaks <- c(-0.4, -0.2, 0, 0.2, 0.4)

# Helper function to create one love plot
make_love_panel_v2 <- function(dat, bvars, var_labels, title_txt,
                               show_y_labels = TRUE) {
  
  fmla <- as.formula(paste("trt ~", paste(bvars, collapse = " + ")))
  
  prep <- bal.tab(
    fmla,
    data        = dat,
    weights     = dat$ow_wt,
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

## MET interruption of 3-4 days #############

cat("\nTrial 2 (3-4 days)...\n")


### First-interruption-in-any-episode #########
cat("  Creating SMD plots: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds")) 

dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = "≥50% reduction/≤40mg",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0),
    prior_episode = as.numeric(prior_episode)
  )

fmla <- as.formula(paste("trt ~", paste(c(baselinevars, "prior_episode"), collapse = " + ")))

# Panel 1
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$ow_wt,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_1 <- love.plot(prep_1, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "1-<50% reduction\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2
d_2 <- subset(dat, trt %in% c("No dose change", "≥50% reduction/≤40mg"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$ow_wt,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_2 <- love.plot(prep_2, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "≥50% reduction/≤40mg\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$ow_wt,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_3 <- love.plot(prep_3, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "Dose increase\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

row_tr2_all <- p_1 + plot_spacer() + p_2 + plot_spacer() + p_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

row_tr2_all <- plot_spacer() + row_tr2_all + plot_layout(heights = c(0.05, 1))

plot_tr2_all_final <- ggdraw(row_tr2_all) +
  draw_label("A. Methadone interruption of 3-4 days",
             x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  draw_label("Note: SMD ≥ 0.1 included in outcome models.",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/June09/", "smd_plots_tr2_ow.jpg"),
     width = 4200, height = 2200, res = 300)
print(plot_tr2_all_final)
dev.off()

cat("  ✓ Trial 2, First-interruption-in-any-episode SMD figure saved.\n")


## MET interruption of 5-14 days #############

cat("\nTrial 3 (5-14 days)...\n")


### First-interruption-in-any-episode #########
cat("  Creating SMD plots: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds")) 

dat <- ps_fit$data_weights |>
  mutate(
    trt = recode(trt,
                 "no reversion"       = "No dose change",
                 "1-<50%"             = "1-<50% reduction",
                 ">=50%/<=40mg"       = "≥50% reduction/≤40mg",
                 "increase"           = "Dose increase"),
    incar_grp_bl = ifelse(incar_grp_bl == "1_pastyear", 1, 0),
    prior_episode = as.numeric(prior_episode)
  )

fmla <- as.formula(paste("trt ~", paste(c(baselinevars, "prior_episode"), collapse = " + ")))

# Panel 1
d_1 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
d_1$trt <- droplevels(d_1$trt)
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$ow_wt,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_1 <- love.plot(prep_1, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "1-<50% reduction\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none")

# Panel 2
d_2 <- subset(dat, trt %in% c("No dose change", "≥50% reduction/≤40mg"))
d_2$trt <- droplevels(d_2$trt)
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$ow_wt,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_2 <- love.plot(prep_2, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "≥50% reduction/≤40mg\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

# Panel 3
d_3 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
d_3$trt <- droplevels(d_3$trt)
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$ow_wt,
                  estimand = "ATE", m.threshold = 0.1, s.d.denom = "pooled",
                  un = TRUE, stats = "mean.diffs")
p_3 <- love.plot(prep_3, stat = "mean.diffs", stars = "none",
                 thresholds = c(m = 0.1), abs = FALSE,
                 var.names = VAR_LABELS[names(VAR_LABELS) %in% c(baselinevars, "prior_episode")],
                 sample.names = c("Before", "After"), position = "none", line = FALSE,
                 colors = c("black", "red"), shapes = c(16, 16)) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = c(-0.4, -0.2, 0, 0.2, 0.4), expand = c(0, 0)) +
  labs(x = "", y = "", color = "", shape = "", title = "Dose increase\nvs\nno change") +
  theme(title = element_text(color = "black", size = 12),
        panel.grid = element_blank(), axis.line = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_blank(),
        legend.position = "none")

row_tr3_all <- p_1 + plot_spacer() + p_2 + plot_spacer() + p_3 +
  plot_layout(widths = c(1, 0.05, 1, 0.05, 1))

row_tr3_all <- plot_spacer() + row_tr3_all + plot_layout(heights = c(0.05, 1))

plot_tr3_all_final <- ggdraw(row_tr3_all) +
  draw_label("B. Methadone interruption of 5-14 days",
             x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  draw_label("Note: SMD ≥ 0.1 included in outcome models",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/June09/", "smd_plots_tr3_ow.jpg"),
     width = 4200, height = 2200, res = 300)
print(plot_tr3_all_final)
dev.off()

cat("  ✓ Trial 2, First-interruption-in-any-episode SMD figure saved.\n")






# 6. CUMULATIVE INCIDENCE ####
# Risk curves and confidence intervals for overdose and discontinuation
# Two trials × Two cohorts

cat("\n=== Computing Cumulative Incidence ===\n")

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- c(30, 60, 90, 180, 365)

# ── INITIALIZE RESULTS STRUCTURE ──

risk <- list(
  tr2 = list(
    first = list(od = list(), dis = list()),
    all = list(od = list(), dis = list())
  ),
  tr3 = list(
    first = list(od = list(), dis = list()),
    all = list(od = list(), dis = list())
  )
)

## MET interruption of 3-4 days #############

cat("\nTrial 2 (3-4 days)...\n")

### FIRST-INTERRUPTION-IN-FIRST-EPISODE ####
cat("  Computing: First-interruption-in-first-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_first.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)


risk$tr2$first$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

risk$tr2$first$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr2$first$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

risk$tr2$first$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

rm(long_dat); gc()

### FIRST-INTERRUPTION-IN-ANY-EPISODE ####
cat("  Computing: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$all$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "ow_wt"
)

risk$tr2$all$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "ow_wt"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)

risk$tr2$all$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs =c("ddose_prev", "prior_episode"),  
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "ow_wt"
)

risk$tr2$all$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "ow_wt"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/OW/", "risk_tr2_34d.rds"))

cat("  ✓ Trial 2 risk estimates saved.\n")

## MET interruption of 5-14 days #############

cat("\nTrial 3 (5-14 days)...\n")

### FIRST-INTERRUPTION-IN-FIRST-EPISODE ####
cat("  Computing: First-interruption-in-first-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_first.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$first$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

risk$tr3$first$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$first$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt"
)


risk$tr3$first$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

rm(long_dat); gc()

### FIRST-INTERRUPTION-IN-ANY-EPISODE ####
cat("  Computing: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$all$od$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

risk$tr3$all$od$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$all$dis$curve <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

risk$tr3$all$dis$ci <- predict_cuminc(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "ow_wt"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/OW/", "risk_tr3_514d.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")

