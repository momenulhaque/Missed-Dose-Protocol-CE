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
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/GT30D/"
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
  prior_gap     = "Total OAT interruptions since initiation",
  prior_episode = "Number of prior episode"
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
    mutate(trt = factor(trt, levels = TRT_LEVELS),
           prior_episode = as.numeric(episode - 1)) |>
    group_by(studyid_rev, episode) |>
    mutate(
      prior_gap     = cumsum(gap_days),
      prior_gap_num = row_number()
    ) %>%
    ungroup()
  
}

# 1. LOAD AND FILTER RAW SAS DATA #########

# First-interruption-in-first-episode
raw_first <- read_sas(file.path(DATA_DIR, "TT_analytic_incid_sens_nw_m.sas7bdat"))

dat_tr2_first <- raw_first |>
  filter(trial_gap %in% 3:4) |>
  clean_trt() |>
  filter(ret_length > 30)

dat_tr3_first <- raw_first |>
  filter(trial_cat == "5-14") |>
  clean_trt() |>
  filter(ret_length > 30)

rm(raw_first); gc()

# First-interruption-in-any-episode 
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_all <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all() |>
  filter(ret_length > 30)

dat_tr3_all <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all() |>
  filter(ret_length > 30)


rm(raw_all); gc()


# 2. PROPENSITY SCORE WEIGHTS ################
## MET 3-4 days: first-interruption-in-any-episode #####

ps_tr2_all <- ps_model_wt(
  dat             = dat_tr2_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_episode"),
  continuous_vars = c(continuous_vars, "prior_episode"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr2_all, paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))

## MET 5-14 days: first-interruption-in-any-episode #####
ps_tr3_all <- ps_model_wt(
  dat             = dat_tr3_all,
  exposure        = "trt",
  baselinevars    = c(baselinevars, "prior_episode"),
  continuous_vars = c(continuous_vars, "prior_episode"),
  out_vars        = out_vars,
  clust_vars      = clust_vars
)

saveRDS(ps_tr3_all, paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))



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
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$wt_trunc,
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
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$wt_trunc,
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
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$wt_trunc,
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
  draw_label("Note: SMD ≥ 0.1 included in outcome models.",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/", "smd_love_plots_tr2_all.jpg"),
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
prep_1 <- bal.tab(fmla, data = d_1, weights = d_1$wt_trunc,
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
prep_2 <- bal.tab(fmla, data = d_2, weights = d_2$wt_trunc,
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
prep_3 <- bal.tab(fmla, data = d_3, weights = d_3$wt_trunc,
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
  draw_label("Note: SMD ≥ 0.1 included in outcome models",
             x = 0.5, y = 0.01, hjust = 0.5, vjust = 0, size = 11)

jpeg(paste0(RES_DIR, "Figures/", "smd_love_plots_tr3_all.jpg"),
     width = 4200, height = 2200, res = 300)
print(plot_tr3_all_final)
dev.off()

cat("  ✓ Trial 2, First-interruption-in-any-episode SMD figure saved.\n")






# 7. RMST & CUMULATIVE INCIDENCE ####
# Risk curves and confidence intervals for overdose and discontinuation
predict_cuminc_rmst_fast <- function(Dat_expanded, trt_var, covs, horizons,
                                     cluster_var, samples = 999, conf_int = TRUE,
                                     ref_group = "no reversion", wt_type = "wt_trunc",
                                     rmst_tau = 365) {
  
  # 1. Fit pooled logistic regression
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  model <- glm(plr_formula, data = Dat_expanded,
               family = quasibinomial(), weights = Dat_expanded[[wt_type]])
  
  V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                        sandwich = TRUE, fix = FALSE)
  coefs <- coef(model)
  
  # 2. Reference group & Factor Level Preservation
  is_trt_factor <- is.factor(Dat_expanded[[trt_var]])
  if (is_trt_factor) {
    trt_levels <- levels(Dat_expanded[[trt_var]])
  } else {
    trt_levels <- as.character(unique(Dat_expanded[[trt_var]]))
  }
  
  if (is.null(ref_group)) ref_group <- trt_levels[1]
  
  # 3. Pre-compute Design Matrices using data.table for extreme speed
  dt <- as.data.table(Dat_expanded)
  setorderv(dt, c(cluster_var, "stop")) 
  
  # Build a list of model matrices (one for each treatment level) ONCE.
  mm_list <- lapply(trt_levels, function(lvl) {
    dt_temp <- copy(dt)
    
    # Explicitly preserve factor levels so we don't trigger the contrasts error
    if (is_trt_factor) {
      dt_temp[[trt_var]] <- factor(lvl, levels = trt_levels)
    } else {
      dt_temp[[trt_var]] <- lvl
    }
    
    # Use xlev = model$xlevels to strictly enforce the model's known factor levels
    model.matrix(delete.response(terms(model)), data = dt_temp, xlev = model$xlevels)
  })
  names(mm_list) <- trt_levels
  
  # 4. Vectorized Risk Function (Runs strictly on pre-computed matrices)
  fast_risk_fun <- function(coef_draw) {
    res_list <- vector("list", length(trt_levels))
    
    for (i in seq_along(trt_levels)) {
      lvl <- trt_levels[i]
      
      # Matrix multiplication for all observations at once
      dt[, p_event := plogis(mm_list[[lvl]] %*% coef_draw)]
      
      # Fast cumulative product of survival probabilities grouped by cluster
      dt[, cum_surv := cumprod(1 - p_event), by = cluster_var]
      
      # Extract mean risk at exact specified horizons
      risks_h <- numeric(length(horizons))
      for (j in seq_along(horizons)) {
        h <- horizons[j]
        risk_val <- 1 - dt[stop <= h, .(last_surv = cum_surv[.N]), by = cluster_var][, mean(last_surv)]
        risks_h[j] <- risk_val
      }
      
      res_list[[i]] <- data.table(trt = lvl, horizon = horizons, risk = risks_h)
    }
    rbindlist(res_list)
  }
  
  # 5. Point Estimates
  point <- fast_risk_fun(coefs)
  
  # Return early if no CI requested
  if (!conf_int) {
    res <- list(risks = as_tibble(point))
    
    if (!is.null(rmst_tau)) {
      if (max(horizons) < rmst_tau) stop("max(horizons) must be >= rmst_tau to compute RMST.")
      
      point_rmst <- point %>%
        filter(horizon <= rmst_tau) %>%
        group_by(trt) %>%
        arrange(horizon) %>%
        mutate(
          survival = 1 - risk,
          time_diff = horizon - dplyr::lag(horizon, default = 0),
          area = dplyr::lag(survival, default = 1) * time_diff
        ) %>%
        summarise(RMST = sum(area), .groups = "drop") %>%
        mutate(
          ref_RMST = RMST[trt == ref_group],
          RMST_diff = ifelse(trt == ref_group, 0, RMST - ref_RMST)
        ) %>% select(-ref_RMST)
      
      res$rmst <- point_rmst
    }
    return(res)
  }
  
  # 6. Sequential Simulation for Confidence Intervals (Safe Mode with Fallback)
  sim_coefs <- tryCatch({
    # Primary Attempt: Standard multivariate normal draw
    mvtnorm::rmvnorm(n = samples, mean = coefs, sigma = V)
    
  }, error = function(e) {
    # Fallback: Matrix correction for asymmetry or non-positive-definiteness
    warning(paste("Standard mvtnorm failed:", e$message, 
                  "\nApplying nearest positive-definite projection to variance matrix."))
    
    # Force perfect symmetry
    V_sym <- (V + t(V)) / 2
    
    # Project to nearest positive-definite matrix
    V_pd <- as.matrix(Matrix::nearPD(V_sym)$mat)
    
    # Draw from the corrected matrix
    mvtnorm::rmvnorm(n = samples, mean = coefs, sigma = V_pd)
  })
  
  # Standard lapply - completely avoids parallel processing instability
  sim_list <- lapply(seq_len(nrow(sim_coefs)), function(i) {
    fast_risk_fun(sim_coefs[i, ])
  })
  
  # ---> RESTORED: Combine the list into the required dataframe <---
  sim_df <- rbindlist(sim_list, idcol = "sim")
  
  
  
  # 7. Compute Final Statistics (Standard dplyr format)
  risk_ci <- sim_df %>%
    group_by(trt, horizon) %>%
    summarise(
      risk_lower = quantile(risk, 0.025),
      risk_upper = quantile(risk, 0.975),
      risk = mean(risk), 
      .groups = "drop"
    )
  
  sim_contrasts <- sim_df %>%
    group_by(sim, horizon) %>%
    mutate(ref_risk = risk[trt == ref_group]) %>%
    filter(trt != ref_group) %>%
    mutate(RD = risk - ref_risk, RR = risk / ref_risk)
  
  contrast_ci <- sim_contrasts %>%
    group_by(trt, horizon) %>%
    summarise(
      RD_lower = quantile(RD, 0.025),
      RD_upper = quantile(RD, 0.975),
      RD = mean(RD),
      RR_lower = quantile(RR, 0.025),
      RR_upper = quantile(RR, 0.975),
      RR = mean(RR),
      .groups = "drop"
    )
  
  res <- list(
    risks = risk_ci, 
    contrasts = contrast_ci, 
    ref_group = ref_group,
    sim_data = sim_df 
  )
  
  
  # 8. Compute RMST distributions if requested
  if (!is.null(rmst_tau)) {
    if (max(horizons) < rmst_tau) stop("max(horizons) must be >= rmst_tau to compute RMST.")
    
    sim_rmst <- sim_df %>%
      filter(horizon <= rmst_tau) %>%
      group_by(sim, trt) %>%
      arrange(horizon) %>%
      mutate(
        survival = 1 - risk,
        time_diff = horizon - dplyr::lag(horizon, default = 0),
        area = dplyr::lag(survival, default = 1) * time_diff
      ) %>%
      summarise(RMST = sum(area), .groups = "drop")
    
    rmst_ci <- sim_rmst %>%
      group_by(trt) %>%
      summarise(
        RMST_lower = quantile(RMST, 0.025),
        RMST_upper = quantile(RMST, 0.975),
        RMST = mean(RMST),
        .groups = "drop"
      )
    
    rmst_contrasts <- sim_rmst %>%
      group_by(sim) %>%
      mutate(ref_RMST = RMST[trt == ref_group]) %>%
      filter(trt != ref_group) %>%
      mutate(RMST_diff = RMST - ref_RMST) %>%
      group_by(trt) %>%
      summarise(
        RMST_diff_lower = quantile(RMST_diff, 0.025),
        RMST_diff_upper = quantile(RMST_diff, 0.975),
        RMST_diff = mean(RMST_diff),
        .groups = "drop"
      )
    
    res$rmst <- rmst_ci
    res$rmst_contrasts <- rmst_contrasts
  }
  
  return(res)
}





cat("\n=== Computing Cumulative Incidence ===\n")

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))

# ── INITIALIZE RESULTS STRUCTURE ──

risk <- list(
  tr2 = list(od = list(), dis = list()),
  tr3 = list(od = list(), dis = list())
)

## MET interruption of 3-4 days #############
### FIRST-INTERRUPTION-IN-ANY-EPISODE ####
cat("  Computing: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

risk$tr2$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ret_length", "ddose_prev"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 999, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ret_length", "ddose_prev"), 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 999, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)

risk$tr2$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs = c("year", "ret_length", "ddose_prev"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 999, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("year", "ret_length", "ddose_prev"), 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 999, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2.rds"))


cat("  ✓ Trial 2 risk estimates saved.\n")

## MET interruption of 5-14 days #############

### FIRST-INTERRUPTION-IN-ANY-EPISODE ####
cat("  Computing: First-interruption-in-any-episode\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

# Overdose
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_od_acute",
                                 event_var = "sens_od_acute",
                                 interval_length = 14)

risk$tr3$od$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 999, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$od$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 999, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$dis$curve <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 999, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$dis$ci <- predict_cuminc_rmst_fast(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 999, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")
