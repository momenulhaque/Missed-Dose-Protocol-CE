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
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

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

# Shared pre-processing: collapse all reduction categories into one
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

# 1. First-interruption-in-any-episode ############
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_first <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_all() 

dat_tr3_first <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_all()


# Combine datasets and define clinical epochs (UPDATED TO TWO STRATA)
df_explore_first <- bind_rows(
  dat_tr2_first %>% mutate(Trial = "Methadone interruption of 3-4 days"),
  dat_tr3_first %>% mutate(Trial = "Methadone interruption of 5-14 days")
) %>%
  mutate(
    # Create the categorical epochs (<= 30 vs > 30)
    ret_length_cat = case_when(
      ret_length <= 30 ~ "\u2264 30 days",
      ret_length > 30 ~ "> 30 days"
    ),
    ret_length_cat = factor(ret_length_cat, 
                            levels = c("\u2264 30 days", "> 30 days")),
    
    # Calculate the amount of increase (since ddose_red is before - after, an increase is negative)
    increase_amount = ifelse(trt == "increase", abs(ddose_red), NA),
    
    # Format the dose adjustment strategies into 5 detailed categories
    trt_detailed = case_when(
      trt == "no reversion" ~ "No dose change",
      trt == "1-<50%" ~ "<50% reduction",
      trt == ">=50%/<=40mg" ~ "\u226550% reduction/\u226440mg",
      trt == "increase" & increase_amount <= 10 ~ "Dose increase \u2264 10mg",
      trt == "increase" & increase_amount > 10 ~ "Dose increase > 10mg",
      TRUE ~ NA_character_
    ),
    trt_detailed = factor(trt_detailed, levels = c("No dose change", 
                                                   "<50% reduction", 
                                                   "\u226550% reduction/\u226440mg", 
                                                   "Dose increase \u2264 10mg", 
                                                   "Dose increase > 10mg"))
  )

# Define the new 5-color palette (Using two distinct reds for the increases)
trt_colors_detailed <- c(
  "No dose change" = "#000000",                   # Black
  "<50% reduction" = "#2CA02C",                   # Green
  "\u226550% reduction/\u226440mg" = "#FF7F0E",  # Orange
  "Dose increase \u2264 10mg" = "#9467BD",        # Muted Purple
  "Dose increase > 10mg" = "#D62728"              # Dark Red
)


plot_ecdf_first <- ggplot(df_explore_first, aes(x = ret_length, color = Trial)) +
  stat_ecdf(geom = "step", size = 1.2) +
  
  # Clinical anchor line matching the bottom panel epochs
  geom_vline(xintercept = 30, linetype = "dashed", color = "gray60") +
  
  scale_x_continuous(limits = c(0, 365), breaks = c(0, 30, 90, 180, 270, 365)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c("Methadone interruption of 3-4 days" = "#1f78b4", "Methadone interruption of 5-14 days" = "#33a02c")) +
  
  labs(
    title = "A. Cumulative proportion of first interruption occurred in OAT episode",
    x = "Duration of OAT episode prior to the index interruption (in days)",
    y = "Cumulative proportion first interruptions (%)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.y      = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.x      = element_line(linewidth = 1.2, color = "black"),
    axis.text.x       = element_text(color = "black", size = 14),
    axis.text.y       = element_text(color = "black", size = 14),
    axis.title.y      = element_text(size = 14),
    legend.position = c(0.8, 0.2), # Inset legend
    legend.title      = element_blank(),
    legend.text       = element_text(size = 14, color = "black"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key        = element_blank(),
    strip.background  = element_rect(fill = "lightgray", color = NA),
    strip.text        = element_text(face = "bold", size = 13)
    
  )


plot_prop_first <- ggplot(df_explore_first, aes(x = ret_length_cat, fill = trt_detailed)) +
  geom_bar(position = "fill", color = "black", alpha = 0.85, width = 0.7) +
  
  facet_wrap(~ Trial) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = trt_colors_detailed) +
  
  labs(
    title = "B. Proportion of dose adjustment after first interruption in any OAT episode",
    x = "\nDuration of OAT episode prior to the index interruption",
    y = "Proportion of dose adjustment strategy (%)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.y      = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.x      = element_line(linewidth = 1.2, color = "black"),
    axis.text.x       = element_text(color = "black", size = 14),
    axis.text.y       = element_text(color = "black", size = 14),
    axis.title.y      = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text       = element_text(size = 14, color = "black"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key        = element_blank(),
    strip.background  = element_rect(fill = "lightgray", color = NA),
    strip.text        = element_text(face = "bold", size = 13)
    
  )


# Combine the plots with patchwork
final_composite_first <- plot_ecdf_first / plot_spacer() / plot_prop_first +
  plot_layout(heights = c(1, 0.05, 1.2))


# Save the figure
jpeg(paste0(RES_DIR, "Figures/June09/Dist_of_interruptions_and_doses.jpg"), 
     width = 3900, height = 4600, res = 300)
print(final_composite_first)
dev.off()


df_explore_first %>% count(trt)


## % ##############
df_explore_first %>% 
  filter(Trial == "Methadone interruption of 5-14 days") %>% 
  count(trt_detailed)


t = c(620, 199, 68, 37)
gt10=c(199, 37)
sum(gt10)/sum(t)
