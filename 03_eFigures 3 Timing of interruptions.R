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

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")   

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

# ── 5-Category Pre-processing ──────────────────────────────────────────────────
TRT_LEVELS_5 <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase <50%", "increase >=50%")

clean_trt_5cat <- function(dat) {
  dat |>
    mutate(
      trt_s456 = case_when(
        ddose_rev == 0 ~ "no reversion",
        ddose_rev > 0 & ddose_rev < 50  ~ "1-<50%",
        (!is.na(ddose) & ddose < ddose_prev & ddose <= 40) | ddose_rev >= 50 ~ ">=50%/<=40mg",
        ddose_rev < 0 & ddose_rev > -50  ~ "increase <50%",
        ddose_rev <= -50  ~ "increase >=50%", 
        TRUE ~ NA_character_
      ),
      trt_s456 = factor(trt_s456, levels = TRT_LEVELS_5),
      prior_episode = as.numeric(episode - 1)
    ) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS_5) |>
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
  clean_trt_5cat() 

dat_tr3_first <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_5cat()

# Combine datasets and define clinical epochs & 5-category details
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
    
    # Format the dose adjustment strategies into 5 detailed categories
    trt_detailed = case_when(
      trt == "no reversion" ~ "No dose change",
      trt == "1-<50%" ~ "<50% reduction",
      trt == ">=50%/<=40mg" ~ "\u226550% reduction/\u226440mg",
      trt == "increase <50%" ~ "<50% increase",
      trt == "increase >=50%" ~ "\u226550% increase",
      TRUE ~ NA_character_
    ),
    trt_detailed = factor(trt_detailed, levels = c("No dose change", 
                                                   "<50% reduction", 
                                                   "\u226550% reduction/\u226440mg", 
                                                   "<50% increase", 
                                                   "\u226550% increase"))
  )

# Define the new 5-color palette (Matches the forest plots)
trt_colors_detailed <- c(
  "No dose change" = "black",                    
  "<50% reduction" = "#2CA02C",                  
  "\u226550% reduction/\u226440mg" = "#FF7F0E",  
  "<50% increase" = "#9467BD",                   
  "\u226550% increase" = "#A00000"               
)


## A. TIMING OF FIRST INTERRUPTION (ECDF) ######

plot_ecdf_first <- ggplot(df_explore_first, aes(x = ret_length, color = Trial)) +
  # Set pad = FALSE to prevent the line from anchoring to y=0 at the left edge
  stat_ecdf(geom = "step", size = 1.2, pad = FALSE) +
  
  # Clinical anchor line
  geom_vline(xintercept = 30, linetype = "dashed", color = "gray60") +
  
  scale_x_continuous(limits = c(0, 365), breaks = c(0, 30, 90, 180, 270, 365)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c("Methadone interruption of 3-4 days" = "#1f78b4", 
                                "Methadone interruption of 5-14 days" = "#33a02c")) +
  
  labs(
    title = NULL,
    x = "Time from OAT initiation to interruption (in days)",
    y = "Cumulative proportion of interruptions (%)"
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
    legend.position   = c(0.75, 0.2), # Inset legend
    legend.title      = element_blank(),
    legend.text       = element_text(size = 14, color = "black"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key        = element_blank(),
    strip.background  = element_rect(fill = "lightgray", color = NA),
    strip.text        = element_text(face = "bold", size = 13)
  )

jpeg(paste0(RES_DIR, "Figures/June09/eFigure 3A Timing_of_interruptions_ECDF.jpg"), 
     width = 3800, height = 2500, res = 400)
print(plot_ecdf_first)
dev.off()


## B. DOSE ADJUSTMENT PROPORTIONS (BAR CHART) ######

plot_prop_first <- ggplot(df_explore_first, aes(x = ret_length_cat, fill = trt_detailed)) +
  geom_bar(position = "fill", color = "black", alpha = 0.85, width = 0.7) +
  facet_wrap(~ Trial) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = trt_colors_detailed) +
  labs(
    title = NULL,
    x = "Time from OAT initiation to interruption (in days)",
    y = "Proportion of episodes (%)"
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
    axis.title.x      = element_text(margin = ggplot2::margin(t = 15)), 
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 14, color = "black"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key        = element_blank(),
    strip.background  = element_rect(fill = "lightgray", color = "black"), 
    strip.text        = element_text(face = "bold", size = 13)
  )

jpeg(paste0(RES_DIR, "Figures/June09/eFigure 3B Dose_adjustment_proportions.jpg"), 
     width = 4800, height = 2600, res = 400)
print(plot_prop_first)
dev.off()


df_explore_first %>% count(trt_s456)
df_explore_first %>% group_by(Trial) %>% filter(trt_s456=="increase") %>% count(trt_detailed)
