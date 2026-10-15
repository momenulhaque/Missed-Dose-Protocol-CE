# R2R version #####


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

source("H:/GitHub/Missed-Dose-Protocols-PRI/00_functions.R")   

DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/LE30D/"

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

raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_first <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_5cat() 

dat_tr3_first <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_5cat()

# Combine datasets and define clinical epochs & 5-category details
df_explore_first <- bind_rows(
  dat_tr2_first %>% mutate(Trial = "Interruptions of 3-4 days"),
  dat_tr3_first %>% mutate(Trial = "Interruptions of 5-14 days")
) %>%
  mutate(
    # Create the categorical epochs (keeping all data)
    ret_length_cat = case_when(
      ret_length <= 15 ~ "\u2264 15 days",
      ret_length > 15 & ret_length <= 30 ~ "16-30 days",
      ret_length > 30 & ret_length <= 60 ~ "31-60 days",
      ret_length > 60 & ret_length <= 120 ~ "61-120 days",
      ret_length > 120 & ret_length <= 180 ~ "121-180 days",
      ret_length > 180 ~ "> 180 days"
    ),
    ret_length_cat = factor(ret_length_cat, 
                            levels = c("\u2264 15 days", "16-30 days", 
                                       "31-60 days", "61-120 days", 
                                       "121-180 days", "> 180 days")),
    
    # Format the dose adjustment strategies into 5 detailed categories
    trt_detailed = case_when(
      trt == "no reversion" ~ "No dose change",
      trt == "1-<50%" ~ "<50% reduction",
      trt == ">=50%/<=40mg" ~ "\u226550% reduction or \u226440mg",
      trt == "increase <50%" ~ "<50% increase",
      trt == "increase >=50%" ~ "\u226550% increase",
      TRUE ~ NA_character_
    ),
    trt_detailed = factor(trt_detailed, levels = c("No dose change", 
                                                   "<50% reduction", 
                                                   "\u226550% reduction or \u226440mg", 
                                                   "<50% increase", 
                                                   "\u226550% increase"))
  )



## A. TIMING OF FIRST INTERRUPTION (ECDF) WITH ALIGNED TABLE #####

# 1. Create the base ECDF plot 
plot_ecdf_base <- ggplot(df_explore_first, aes(x = ret_length, color = Trial, linetype = Trial)) +
  stat_ecdf(geom = "step", size = 1.2, pad = FALSE) +
  
  # Clinical anchor lines
  geom_vline(xintercept = c(15, 30, 60, 120), linetype = "dashed", color = "gray60", alpha = 0.6) +
  
  coord_cartesian(xlim = c(0, 180)) +
  scale_x_continuous(breaks = c(0, 15, 30, 60, 120, 180)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, 0.25)) +
  
  # Set colors and linetypes 
  scale_color_manual(values = c("Interruptions of 3-4 days" = "#1f78b4", 
                                "Interruptions of 5-14 days" = "#33a02c")) +
  scale_linetype_manual(values = c("Interruptions of 3-4 days" = "solid", 
                                   "Interruptions of 5-14 days" = "dashed")) +
  
  labs(
    title = NULL,
    x = NULL, 
    y = "Cumulative incidence of first methadone\ninterruption (%)" 
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
    axis.title.y      = element_text(size = 14, margin = ggplot2::margin(r = 10)),
    legend.position   = c(0.75, 0.25), 
    legend.title      = element_blank(),
    legend.text       = element_text(size = 14, color = "black"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key        = element_blank(),
    legend.key.width  = unit(1.5, "cm") 
  )

# 2. Prepare the data for the Cumulative Proportion table
# We evaluate the cumulative proportion precisely at the x-axis tick marks
days_to_evaluate <- c(0, 15, 30, 60, 120, 180)

df_table_text <- purrr::map_dfr(days_to_evaluate, function(day) {
  df_explore_first %>%
    group_by(Trial) %>%
    summarise(
      Total_N = n(),
      cum_count = sum(ret_length <= day),
      .groups = "drop"
    ) %>%
    mutate(
      x_pos = day,
      cum_prop = cum_count / Total_N,
      # Format as a clean percentage (e.g., "45%")
      lbl = scales::percent(cum_prop, accuracy = 1)
    )
}) %>%
  # Append the Total N to the Trial name for the y-axis label
  mutate(Trial_Label = paste0(Trial, "\n(total episodes ", Total_N, ")"))

# 3. Create the text-only table plot 
# Using Trial_Label for the y-axis, but Trial for color to maintain mapping
plot_table <- ggplot(df_table_text, aes(x = x_pos, y = fct_rev(Trial_Label), label = lbl, color = Trial)) +
  geom_text(size = 4.2, fontface = "bold") +
  
  # Ensure the x-axis limits exactly match the ECDF plot for perfect alignment
  coord_cartesian(xlim = c(0, 180)) +
  scale_color_manual(values = c("Interruptions of 3-4 days" = "#1f78b4", 
                                "Interruptions of 5-14 days" = "#33a02c")) +
  
  labs(x = "Days from methadone episode initiation to first interruption", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks      = element_blank(),
    axis.title.x    = element_text(margin = ggplot2::margin(t = 15)),
    axis.text.y     = element_text(color = "black", size = 12, hjust = 1), 
    legend.position = "none"
  )

# 4. Stitch them together using patchwork
plot_ecdf_combined <- plot_ecdf_base / plot_table + 
  plot_layout(heights = c(4, 1)) 

jpeg(paste0(RES_DIR, "Figures/eFigure 2A Timing_of_interruptions_ECDF_with_Table.jpg"), 
     width = 3800, height = 2800, res = 400) 
print(plot_ecdf_combined)
dev.off()




## B. DOSE ADJUSTMENT PROPORTIONS (CLUSTERED BAR CHART) ######

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
library(ggpattern)

source("H:/GitHub/Missed-Dose-Protocols-PRI/00_functions.R")   

DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/LE30D/"

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

raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_tr2_first <- raw_all |>
  filter(trial_gap %in% 3:4) |>
  clean_trt_5cat() 

dat_tr3_first <- raw_all |>
  filter(trial_cat == "5-14") |>
  clean_trt_5cat()

# Combine datasets and define clinical epochs & 5-category details
df_explore_first <- bind_rows(
  dat_tr2_first %>% mutate(Trial = "Interruption of 3-4 days"),
  dat_tr3_first %>% mutate(Trial = "Interruption of 5-14 days")
) %>%
  mutate(
    # Create the categorical epochs (keeping all data)
    ret_length_cat = case_when(
      ret_length <= 15 ~ "\u2264 15 days",
      ret_length > 15 & ret_length <= 30 ~ "16-30 days",
      ret_length > 30 & ret_length <= 60 ~ "31-60 days",
      ret_length > 60 & ret_length <= 120 ~ "61-120 days",
      ret_length > 120 & ret_length <= 180 ~ "121-180 days",
      ret_length > 180 ~ "> 180 days"
    ),
    ret_length_cat = factor(ret_length_cat, 
                            levels = c("\u2264 15 days", "16-30 days", 
                                       "31-60 days", "61-120 days", 
                                       "121-180 days", "> 180 days")),
    
    # Format the dose adjustment strategies into 5 detailed categories
    trt_detailed = case_when(
      trt == "no reversion" ~ "No dose change",
      trt == "1-<50%" ~ "<50% reduction",
      trt == ">=50%/<=40mg" ~ "\u226550% reduction or \u226440mg",
      trt == "increase <50%" ~ "<50% increase",
      trt == "increase >=50%" ~ "\u226550% increase",
      TRUE ~ NA_character_
    ),
    trt_detailed = factor(trt_detailed, levels = c("No dose change", 
                                                   "<50% reduction", 
                                                   "\u226550% reduction or \u226440mg", 
                                                   "<50% increase", 
                                                   "\u226550% increase"))
  )


# 1. Prepare conditional data matching the text's exact cohorts and epochs
plot_data_cond <- df_explore_first %>% 
  # Drop "No dose change" as it does not have a magnitude sub-category
  filter(trt_detailed != "No dose change") %>%
  mutate(
    # Collapse time bins to match the manuscript text ("within 30 days" vs "after")
    epoch_text = ifelse(ret_length <= 30, "\u2264 30 days", "> 30 days"),
    epoch_text = factor(epoch_text, levels = c("\u2264 30 days", "> 30 days")),
    
    # Define adjustment direction to create distinct, accurate denominators
    direction = case_when(
      grepl("reduction", trt_detailed) ~ "Dose reduction after interruption",
      grepl("increase", trt_detailed) ~ "Dose increase after interruption"
    ),
    direction = factor(direction, levels = c("Dose reduction after interruption", 
                                             "Dose increase after interruption"))
  ) %>%
  group_by(Trial, direction, epoch_text, trt_detailed) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(
    # Conditional probability: divided strictly by the number of people in that direction
    total_in_direction = sum(count),
    prop = count / total_in_direction,
    lbl = scales::percent(prop, accuracy = 1)
  ) %>%
  ungroup()




# 2. Create the dual-faceted conditional plot with textures (No color)
plot_prop_cond <- ggplot(plot_data_cond, 
                         aes(x = epoch_text, y = prop, 
                             pattern = trt_detailed)) + 
  
  geom_col_pattern(
    position = position_dodge(width = 0.85), 
    fill = "white",              # Hardcoded to white to strip color but keep contrast
    color = "black",
    width = 0.75,
    pattern_color = "black",     # Color of the pattern lines
    pattern_fill = "black",      # Color of pattern shapes (like circles)
    pattern_angle = 45,          # Angle for stripes
    pattern_density = 0.06,      # Density of the pattern
    pattern_spacing = 0.05,      # Distance between pattern lines
    pattern_key_scale_factor = 0.6 
  ) +
  
  geom_text(aes(label = lbl, group = trt_detailed), 
            position = position_dodge(width = 0.85), 
            vjust = -0.5, size = 4, color = "black", fontface = "bold") +
  
  facet_grid(direction ~ Trial) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     expand = expansion(mult = c(0, 0.15))) + 
  
  scale_pattern_manual(name = NULL, values = c("none", "stripe", "crosshatch", "circle")) +
  
  # Update guide target from fill to pattern to maintain legend styling
  guides(pattern = guide_legend(override.aes = list(fill = "white", color = "black"))) +
  
  labs(
    title = NULL,
    x = "Days from methadone initiation to first interruption",
    y = "Proportion of Episodes (%)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(color = "black", size = 13),
    axis.text.y        = element_text(color = "black", size = 14),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    strip.background   = element_rect(fill = "lightgray", color = "black"), 
    strip.text         = element_text(face = "bold", size = 12)
  )

jpeg(paste0(RES_DIR, "Figures/eFigure 2B Dose_adjustment_conditional.jpg"), 
     width = 5200, height = 3400, res = 400) 
print(plot_prop_cond)
dev.off()



write.csv(plot_data_cond, paste0(RES_DIR, "Tables/eTable_con.csv")) 
