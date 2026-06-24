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
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/BNX/"

# ── 5-Category Pre-processing ──────────────────────────────────────────────────
TRT_LEVELS_5 <- c("no reversion", "1-<50%", ">=50%/<=4mg", "increase <50%", "increase >=50%")

clean_trt_5cat <- function(dat) {
  dat |>
    mutate(
      trt_s5 = case_when(
        ddose_rev == 0 ~ "no reversion",
        ddose_rev > 0 & ddose_rev < 50  ~ "1-<50%",
        (!is.na(ddose) & ddose < ddose_prev & ddose <= 40) | ddose_rev >= 50 ~ ">=50%/<=4mg",
        ddose_rev < 0 & ddose_rev > -50  ~ "increase <50%",
        ddose_rev <= -50  ~ "increase >=50%", 
        TRUE ~ NA_character_
      ),
      trt_s5 = factor(trt_s5, levels = TRT_LEVELS_5),
      prior_episode = as.numeric(episode - 1)
    ) |>
    rename(trt = trt_s5) |>
    filter(trt %in% TRT_LEVELS_5) |>
    group_by(studyid_rev, episode) |>
    mutate(
      prior_gap     = cumsum(gap_days),
      prior_gap_num = row_number()
    ) %>%
    ungroup()
}

# 1. First-interruption-in-any-episode ############
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_b.sas7bdat"))

dat_tr4_first <- raw_all |>
  filter(trial_gap %in% 3:5) |>
  clean_trt_5cat() 

dat_tr5_first <- raw_all |>
  filter(trial_cat == "6-14") |>
  clean_trt_5cat()

# Combine datasets and define clinical epochs & 5-category details
df_explore_first <- bind_rows(
  dat_tr4_first %>% mutate(Trial = "Interruption of 3-5 days"),
  dat_tr5_first %>% mutate(Trial = "Interruption of 6-14 days")
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
      trt == ">=50%/<=4mg" ~ "\u226550% reduction/\u22644mg",
      trt == "increase <50%" ~ "<50% increase",
      trt == "increase >=50%" ~ "\u226550% increase",
      TRUE ~ NA_character_
    ),
    trt_detailed = factor(trt_detailed, levels = c("No dose change", 
                                                   "<50% reduction", 
                                                   "\u226550% reduction/\u22644mg", 
                                                   "<50% increase", 
                                                   "\u226550% increase"))
  )



## A. TIMING OF FIRST INTERRUPTION (ECDF) WITH ALIGNED TABLE ######

# 1. Create the base ECDF plot (Removing the x-axis title so it can move to the bottom table)
plot_ecdf_base <- ggplot(df_explore_first, aes(x = ret_length, color = Trial)) +
  stat_ecdf(geom = "step", size = 1.2, pad = FALSE) +
  
  # Clinical anchor lines
  geom_vline(xintercept = c(15, 30, 60, 120), linetype = "dashed", color = "gray60", alpha = 0.6) +
  
  coord_cartesian(xlim = c(0, 180)) +
  scale_x_continuous(breaks = c(0, 15, 30, 60, 120, 180)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c("Interruption of 3-5 days" = "#1f78b4", 
                                "Interruption of 6-14 days" = "#33a02c")) +
  
  labs(
    title = NULL,
    x = NULL, # X-axis title removed here, will be added below the table
    y = "Cumulative proportion\nof interruptions (%)" # Added newline to prevent squeezing
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
    legend.key        = element_blank()
  )

# 2. Prepare the data for the interval table
df_table_text <- df_explore_first %>%
  group_by(Trial) %>%
  mutate(Total_N = n()) %>% # Calculate using the full 365-day denominator
  group_by(Trial, ret_length_cat) %>%
  summarise(
    prop = n() / first(Total_N),
    lbl = scales::percent(prop, accuracy = 1),
    .groups = "drop"
  ) %>%
  mutate(
    # Set the x-coordinates to the midpoint of each vertical line bin
    x_pos = case_when(
      ret_length_cat == "\u2264 15 days" ~ 7.5,
      ret_length_cat == "16-30 days" ~ 22.5,
      ret_length_cat == "31-60 days" ~ 45,
      ret_length_cat == "61-120 days" ~ 90,
      ret_length_cat == "121-180 days" ~ 150,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(x_pos)) # Drop the >180 bin so it doesn't plot off-screen

# 3. Create the text-only table plot
plot_table <- ggplot(df_table_text, aes(x = x_pos, y = fct_rev(Trial), label = lbl, color = Trial)) +
  geom_text(size = 4.5, fontface = "bold") +
  
  # Ensure the x-axis limits exactly match the ECDF plot for perfect alignment
  coord_cartesian(xlim = c(0, 180)) +
  scale_color_manual(values = c("Interruption of 3-5 days" = "#1f78b4", 
                                "Interruption of 6-14 days" = "#33a02c")) +
  
  labs(x = "Days from OAT initiation to interruption", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks      = element_blank(),
    axis.title.x    = element_text(margin = ggplot2::margin(t = 15)),
    # Right-align the trial names so they look like table row headers
    axis.text.y     = element_text(color = "black", size = 12, hjust = 1), 
    legend.position = "none"
  )

# 4. Stitch them together using patchwork
# The heights argument ensures the graph is taller than the table row
plot_ecdf_combined <- plot_ecdf_base / plot_table + 
  plot_layout(heights = c(4, 1)) 

jpeg(paste0(RES_DIR, "Figures/eFigure 3A Timing_of_interruptions_ECDF_with_Table.jpg"), 
     width = 3800, height = 2800, res = 400) # Slightly increased height to fit the table
print(plot_ecdf_combined)
dev.off()


## B. DOSE ADJUSTMENT PROPORTIONS (BAR CHART) ######
# Define the new 5-color palette (Matches the forest plots)
trt_colors_detailed <- c(
  "No dose change" = "black",                    
  "<50% reduction" = "#2CA02C",                  
  "\u226550% reduction/\u22644mg" = "#FF7F0E",  
  "<50% increase" = "#9467BD",                   
  "\u226550% increase" = "#A00000"               
)


# Filter out the > 180 days category just for the plot
plot_prop_first <- ggplot(df_explore_first %>% filter(ret_length <= 180), 
                          aes(x = ret_length_cat, fill = trt_detailed)) +
  geom_bar(position = "fill", color = "black", alpha = 0.85, width = 0.7) +
  facet_wrap(~ Trial) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = trt_colors_detailed) +
  labs(
    title = NULL,
    x = "Days from OAT initiation to interruption",
    y = "Proportion of episodes (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks.y      = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.x      = element_line(linewidth = 1.2, color = "black"),
    # Angled x-axis text to prevent overlap of the 5 new bins
    axis.text.x       = element_text(color = "black", size = 13, angle = 45, hjust = 1),
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

jpeg(paste0(RES_DIR, "Figures/eFigure 3B Dose_adjustment_proportions.jpg"), 
     width = 5200, height = 2800, res = 400)
print(plot_prop_first)
dev.off()