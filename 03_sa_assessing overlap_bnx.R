# 1. PS plots ############

library(tidyverse)
library(nnet)

# 1. Define colors
COMPARISON_COLORS <- c(
  "No dose change"       = "#000000",
  "Dose increase"        = "#D62728",
  "<50% reduction"       = "#2CA02C",
  "≥50% reduction/≤4mg" = "#FF7F0E"
)

# 2. Load PS model objects
ps_tr4_all <- readRDS(paste0(RES_DIR, "PScore/ps_weight_tr4_all.rds"))$data_weights
ps_tr5_all <- readRDS(paste0(RES_DIR, "PScore/ps_weight_tr5_all.rds"))$data_weights

# 3. Recode and label
recode_ps <- function(dat, trial_label) {
  dat %>%
    mutate(
      group = recode(trt,
                     "no reversion" = "No dose change",
                     "1-<50%"       = "<50% reduction",
                     ">=50%/<=4mg" = "≥50% reduction/≤4mg",
                     "increase"     = "Dose increase"
      ),
      group = factor(group, levels = c(
        "No dose change", "<50% reduction",
        "≥50% reduction/≤4mg", "Dose increase"
      )),
      Trial = trial_label
    ) %>%
    select(group, ps, Trial)
}

tr4 <- recode_ps(ps_tr4_all, "Buprenorphine/naloxone interruption of 3-5 days")
tr5 <- recode_ps(ps_tr5_all, "Buprenorphine/naloxone interruption of 6-14 days")
dat_all <- bind_rows(tr4, tr5)

# 4. Plot function
ps_score_plot_c3 <- function(dat,
                             x_lab = "Propensity Score",
                             y_lab = "Density",
                             legend_position = c(0.02, 0.02),
                             legend_justification = c(0, 0),
                             x_points = c(0, .25, .5, .75, 1)) {
  ggplot(dat, aes(x = ps, colour = group, fill = group)) +
    geom_histogram(aes(y = after_stat(density)),
                   position = "identity", alpha = 0.2, bins = 30) +
    geom_density(alpha = 0.05, linewidth = 0.5) +
    scale_colour_manual(values = COMPARISON_COLORS) +
    scale_fill_manual(values = COMPARISON_COLORS) +
    facet_wrap(~Trial, nrow = 2, scales = "free") +
    labs(title = NULL, x = x_lab, y = y_lab) +
    scale_x_continuous(breaks = x_points) +
    theme_minimal(base_size = 13) +
    theme(
      title            = element_text(color = "black", size = 14),
      panel.grid       = element_blank(),
      axis.line        = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks.y     = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x     = element_line(linewidth = 1.2, color = "black"),
      axis.text        = element_text(color = "black", size = 14),
      axis.title.y     = element_text(size = 14),
      legend.position      = legend_position,
      legend.justification = legend_justification,
      legend.direction     = "horizontal",
      legend.title         = element_blank(),
      legend.text          = element_text(size = 14, color = "black"),
      legend.background    = element_rect(fill = "white", color = "black"),
      legend.key           = element_blank(),
      strip.background     = element_rect(fill = "lightgray", color = NA),
      strip.text           = element_text(face = "bold", size = 13)
    )
}

# 5. Generate and save
plot2 <- ps_score_plot_c3(dat_all, legend_position = "bottom") +
  guides(color = guide_legend(nrow = 1))

jpeg(paste0(RES_DIR, "Figures/PS_plots.jpg"),
     width = 2600, height = 2500, res = 300)
print(plot2)
dev.off()




# 2. PS MIRROR PLOT ########

library(tidyverse)

# 1. Define colors 
COMPARISON_COLORS <- c(
  "No dose change"       = "#000000", 
  "Dose increase"        = "#D62728",  
  "<50% reduction"     = "#2CA02C",  
  "≥50% reduction/≤4mg" = "#FF7F0E"  
)

# 2. Helper function to generate the plot
build_mirror_ess_plot <- function(ps_obj) {
  
  # Format dataset labels
  dat_plot <- ps_obj$data_weights %>%
    mutate(
      trt = recode(trt,
                   "no reversion" = "No dose change",
                   "1-<50%" = "<50% reduction",
                   ">=50%/<=4mg" = "≥50% reduction/≤4mg",
                   "increase" = "Dose increase"),
      trt = factor(trt, levels = names(COMPARISON_COLORS))
    )
  
  # Calculate ESS
  ess_summary <- dat_plot %>%
    group_by(trt) %>%
    summarise(
      raw_n = n(),
      ess = (sum(wt_trunc))^2 / sum(wt_trunc^2),
      ess_ow = (sum(ow_wt))^2 / sum(ow_wt^2),
      .groups = "drop"
    ) %>%
    mutate(
      label_text = paste0("Original N = ", formatC(raw_n, format="f", digits=0, big.mark=","),
                          "\nIPTW-ESS = ", formatC(ess, format="f", digits=0, big.mark=","),
                          "\nOW-ESS = ", formatC(ess_ow, format="f", digits=0, big.mark=","))
    )
  
  # Build Plot
  p_mirror <- ggplot(dat_plot, aes(x = ps)) +
    
    # TOP HALF: Unweighted raw count
    geom_histogram(aes(y = after_stat(count)), bins = 30, 
                   fill = "grey70", color = "grey50", alpha = 0.8) +
    
    # BOTTOM HALF: Truncated weighted count
    geom_histogram(aes(y = -after_stat(count), weight = wt_trunc, fill = trt), 
                   bins = 30, color = "black", alpha = 0.8) +
    
    # Facet into 4 panels
    facet_wrap(~ trt, scales = "free_y") +
    
    # Inject ESS text annotations
    geom_text(data = ess_summary, 
              aes(x = Inf, y = Inf, label = label_text),
              hjust = 1.05, vjust = 1.5, size = 4.5, fontface = "bold", color = "black", inherit.aes = FALSE) +
    
    scale_y_continuous(labels = abs) +
    scale_fill_manual(values = COMPARISON_COLORS) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    
    # Titles removed, caption added for the note at the bottom
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Estimated Propensity Score",
      y = "Count / Weighted Count",
      caption = "Note: Top (Grey) represents the unweighted sample.\nBottom (Colored) represents the population after using inverse propbability of treatment weight (IPTW).\nOW-ESS, effective sample size under overlap weights.\nIPTW-ESS, effective sample size under IPTW"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(), 
      axis.line = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm"), 
      axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
      axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
      axis.text.x = element_text(color = "black", size = 11),
      axis.text.y = element_text(color = "black", size = 10),
      # THE FIX: Explicitly call ggplot2::margin
      plot.caption = element_text(hjust = 0, size = 12, color = "black", margin = ggplot2::margin(t = 15)) 
    )
  
  return(p_mirror)
}

# 3. Load the model objects (if not already in environment)
ps_tr4_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr4_all.rds"))
ps_tr5_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr5_all.rds"))

# 4. Generate and save Trial 4 Plot (3-5 Days)
p_tr4 <- build_mirror_ess_plot(ps_tr4_all)

jpeg(paste0(RES_DIR, "Figures/Diagnostic_TR4_Mirrored_with_ESS.jpg"), 
     width = 3600, height = 2600, res = 300)
print(p_tr4)
dev.off()

# 5. Generate and save Trial 3 Plot (5-14 Days)
p_tr5 <- build_mirror_ess_plot(ps_tr5_all)

jpeg(paste0(RES_DIR, "Figures/Diagnostic_TR5_Mirrored_with_ESS.jpg"), 
     width = 3600, height = 2600, res = 300)
print(p_tr5)
dev.off()







