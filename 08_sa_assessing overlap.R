# PS MIRROR PLOT ########


library(tidyverse)

# 1. Define colors 
COMPARISON_COLORS <- c(
  "No dose change"       = "#000000", 
  "Dose increase"        = "#D62728",  
  "1-<50% reduction"     = "#2CA02C",  
  "≥50% reduction/≤40mg" = "#FF7F0E"  
)

# 2. Helper function to generate the plot
build_mirror_ess_plot <- function(ps_obj) {
  
  # Format dataset labels
  dat_plot <- ps_obj$data_weights %>%
    mutate(
      trt = recode(trt,
                   "no reversion" = "No dose change",
                   "1-<50%" = "1-<50% reduction",
                   ">=50%/<=40mg" = "≥50% reduction/≤40mg",
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
                          "\nTruncated ESS = ", formatC(ess, format="f", digits=0, big.mark=","),
                          "\nOW ESS = ", formatC(ess_ow, format="f", digits=0, big.mark=","))
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
      caption = "Note: Top (Grey) represents the unweighted sample size. Bottom (Colored) represents the weighted pseudo-population using truncated IPTW."
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
ps_tr2_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
ps_tr3_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

# 4. Generate and save Trial 2 Plot (3-4 Days)
p_tr2 <- build_mirror_ess_plot(ps_tr2_all)

jpeg(paste0(RES_DIR, "Figures/June09/Diagnostic_TR2_Mirrored_with_ESS.jpg"), 
     width = 3600, height = 2600, res = 300)
print(p_tr2)
dev.off()

# 5. Generate and save Trial 3 Plot (5-14 Days)
p_tr3 <- build_mirror_ess_plot(ps_tr3_all)

jpeg(paste0(RES_DIR, "Figures/June09/Diagnostic_TR3_Mirrored_with_ESS.jpg"), 
     width = 3600, height = 2600, res = 300)
print(p_tr3)
dev.off()





# PS OVERLAP PLOT ###########
library(tidyverse)
library(nnet)

# 1. Define colors 
COMPARISON_COLORS <- c(
  "No dose change"       = "#000000", 
  "Dose increase"        = "#D62728",  
  "1-<50% reduction"     = "#2CA02C",  
  "≥50% reduction/≤40mg" = "#FF7F0E"  
)

# 2. Helper function to generate Pairwise Overlap Density Plots with ALIGNED axes
build_pairwise_overlap_plot <- function(ps_obj) {
  
  # A. Format dataset labels and ensure row indexing is intact
  dat <- ps_obj$data_weights %>%
    mutate(
      trt = recode(trt,
                   "no reversion" = "No dose change",
                   "1-<50%" = "1-<50% reduction",
                   ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                   "increase" = "Dose increase"),
      trt = factor(trt, levels = names(COMPARISON_COLORS)),
      row_id = row_number()
    )
  
  # B. Extract the full Nx4 matrix of probabilities
  ps_matrix <- predict(ps_obj$model, type = "probs")
  
  # Rename columns safely to match our clean labels
  colnames(ps_matrix) <- recode(colnames(ps_matrix),
                                "no reversion" = "No dose change",
                                "1-<50%" = "1-<50% reduction",
                                ">=50%/<=40mg" = "≥50% reduction/≤40mg",
                                "increase" = "Dose increase")
  
  # C. Create pairwise comparison datasets
  active_strats <- c("1-<50% reduction", "≥50% reduction/≤40mg", "Dose increase")
  
  plot_data <- bind_rows(
    lapply(active_strats, function(active_trt) {
      
      # Filter to only the specific active strategy and the reference group
      d_sub <- dat %>% filter(trt %in% c("No dose change", active_trt))
      
      # Extract the propensity score for the ACTIVE treatment for these patients
      d_sub$ps_score <- ps_matrix[d_sub$row_id, active_trt]
      
      # Label the facet
      d_sub$comparison <- paste(active_trt, "vs No dose change")
      return(d_sub)
    })
  )
  
  # Lock the facet order
  plot_data$comparison <- factor(plot_data$comparison, 
                                 levels = paste(active_strats, "vs No dose change"))
  
  # D. Build the Mirrored Density Plot
  p_pairwise <- ggplot(plot_data, aes(x = ps_score, fill = trt)) +
    
    # TOP HALF: Active Strategy Density
    geom_density(data = filter(plot_data, trt != "No dose change"),
                 aes(y = after_stat(density)), alpha = 0.8, color = "black") +
    
    # BOTTOM HALF: No dose change Density
    geom_density(data = filter(plot_data, trt == "No dose change"),
                 aes(y = -after_stat(density)), alpha = 0.8, color = "black") +
    
    # ALIGNED AXES FIX: Removed scales = "free_y" to lock the y-axis across all panels
    facet_wrap(~ comparison, ncol = 3) +
    
    # ALIGNED AXES FIX: Locked x-axis strictly from 0 to 1
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
    scale_y_continuous(labels = abs) +
    scale_fill_manual(values = COMPARISON_COLORS) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Estimated Propensity Score for Adjustment Strategy",
      y = "Density",
      caption = "Note: Top (Colored) displays the density of the active dose adjustment strategy. Bottom (Black) displays the reference group."
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray85"),
      axis.line = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm"), 
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      axis.text = element_text(color = "black", size = 11),
      plot.caption = element_text(hjust = 0, size = 12, color = "black", margin = ggplot2::margin(t = 15))
    )
  
  return(p_pairwise)
}

# 3. Load the model objects
ps_tr2_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
ps_tr3_all <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

# 4. Generate and save Trial 2 Plot
p_pair_tr2 <- build_pairwise_overlap_plot(ps_tr2_all)
jpeg(paste0(RES_DIR, "Figures/June09/Diagnostic_TR2_Pairwise_Overlap_Aligned.jpg"), 
     width = 3800, height = 2000, res = 300)
print(p_pair_tr2)
dev.off()

# 5. Generate and save Trial 3 Plot
p_pair_tr3 <- build_pairwise_overlap_plot(ps_tr3_all)
jpeg(paste0(RES_DIR, "Figures/June09/Diagnostic_TR3_Pairwise_Overlap_Aligned.jpg"), 
     width = 3800, height = 2000, res = 300)
print(p_pair_tr3)
dev.off()
