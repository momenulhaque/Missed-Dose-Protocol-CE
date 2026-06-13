# Install if necessary: install.packages("riskRegression")
library(survival)
library(riskRegression)
library(dplyr)
library(tibble)

library(survival)
library(dplyr)
library(tibble)

predict_cuminc_bellavia <- function(Dat, trt_var, mod_var, mod_seq, covs = NULL, 
                                    horizons, time_var, event_var, cluster_var, 
                                    samples = 499, conf_int = TRUE,
                                    ref_group = "no reversion", wt_type = "wt_trunc",
                                    seed = 12345) {
  
  original_levels <- levels(factor(Dat[[trt_var]]))
  if(is.null(ref_group)) ref_group <- original_levels[1]
  
  # 1. Build Cox Formula
  if (is.null(covs) || length(covs) == 0) {
    cox_formula <- as.formula(paste("Surv(", time_var, ",", event_var, ") ~", trt_var, "* ns(", mod_var, ", df=3)"))
  } else {
    covs_str <- paste0(covs, collapse = " + ")
    cox_formula <- as.formula(paste("Surv(", time_var, ",", event_var, ") ~", trt_var, "* ns(", mod_var, ", df=3) +", covs_str))
  }
  
  cat("Fitting Cox model:", deparse(cox_formula), "\n")
  
  # 2. Fit Weighted Cox Model (Suppressing sparse-data infinite coefficient warnings for the splines)
  model <- suppressWarnings(coxph(cox_formula, data = Dat, weights = Dat[[wt_type]], 
                                  id = Dat[[cluster_var]], robust = TRUE, x = TRUE, y = TRUE))
  
  # Helper: Calculate standardized risks using native survival::survfit
  calc_risks <- function(fit, data_to_predict) {
    grid_combos <- expand.grid(trt_level = original_levels, mod_val = mod_seq, stringsAsFactors = FALSE)
    
    res_list <- lapply(1:nrow(grid_combos), function(i) {
      trt_l <- grid_combos$trt_level[i]
      mod_v <- grid_combos$mod_val[i]
      
      nd <- data_to_predict
      nd[[trt_var]] <- factor(trt_l, levels = original_levels)
      nd[[mod_var]] <- mod_v
      
      # Use base survival package to predict expected survival for all individuals
      sf <- survfit(fit, newdata = nd, se.fit = FALSE)
      
      # Extract survival exactly at the target horizons
      sf_sum <- summary(sf, times = horizons, extend = TRUE)
      
      # survfit returns a matrix: rows are horizons, columns are individuals
      surv_mat <- as.matrix(sf_sum$surv)
      
      # Standardize: Average the survival probability across all individuals, then convert to Risk
      avg_surv <- rowMeans(surv_mat, na.rm = TRUE)
      avg_risks <- 1 - avg_surv
      
      tibble(trt = trt_l, !!mod_var := mod_v, horizon = horizons, risk = avg_risks)
    })
    bind_rows(res_list)
  }
  
  cat("Calculating point estimates...\n")
  point_risks <- calc_risks(model, Dat)
  
  # Helper: Calculate Risk Differences and Ratios
  calc_contrasts <- function(risks_df) {
    risks_df %>%
      group_by(horizon, !!sym(mod_var)) %>%
      mutate(ref_risk = risk[trt == ref_group]) %>%
      filter(trt != ref_group) %>%
      mutate(RD = risk - ref_risk, RR = risk / ref_risk) %>%
      ungroup()
  }
  
  point_contrasts <- calc_contrasts(point_risks)
  
  # Early exit if CIs are not requested (Fast Mode)
  if(!conf_int) {
    point_contrasts <- point_contrasts %>%
      mutate(RD_lower = NA, RD_upper = NA, RR_lower = NA, RR_upper = NA, risk_lower = NA, risk_upper = NA)
    point_risks <- point_risks %>% mutate(risk_lower = NA, risk_upper = NA)
    return(list(risks = point_risks, contrasts = point_contrasts, ref_group = ref_group))
  }
  
  # 3. Non-Parametric Cluster Bootstrap
  cat("Running non-parametric cluster bootstrap (", samples, " samples)...\n", sep="")
  set.seed(seed)
  unique_clusters <- unique(Dat[[cluster_var]])
  n_clusters <- length(unique_clusters)
  
  sim_list <- list()
  for(b in 1:samples) {
    if (b %% 10 == 0) cat("  Bootstrap sample", b, "/", samples, "\n")
    
    boot_clusters <- sample(unique_clusters, size = n_clusters, replace = TRUE)
    
    boot_dat <- bind_rows(lapply(1:n_clusters, function(j) {
      c_data <- Dat[Dat[[cluster_var]] == boot_clusters[j], ]
      c_data$boot_id <- j 
      return(c_data)
    }))
    
    boot_fit <- suppressWarnings(coxph(cox_formula, data = boot_dat, weights = boot_dat[[wt_type]], 
                                       id = boot_dat$boot_id, robust = TRUE, x = TRUE, y = TRUE))
    
    boot_risks <- suppressWarnings(calc_risks(boot_fit, Dat))
    sim_list[[b]] <- boot_risks
  }
  
  # 4. Extract Confidence Intervals
  sim_df <- bind_rows(sim_list, .id = "sim")
  
  risk_ci <- sim_df %>%
    group_by(trt, horizon, !!sym(mod_var)) %>%
    summarise(risk_lower = quantile(risk, 0.025, na.rm=TRUE), 
              risk_upper = quantile(risk, 0.975, na.rm=TRUE), .groups = "drop") %>%
    left_join(point_risks, by = c("trt", "horizon", mod_var))
  
  sim_contrasts <- bind_rows(lapply(sim_list, calc_contrasts), .id = "sim")
  
  contrast_ci <- sim_contrasts %>%
    group_by(trt, horizon, !!sym(mod_var)) %>%
    summarise(
      RD_lower = quantile(RD, 0.025, na.rm=TRUE), RD_upper = quantile(RD, 0.975, na.rm=TRUE),
      RR_lower = quantile(RR, 0.025, na.rm=TRUE), RR_upper = quantile(RR, 0.975, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    left_join(point_contrasts %>% select(-ref_risk), by = c("trt", "horizon", mod_var))
  
  list(risks = risk_ci, contrasts = contrast_ci, ref_group = ref_group)
}





# MASTER CONTROL SWITCH
RUN_BOOTSTRAP_CI <- FALSE    
N_BOOT_SAMPLES   <- 499      

## overdose ##############################
### MET interruption of 3-4 days ####

cat("  Computing: Overdose - Trial 2 (3-4 days)\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
dat_unexpanded <- ps_fit$data_weights # No expansion needed!

max_ret <- quantile(dat_unexpanded$ret_length, 0.95, na.rm = TRUE)
ret_seq_grid <- seq(0, max_ret, by = 5)

risk_em_ci_tr2_od <- predict_cuminc_bellavia(
  Dat           = dat_unexpanded, 
  trt_var       = "trt", 
  mod_var       = "ret_length",     
  mod_seq       = ret_seq_grid,     
  covs          = c("ddose_prev", "prior_episode"),
  horizons      = c(30, 90, 180, 365), 
  time_var      = "time_od_acute",  # <-- Cox Needs Time
  event_var     = "sens_od_acute",  # <-- Cox Needs Event
  cluster_var   = "studyid_rev", 
  samples       = N_BOOT_SAMPLES,     
  conf_int      = RUN_BOOTSTRAP_CI,   
  ref_group     = "no reversion", 
  wt_type       = "wt_trunc"
)

saveRDS(risk_em_ci_tr2_od, paste0(RES_DIR, "Risk/", "risk_tr2_34d_od_em_bel.rds"))

### MET interruption of 5-14 days ####

cat("  Computing: Overdose - Trial 3 (5-14 days)\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))
dat_unexpanded <- ps_fit$data_weights

max_ret <- quantile(dat_unexpanded$ret_length, 0.95, na.rm = TRUE)
ret_seq_grid <- seq(0, max_ret, by = 5)

risk_em_ci_tr3_od <- predict_cuminc_bellavia(
  Dat           = dat_unexpanded, 
  trt_var       = "trt", 
  mod_var       = "ret_length",     
  mod_seq       = ret_seq_grid,     
  covs          = c("ddose_prev", "prior_episode"),
  horizons      = c(30, 90, 180, 365), 
  time_var      = "time_od_acute",  
  event_var     = "sens_od_acute",  
  cluster_var   = "studyid_rev", 
  samples       = N_BOOT_SAMPLES,     
  conf_int      = RUN_BOOTSTRAP_CI,   
  ref_group     = "no reversion", 
  wt_type       = "wt_trunc"
)

saveRDS(risk_em_ci_tr3_od, paste0(RES_DIR, "Risk/", "risk_tr3_od_em_bel.rds"))

## treatment discontinuation ##############################
### MET interruption of 3-4 days ####

cat("  Computing: Discontinuation - Trial 2 (3-4 days)\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))
dat_unexpanded <- ps_fit$data_weights

max_ret <- quantile(dat_unexpanded$ret_length, 0.95, na.rm = TRUE)
ret_seq_grid <- seq(0, max_ret, by = 5)

risk_em_ci_tr2_dis <- predict_cuminc_bellavia(
  Dat           = dat_unexpanded, 
  trt_var       = "trt", 
  mod_var       = "ret_length",     
  mod_seq       = ret_seq_grid,     
  covs          = c("ddose_prev", "prior_episode"),
  horizons      = c(30, 90, 180, 365), 
  time_var      = "time_discont",  # <-- Discontinuation Time
  event_var     = "sens_discont",  # <-- Discontinuation Event
  cluster_var   = "studyid_rev", 
  samples       = N_BOOT_SAMPLES,     
  conf_int      = RUN_BOOTSTRAP_CI,   
  ref_group     = "no reversion", 
  wt_type       = "wt_trunc"
)

saveRDS(risk_em_ci_tr2_dis, paste0(RES_DIR, "Risk/", "risk_tr2_34d_dis_em_bel.rds"))

### MET interruption of 5-14 days ####

cat("  Computing: Discontinuation - Trial 3 (5-14 days)\n")
ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))
dat_unexpanded <- ps_fit$data_weights

max_ret <- quantile(dat_unexpanded$ret_length, 0.95, na.rm = TRUE)
ret_seq_grid <- seq(0, max_ret, by = 5)

risk_em_ci_tr3_dis <- predict_cuminc_bellavia(
  Dat           = dat_unexpanded, 
  trt_var       = "trt", 
  mod_var       = "ret_length",     
  mod_seq       = ret_seq_grid,     
  covs          = c("ddose_prev", "prior_episode"),
  horizons      = c(30, 90, 180, 365), 
  time_var      = "time_discont",  
  event_var     = "sens_discont",  
  cluster_var   = "studyid_rev", 
  samples       = N_BOOT_SAMPLES,     
  conf_int      = RUN_BOOTSTRAP_CI,   
  ref_group     = "no reversion", 
  wt_type       = "wt_trunc"
)

saveRDS(risk_em_ci_tr3_dis, paste0(RES_DIR, "Risk/", "risk_tr3_dis_em_bel.rds"))






## plotting ###########################

library(dplyr)
library(tibble)
library(sandwich)
library(mvtnorm)
library(dplyr)
library(ggplot2)


plot_rd_em <- function(risk_em, 
                       risk_at=90, 
                       ret_limit=180, 
                       y_min=-5, y_max=5,
                       title_name = "Adjusted risk difference varied with the duration of episode before index interruption",
                       legend_option=FALSE){
  
  plot_dat_rd_ci <- risk_em$contrasts %>%
    filter(horizon == risk_at, ret_length <= ret_limit) %>%
    mutate(
      # Convert estimates and bounds to percentages
      RD_pct = RD * 100,
      RD_lower_pct = RD_lower * 100,
      RD_upper_pct = RD_upper * 100,
      
      # Format labels to match your legend
      contrast = case_when(
        trt == "1-<50%" ~ "<50% reduction\nvs\nno change",
        trt == ">=50%/<=40mg" ~ "≥50% reduction/≤40mg\nvs\nno change",
        trt == "increase" ~ "Dose increase\nvs\nno change",
        TRUE ~ paste0(trt, "\nvs\nno change")
      )
    )
  
  
  
  
  # Define your colors again (make sure they match your original plot)
  comparison_colors_v1 = c(
    "<50% reduction\nvs\nno change" = "#2CA02C",                  # Green
    "≥50% reduction/≤40mg\nvs\nno change" = "#FF7F0E", # Orange
    "Dose increase\nvs\nno change" = "#D62728"                            # Red
  )
  
  p_curve_with_ci <- ggplot(plot_dat_rd_ci, aes(x = ret_length, y = RD_pct, color = contrast, fill = contrast)) +
    
    # 1. The Null Line (Must be drawn first so it sits behind the data)
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "red4") +
    
    # 2. The 95% CI Ribbons (Highly transparent, no borders)
    #geom_ribbon(aes(ymin = RD_lower_pct, ymax = RD_upper_pct), alpha = 0.15, color = NA) +
    
    # 3. The Central Estimate Lines (Drawn last so they sit on top)
    geom_line(size = 1.2) +
    
    # Map both color (for lines) and fill (for ribbons) to your palette
    scale_color_manual(values = comparison_colors_v1) +
    scale_fill_manual(values = comparison_colors_v1) +
    
    # Axes adjustments
    scale_x_continuous(breaks = seq(0, ret_limit, by=30)) + 
    scale_y_continuous(limits = c(y_min, y_max)) + # CRITICAL: Adjust based on your actual RD bounds!
    
    labs(
      title = title_name,
      x = "Duration of OAT episode prior to the index interruption (in days)",
      y = paste0("% difference in risk at ", risk_at, " days")
    ) +
    
    theme_minimal() +
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.x = element_line(linewidth = 1.1, color = "black"),
      axis.ticks.y = element_line(linewidth = 1.1, color = "black"),
      axis.ticks.length  = unit(0.15, "cm"),
      axis.text = element_text(color = "black", size = 14),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      legend.position = if(legend_option) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 14)
    )
  print(p_curve_with_ci)
}



risk_2_od <- readRDS(paste0(RES_DIR, "Risk/", "risk_tr2_34d_od_em_bel.rds"))
risk_3_od <- readRDS(paste0(RES_DIR, "Risk/", "risk_tr3_dis_em_bel.rds"))
risk_2_dis <- readRDS(paste0(RES_DIR, "Risk/", "risk_tr2_34d_dis_em_bel.rds"))
risk_3_dis <- readRDS(paste0(RES_DIR, "Risk/", "risk_tr3_dis_em_bel.rds"))



plot1 = plot_rd_em(risk_em=risk_2_od, 
                   risk_at=90,
                   ret_limit=180, y_min=-10, y_max=10,
                   title_name = NULL,
                   legend_option = FALSE)                   


plot2 = plot_rd_em(risk_em=risk_3_od, 
                   risk_at=90,
                   ret_limit=180, y_min=-10, y_max=15,
                   title_name = NULL,
                   legend_option = T)                   


plot3 = plot_rd_em(risk_em=risk_2_dis, 
                   risk_at=365,
                   ret_limit=180, y_min=-20, y_max=20,
                   title_name = NULL,
                   legend_option = FALSE)                   



plot4 = plot_rd_em(risk_em=risk_3_dis, 
                   risk_at=365,
                   ret_limit=180, y_min=-20, y_max=20,
                   title_name = NULL,
                   legend_option = T)                   










row_od_all <-  plot_spacer()/plot1/plot_spacer()/plot2 +
  plot_layout(heights = c(0.05, 1, 0.08, 1))


row_od_all_final <- ggdraw(row_od_all) +
  draw_label("A. Methadone interruption of 3-4 days: First-interruption-in-any-episode\n",
             x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  draw_label("B. Methadone interruption of 5-14 days: First-interruption-in-any-episode\n",
             x = 0.02, y = 0.53, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  
  
jpeg(paste0(RES_DIR, "Figures/", "risk_od_em_bel.jpg"),
       width = 3200, height = 3500, res = 300)
print(row_od_all_final)
dev.off()








row_dis_all <-  plot_spacer()/plot3/plot_spacer()/plot4 +
  plot_layout(heights = c(0.05, 1, 0.08, 1))


row_dis_all_final <- ggdraw(row_dis_all) +
  draw_label("A. Methadone interruption of 3-4 days: First-interruption-in-any-episode\n",
             x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  draw_label("B. Methadone interruption of 5-14 days: First-interruption-in-any-episode\n",
             x = 0.02, y = 0.53, hjust = 0, vjust = 1, size = 16, fontface = "bold") +
  
  
jpeg(paste0(RES_DIR, "Figures/", "risk_dis_em_bel.jpg"),
       width = 3200, height = 3500, res = 300)
print(row_dis_all_final)
dev.off()



