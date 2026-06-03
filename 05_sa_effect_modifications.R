install.packages("interactionRCS")
library(splines)
library(interactionRCS)
library(survival)
library(ggplot2)
library(scales)
library(cobalt)
library(glmtoolbox)

library(dplyr)
library(tibble)
library(sandwich)
library(mvtnorm)


source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")   # ps_model_wt(), predict_cuminc(), expand_person_period()

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"



predict_cuminc_em <- function(Dat_expanded, trt_var, mod_var, mod_seq, covs = NULL, 
                              horizons, cluster_var, time_var_name = "stop", 
                              samples = 1000, conf_int = TRUE,
                              ref_group = "no reversion", wt_type = "wt_trunc") {
  
  if (is.null(covs) || length(covs) == 0) {
    plr_formula <- as.formula(paste("event ~", trt_var, "* ns(", mod_var, ", df=3) + ns(interval, df=3)"))
  } else {
    covs_str <- paste0(covs, collapse = " + ")
    plr_formula <- as.formula(paste("event ~", trt_var, "* ns(", mod_var, ", df=3) + ns(interval, df=3) +", covs_str))
  }
  
  cat("Fitting model:", deparse(plr_formula), "\n")
  model <- glm(plr_formula, data = Dat_expanded, family = quasibinomial(), weights = Dat_expanded[[wt_type]])
  
  V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL, sandwich = TRUE, fix = FALSE)
  coefs <- coef(model)
  original_levels <- levels(factor(Dat_expanded[[trt_var]]))
  
  if(is.null(ref_group)) ref_group <- original_levels[1]
  
  risk_fun <- function(coefs, horizon) {
    nd <- Dat_expanded %>% filter(.data[[time_var_name]] <= horizon)
    grid_combos <- expand.grid(trt_level = original_levels, mod_val = mod_seq, stringsAsFactors = FALSE)
    
    risks <- lapply(1:nrow(grid_combos), function(i) {
      trt_l <- grid_combos$trt_level[i]
      mod_v <- grid_combos$mod_val[i]
      
      nd_mod <- nd
      nd_mod[[trt_var]] <- factor(trt_l, levels = original_levels)
      nd_mod[[mod_var]] <- mod_v
      
      mm_mod <- model.matrix(delete.response(terms(model)), nd_mod)
      nd_mod$p_event <- plogis(mm_mod %*% coefs)
      
      df_pred <- nd_mod %>%
        group_by(.data[[cluster_var]]) %>%
        summarise(surv = prod(1 - p_event), .groups = "drop") %>%
        mutate(risk = 1 - surv)
      
      tibble(trt = trt_l, !!mod_var := mod_v, horizon = horizon, risk = mean(df_pred$risk))
    })
    bind_rows(risks)
  }
  
  cat("Calculating point estimates...\n")
  point_risks <- bind_rows(lapply(horizons, function(h) risk_fun(coefs, h)))
  
  # --- NEW: Compute point contrasts when conf_int is FALSE ---
  if(!conf_int) {
    point_contrasts <- point_risks %>%
      group_by(horizon, !!sym(mod_var)) %>%
      mutate(ref_risk = risk[trt == ref_group]) %>%
      filter(trt != ref_group) %>%
      mutate(
        RD = risk - ref_risk,
        RR = risk / ref_risk,
        # Fill dummy CI columns so downstream plotting code doesn't break
        RD_lower = NA, RD_upper = NA, 
        RR_lower = NA, RR_upper = NA,
        risk_lower = NA, risk_upper = NA
      ) %>%
      ungroup()
    
    # Add dummy CIs to risks as well
    point_risks <- point_risks %>% mutate(risk_lower = NA, risk_upper = NA)
    
    return(list(risks = point_risks, contrasts = point_contrasts, ref_group = ref_group))
  }
  
  cat("Running parametric bootstrap (", samples, " samples)...\n", sep="")
  sim_coefs <- mvtnorm::rmvnorm(n = samples, mean = coefs, sigma = V)
  
  sim_list <- lapply(seq_len(nrow(sim_coefs)), function(i) {
    bind_rows(lapply(horizons, function(h) risk_fun(sim_coefs[i,], h)))
  })
  sim_df <- bind_rows(sim_list, .id = "sim")
  
  risk_ci <- sim_df %>%
    group_by(trt, horizon, !!sym(mod_var)) %>%
    summarise(risk_lower = quantile(risk, 0.025), risk_upper = quantile(risk, 0.975), risk = mean(risk), .groups = "drop")
  
  contrast_ci <- sim_df %>%
    group_by(sim, horizon, !!sym(mod_var)) %>%
    mutate(ref_risk = risk[trt == ref_group]) %>%
    filter(trt != ref_group) %>%
    mutate(RD = risk - ref_risk, RR = risk / ref_risk) %>%
    group_by(trt, horizon, !!sym(mod_var)) %>%
    summarise(
      RD_lower = quantile(RD, 0.025), RD_upper = quantile(RD, 0.975), RD = mean(RD),
      RR_lower = quantile(RR, 0.025), RR_upper = quantile(RR, 0.975), RR = mean(RR),
      .groups = "drop"
    )
  
  list(risks = risk_ci, contrasts = contrast_ci, ref_group = ref_group)
}






## overdose ##############################
### MET interruption of 3-4 days ####

cat("  Computing: Overdose - Trial 2 (3-4 days)\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))

long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

max_ret <- quantile(long_dat$ret_length, 0.95, na.rm = TRUE)
ret_seq_grid <- seq(0, max_ret, by = 3)

# MASTER CONTROL SWITCH

RUN_BOOTSTRAP_CI <- FALSE    
N_BOOT_SAMPLES   <- 499   

# Example function call:
risk_em_ci_tr2_od <- predict_cuminc_em(
  Dat_expanded  = long_dat, 
  trt_var       = "trt", 
  mod_var       = "ret_length",     
  mod_seq       = ret_seq_grid,     
  covs          = c("ddose_prev", "prior_episode"),
  horizons      = c(30, 90, 180, 365), 
  cluster_var   = "studyid_rev", 
  time_var_name = "stop",          
  samples       = N_BOOT_SAMPLES,     # Linked to global
  conf_int      = RUN_BOOTSTRAP_CI,   # Linked to global switch
  ref_group     = "no reversion", 
  wt_type       = "wt_trunc"
)


saveRDS(risk_em_ci_tr2_od, paste0(RES_DIR, "Risk/", "risk_tr2_34d_od_em.rds"))


### MET interruption of 5-14 days ####

cat("  Computing: Overdose - Trial 3 (5-14 days)\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_od_acute",
                                 event_var = "sens_od_acute", 
                                 interval_length = 14)

max_ret <- quantile(long_dat$ret_length, 0.95, na.rm = TRUE)
ret_seq_grid <- seq(0, max_ret, by = 3)

risk_em_ci_tr3_od <- predict_cuminc_em(
  Dat_expanded  = long_dat, 
  trt_var       = "trt", 
  mod_var       = "ret_length",     
  mod_seq       = ret_seq_grid,     
  covs          = c("ddose_prev", "prior_episode"),
  horizons      = c(30, 90, 180, 365), 
  cluster_var   = "studyid_rev", 
  time_var_name = "stop",          
  samples       = N_BOOT_SAMPLES,     # Linked to global
  conf_int      = RUN_BOOTSTRAP_CI,   # Linked to global switch
  ref_group     = "no reversion", 
  wt_type       = "wt_trunc"
)

saveRDS(risk_em_ci_tr3_od, paste0(RES_DIR, "Risk/", "risk_tr3_od_em.rds"))


## treatment discontinuation ##############################
### MET interruption of 3-4 days ####

cat("  Computing: Discontinuation - Trial 2 (3-4 days)\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr2_34d_all.rds"))

long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)

max_ret <- quantile(long_dat$ret_length, 0.95, na.rm = TRUE)
ret_seq_grid <- seq(0, max_ret, by = 3)

risk_em_ci_tr2_dis <- predict_cuminc_em(
  Dat_expanded  = long_dat, 
  trt_var       = "trt", 
  mod_var       = "ret_length",     
  mod_seq       = ret_seq_grid,     
  covs          = c("ddose_prev", "prior_episode"),
  horizons      = c(30, 90, 180, 365), 
  cluster_var   = "studyid_rev", 
  time_var_name = "stop",          
  samples       = N_BOOT_SAMPLES,     # Linked to global
  conf_int      = RUN_BOOTSTRAP_CI,   # Linked to global switch
  ref_group     = "no reversion", 
  wt_type       = "wt_trunc"
)

saveRDS(risk_em_ci_tr2_dis, paste0(RES_DIR, "Risk/", "risk_tr2_34d_dis_em.rds"))


### MET interruption of 5-14 days ####

cat("  Computing: Discontinuation - Trial 3 (5-14 days)\n")

ps_fit <- readRDS(paste0(RES_DIR, "PScore/", "ps_weight_tr3_all.rds"))

# FIXED: Changed from od_acute to discont
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)

max_ret <- quantile(long_dat$ret_length, 0.95, na.rm = TRUE)
ret_seq_grid <- seq(0, max_ret, by = 3)

risk_em_ci_tr3_dis <- predict_cuminc_em(
  Dat_expanded  = long_dat, 
  trt_var       = "trt", 
  mod_var       = "ret_length",     
  mod_seq       = ret_seq_grid,     
  covs          = c("ddose_prev", "prior_episode"),
  horizons      = c(30, 90, 180, 365), 
  cluster_var   = "studyid_rev", 
  time_var_name = "stop",          
  samples       = N_BOOT_SAMPLES,     # Linked to global
  conf_int      = RUN_BOOTSTRAP_CI,   # Linked to global switch
  ref_group     = "no reversion", 
  wt_type       = "wt_trunc"
)


saveRDS(risk_em_ci_tr3_dis, paste0(RES_DIR, "Risk/", "risk_tr3_dis_em.rds"))




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



risk_2_od <- readRDS(paste0(RES_DIR, "Risk/", "risk_tr2_34d_od_em.rds"))
risk_3_od <- readRDS(paste0(RES_DIR, "Risk/", "risk_tr3_dis_em.rds"))
risk_2_dis <- readRDS(paste0(RES_DIR, "Risk/", "risk_tr2_34d_dis_em.rds"))
risk_3_dis <- readRDS(paste0(RES_DIR, "Risk/", "risk_tr3_dis_em.rds"))



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
  

jpeg(paste0(RES_DIR, "Figures/", "risk_od_em.jpg"),
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
  
  
  jpeg(paste0(RES_DIR, "Figures/", "risk_dis_em.jpg"),
       width = 3200, height = 3500, res = 300)
print(row_dis_all_final)
dev.off()



