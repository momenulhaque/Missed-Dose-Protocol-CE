predict_cuminc_rmst <- function(Dat_expanded, trt_var, covs, horizons,
                                cluster_var, samples = 1000, conf_int = TRUE,
                                ref_group = "no reversion", wt_type="wt_trunc",
                                rmst_tau = 365) {
  
  # Formula for pooled logistic regression
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  # Fit GLM using quasibinomial to suppress fractional weight warnings
  model <- glm(plr_formula, data = Dat_expanded,
               family = quasibinomial(), weights = Dat_expanded[[wt_type]])
  
  # Cluster-robust covariance
  V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                        sandwich = TRUE, fix = FALSE)
  
  coefs <- coef(model)
  
  # Function to calculate risk at a given horizon
  risk_fun <- function(coef_draw, h) {
    nd <- Dat_expanded %>% filter(stop <= h)
    
    risks <- lapply(unique(Dat_expanded[[trt_var]]), function(trt_level) {
      nd_trt <- nd
      nd_trt[[trt_var]] <- trt_level
      
      mm_trt <- model.matrix(delete.response(terms(model)), nd_trt)
      p_trt <- plogis(mm_trt %*% coef_draw)
      nd_trt$p_event <- p_trt
      
      df_pred <- nd_trt %>%
        group_by(.data[[cluster_var]]) %>%
        summarise(surv = prod(1 - p_event), .groups = "drop") %>%
        mutate(risk = 1 - surv)
      
      tibble(trt = as.character(trt_level), horizon = h, risk = mean(df_pred$risk))
    })
    
    bind_rows(risks)
  }
  
  # Reference group check
  if(is.null(ref_group)) {
    ref_group <- as.character(unique(Dat_expanded[[trt_var]]))[1]
  }
  
  # Point estimates
  point <- bind_rows(lapply(horizons, function(h) risk_fun(coefs, h)))
  
  if(!conf_int) {
    res <- list(risks = point)
    
    # Point estimate for RMST if requested
    if(!is.null(rmst_tau)) {
      if(max(horizons) < rmst_tau) stop("max(horizons) must be >= rmst_tau to compute RMST.")
      
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
        ) %>% dplyr::select(-ref_RMST) # <-- Explicitly namespaced to avoid MASS collision
      
      res$rmst <- point_rmst
    }
    return(res)
  }
  
  # Simulate coefficient draws for CI
  sim_coefs <- mvtnorm::rmvnorm(n = samples, mean = coefs, sigma = V)
  
  sim_list <- lapply(seq_len(nrow(sim_coefs)), function(i) {
    bind_rows(lapply(horizons, function(h) risk_fun(sim_coefs[i,], h)))
  })
  sim_df <- bind_rows(sim_list, .id = "sim")
  
  # Risk CI
  risk_ci <- sim_df %>%
    group_by(trt, horizon) %>%
    summarise(
      risk_lower = quantile(risk, 0.025),
      risk_upper = quantile(risk, 0.975),
      risk = mean(risk),  # point estimate from simulation
      .groups = "drop"
    )
  
  # Compute RD and RR in each simulation
  sim_contrasts <- sim_df %>%
    group_by(sim, horizon) %>%
    mutate(ref_risk = risk[trt == ref_group]) %>%
    filter(trt != ref_group) %>%
    mutate(
      RD = risk - ref_risk,
      RR = risk / ref_risk
    )
  
  # Summarise RD and RR
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
    ref_group = ref_group
  )
  
  # Compute RMST distribution if requested
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


# 7. RMST & CUMULATIVE INCIDENCE ####
# Risk curves and confidence intervals for overdose and discontinuation
# Two trials × Two cohorts

cat("\n=== Computing Cumulative Incidence ===\n")

HORIZONS_CURVE <- c(seq(14, 357, by = 14), 365)
HORIZONS_CI    <- sort(c(HORIZONS_CURVE, c(30, 60, 90, 180)))

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


risk$tr2$first$od$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc", rmst_tau = 365
)

risk$tr2$first$od$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc", rmst_tau = 365
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr2$first$dis$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr2$first$dis$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("po_preOUD", "phy_att_gr2", "ret_length", "sud", "ddose_prev"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
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

risk$tr2$all$od$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$all$od$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"), 
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights, 
                                 time_var = "time_discont",
                                 event_var = "sens_discont", 
                                 interval_length = 14)

risk$tr2$all$dis$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat,
  trt_var = "trt", 
  covs =c("ddose_prev", "prior_episode"),  
  horizons = HORIZONS_CURVE,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = FALSE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

risk$tr2$all$dis$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, 
  trt_var = "trt", 
  covs = c("ddose_prev", "prior_episode"),
  horizons = HORIZONS_CI,
  cluster_var = "studyid_rev", 
  samples = 499, 
  conf_int = TRUE,
  ref_group = "no reversion", 
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr2, paste0(RES_DIR, "Risk/", "risk_tr2_34d_rmst.rds"))


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

risk$tr3$first$od$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$first$od$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$first$dis$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)


risk$tr3$first$dis$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "ret_length"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
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

risk$tr3$all$od$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$all$od$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

# Discontinuation
long_dat <- expand_person_period(ps_fit$data_weights,
                                 time_var  = "time_discont",
                                 event_var = "sens_discont",
                                 interval_length = 14)

risk$tr3$all$dis$curve <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CURVE, cluster_var = "studyid_rev",
  samples = 499, conf_int = FALSE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

risk$tr3$all$dis$ci <- predict_cuminc_rmst(
  Dat_expanded = long_dat, trt_var = "trt",
  covs = c("ddose_prev", "year"),
  horizons = HORIZONS_CI, cluster_var = "studyid_rev",
  samples = 499, conf_int = TRUE, ref_group = "no reversion",
  wt_type = "wt_trunc"
)

rm(long_dat); gc()

saveRDS(risk$tr3, paste0(RES_DIR, "Risk/", "risk_tr3_514d_rmst.rds"))

cat("  ✓ Trial 3 risk estimates saved.\n")
