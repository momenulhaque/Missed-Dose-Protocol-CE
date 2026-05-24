library(tidyverse)
library(survival)
library(haven)
library(Publish)
library(ggplot2)
library(survminer)
library(twang)
require(tableone)
require(survey)
library(lubridate)
library(naniar) # For missing value
library(survival)
library(cobalt)
library(autoCovariateSelection)
library(glmnet)
library(parallel)
library(doParallel)
library(ranger)
library(splines)
library(marginaleffects)
library(sandwich)
library(clubSandwich)
library(CFsurvival)
library(survSuperLearner)
library(earth)
library(randomForest)


expand_person_period <- function(Dat, time_var, event_var, interval_length = 7){
  Dat <- Dat %>%
    mutate(max_int = ceiling(.data[[time_var]] / interval_length))
  
  Dat_expanded <- Dat %>%
    tidyr::uncount(weights = max_int, .id = "interval", .remove = FALSE) %>%
    mutate(
      start = (interval - 1) * interval_length,
      stop  = interval * interval_length,
      event = ifelse(interval == max_int, .data[[event_var]], 0)
    )
  return(Dat_expanded)
}

ps_model_ed <- function(dat, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99, 
                        out_vars, clust_vars,
                        baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                                         "age_40_49", "age_ge50", "rural",
                                         "year", "po_preOUD", "phy_att_gr2", 
                                         "social_assist", "ret_length", "home_grp_b", 
                                         "incar_grp_bl", "drug_acute", "cci_gt0", 
                                         "copd_asthma", "sed_psy_med", "po", "sud", 
                                         "aud", "smd", "hcv", "tob", "udt", "ddose_prev", 
                                         "clin_size", "virtual_care", "pss"),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length")) {
  
  library(dplyr)
  library(nnet)
  
  dat1 <- dat %>% select(all_of(out_vars), all_of(clust_vars), all_of(exposure), all_of(baselinevars))
  
  cont_vars <- continuous_vars[continuous_vars %in% baselinevars]
  
  dat1 <- dat1 %>%
    mutate(across(-c(all_of(clust_vars), all_of(cont_vars), all_of(out_vars)), as.factor))
  
  model_data <- dat1 %>%
    filter(if_all(c(all_of(baselinevars), exposure), ~ !is.na(.)))
  
  # Drop baseline variables with only one level
  valid_baselinevars <- baselinevars[sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)]
  
  # Define formula
  formula_ps <- as.formula(paste(exposure, "~", paste(valid_baselinevars, collapse = " + ")))
  
  # Ensure exposure is a factor
  model_data[[exposure]] <- as.factor(model_data[[exposure]])
  levels_exp <- levels(model_data[[exposure]])
  
  # Fit multinomial PS model
  if (nlevels(model_data[[exposure]]) > 2) {
    ps_fit <- multinom(formula_ps, data = model_data, trace = FALSE)
    ps_probs <- predict(ps_fit, type = "probs")
    
    # Ensure matrix structure
    if (is.vector(ps_probs)) {
      ps_probs <- matrix(ps_probs, ncol = length(levels_exp))
      colnames(ps_probs) <- levels_exp
    }
    
    # Get marginal probabilities (stabilizer)
    marginal_probs <- prop.table(table(model_data[[exposure]]))
    
    # Compute weights
    obs_exposure <- model_data[[exposure]]
    model_data$ps <- mapply(function(i, trt) ps_probs[i, trt], 
                            i = seq_along(obs_exposure), 
                            trt = obs_exposure)
    model_data$stab_numerator <- marginal_probs[as.character(obs_exposure)]
    model_data$wt_stab <- as.numeric(as.character(model_data$stab_numerator / model_data$ps))
    
  } else {
    # Binary logistic regression
    ps_fit <- glm(formula_ps, data = model_data, family = binomial(link = "logit"))
    ps_prob <- predict(ps_fit, type = "response")
    y <- model_data[[exposure]]
    
    model_data$ps <- ifelse(y == 1, ps_prob, 1 - ps_prob)
    stab_num <- rep(mean(y == 1), length(y))
    stab_num <- ifelse(y == 1, stab_num, 1 - stab_num)
    model_data$wt_stab <- stab_num / model_data$ps
  }
  
  # Truncate weights
  p_low <- quantile(model_data$wt_stab, trunc_low)
  p_high <- quantile(model_data$wt_stab, trunc_high)
  model_data$wt_trunc <- pmin(pmax(model_data$wt_stab, p_low), p_high)
  
  return(list(model = ps_fit, data_weights = model_data))
}

ps_model_cbps <- function(dat, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99, 
                          out_vars, clust_vars,
                          baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                                           "age_40_49", "age_ge50", "rural",
                                           "year", "po_preOUD", "phy_att_gr2", 
                                           "social_assist", "ret_length", "home_grp_b", 
                                           "incar_grp_bl", "drug_acute", "cci_gt0", 
                                           "copd_asthma", "sed_psy_med", "po", "sud", 
                                           "aud", "smd", "hcv", "tob", "udt", "ddose_prev", 
                                           "clin_size", "virtual_care", "pss"),
                          continuous_vars = c("age_st", "ddose_prev", "year", "ret_length")) {
  
  library(dplyr)
  library(CBPS)
  
  dat1 <- dat[ , c(out_vars, clust_vars, exposure, baselinevars)]
  cont_vars <- continuous_vars[continuous_vars %in% baselinevars]
  
  dat1 <- dat1 %>%
    mutate(across(-c(all_of(clust_vars), all_of(cont_vars), all_of(out_vars)), as.factor))
  
  model_data <- dat1 %>%
    filter(if_all(all_of(c(baselinevars, exposure)), ~ !is.na(.)))
  
  
  # Drop baseline variables with only one level
  valid_baselinevars <- baselinevars[sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)]
  
  # Define formula
  formula_ps <- as.formula(paste(exposure, "~", paste(valid_baselinevars, collapse = " + ")))
  
  # Ensure exposure is a factor
  model_data[[exposure]] <- as.factor(model_data[[exposure]])
  levels_exp <- levels(model_data[[exposure]])
  
  # Fit CBPS model
  ps_fit <- CBPS::CBPS(formula_ps, data = model_data, ATT = 0, method = "over")
  
  ps_probs <- ps_fit$fitted.values	
  
  obs_exposure <- model_data[[exposure]]
  model_data$ps <- mapply(function(i, trt) ps_probs[i, trt], 
                          i = seq_along(obs_exposure), 
                          trt = obs_exposure)
  # Extract weights
  model_data$wt <- ps_fit$weights
  
  # Truncate weights (if desired)
  p_low <- quantile(model_data$wt, trunc_low)
  p_high <- quantile(model_data$wt, trunc_high)
  model_data$wt_trunc <- pmin(pmax(model_data$wt, p_low), p_high)
  
  return(list(model = ps_fit, data_weights = model_data))
}


ps_model_weight <- function(dat, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99, 
                            out_vars, clust_vars,
                            baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                                             "age_40_49", "age_ge50", "rural",
                                             "year", "po_preOUD", "phy_att_gr2", 
                                             "social_assist", "ret_length", "home_grp_b", 
                                             "incar_grp_bl", "drug_acute", "cci_gt0", 
                                             "copd_asthma", "sed_psy_med", "po", "sud", 
                                             "aud", "smd", "hcv", "tob", "udt", "ddose_prev", 
                                             "clin_size", "virtual_care", "pss"),
                            continuous_vars = c("age_st", "ddose_prev", "year", "ret_length")) {
  
  library(dplyr)
  
  dat1 <- dat[ , c(out_vars, clust_vars, exposure, baselinevars)]
  cont_vars <- continuous_vars[continuous_vars %in% baselinevars]
  
  dat1 <- dat1 %>%
    mutate(across(-c(all_of(clust_vars), all_of(cont_vars), all_of(out_vars)), as.factor))
  
  model_data <- dat1 %>%
    filter(if_all(all_of(c(baselinevars, exposure)), ~ !is.na(.)))
  
  
  # Drop baseline variables with only one level
  valid_baselinevars <- baselinevars[sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)]
  
  # Define formula
  formula_ps <- as.formula(paste(exposure, "~", paste(valid_baselinevars, collapse = " + ")))
  
  # Ensure exposure is a factor
  model_data[[exposure]] <- as.factor(model_data[[exposure]])
  levels_exp <- levels(model_data[[exposure]])
  
  # Fit CBPS model
  ps_fit <- WeightIt::weightit(formula_ps,
                               data = model_data,
                               method = "cbps",
                               estimand = "ATE")
  # Extract weights
  model_data$wt <- ps_fit$weights
  
  # Truncate weights (if desired)
  p_low <- quantile(model_data$wt, trunc_low)
  p_high <- quantile(model_data$wt, trunc_high)
  model_data$wt_trunc <- pmin(pmax(model_data$wt, p_low), p_high)
  
  return(list(model = ps_fit, data_weights = model_data))
}

ps_model_psweight <- function(dat,
                              exposure = "trt",
                              out_vars = NULL,
                              outcome,
                              clust_vars = NULL,
                              baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                                               "age_40_49", "age_ge50", "rural",
                                               "year", "po_preOUD", "phy_att_gr2", 
                                               "social_assist", "ret_length", "home_grp_b", 
                                               "incar_grp_bl", "drug_acute", "cci_gt0", 
                                               "copd_asthma", "sed_psy_med", "po", "sud", 
                                               "aud", "smd", "hcv", "tob", "udt", "ddose_prev", 
                                               "clin_size", "virtual_care", "pss"),
                              continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"),
                              trunc_low = 0, trunc_high = 0.99) {
  
  library(dplyr)
  library(PSweight)
  
  model_data <- dat %>%
    dplyr::select(all_of(c(out_vars, clust_vars, exposure, baselinevars))) %>%
    mutate(across(all_of(setdiff(baselinevars, continuous_vars)), as.factor)) %>%
    filter(if_all(all_of(c(exposure, baselinevars)), ~ !is.na(.)))
  valid_baselinevars <- baselinevars[sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)]
  f_ps <- as.formula(paste(exposure, "~", paste(valid_baselinevars, collapse = " + ")))
  
  # IPW 
  ps_out <- PSweight::PSweight(
    ps.formula = f_ps,
    data = as.data.frame(model_data),
    zname = c("no reversion", "reversion", "increase"),
    yname = outcome,
    weight ="IPW",
    family = binomial(link = "logit"),
    ps.method = "glm"
  )
  
  
  ps_probs <- ps_out$propensity
  marginal_probs <- prop.table(table(model_data[[exposure]]))
  
  obs_exposure <- model_data[[exposure]]
  model_data$ipw_ps <- mapply(function(i, trt) ps_probs[i, trt], 
                              i = seq_along(obs_exposure), 
                              trt = obs_exposure)
  model_data$ipw_wt <- as.numeric(1/ model_data$ipw_ps)
  
  model_data$num_prob <- marginal_probs[as.character(obs_exposure)]
  model_data$ipw_wt_stab <- as.numeric(as.character(model_data$num_prob / model_data$ipw_ps))
  
  p_low <- quantile(model_data$ipw_wt_stab, trunc_low)
  p_high <- quantile(model_data$ipw_wt_stab, trunc_high)
  model_data$ipw_wt_trunc <- pmin(pmax(model_data$ipw_wt_stab, p_low), p_high)
  
  
  
  # overlap weights
  tilting_fn <- apply(ps_probs, 1, function(x) sum(1/x)) 
  
  model_data$ow_wt <- as.numeric(1/ (tilting_fn*model_data$ipw_ps))
  
  
  return(list(
    model = ps_out,
    data_weights = model_data
  ))
}


ps_model_ow <- function(dat, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99, 
                        out_vars, clust_vars,
                        baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                                         "age_40_49", "age_ge50", "rural",
                                         "year", "po_preOUD", "phy_att_gr2", 
                                         "social_assist", "ret_length", "home_grp_b", 
                                         "incar_grp_bl", "drug_acute", "cci_gt0", 
                                         "copd_asthma", "sed_psy_med", "po", "sud", 
                                         "aud", "smd", "hcv", "tob", "udt", "ddose_prev", 
                                         "clin_size", "virtual_care", "pss"),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length")) {
  
  library(dplyr)
  library(nnet)
  
  model_data <- dat %>%
    dplyr::select(all_of(c(out_vars, clust_vars, exposure, baselinevars))) %>%
    mutate(across(all_of(setdiff(baselinevars, continuous_vars)), as.factor)) %>%
    filter(if_all(all_of(c(exposure, baselinevars)), ~ !is.na(.)))
  
  valid_baselinevars <- baselinevars[sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)]
  formula_ps <- as.formula(paste(exposure, "~", paste(valid_baselinevars, collapse = " + ")))
  
  # Ensure exposure is a factor
  model_data[[exposure]] <- as.factor(model_data[[exposure]])
  levels_exp <- levels(model_data[[exposure]])
  
  # Fit multinomial PS model
  ps_fit <- multinom(formula_ps, data = model_data, trace = FALSE)
  ps_probs <- predict(ps_fit, type = "probs")
  
  ps_probs <- matrix(ps_probs, ncol = length(levels_exp))
  colnames(ps_probs) <- levels_exp
  
  # Get marginal probabilities (stabilizer)
  marginal_probs <- prop.table(table(model_data[[exposure]]))
  
  # Compute weights
  obs_exposure <- model_data[[exposure]]
  model_data$ps <- mapply(function(i, trt) ps_probs[i, trt], 
                          i = seq_along(obs_exposure), 
                          trt = obs_exposure)
  model_data$stab_numerator <- marginal_probs[as.character(obs_exposure)]
  model_data$wt_stab <- as.numeric(as.character(model_data$stab_numerator / model_data$ps))
  
  # Truncate weights
  p_low <- quantile(model_data$wt_stab, trunc_low)
  p_high <- quantile(model_data$wt_stab, trunc_high)
  model_data$wt_trunc <- pmin(pmax(model_data$wt_stab, p_low), p_high)
  
  ## ow
  tilting_fn <- apply(ps_probs, 1, function(x) sum(1/x)) 
  
  model_data$ow_wt <- as.numeric(1/ (tilting_fn*model_data$ps))
  
  return(list(model = ps_fit, data_weights = model_data))
}



fit_PLR <- function(Dat_expanded, trt_var = "trt", covs = NULL){
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  glm(plr_formula, data = Dat_expanded,
      family = binomial(), weights = wt_trunc)
}

fit_PLR_HR <- function(Dat_expanded, trt_var = "trt", covs = NULL, cluster_var = "studyid_rev"){
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  model <- glm(plr_formula, data = Dat_expanded,
               family = binomial(), weights = wt_trunc)
  
  var_matrix <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                                 sandwich = TRUE,
                                 fix = FALSE)
  se <- sqrt(diag(var_matrix))
  output <- data.frame(names = names(model$coefficients))
  output$estimate <- model$coefficients
  output$robust_se <- se[names(model$coefficients)]
  output$`2.5%` <- output$estimate - (1.96 * output$robust_se)
  output$`97.5%` <- output$estimate + (1.96 * output$robust_se)
  output$OR <- exp(output$estimate)
  output$lower = exp(output$`2.5%`)
  output$upper = exp(output$`97.5%`)
  output$z <- output$estimate / output$robust_se
  output$p_value <- 2 * (1 - pnorm(abs(output$z)))
  output <- output %>% 
    filter(grepl("^trt", names))
  return(output)
}

fit_PLR_HR_v1 <- function(Dat_expanded, trt_var = "trt", covs = NULL, cluster_var = "studyid_rev"){
  covs1 <- paste0(c(trt_var, "ns(interval, df = 2)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  model <- glm(plr_formula, data = Dat_expanded,
               family = binomial(), weights = wt_trunc)
  
  var_matrix <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                                 sandwich = TRUE,
                                 fix = FALSE)
  se <- sqrt(diag(var_matrix))
  output <- data.frame(names = names(model$coefficients))
  output$estimate <- model$coefficients
  output$robust_se <- se[names(model$coefficients)]
  output$`2.5%` <- output$estimate - (1.96 * output$robust_se)
  output$`97.5%` <- output$estimate + (1.96 * output$robust_se)
  output$OR <- exp(output$estimate)
  output$lower = exp(output$`2.5%`)
  output$upper = exp(output$`97.5%`)
  output$z <- output$estimate / output$robust_se
  output$p_value <- 2 * (1 - pnorm(abs(output$z)))
  output <- output %>% 
    filter(grepl("^trt", names))
  return(output)
}

iv_fit_PLR_HR <- function(Dat_expanded, trt_var = "trt", covs = NULL, cluster_var = "studyid_rev"){
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  model <- glm(plr_formula, data = Dat_expanded,
               family = binomial())
  
  var_matrix <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                                 sandwich = TRUE,
                                 fix = FALSE)
  se <- sqrt(diag(var_matrix))
  output <- data.frame(names = names(model$coefficients))
  output$estimate <- model$coefficients
  output$robust_se <- se[names(model$coefficients)]
  output$`2.5%` <- output$estimate - (1.96 * output$robust_se)
  output$`97.5%` <- output$estimate + (1.96 * output$robust_se)
  output$OR <- exp(output$estimate)
  output$lower = exp(output$`2.5%`)
  output$upper = exp(output$`97.5%`)
  output$z <- output$estimate / output$robust_se
  output$p_value <- 2 * (1 - pnorm(abs(output$z)))
  output <- output %>% 
    filter(grepl("^trt", names))
  return(output)
}


predict_cuminc <- function(Dat_expanded, trt_var, covs, horizons,
                           cluster_var, samples = 1000, conf_int = TRUE,
                           ref_group = "no reversion") {
  
  # Formula for pooled logistic regression
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  # Fit GLM
  model <- glm(plr_formula, data = Dat_expanded,
               family = binomial(), weights = wt_trunc)
  
  # Cluster-robust covariance
  V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                        sandwich = TRUE, fix = FALSE)
  
  coefs <- coef(model)
  
  # Function to calculate risk at a given horizon
  risk_fun <- function(coefs, horizon) {
    nd <- Dat_expanded %>% filter(stop <= horizon)
    
    risks <- lapply(unique(Dat_expanded[[trt_var]]), function(trt_level) {
      nd_trt <- nd
      nd_trt[[trt_var]] <- trt_level
      
      mm_trt <- model.matrix(delete.response(terms(model)), nd_trt)
      p_trt <- plogis(mm_trt %*% coefs)
      nd_trt$p_event <- p_trt
      
      df_pred <- nd_trt %>%
        group_by(.data[[cluster_var]]) %>%
        summarise(surv = prod(1 - p_event), .groups = "drop") %>%
        mutate(risk = 1 - surv)
      
      tibble(trt = trt_level, horizon = horizon, risk = mean(df_pred$risk))
    })
    
    bind_rows(risks)
  }
  
  # Point estimates
  point <- bind_rows(lapply(horizons, function(h) risk_fun(coefs, h)))
  
  if(!conf_int) return(point)
  
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
  
  # Reference group
  if(is.null(ref_group)) {
    ref_group <- levels(factor(Dat_expanded[[trt_var]]))[1]
  }
  
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
  
  list(
    risks = risk_ci,
    contrasts = contrast_ci,
    ref_group = ref_group
  )
}



predict_cuminc_ow <- function(Dat_expanded, trt_var, covs, horizons,
                           cluster_var, samples = 1000, conf_int = TRUE,
                           ref_group = "no reversion", wt_type="wt_trunc") {
  
  # Formula for pooled logistic regression
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  # Fit GLM
  model <- glm(plr_formula, data = Dat_expanded,
               family = binomial(), weights = Dat_expanded[[wt_type]])
  
  # Cluster-robust covariance
  V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
                        sandwich = TRUE, fix = FALSE)
  
  coefs <- coef(model)
  
  # Function to calculate risk at a given horizon
  risk_fun <- function(coefs, horizon) {
    nd <- Dat_expanded %>% filter(stop <= horizon)
    
    risks <- lapply(unique(Dat_expanded[[trt_var]]), function(trt_level) {
      nd_trt <- nd
      nd_trt[[trt_var]] <- trt_level
      
      mm_trt <- model.matrix(delete.response(terms(model)), nd_trt)
      p_trt <- plogis(mm_trt %*% coefs)
      nd_trt$p_event <- p_trt
      
      df_pred <- nd_trt %>%
        group_by(.data[[cluster_var]]) %>%
        summarise(surv = prod(1 - p_event), .groups = "drop") %>%
        mutate(risk = 1 - surv)
      
      tibble(trt = trt_level, horizon = horizon, risk = mean(df_pred$risk))
    })
    
    bind_rows(risks)
  }
  
  # Point estimates
  point <- bind_rows(lapply(horizons, function(h) risk_fun(coefs, h)))
  
  if(!conf_int) return(point)
  
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
  
  # Reference group
  if(is.null(ref_group)) {
    ref_group <- levels(factor(Dat_expanded[[trt_var]]))[1]
  }
  
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
  
  list(
    risks = risk_ci,
    contrasts = contrast_ci,
    ref_group = ref_group
  )
}



# predict_cuminc <- function(Dat_expanded, trt_var, covs, horizons,
#                            cluster_var, samples = 100, conf_int = TRUE,
#                            ref_group = "no reversion"){
#   
#   covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
#   plr_formula <- as.formula(paste("event ~", covs1))
#   
#   model <- glm(plr_formula, data = Dat_expanded,
#                family = binomial(), weights = wt_trunc)
#   
#   V <- sandwich::vcovCL(model, cluster = Dat_expanded[[cluster_var]], type = NULL,
#                                  sandwich = TRUE,
#                                  fix = FALSE)
# 
#   coefs <- coef(model)
#   
#   
#   risk_fun <- function(coefs, horizon){
#     nd <- Dat_expanded %>% filter(stop <= horizon)
#     
#     risks <- lapply(unique(Dat_expanded[[trt_var]]), function(trt_level){
#       nd_trt <- nd
#       nd_trt[[trt_var]] <- trt_level
#       
#       mm_trt <- model.matrix(delete.response(terms(model)), nd_trt)
#       p_trt <- plogis(mm_trt %*% coefs)
#       nd_trt$p_event <- p_trt
#       
#       df_pred <- nd_trt %>%
#         group_by(.data[[cluster_var]]) %>%
#         summarise(surv = prod(1 - p_event), .groups = "drop") %>%
#         mutate(risk = 1 - surv)
#       
#       tibble(trt = trt_level, horizon = horizon, risk = mean(df_pred$risk))
#     })
#     
#     bind_rows(risks)
#   }
# 
#   
#   point <- bind_rows(lapply(horizons, function(h) risk_fun(coefs, h)))
#   
#   if(!conf_int){
#     return(point)
#   }
#   
#   sim_coefs <- mvtnorm::rmvnorm(n = samples, mean = coefs, sigma = V)
#   
#   sim_list <- lapply(seq_len(nrow(sim_coefs)), function(i){
#     bind_rows(lapply(horizons, function(h) risk_fun(sim_coefs[i,], h)))
#   })
#   sim_df <- bind_rows(sim_list, .id = "sim")
#   
#   risk_ci <- sim_df %>%
#     group_by(trt, horizon) %>%
#     summarise(
#       lower = quantile(risk, 0.025),
#       upper = quantile(risk, 0.975),
#       .groups = "drop"
#     ) %>%
#     left_join(point, by = c("trt","horizon"))
#   
#   if(is.null(ref_group)){
#     ref_group <- levels(factor(Dat_expanded[[trt_var]]))[1]  # default reference
#   }
#   
#   # compute contrasts in each simulation draw
#   sim_contrasts <- sim_df %>%
#     group_by(sim, horizon) %>%
#     mutate(ref_risk = risk[trt == ref_group]) %>%
#     filter(trt != ref_group) %>%
#     mutate(
#       RD = risk - ref_risk,
#       RR = risk / ref_risk
#     )
#   
#   # summarise contrasts
#   contrast_ci <- sim_contrasts %>%
#     group_by(trt, horizon) %>%
#     summarise(
#       RD_lower = quantile(RD, 0.025),
#       RD_upper = quantile(RD, 0.975),
#       RR_lower = quantile(RR, 0.025),
#       RR_upper = quantile(RR, 0.975),
#       .groups = "drop"
#     ) %>%
#     # add point estimates
#     left_join(
#       point %>%
#         group_by(horizon) %>%
#         mutate(ref_risk = risk[trt == ref_group]) %>%
#         filter(trt != ref_group) %>%
#         mutate(
#           RD = risk - ref_risk,
#           RR = risk / ref_risk
#         ) %>%
#         select(trt, horizon, RD, RR),
#       by = c("trt","horizon")
#     )
#   
#   list(risks = risk_ci, contrasts = contrast_ci, ref_group = ref_group)
# }



predict_cuminc_coxPH <- function(Dat, time_var, event_var, covs = c("ret_length", "ddose_prev"),
                                 weight_var = "wt_trunc", cluster_var = "studyid_rev",
                                 horizons_risk = c(seq(14, 357, by=14), 365),
                                 horizons_ci = c(30, 60, 90, 180, 365)) {
  
  horizons <- sort(unique(c(horizons_risk, horizons_ci)))
  trt_levels <- unique(Dat$trt)
  ref_trt <- "no reversion"
  other_trts <- setdiff(trt_levels, ref_trt)
  
  points_list <- list()
  contrasts_list <- list()
  
  for(group in other_trts){
    dat <- Dat %>% filter(trt %in% c(ref_trt, group))
    dat$trt <- droplevels(dat$trt)
    dat$trt <- as.numeric(dat$trt == group)
    
    # Fit Cox model
    f <- as.formula(
      paste0("Surv(", time_var, ", ", event_var, ") ~ trt + ", paste(covs, collapse = "+"))
    )
    cox_model <- coxph(f, data = dat, weights = dat[[weight_var]], method = "breslow")
    
    # Standardized estimates
    res <- stdCoxph(fit = cox_model, data = dat, X = "trt",
                    x = c(0, 1), t = horizons,
                    clusterid = cluster_var)
    
    # functions 
    ref_risk_fun <- function(est) 1 - est[, 1]
    risk_fun <- function(est) 1 - est[, 2]
    rd_fun <- function(est) apply(1 - est, 1, function(x) x[2] - x[1])
    rr_fun <- function(est) apply(1 - est, 1, function(x) x[2] / x[1])
    
    
    # estimation
    ref_risk <- ref_risk_fun(res$est) 
    ref_risk <- data.frame(risk=ref_risk) %>% mutate(horizon = horizons,
                                                     trt = ref_trt)
    
    
    risk <-  risk_fun(res$est)
    cum_inc <- data.frame(risk=risk) %>% mutate(horizon = horizons,
                                                trt = group)
    
    
    rd_point <- rd_fun(res$est)
    rr_point <- rr_fun(res$est)
    
    ref_risk_ci <- data.frame(confint(res, fun = ref_risk_fun))
    risk_ci <-  data.frame(confint(res, fun = risk_fun))
    
    rd_ci <- data.frame(confint(res, fun = rd_fun))
    rr_ci <- data.frame(confint(res, fun = rr_fun))
    
    contrast <- tibble(
      trt = group,
      horizon = horizons,
      RD = rd_point,
      RD_lower = rd_ci$lower,
      RD_upper = rd_ci$upper,
      RR = rr_point,
      RR_lower = rr_ci$lower,
      RR_upper = rr_ci$upper
    ) %>% filter(horizon %in% horizons_ci)
    
    points_list[[ref_trt]] <- bind_cols(ref_risk, ref_risk_ci)
    points_list[[group]] <- bind_cols(cum_inc, risk_ci)
    contrasts_list[[group]] <- contrast
  }
  
  points <- bind_rows(points_list) %>%  dplyr::arrange(trt, horizon) %>% 
    filter(horizon %in% horizons_risk) 
  
  
  risk <- bind_rows(points_list) %>%  dplyr::arrange(trt, horizon) %>% 
    filter(horizon %in% horizons_ci) %>% 
    rename(
      risk_lower = lower,
      risk_upper = upper
    ) %>% select(trt, horizon, risk_lower, risk_upper, risk)
  
  contrasts <- bind_rows(contrasts_list) %>% select(trt, horizon, RD_lower, RD_upper, RD,  
                                                    RR_lower, RR_upper, RR)
  
  return(list(points = points, risk=risk, contrasts = contrasts, ref_group="no reversion"))
}



### hdPS related ###########
ps_model_lasso <- function(dat, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99, 
                           out_vars, clust_vars,
                           baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                                            "age_40_49", "age_ge50", "rural",
                                            "year", "po_preOUD", "phy_att_gr2", 
                                            "social_assist", "ret_length", "home_grp_b", 
                                            "incar_grp_bl", "drug_acute", "cci_gt0", 
                                            "copd_asthma", "sed_psy_med", "po", "sud", 
                                            "aud", "smd", "hcv", "tob", "udt", "ddose_prev", 
                                            "clin_size", "virtual_care", "pss"),
                           continuous_vars = c("age_st", "ddose_prev", "year", "ret_length")) {
  
  library(dplyr)
  library(glmnet)

  dat1 <- dat[, c(out_vars, clust_vars, exposure, baselinevars)]
  
  cont_vars <- continuous_vars[continuous_vars %in% baselinevars]
  
  dat1 <- dat1 %>%
    mutate(across(-c(all_of(clust_vars), all_of(cont_vars), all_of(out_vars)), as.factor))
  
  model_data <- dat1 %>% 
    filter(if_all(c(all_of(baselinevars), exposure), ~ !is.na(.)))
  
  # Drop variables with only one level (after filtering)
  valid_baselinevars <- baselinevars[sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)]
  
  # Prepare design matrix and response
  X <- model.matrix(~ ., data = model_data[, valid_baselinevars])[,-1]
  y <- as.factor(model_data[[exposure]])
  levels_y <- levels(y)
  model_data$exposure <- y
  
  # Fit LASSO model
  set.seed(123)
  ps_fit <- cv.glmnet(x = X, y = y, alpha = 1, family = if (nlevels(y) > 2) "multinomial" else "binomial")
  
  # Predict probabilities at lambda.min
  ps_probs <- predict(ps_fit, newx = X, s = "lambda.min", type = "response")
  
  # Calculate stabilized weights
  if (nlevels(y) > 2) {
    exposure_index <- match(as.character(y), levels_y)
    ps_vec <- sapply(1:length(y), function(i) ps_probs[i, exposure_index[i], 1])
    marginal_probs <- prop.table(table(y))
    stab_num <- marginal_probs[as.character(y)]
  } else {
    ps_vec <- as.vector(ps_probs)
    exposure_num <- as.numeric(as.character(y))
    stab_num <- rep(mean(exposure_num), length(exposure_num))
    ps_vec <- ifelse(exposure_num == 1, ps_vec, 1 - ps_vec)
  }
  
  wt <- as.numeric(as.character(stab_num / ps_vec))
  
  # Truncate weights
  p_low <- quantile(wt, trunc_low)
  p_high <- quantile(wt, trunc_high)
  wt_trunc <- pmin(pmax(wt, p_low), p_high)
  
  
  model_data$ps <- ps_vec
  model_data$wt_stab <- wt
  model_data$wt_trunc <- wt_trunc
  
  return(list(model = ps_fit, data_weights = model_data))
}


ps_model_lasso_ext <- function(dat, 
                               exposure = "trt", 
                               trunc_low = 0.01, 
                               trunc_high = 0.99, 
                               out_vars, 
                               clust_vars,
                               baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                                                "age_40_49", "age_ge50", "rural",
                                                "year", "po_preOUD", "phy_att_gr2", 
                                                "social_assist", "ret_length", "home_grp_b", 
                                                "incar_grp_bl", "drug_acute", "cci_gt0", 
                                                "copd_asthma", "sed_psy_med", "po", "sud", 
                                                "aud", "smd", "hcv", "tob", "udt", "ddose_prev", 
                                                "clin_size", "virtual_care", "pss"),
                               continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"),
                               interact_vars = NULL,    # user-specified terms like "x:y"
                               squared_vars = NULL) {
  
  library(dplyr)
  library(glmnet)
  
  dat1 <- dat[, c(out_vars, clust_vars, exposure, baselinevars)]
  
  cont_vars <- continuous_vars[continuous_vars %in% baselinevars]
  
  dat1 <- dat1 %>%
    mutate(across(-c(all_of(clust_vars), all_of(cont_vars), all_of(out_vars)), as.factor))
  
  model_data <- dat1 %>% 
    filter(if_all(c(all_of(baselinevars), exposure), ~ !is.na(.)))
  
  # Drop variables with only one level
  valid_baselinevars <- baselinevars[sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)]
  
  # --- Add squared terms ---
  if (!is.null(squared_vars)) {
    for (v in squared_vars) {
      if (v %in% names(model_data)) {
        sq_name <- paste0(v, "_sq")
        model_data[[sq_name]] <- model_data[[v]]^2
        valid_baselinevars <- c(valid_baselinevars, sq_name)
      }
    }
  }
  
  # --- Add user-specified interaction terms ---
  if (!is.null(interact_vars)) {
    for (term in interact_vars) {
      vars <- unlist(strsplit(term, ":"))
      if (length(vars) == 2 && all(vars %in% names(model_data))) {
        inter_name <- paste(vars, collapse = "_x_")
        model_data[[inter_name]] <- model_data[[vars[1]]] * model_data[[vars[2]]]
        valid_baselinevars <- c(valid_baselinevars, inter_name)
      } else {
        warning(paste("Skipping interaction", term, "- variable(s) not found in data"))
      }
    }
  }
  
  # --- Prepare design matrix and response ---
  X <- model.matrix(~ ., data = model_data[, valid_baselinevars])[,-1]
  y <- as.factor(model_data[[exposure]])
  levels_y <- levels(y)
  model_data$exposure <- y
  
  # --- Fit LASSO model ---
  set.seed(123)
  ps_fit <- cv.glmnet(
    x = X, 
    y = y, 
    alpha = 1, 
    family = if (nlevels(y) > 2) "multinomial" else "binomial"
  )
  
  # --- Predict probabilities ---
  ps_probs <- predict(ps_fit, newx = X, s = "lambda.min", type = "response")
  
  # --- Stabilized weights ---
  if (nlevels(y) > 2) {
    exposure_index <- match(as.character(y), levels_y)
    ps_vec <- sapply(1:length(y), function(i) ps_probs[i, exposure_index[i], 1])
    marginal_probs <- prop.table(table(y))
    stab_num <- marginal_probs[as.character(y)]
  } else {
    ps_vec <- as.vector(ps_probs)
    exposure_num <- as.numeric(as.character(y))
    stab_num <- rep(mean(exposure_num), length(exposure_num))
    ps_vec <- ifelse(exposure_num == 1, ps_vec, 1 - ps_vec)
  }
  
  wt <- as.numeric(as.character(stab_num / ps_vec))
  
  # --- Truncate weights ---
  p_low <- quantile(wt, trunc_low)
  p_high <- quantile(wt, trunc_high)
  wt_trunc <- pmin(pmax(wt, p_low), p_high)
  
  model_data$ps <- ps_vec
  model_data$wt_stab <- wt
  model_data$wt_trunc <- wt_trunc
  
  return(list(model = ps_fit, data_weights = model_data))
}



get_dat_all_incid_m <- function(dat1){
  
  id_set <- dat1 %>% pull(unique_id)
  
  # Helper to read, mutate, and filter SAS data by unique_id
  read_and_filter <- function(file, code_col, dim_name) {
    df <- read_sas(file) %>%
      mutate(unique_id = paste0(studyid_rev, "_", time0)) %>%
      filter(paste0(studyid_rev, "_", time0) %in% id_set) %>%
      transmute(
        studyid_rev, time0,
        unique_id,
        code = na_if(as.character(.data[[code_col]]), ""),
        dim = dim_name
      ) %>%
      na.omit()
    return(df)
  }
  
  dat.proxy <- bind_rows(
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_diagx_incid_m.sas7bdat", "diagx_trunc", "diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_proc_incid_m.sas7bdat", "icode_trunc", "proc"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_diag_incid_m.sas7bdat", "diagcd_trunc", "msp_diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_fee_incid_m.sas7bdat", "fitm", "msp_fee"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Pnet_drug_rev_incid_m.sas7bdat", "TC5_CD_trunc", "pnet")
  )
  
  # # Generate candidate empirical covariates
  out1 <- get_candidate_covariates(
    df = dat.proxy,
    domainVarname = "dim",
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set,
    n = 500,
    min_num_patients = 20
  )$covars_data
  
  
  # Assess recurrence
  out2 <- get_recurrence_covariates(
    df = out1,
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set
  )$recurrence_data
  
  # Merge with input
  dat.all <- dat1 %>%
    left_join(out2, by = "unique_id")
  
  return(dat.all)
}

get_dat_all_prev_m <- function(dat1){
  
  id_set <- dat1 %>% pull(unique_id)
  
  # Helper to read, mutate, and filter SAS data by unique_id
  read_and_filter <- function(file, code_col, dim_name) {
    df <- read_sas(file) %>%
      mutate(unique_id = paste0(studyid_rev, "_", time0)) %>%
      filter(paste0(studyid_rev, "_", time0) %in% id_set) %>%
      transmute(
        studyid_rev, time0,
        unique_id,
        code = na_if(as.character(.data[[code_col]]), ""),
        dim = dim_name
      ) %>%
      na.omit()
    return(df)
  }
  
  dat.proxy <- bind_rows(
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_diagx_prev_m.sas7bdat", "diagx_trunc", "diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_proc_prev_m.sas7bdat", "icode_trunc", "proc"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_diag_prev_m.sas7bdat", "diagcd_trunc", "msp_diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_fee_prev_m.sas7bdat", "fitm", "msp_fee"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Pnet_drug_rev_prev_m.sas7bdat", "TC5_CD_trunc", "pnet")
  )
  
  # # Generate candidate empirical covariates
  out1 <- get_candidate_covariates(
    df = dat.proxy,
    domainVarname = "dim",
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set,
    n = 500,
    min_num_patients = 20
  )$covars_data
  
  
  # Assess recurrence
  out2 <- get_recurrence_covariates(
    df = out1,
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set
  )$recurrence_data
  
  # Merge with input
  dat.all <- dat1 %>%
    left_join(out2, by = "unique_id")
  
  return(dat.all)
}

get_dat_all_rprev_m <- function(dat1){
  
  id_set <- dat1 %>% pull(unique_id)
  
  # Helper to read, mutate, and filter SAS data by unique_id
  read_and_filter <- function(file, code_col, dim_name) {
    df <- read_sas(file) %>%
      mutate(unique_id = paste0(studyid_rev, "_", time0)) %>%
      filter(paste0(studyid_rev, "_", time0) %in% id_set) %>%
      transmute(
        studyid_rev, time0,
        unique_id,
        code = na_if(as.character(.data[[code_col]]), ""),
        dim = dim_name
      ) %>%
      na.omit()
    return(df)
  }
  
  dat.proxy <- bind_rows(
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_diagx_rprev_m.sas7bdat", "diagx_trunc", "diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_proc_rprev_m.sas7bdat", "icode_trunc", "proc"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_diag_rprev_m.sas7bdat", "diagcd_trunc", "msp_diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_fee_rprev_m.sas7bdat", "fitm", "msp_fee"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Pnet_drug_rev_rprev_m.sas7bdat", "TC5_CD_trunc", "pnet")
  )
  
  # # Generate candidate empirical covariates
  out1 <- get_candidate_covariates(
    df = dat.proxy,
    domainVarname = "dim",
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set,
    n = 100,
    min_num_patients = 150
  )$covars_data
  
  
  # Assess recurrence
  out2 <- get_recurrence_covariates(
    df = out1,
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set
  )$recurrence_data
  
  # Merge with input
  dat.all <- dat1 %>%
    left_join(out2, by = "unique_id")
  
  return(dat.all)
}

get_dat_all_incid_b <- function(dat1){
  
  id_set <- dat1 %>% pull(unique_id)
  
  # Helper to read, mutate, and filter SAS data by unique_id
  read_and_filter <- function(file, code_col, dim_name) {
    df <- read_sas(file) %>%
      mutate(unique_id = paste0(studyid_rev, "_", time0)) %>%
      filter(paste0(studyid_rev, "_", time0) %in% id_set) %>%
      transmute(
        studyid_rev, time0,
        unique_id,
        code = na_if(as.character(.data[[code_col]]), ""),
        dim = dim_name
      ) %>%
      na.omit()
    return(df)
  }
  
  dat.proxy <- bind_rows(
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_diagx_incid_b.sas7bdat", "diagx_trunc", "diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_proc_incid_b.sas7bdat", "icode_trunc", "proc"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_diag_incid_b.sas7bdat", "diagcd_trunc", "msp_diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_fee_incid_b.sas7bdat", "fitm", "msp_fee"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Pnet_drug_rev_incid_b.sas7bdat", "TC5_CD_trunc", "pnet")
  )
  
  # # Generate candidate empirical covariates
  out1 <- get_candidate_covariates(
    df = dat.proxy,
    domainVarname = "dim",
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set,
    n = 500,
    min_num_patients = 20
  )$covars_data
  
  
  # Assess recurrence
  out2 <- get_recurrence_covariates(
    df = out1,
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set
  )$recurrence_data
  
  # Merge with input
  dat.all <- dat1 %>%
    left_join(out2, by = "unique_id")
  
  return(dat.all)
}

get_dat_all_prev_b <- function(dat1){
  
  id_set <- dat1 %>% pull(unique_id)
  
  # Helper to read, mutate, and filter SAS data by unique_id
  read_and_filter <- function(file, code_col, dim_name) {
    df <- read_sas(file) %>%
      mutate(unique_id = paste0(studyid_rev, "_", time0)) %>%
      filter(paste0(studyid_rev, "_", time0) %in% id_set) %>%
      transmute(
        studyid_rev, time0,
        unique_id,
        code = na_if(as.character(.data[[code_col]]), ""),
        dim = dim_name
      ) %>%
      na.omit()
    return(df)
  }
  
  dat.proxy <- bind_rows(
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_diagx_prev_b.sas7bdat", "diagx_trunc", "diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_proc_prev_b.sas7bdat", "icode_trunc", "proc"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_diag_prev_b.sas7bdat", "diagcd_trunc", "msp_diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_fee_prev_b.sas7bdat", "fitm", "msp_fee"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Pnet_drug_rev_prev_b.sas7bdat", "TC5_CD_trunc", "pnet")
  )
  
  # # Generate candidate empirical covariates
  out1 <- get_candidate_covariates(
    df = dat.proxy,
    domainVarname = "dim",
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set,
    n = 500,
    min_num_patients = 20
  )$covars_data
  
  
  # Assess recurrence
  out2 <- get_recurrence_covariates(
    df = out1,
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set
  )$recurrence_data
  
  # Merge with input
  dat.all <- dat1 %>%
    left_join(out2, by = "unique_id")
  
  return(dat.all)
}

get_dat_all_rprev_b <- function(dat1){
  
  id_set <- dat1 %>% pull(unique_id)
  
  # Helper to read, mutate, and filter SAS data by unique_id
  read_and_filter <- function(file, code_col, dim_name) {
    df <- read_sas(file) %>%
      mutate(unique_id = paste0(studyid_rev, "_", time0)) %>%
      filter(paste0(studyid_rev, "_", time0) %in% id_set) %>%
      transmute(
        studyid_rev, time0,
        unique_id,
        code = na_if(as.character(.data[[code_col]]), ""),
        dim = dim_name
      ) %>%
      na.omit()
    return(df)
  }
  
  dat.proxy <- bind_rows(
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_diagx_rprev_b.sas7bdat", "diagx_trunc", "diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_proc_rprev_b.sas7bdat", "icode_trunc", "proc"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_diag_rprev_b.sas7bdat", "diagcd_trunc", "msp_diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_fee_rprev_b.sas7bdat", "fitm", "msp_fee"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Pnet_drug_rev_rprev_b.sas7bdat", "TC5_CD_trunc", "pnet")
  )
  
  # # Generate candidate empirical covariates
  out1 <- get_candidate_covariates(
    df = dat.proxy,
    domainVarname = "dim",
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set,
    n = 200,
    min_num_patients = 100
  )$covars_data
  
  
  # Assess recurrence
  out2 <- get_recurrence_covariates(
    df = out1,
    eventCodeVarname = "code",
    patientIdVarname = "unique_id",
    patientIdVector = id_set
  )$recurrence_data
  
  # Merge with input
  dat.all <- dat1 %>%
    left_join(out2, by = "unique_id")
  
  return(dat.all)
}


select_top200_empirical_covariates <- function(dat.all,
                                               time_var = "time_discont",
                                               event_var = "sens_discont",
                                               max_vars = 200,
                                               nfolds = 10,
                                               parallel = TRUE,
                                               vars.empirical) {
  library(survival)
  library(glmnet)
  library(doParallel)
  library(dplyr)
  
  valid_vars.empirical <- vars.empirical[sapply(dat.all[vars.empirical], function(x) length(unique(x)) > 1)]
  
  
  # If ≤ max_vars, return all empirical variables directly
  if (length(valid_vars.empirical) <= max_vars) {
    return(valid_vars.empirical)
  }
  
  # Prepare design matrix and response
  formula.out <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", paste(valid_vars.empirical, collapse = " + ")))
  X <- model.matrix(formula.out, data = dat.all)[, -1]
  Y <- Surv(dat.all[[time_var]], dat.all[[event_var]])
  
  if(parallel){
    # Set up parallel LASSO fitting
    n_cores <- parallel::detectCores()
    cl <- makeCluster(n_cores - 1)
    registerDoParallel(cl)
    
    set.seed(123)
    fit.lasso <- cv.glmnet(x = X, y = Y, nfolds = nfolds, parallel = parallel, alpha = 1, family = "cox")
    stopCluster(cl)
  }
  
  set.seed(123)
  fit.lasso <- cv.glmnet(x = X, y = Y, nfolds = nfolds, parallel = FALSE, alpha = 1, family = "cox")
  
  
  # Extract non-zero coefficients at lambda.min
  empvars.lasso <- coef(fit.lasso, s = fit.lasso$lambda.min)
  empvars.lasso <- as.matrix(empvars.lasso)
  
  # Convert to data frame and rank by absolute coefficient magnitude
  empvars.lasso <- data.frame(
    vars = rownames(empvars.lasso),
    coef = empvars.lasso[, 1],
    stringsAsFactors = FALSE
  ) %>%
    mutate(coef.abs = abs(coef)) %>%
    arrange(desc(coef.abs))
  
  # Select top N variables
  top_vars <- empvars.lasso$vars[1:max_vars]
  
  return(top_vars)
}

ps_model_xgb <- function(dat=ps_fit$data_weights,
                         exposure = "trt", 
                         trunc_low = 0, 
                         trunc_high = 0.99, 
                         out_vars, clust_vars,
                         baselinevars,
                         continuous_vars = continuous_vars,
                         emp_vars=vars.emp,
                         num_class = NULL) {
  
  library(dplyr)
  library(xgboost)
  
  dat1 <- dat %>%
    select(all_of(out_vars), all_of(clust_vars), all_of(exposure), all_of(baselinevars), all_of(emp_vars))
  
  cont_vars <- continuous_vars[continuous_vars %in% baselinevars]
  
  
  
  model_data <- dat1 %>%
    filter(if_all(c(all_of(baselinevars), all_of(emp_vars), exposure), ~ !is.na(.)))
  
  valid_baselinevars <- c(baselinevars, emp_vars)[sapply(model_data[c(baselinevars, emp_vars)], function(x) length(unique(x)) > 1)]
  
  # Exposure as factor → numeric labels for xgboost
  model_data[[exposure]] <- as.factor(model_data[[exposure]])
  levels_exp <- levels(model_data[[exposure]])
  n_exp <- length(levels_exp)
  label_y <- as.integer(model_data[[exposure]]) - 1   # required: 0,1,...,K-1
  
  # Build model matrix
  X <- model.matrix(~ . -1, data = model_data[, valid_baselinevars, drop = FALSE])
  dtrain <- xgb.DMatrix(data = X, label = label_y)
  
  # Choose objective
  obj <- if (n_exp > 2) "multi:softprob" else "binary:logistic"
  if (is.null(num_class)) num_class <- if (n_exp > 2) n_exp else 1
  
  # Fit xgboost
  ps_fit <- xgboost(data = dtrain,
                    objective = obj,
                    num_class = num_class,
                    nthread = parallel::detectCores(),
                    verbose = 0, 
                    nrounds = 100)
  
  # Predict probabilities
  ps_probs <- predict(ps_fit, newdata = dtrain)
  
  if (n_exp > 2) {
    ps_probs <- matrix(ps_probs, ncol = n_exp, byrow = TRUE)
    colnames(ps_probs) <- levels_exp
    obs_exposure <- model_data[[exposure]]
    
    # Extract each subject’s probability of observed treatment
    model_data$ps <- mapply(function(i, trt) ps_probs[i, trt],
                            i = seq_along(obs_exposure),
                            trt = obs_exposure)
    
    # Stabilized weights
    marginal_probs <- prop.table(table(obs_exposure))
    model_data$stab_numerator <- marginal_probs[as.character(obs_exposure)]
    model_data$wt_stab <- model_data$stab_numerator / model_data$ps
    
  } else {
    # Binary exposure
    ps_prob <- ps_probs
    y <- as.numeric(as.character(model_data[[exposure]]))
    
    model_data$ps <- ifelse(y == 1, ps_prob, 1 - ps_prob)
    stab_num <- ifelse(y == 1, mean(y == 1), 1 - mean(y == 1))
    model_data$wt_stab <- as.numeric(stab_num / model_data$ps)
  }
  
  # Truncate weights
  p_low <- quantile(model_data$wt_stab, trunc_low)
  p_high <- quantile(model_data$wt_stab, trunc_high)
  model_data$wt_trunc <- as.numeric(pmin(pmax(model_data$wt_stab, p_low), p_high))
  
  return(list(model = ps_fit, data_weights = model_data))
}






library(data.table)
library(lubridate)

compute_prescriber_pct_dt <- function(dat, oat, trt_name, cond = NULL) {
  
  # Convert to data.table
  setDT(dat)
  setDT(oat)
  
  # Apply condition if provided
  if (!is.null(cond)) {
    dat <- dat[eval(parse(text=cond))]
    oat <- oat[eval(parse(text=cond))]
  }
  
  
  # Step 1: Pre-aggregate oat by prescriber & episode
  oat_counts <- oat[, .(
    grp1_eps = sum(get(trt_name) == "no reversion", na.rm=TRUE),
    grp2_eps = sum(get(trt_name) %in% c("1-<50%", ">=50%/<=40mg"), na.rm=TRUE),
    grp3_eps = sum(get(trt_name) == "1-<50%", na.rm=TRUE),
    grp4_eps = sum(get(trt_name) == ">=50%/<=40mg", na.rm=TRUE),
    grp5_eps = sum(get(trt_name) == "increase", na.rm=TRUE),
    n_eps    = .N
  ), by = .(prac_studyid_rev, time0)]
  
  # Step 2: Set keys for interval join
  setkey(dat, prac_studyid_rev, time0_12m_prior, time0)
  setkey(oat_counts, prac_studyid_rev, time0)
  
  # Step 3: Interval join: oat episodes in past 12 months
  joined <- foverlaps(
    oat_counts[, .(prac_studyid_rev, oat_start = time0, oat_end = time0, grp1_eps, grp2_eps, grp3_eps, grp4_eps, grp5_eps, n_eps)],
    dat[, .(prac_studyid_rev, time0_12m_prior, time0, studyid_rev)],
    by.x = c("prac_studyid_rev", "oat_start", "oat_end"),
    by.y = c("prac_studyid_rev", "time0_12m_prior", "time0"),
    type = "within",  # select oat rows within dat 12-month interval
    nomatch = 0L
  )
  
  # Step 4: Sum counts per patient episode
  result <- joined[, .(
    grp1_eps = sum(grp1_eps),
    grp2_eps = sum(grp2_eps),
    grp3_eps = sum(grp3_eps),
    grp4_eps = sum(grp4_eps),
    grp5_eps = sum(grp5_eps),
    n_eps    = sum(n_eps)
  ), by = .(studyid_rev, time0, prac_studyid_rev)]
  
  # Step 5: Calculate percentages
  result[, `:=` (
    grp1_eps_pct = fifelse(n_eps > 0, grp1_eps / n_eps * 100, NA_real_),
    grp2_eps_pct = fifelse(n_eps > 0, grp2_eps / n_eps * 100, NA_real_),
    grp3_eps_pct = fifelse(n_eps > 0, grp3_eps / n_eps * 100, NA_real_),
    grp4_eps_pct = fifelse(n_eps > 0, grp4_eps / n_eps * 100, NA_real_),
    grp5_eps_pct = fifelse(n_eps > 0, grp5_eps / n_eps * 100, NA_real_)
  )]
  
  # Step 6: Merge back with original dat to keep all rows
  final <- merge(dat, result, by = c("studyid_rev", "time0", "prac_studyid_rev"), all.x = TRUE)
  
  return(final[])
}


compute_facility_pct_dt <- function(dat, oat, trt_name, cond = NULL) {
  
  # Convert to data.table
  setDT(dat)
  setDT(oat)
  
  # Apply condition if provided
  if (!is.null(cond)) {
    dat <- dat[eval(parse(text = cond))]
    oat <- oat[eval(parse(text = cond))]
  }
  
  # Step 1: Pre-aggregate oat by facility & episode
  oat_counts <- oat[, .(
    grp1_eps = sum(get(trt_name) == "no reversion", na.rm=TRUE),
    grp2_eps = sum(get(trt_name) %in% c("1-<50%", ">=50%/<=40mg"), na.rm=TRUE),
    grp3_eps = sum(get(trt_name) == "1-<50%", na.rm=TRUE),
    grp4_eps = sum(get(trt_name) == ">=50%/<=40mg", na.rm=TRUE),
    grp5_eps = sum(get(trt_name) == "increase", na.rm=TRUE),
    n_eps    = .N
  ), by = .(com, time0)]
  
  # Step 2: Set keys for interval join
  setkey(dat, com, time0_12m_prior, time0)
  setkey(oat_counts, com, time0)
  
  # Step 3: Interval join: oat episodes in past 12 months
  joined <- foverlaps(
    oat_counts[, .(com, oat_start = time0, oat_end = time0, grp1_eps, grp2_eps, grp3_eps, grp4_eps, grp5_eps, n_eps)],
    dat[, .(com, time0_12m_prior, time0, time0_id = time0)],
    by.x = c("com", "oat_start", "oat_end"),
    by.y = c("com", "time0_12m_prior", "time0"),
    type = "within",
    nomatch = 0L
  )
  
  # Step 4: Sum counts per facility × episode
  result <- joined[, .(
    grp1_eps = sum(grp1_eps),
    grp2_eps = sum(grp2_eps),
    grp3_eps = sum(grp3_eps),
    grp4_eps = sum(grp4_eps),
    grp5_eps = sum(grp5_eps),
    n_eps    = sum(n_eps)
  ), by = .(com, time0_id)]
  
  # Step 5: Calculate percentages
  result[, `:=` (
    grp1_eps_pct = fifelse(n_eps > 0, grp1_eps / n_eps * 100, NA_real_),
    grp2_eps_pct = fifelse(n_eps > 0, grp2_eps / n_eps * 100, NA_real_),
    grp3_eps_pct = fifelse(n_eps > 0, grp3_eps / n_eps * 100, NA_real_),
    grp4_eps_pct = fifelse(n_eps > 0, grp4_eps / n_eps * 100, NA_real_),
    grp5_eps_pct = fifelse(n_eps > 0, grp5_eps / n_eps * 100, NA_real_)
  )]
  
  # Step 6: Merge back with original dat to keep all rows
  final <- merge(dat, result, by.x = c("com", "time0"), by.y = c("com", "time0_id"), all.x = TRUE)
  
  return(final[])
}




compute_prescriber_pct_dt_b <- function(dat, oat, trt_name, cond = NULL) {
  
  # Convert to data.table
  setDT(dat)
  setDT(oat)
  
  # Apply condition if provided
  if (!is.null(cond)) {
    dat <- dat[eval(parse(text=cond))]
    oat <- oat[eval(parse(text=cond))]
  }
  
  
  # Step 1: Pre-aggregate oat by prescriber & episode
  oat_counts <- oat[, .(
    grp1_eps = sum(get(trt_name) == "no reversion", na.rm=TRUE),
    grp2_eps = sum(get(trt_name) %in% c("1-<50%", ">=50%/<=4mg"), na.rm=TRUE),
    grp3_eps = sum(get(trt_name) == "1-<50%", na.rm=TRUE),
    grp4_eps = sum(get(trt_name) == ">=50%/<=4mg", na.rm=TRUE),
    grp5_eps = sum(get(trt_name) == "increase", na.rm=TRUE),
    n_eps    = .N
  ), by = .(prac_studyid_rev, time0)]
  
  # Step 2: Set keys for interval join
  setkey(dat, prac_studyid_rev, time0_12m_prior, time0)
  setkey(oat_counts, prac_studyid_rev, time0)
  
  # Step 3: Interval join: oat episodes in past 12 months
  joined <- foverlaps(
    oat_counts[, .(prac_studyid_rev, oat_start = time0, oat_end = time0, grp1_eps, grp2_eps, grp3_eps, grp4_eps, grp5_eps, n_eps)],
    dat[, .(prac_studyid_rev, time0_12m_prior, time0, studyid_rev)],
    by.x = c("prac_studyid_rev", "oat_start", "oat_end"),
    by.y = c("prac_studyid_rev", "time0_12m_prior", "time0"),
    type = "within",  # select oat rows within dat 12-month interval
    nomatch = 0L
  )
  
  # Step 4: Sum counts per patient episode
  result <- joined[, .(
    grp1_eps = sum(grp1_eps),
    grp2_eps = sum(grp2_eps),
    grp3_eps = sum(grp3_eps),
    grp4_eps = sum(grp4_eps),
    grp5_eps = sum(grp5_eps),
    n_eps    = sum(n_eps)
  ), by = .(studyid_rev, time0, prac_studyid_rev)]
  
  # Step 5: Calculate percentages
  result[, `:=` (
    grp1_eps_pct = fifelse(n_eps > 0, grp1_eps / n_eps * 100, NA_real_),
    grp2_eps_pct = fifelse(n_eps > 0, grp2_eps / n_eps * 100, NA_real_),
    grp3_eps_pct = fifelse(n_eps > 0, grp3_eps / n_eps * 100, NA_real_),
    grp4_eps_pct = fifelse(n_eps > 0, grp4_eps / n_eps * 100, NA_real_),
    grp5_eps_pct = fifelse(n_eps > 0, grp5_eps / n_eps * 100, NA_real_)
  )]
  
  # Step 6: Merge back with original dat to keep all rows
  final <- merge(dat, result, by = c("studyid_rev", "time0", "prac_studyid_rev"), all.x = TRUE)
  
  return(final[])
}


compute_facility_pct_dt_b <- function(dat, oat, trt_name, cond = NULL) {
  
  # Convert to data.table
  setDT(dat)
  setDT(oat)
  
  # Apply condition if provided
  if (!is.null(cond)) {
    dat <- dat[eval(parse(text = cond))]
    oat <- oat[eval(parse(text = cond))]
  }
  
  # Step 1: Pre-aggregate oat by facility & episode
  oat_counts <- oat[, .(
    grp1_eps = sum(get(trt_name) == "no reversion", na.rm=TRUE),
    grp2_eps = sum(get(trt_name) %in% c("1-<50%", ">=50%/<=4mg"), na.rm=TRUE),
    grp3_eps = sum(get(trt_name) == "1-<50%", na.rm=TRUE),
    grp4_eps = sum(get(trt_name) == ">=50%/<=4mg", na.rm=TRUE),
    grp5_eps = sum(get(trt_name) == "increase", na.rm=TRUE),
    n_eps    = .N
  ), by = .(com, time0)]
  
  # Step 2: Set keys for interval join
  setkey(dat, com, time0_12m_prior, time0)
  setkey(oat_counts, com, time0)
  
  # Step 3: Interval join: oat episodes in past 12 months
  joined <- foverlaps(
    oat_counts[, .(com, oat_start = time0, oat_end = time0, grp1_eps, grp2_eps, grp3_eps, grp4_eps, grp5_eps, n_eps)],
    dat[, .(com, time0_12m_prior, time0, time0_id = time0)],
    by.x = c("com", "oat_start", "oat_end"),
    by.y = c("com", "time0_12m_prior", "time0"),
    type = "within",
    nomatch = 0L
  )
  
  # Step 4: Sum counts per facility × episode
  result <- joined[, .(
    grp1_eps = sum(grp1_eps),
    grp2_eps = sum(grp2_eps),
    grp3_eps = sum(grp3_eps),
    grp4_eps = sum(grp4_eps),
    grp5_eps = sum(grp5_eps),
    n_eps    = sum(n_eps)
  ), by = .(com, time0_id)]
  
  # Step 5: Calculate percentages
  result[, `:=` (
    grp1_eps_pct = fifelse(n_eps > 0, grp1_eps / n_eps * 100, NA_real_),
    grp2_eps_pct = fifelse(n_eps > 0, grp2_eps / n_eps * 100, NA_real_),
    grp3_eps_pct = fifelse(n_eps > 0, grp3_eps / n_eps * 100, NA_real_),
    grp4_eps_pct = fifelse(n_eps > 0, grp4_eps / n_eps * 100, NA_real_),
    grp5_eps_pct = fifelse(n_eps > 0, grp5_eps / n_eps * 100, NA_real_)
  )]
  
  # Step 6: Merge back with original dat to keep all rows
  final <- merge(dat, result, by.x = c("com", "time0"), by.y = c("com", "time0_id"), all.x = TRUE)
  
  return(final[])
}



library(splines)
library(lme4)
library(lmerTest)

fit_iv_F_lmm <- function(data, trt_var, iv_var, baselinevars,
                         prescriber_id = "prac_studyid_rev",
                         spline_df = 3, center_iv = TRUE,
                         verbose = TRUE) {
  if (center_iv) {
    data[[iv_var]] <- scale(data[[iv_var]], center = TRUE, scale = FALSE)
  }
  
  spline_basis <- as.data.frame(ns(data[[iv_var]], df = spline_df))
  colnames(spline_basis) <- paste0("iv_ns", seq_len(ncol(spline_basis)))
  dat2 <- cbind(data, spline_basis)
  
  fixed_terms <- c(colnames(spline_basis), baselinevars)
  lmer_formula <- as.formula(
    paste0(trt_var, " ~ ", paste(fixed_terms, collapse = " + "),
           " + (1 | ", prescriber_id, ")")
  )
  
  # Try fitting mixed model
  lmm_fit <- NULL
  tryCatch({
    if (verbose) message("Fitting LMM (Gaussian) with lmerTest...")
    lmm_fit <- lmer(lmer_formula, data = dat2, REML = FALSE)
  }, error = function(e) {
    message("LMM failed, falling back to OLS only. Reason: ", e$message)
  })
  
  # Partial F-test via OLS
  lm_reduced <- lm(as.formula(
    paste0(trt_var, " ~ ", paste(baselinevars, collapse = " + "))
  ), data = dat2)
  
  lm_full <- lm(as.formula(
    paste0(trt_var, " ~ ", paste(c(baselinevars, colnames(spline_basis)), collapse = " + "))
  ), data = dat2)
  
  anova_res <- anova(lm_reduced, lm_full, test = "F")
  F_value <- anova_res$`F`[2]
  
  # Return results
  list(
    anova = anova_res,
    F_value = F_value
  )
}



library(data.table)
library(data.table)

s_compute_prescriber_pct_dt <- function(dat, oat, cond, trt_name, file) {
  # Ensure data.table
  setDT(dat)
  setDT(oat)
  
  # Apply condition filter
  dat <- dat[eval(parse(text = cond))]
  oat <- oat[eval(parse(text = cond))]
  
  # Keep only needed columns
  dat <- dat[, .(studyid_rev, prac_studyid_rev, episode, time0, time0_12m_prior, trt = get(trt_name))]
  oat <- oat[, .(studyid_rev, prac_studyid_rev, episode, time0, trt = get(trt_name))]
  
  # Get all prescribers
  pracs <- unique(dat$prac_studyid_rev)
  
  # Prepare empty list to collect chunks
  results_list <- vector("list", length(pracs))
  
  message("Processing ", length(pracs), " prescribers in chunks...")
  
  for (i in seq_along(pracs)) {
    p <- pracs[i]
    
    dat_chunk <- dat[prac_studyid_rev == p]
    oat_chunk <- oat[prac_studyid_rev == p]
    
    if (nrow(dat_chunk) == 0 | nrow(oat_chunk) == 0) next
    
    # Join with condition: time0_12m_prior <= oat.time0 < dat.time0
    out_chunk <- oat_chunk[dat_chunk, 
                           on = .(prac_studyid_rev, time0 < time0, time0 >= time0_12m_prior), 
                           allow.cartesian = TRUE,
                           nomatch = 0]
    
    if (nrow(out_chunk) == 0) next
    
    # Compute counts per dat row
    res_chunk <- out_chunk[, .(
      n_eps = uniqueN(paste(studyid_rev, time0)),
      grp1_eps = uniqueN(paste(studyid_rev, time0)[trt == "no reversion"]),
      grp2_eps = uniqueN(paste(studyid_rev, time0)[trt %in% c("1-<50%", ">=50%/<=40mg")]),
      grp3_eps = uniqueN(paste(studyid_rev, time0)[trt == "1-<50%"]),
      grp4_eps = uniqueN(paste(studyid_rev, time0)[trt == ">=50%/<=40mg"]),
      grp5_eps = uniqueN(paste(studyid_rev, time0)[trt == "increase"])
    ), by = .(prac_studyid_rev, studyid_rev, episode, time0)]
    
    # Add percentages
    res_chunk[, `:=`(
      grp1_eps_pct = 100 * grp1_eps / n_eps,
      grp2_eps_pct = 100 * grp2_eps / n_eps,
      grp3_eps_pct = 100 * grp3_eps / n_eps,
      grp4_eps_pct = 100 * grp4_eps / n_eps,
      grp5_eps_pct = 100 * grp5_eps / n_eps
    )]
    
    results_list[[i]] <- res_chunk
    
    if (i %% 50 == 0) message("Processed ", i, "/", length(pracs), " prescribers...")
  }
  
  # Combine all chunks
  results <- rbindlist(results_list, fill = TRUE)
  
  # Save to RDS
  saveRDS(results, file)
  message("Saved results to: ", file)
}



fit_iv_logistic <- function(data, trt_var, iv_var, spline_df = 3, n_points = 200) {
  library(splines)
  
  # Build formula
  form <- as.formula(paste0(trt_var, " ~ ns(", iv_var, ", df = ", spline_df, ")"))
  
  # Fit logistic regression
  model <- glm(form, data = data, family = binomial)
  
  # Prediction grid
  preddata <- data.frame(
    iv_var = seq(min(data[[iv_var]], na.rm = TRUE),
                 max(data[[iv_var]], na.rm = TRUE),
                 length.out = n_points)
  )
  names(preddata) <- iv_var
  
  preds <- predict(model, newdata = preddata, type = "link", se.fit = TRUE)
  
  # 95% CI
  critval <- 1.96
  lwr <- preds$fit - critval * preds$se.fit
  upr <- preds$fit + critval * preds$se.fit
  
  # Back-transform
  preddata$fit <- model$family$linkinv(preds$fit)
  preddata$lwr <- model$family$linkinv(lwr)
  preddata$upr <- model$family$linkinv(upr)
  
  list(
    model = model,
    predictions = preddata
  )
}



plot_iv_logistic <- function(fit, trt_var="trt", iv_var="iv_var", 
                             binwidth = 10, 
                             xlab_text = NULL, 
                             ylab_text = "Proportion of sample",
                             sec_y_text = "Probability of assignment to treatment",
                             y_points1, 
                             y_points2,
                             y_min=0, y_max=1,
                             legend_option=TRUE,
                             plot_title=NULL) {
  
  library(ggplot2)
  
  data = fit$model$data
  preds <- fit$predictions
  
  # Histogram (scaled to proportion)
  p <- ggplot(data, aes_string(x = iv_var)) +
    geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))),
                   binwidth = binwidth, boundary = 0,
                   colour = "black", fill = "white") +
    
    # Add CI ribbon + fit line
    geom_ribbon(data = preds,
                aes_string(x = iv_var, ymin = "lwr", ymax = "upr"),
                fill = "grey70", alpha = 0.5, inherit.aes = FALSE) +
    geom_line(data = preds,
              aes_string(x = iv_var, y = "fit"),
              color = "red", size = 0.7, inherit.aes = FALSE) +
    
    # Axes
    scale_y_continuous(
      name = ylab_text,
      expand = c(0, 0),
      limits = c(y_min, y_max),
      breaks = y_points1,
      sec.axis = sec_axis(~., name = sec_y_text, breaks = y_points2)
    ) +
    scale_x_continuous(name = ifelse(is.null(xlab_text), iv_var, xlab_text),
                       breaks = seq(0, 100, 10)) +
    theme_bw()+
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.text = element_text(color = "black", size = 14),
      axis.title.y = element_text(size = 14),
      legend.position = if(legend_option) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 13)
    )
  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title)
  }
  
  return(p)
}

fit1st_iv_logistic <- function(data, covs, trt_var, iv_var, n_points = 200) {
  # Collapse covariates into formula string
  covs_str <- if (!is.null(covs) && length(covs) > 0) {
    paste(c(covs, iv_var), collapse = " + ")
  } else {
    iv_var
  }
  
  # Build formula
  form <- as.formula(paste(trt_var, "~", covs_str))
  
  # Fit logistic regression
  model <- glm(form, data = data, family = binomial)
  
  # Residuals and predictions merged with original data
  residuals_df <- data.frame(
    data,
    xbeta    = predict(model, type = "link"),     # linear predictor
    predprob = predict(model, type = "response"), # fitted probability
    resid    = residuals(model, type = "response"), # observed - fitted
    pearson  = residuals(model, type = "pearson"),
    deviance = residuals(model, type = "deviance"),
    student  = rstudent(model)                    # studentized
  )
  
  list(
    model = model,
    residuals = residuals_df
  )
}




inequality_test_binaryIV <- function(data=dat1, iv_cont = "grp2_eps_pct",
                                     exposure = "trt",
                                     outcome = "sens_discont"){
  df <- data[, c(iv_cont, exposure, outcome)]
  
  med_val <- median(df[[iv_cont]], na.rm = TRUE)
  df$IV_bin <- ifelse(df[[iv_cont]] <= med_val, 0, 1)
  
  df$X <- df[[exposure]]
  df$Y <- df[[outcome]]
  df$IV <- df$IV_bin
  
  n_uniq_Y <- length(unique(df$Y))
  com <- lapply(1:n_uniq_Y, function(k) unique(df$IV))
  com[[length(com)+1]] <- unique(df$X)
  com <- expand.grid(com, stringsAsFactors = FALSE)
  names(com) <- c(paste0("IV", 1:n_uniq_Y),"X")
  
  prp <- lapply(1:n_uniq_Y, function(i) {
    sapply(1:nrow(com), function(j) {
      num <- sum(df$IV==com[j, i] & df$X==com$X[j] & df$Y==unique(df$Y)[i])
      den <- sum(df$IV==com[j, i])
      ifelse(den == 0, NA, num/den)
    })
  }) %>% do.call(cbind, .) %>% as.data.frame
  names(prp) <- paste0("Y", unique(df$Y), "_IV", 1:n_uniq_Y)
  prp$sum_prop <- rowSums(prp, na.rm = TRUE)
  
  com <- cbind(com, prp)
  
  agg <- aggregate(sum_prop ~ IV1, data=com, FUN=max, na.rm=TRUE)
  max_val <- max(agg$sum_prop, na.rm = TRUE)
  
  message("Median of ", iv_cont, " = ", round(med_val,2))
  message("Max inequality value = ", round(max_val,4))
  if(max_val <= 1){
    message("Result: Inequalities HOLD (no falsification)")
  } else {
    message("Result: Inequalities VIOLATED (IV invalid under this test)")
  }
  
  return(list(
    median_value = med_val,
    max_value = max_val,
    results_table = com,
    agg = agg
  ))
}

library(dplyr)
library(tidyr)
library(sandwich)
library(lmtest)

compute_bias_components_robust <- function(data, iv_bin = "IV_bin", exposure = "trt", covars,
                                           conf = 0.95) {
  alpha <- 1 - conf
  z_val <- qnorm(1 - alpha/2)
  
  # Function to compute PD and robust SE
  compute_once <- function(dat) {
    results <- lapply(covars, function(v) {
      # 1. PD by treatment
      formula_trt <- as.formula(paste0(v, " ~ ", exposure))
      model_trt <- lm(formula_trt, data = dat)
      coefs <- coef(model_trt)
      se <- sqrt(diag(sandwich::vcovHC(model_trt, type = "HC1")))
      mean_trt <- coefs[exposure]
      lower_trt <- mean_trt - z_val * se[exposure]
      upper_trt <- mean_trt + z_val * se[exposure]
      df_trt <- tibble(
        variable = v, biascomponent = "By treatment",
        Mean = mean_trt, lower = lower_trt, upper = upper_trt
      )
      
      # 2. PD by IV (unscaled)
      formula_iv <- as.formula(paste0(v, " ~ ", iv_bin))
      model_iv <- lm(formula_iv, data = dat)
      coefs_iv <- coef(model_iv)
      se_iv <- sqrt(diag(sandwich::vcovHC(model_iv, type = "HC1")))
      mean_iv <- coefs_iv[iv_bin]
      lower_iv <- mean_iv - z_val * se_iv[iv_bin]
      upper_iv <- mean_iv + z_val * se_iv[iv_bin]
      df_iv <- tibble(
        variable = v, biascomponent = "By IV (unscaled)",
        Mean = mean_iv, lower = lower_iv, upper = upper_iv
      )
      
      # 3. PD by IV (scaled)
      # Scale factor = mean(trt | iv=1) - mean(trt | iv=0)
      scale_factor <- dat %>%
        group_by(!!sym(iv_bin)) %>%
        summarise(p = mean(get(exposure), na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = !!sym(iv_bin), values_from = p, names_prefix = "iv_")
      scale_val <- scale_factor$iv_1 - scale_factor$iv_0
      
      mean_iv_scaled <- mean_iv / scale_val
      se_iv_scaled <- se_iv[iv_bin] / scale_val  # delta method
      lower_iv_scaled <- mean_iv_scaled - z_val * se_iv_scaled
      upper_iv_scaled <- mean_iv_scaled + z_val * se_iv_scaled
      df_iv_scaled <- tibble(
        variable = v, biascomponent = "By IV (scaled)",
        Mean = mean_iv_scaled, lower = lower_iv_scaled, upper = upper_iv_scaled
      )
      
      bind_rows(df_trt, df_iv, df_iv_scaled)
    }) %>% bind_rows()
    results
  }
  
  compute_once(data) %>% arrange(biascomponent, variable)
}


compute_bias_components_boot <- function(data, iv_bin = "IV_bin", exposure = "trt", covars,
                                         n_boot = 200, conf = 0.95, seed = 1234) {
  set.seed(seed)
  
  # Function to compute PDs in a dataset (used for both full data and resamples)
  compute_once <- function(dat) {
    # PD by treatment
    trt_pd <- lapply(covars, function(v) {
      tab <- dat %>%
        group_by(!!sym(exposure)) %>%
        summarise(p = mean(get(v), na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = !!sym(exposure), values_from = p, names_prefix = "trt_")
      pd <- tab$trt_1 - tab$trt_0
      tibble(variable = v, Mean = pd, biascomponent = "By treatment")
    }) %>% bind_rows()
    
    # PD by IV (unscaled)
    iv_pd <- lapply(covars, function(v) {
      tab <- dat %>%
        group_by(!!sym(iv_bin)) %>%
        summarise(p = mean(get(v), na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = !!sym(iv_bin), values_from = p, names_prefix = "iv_")
      pd <- tab$iv_1 - tab$iv_0
      tibble(variable = v, Mean = pd, biascomponent = "By IV (unscaled)")
    }) %>% bind_rows()
    
    # Scale factor (difference in treatment by IV)
    bnx_prev <- dat %>%
      group_by(!!sym(iv_bin)) %>%
      summarise(p = mean(get(exposure), na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = !!sym(iv_bin), values_from = p, names_prefix = "iv_")
    scale_factor <- bnx_prev$iv_1 - bnx_prev$iv_0
    
    # PD by IV (scaled)
    iv_scaled <- iv_pd %>%
      mutate(Mean = Mean / scale_factor,
             biascomponent = "By IV (scaled)")
    
    bind_rows(trt_pd, iv_pd, iv_scaled)
  }
  
  # Compute original estimates
  orig <- compute_once(data) %>% arrange(biascomponent, variable)
  
  # Bootstrap
  boot_list <- replicate(n_boot, {
    samp <- data[sample(nrow(data), replace = TRUE), ]
    compute_once(samp)
  }, simplify = FALSE)
  
  # Combine bootstrap samples into a single tibble
  boot_df <- bind_rows(boot_list, .id = "boot_id")
  
  # Summarize bootstrap distribution
  res <- boot_df %>%
    group_by(variable, biascomponent) %>%
    summarise(
      lower = quantile(Mean, 0.025, na.rm = TRUE),
      upper = quantile(Mean, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>% arrange(biascomponent, variable)
  
  # Merge original and bootstrap summary
  merged <- orig %>%
    left_join(res, by = c("variable", "biascomponent"))
  
  return(merged)
}



extract_max_values <- function(res_obj, trial_label = "Trial", include_median = FALSE) {
  designs <- c(
    incid = "Incident users",
    prev  = "Prevalent new users",
    rzero = "All prevalent users"
  )
  outcomes <- c(
    od  = "Overdose-related acute care visit",
    dis = "OAT discontinuation"
  )
  
  # Prepare an empty result if nothing found
  make_empty_df <- function() {
    if (include_median) {
      data.frame(
        Group = character(),
        Subgroup = character(),
        Outcome = character(),
        Design = character(),
        MaxValue = numeric(),
        MedianValue = numeric(),
        Trial = character(),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Group = character(),
        Subgroup = character(),
        Outcome = character(),
        Design = character(),
        MaxValue = numeric(),
        Trial = character(),
        stringsAsFactors = FALSE
      )
    }
  }
  
  out_list <- list()
  idx <- 1L
  
  if (is.null(res_obj) || length(res_obj) == 0) return(make_empty_df())
  
  for (group_name in names(res_obj)) {
    group_data <- res_obj[[group_name]]
    if (is.null(group_data) || length(group_data) == 0) next
    
    for (subgroup_name in names(group_data)) {
      subgroup_data <- group_data[[subgroup_name]]
      if (is.null(subgroup_data) || length(subgroup_data) == 0) next
      
      for (design_name in names(designs)) {
        result <- subgroup_data[[design_name]]
        if (is.null(result)) next
        
        mv <- result$max_value
        if (is.null(mv) || length(mv) == 0 || is.na(mv) || is.nan(mv)) next
        
        outcome_label <- if (subgroup_name %in% names(outcomes)) outcomes[[subgroup_name]] else subgroup_name
        
        if (include_median) {
          medv <- if (!is.null(result$median_value) && length(result$median_value) > 0) result$median_value else NA_real_
          out_list[[idx]] <- data.frame(
            Group = group_name,
            Subgroup = subgroup_name,
            Outcome = outcome_label,
            Design = designs[[design_name]],
            MaxValue = as.numeric(mv),
            MedianValue = as.numeric(medv),
            Trial = trial_label,
            stringsAsFactors = FALSE
          )
        } else {
          out_list[[idx]] <- data.frame(
            Group = group_name,
            Subgroup = subgroup_name,
            Outcome = outcome_label,
            Design = designs[[design_name]],
            MaxValue = as.numeric(mv),
            Trial = trial_label,
            stringsAsFactors = FALSE
          )
        }
        idx <- idx + 1L
      } # design
    } # subgroup
  } # group
  
  if (length(out_list) == 0) return(make_empty_df())
  do.call(rbind, out_list)
}



extract_bias_results <- function(res_obj) {
  designs <- c(
    incid = "Incident users",
    prev  = "Prevalent new users",
    rzero = "All prevalent users"
  )
  
  all_rows <- lapply(names(res_obj), function(group_name) {
    group_data <- res_obj[[group_name]]
    
    lapply(names(group_data), function(design_name) {
      df <- group_data[[design_name]]
      if (is.null(df)) return(NULL)
      
      df$Group  <- group_name
      df$Design <- designs[[design_name]]
      df
    }) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()
  
  return(all_rows)
}