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
  
  # Truncate weights
  ps_probs_trunc <- as.data.frame(
    lapply(ps_probs, function(x) {
      p_low  <- quantile(x, trunc_low)
      p_high <- quantile(x, trunc_high)
      pmin(pmax(x, p_low), p_high)
    })
  )
  
  p_low <- quantile(model_data$ps, trunc_low)
  p_high <- quantile(model_data$ps, trunc_high)
  model_data$ps_trunc <- pmin(pmax(model_data$ps, p_low), p_high)
  
  ## ow
  tilting_fn <- apply(ps_probs, 1, function(x) sum(1/x)) 
  
  tilting_fn_trunc <- apply(ps_probs_trunc, 1, function(x) sum(1/x)) 
  
  model_data$ow_wt <- as.numeric(1/ (tilting_fn*model_data$ps))
  
  model_data$ow_wt_trunc <- as.numeric(1/ (tilting_fn_trunc*model_data$ps_trunc))
  
  return(list(model = ps_fit, data_weights = model_data))
}


ps_model_ow_cl <- function(dat, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99, 
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
  dat <- dat2_rprev
  model_data <- dat %>%
    dplyr::select(all_of(c(out_vars, clust_vars, exposure, baselinevars))) %>%
    mutate(across(all_of(setdiff(baselinevars, continuous_vars)), as.factor)) %>%
    filter(if_all(all_of(c(exposure, baselinevars)), ~ !is.na(.)))
  
  valid_baselinevars <- baselinevars[sapply(model_data[baselinevars], function(x) length(unique(x)) > 1)]
  model_data$studyid_rev <- factor(model_data$studyid_rev)
  formula_ps <- as.formula(paste(exposure, "~", paste(c(valid_baselinevars, "studyid_rev"), collapse = " + ")))
  
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
  
  # Truncate weights
  ps_probs_trunc <- as.data.frame(
    lapply(ps_probs, function(x) {
      p_low  <- quantile(x, trunc_low)
      p_high <- quantile(x, trunc_high)
      pmin(pmax(x, p_low), p_high)
    })
  )
  
  p_low <- quantile(model_data$ps, trunc_low)
  p_high <- quantile(model_data$ps, trunc_high)
  model_data$ps_trunc <- pmin(pmax(model_data$ps, p_low), p_high)
  
  ## ow
  tilting_fn <- apply(ps_probs, 1, function(x) sum(1/x)) 
  
  tilting_fn_trunc <- apply(ps_probs_trunc, 1, function(x) sum(1/x)) 
  
  model_data$ow_wt <- as.numeric(1/ (tilting_fn*model_data$ps))
  
  model_data$ow_wt_trunc <- as.numeric(1/ (tilting_fn_trunc*model_data$ps_trunc))
  
  return(list(model = ps_fit, data_weights = model_data))
}




create_table_ow = function(data, exposure = "trt", trial=1, 
                           incid=TRUE, baselinevars,
                           continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"),
                           wt){
  dat1 = data
  cont_vars <- continuous_vars[continuous_vars %in% baselinevars]
  
  # Table 1
  tr1_table1_before = CreateTableOne(vars = baselinevars, factorVars = setdiff(baselinevars, cont_vars), strata = "trt", data = dat1, addOverall=T)
  tr1_table1_before <- print(tr1_table1_before, smd = TRUE, showAllLevels = F, exact = "stage", 
                             quote = F, noSpaces = T, printToggle = FALSE, nonnormal = cont_vars, test = FALSE)
  rename_vars = c("Total", 
                  "Sex: Male", 
                  "Age[Median (Q1-Q3)]",  
                  "Age category: 18-29", 
                  "30-39", 
                  "40-49", 
                  "50+",
                  "Region: rural", 
                  "Year of reinitiation date [Median (Q1-Q3)", 
                  "Opioid dispensations other than OAT prior to the first OUD indication",
                  "Attachment to the prescriber who reinitiated OAT", 
                  "Receipt of income assistance in the past year", 
                  "Most recent OAT Retention length (in days) [Median (Q1-Q3)",
                  "Unstable housing in past 5 years", 
                  "Incarcerated in the past year", 
                  "Drug-related acute care visits in the past month",
                  "Charlson Comorbidity Index (>0)",
                  "Asthma or COPD in the past year",
                  "Dispensations of psychiatric or sedative medication in the past month",
                  "Dispensations opioid other than OAT in the past month", 
                  "Other substance use disorder",
                  "Alcohol use disorder",
                  "Serious mental disorder",
                  "Hepatitis C",
                  "Tobacco use disorder (past year)",
                  "Urine drug test in the past moth",
                  "Extent of the last dose dispensed [Median (Q1-Q3)",
                  "Large practice size (>300 clients)",
                  "Virtual care in past month",
                  "pss")
  
  if (incid) {
    design <- svydesign(id = ~1, weights = ~dat1[[wt]], data = dat1)
  } else {
    design <- svydesign(id = ~studyid_rev, weights = ~dat1[[wt]], data = dat1)
  }
  
  tr1_table1_after = svyCreateTableOne(vars = baselinevars, factorVars = setdiff(baselinevars, cont_vars), 
                                       strata = "trt", addOverall = F, data = design, smd = T, test = T)
  tr1_table1_after <- print(tr1_table1_after, smd = TRUE, showAllLevels = F, exact = "stage", 
                            quote = F, noSpaces = T, printToggle = FALSE, nonnormal = cont_vars, test = FALSE)
  
  smd_column <- tr1_table1_after[, "SMD"]
  table <- cbind(tr1_table1_before, SMD_after = smd_column)
  
  return(table)
  
}

fit_PLR_HR <- function(Dat_expanded, trt_var = "trt", 
                       covs = NULL, 
                       cluster_var = "studyid_rev",
                       wt){
  covs1 <- paste0(c(trt_var, "ns(interval, df = 3)", covs), collapse = "+")
  plr_formula <- as.formula(paste("event ~", covs1))
  
  model <- glm(plr_formula, data = Dat_expanded,
               family = binomial(), weights = Dat_expanded[[wt]])
  
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


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"

# Primary ###############
# Trial 1: 1-2 days #################
res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))

## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))


### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ret_length"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")


### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))

### od ##############################


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "ret_length", "ddose_prev"),
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = "ret_length",
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")


### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ret_length",
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")




## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")



### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


saveRDS(res_ipw, file = paste0(folder, "res1_12d_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res1_12d_ow.rds"))


# Trial 2: 3-4 days #################
## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("age_st", "sud",
                                       "phy_att_gr2", "social_assist"),  
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "ret_length", 
                                       "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "sud",
                                        "phy_att_gr2", "social_assist"),   
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("age_st", "year", "ret_length", 
                                       "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))

### od ##############################


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev", 
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = c("age_st", "ret_length", "ddose_prev"), 
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev", 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "ret_length", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs =  "ddose_prev", 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev",
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "ddose_prev", 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev",
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res2_34d_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res2_34d_ow.rds"))


# Trial 3 ########################

res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))

## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_ow.rds"))


### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ret_length", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "po", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "po", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("ret_length", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = c("age_st", "year", "ret_length", "ddose_prev"), 
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "ret_length", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev", "prior_gap"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "ret_length", "ddose_prev"),  
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev", "prior_gap"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ret_length", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res3_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res3_ow.rds"))



# Trial 4 ##################

res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))
## incident #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))
### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = c("ddose_prev"), 
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev", "prior_gap"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev", "prior_gap"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res4_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res4_ow.rds"))

# Trial 5 ##################
res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))
## incident #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ddose_prev"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("age_st", "year", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))
### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = "ddose_prev",
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev",
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev",
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev",
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res5_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res5_ow.rds"))




folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/"

# Less than 30 day ###############
# Trial 1: 1-2 days #################
res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))

## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))


### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")


### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))

### od ##############################


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "ddose_prev"),
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = NULL,
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")


### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ddose_prev"),  
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")




## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")



### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


saveRDS(res_ipw, file = paste0(folder, "res1_12d_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res1_12d_ow.rds"))


# Trial 2: 3-4 days #################
## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("age_st", "sud",
                                       "phy_att_gr2", "social_assist"),  
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year",  
                                      "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "sud",
                                        "phy_att_gr2", "social_assist"),   
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("age_st", "year", 
                                       "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))

### od ##############################


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev", 
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = c("age_st", "ddose_prev"), 
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev", 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")


baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs =  "ddose_prev", 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev",
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "ddose_prev", 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev",
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res2_34d_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res2_34d_ow.rds"))


# Trial 3 ########################

res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))

## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_ow.rds"))


### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("age_st", "year","po_preOUD", 
                                       "sed_psy_med", "ddose_prev", "aud"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year","po_preOUD", 
                                      "sed_psy_med", "ddose_prev", 
                                      "udt", "aud"),
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "year","po_preOUD", 
                                        "sed_psy_med", "ddose_prev", "aud"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("age_st", "year","po_preOUD", 
                                       "sed_psy_med", "ddose_prev", 
                                       "udt", "aud"),
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = c("age_st", "year", "ddose_prev"), 
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "ddose_prev"),  
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "year", 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res3_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res3_ow.rds"))



# Trial 4 ##################

res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))
## incident #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))
### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = NULL,
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res4_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res4_ow.rds"))


# Trial 5 ##################
res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))
## incident #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev",
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev",
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "ddose_prev",
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev",
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))
### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = "ddose_prev",
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev",
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev",
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev",
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res5_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res5_ow.rds"))



# 31 to 90 day ###############
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_30_lt90d/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/"

# Trial 1: 1-2 days #################
res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))

## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))


### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev", 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")


### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev", 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))

### od ##############################


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year"),
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = NULL,
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")


### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"),  
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")




## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")



### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


saveRDS(res_ipw, file = paste0(folder, "res1_12d_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res1_12d_ow.rds"))


# Trial 2: 3-4 days #################
## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

#res_ipw$od$incid = NA
res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                              covs = c("age_st", "sud",
                                       "phy_att_gr2", "social_assist"),
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

#res_ow$od$incid = NA
res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                             covs = c("age_st", "year",
                                      "ddose_prev"),
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights,
                                time_var="time_discont",
                                event_var="sens_discont",
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                               covs = c("age_st", "sud",
                                        "phy_att_gr2", "social_assist"),
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                              covs = c("age_st", "year",
                                       "ddose_prev"),
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


# res_ipw$dis$incid = NA
# res_ow$dis$incid = NA
# 

## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "po_preOUD",
                                      "ddose_prev"),  
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = c("age_st", "year", "ddose_prev"), 
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("age_st", "year", "po_preOUD",
                                       "ddose_prev"),  
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")


baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs =  "ddose_prev", 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "ddose_prev", 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res2_34d_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res2_34d_ow.rds"))



# Trial 3 ########################

res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))

## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_ow.rds"))


### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights,
                                time_var="time_od_acute",
                                event_var="sens_od_acute",
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                              covs = c("age_st", "year","po_preOUD",
                                       "sed_psy_med", "ddose_prev", "aud"),
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                             covs = c("age_st", "year","po_preOUD",
                                      "sed_psy_med", "ddose_prev",
                                      "udt", "aud"),
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights,
                                time_var="time_discont",
                                event_var="sens_discont",
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                               covs = c("age_st", "year","po_preOUD",
                                        "sed_psy_med", "ddose_prev", "aud"),
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                              covs = c("age_st", "year","po_preOUD",
                                       "sed_psy_med", "ddose_prev",
                                       "udt", "aud"),
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

# res_ow$od$incid = NA
# res_ow$dis$incid = NA
# res_ipw$od$incid = NA
# res_ipw$dis$incid = NA

## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = c("age_st", "year", "ddose_prev"), 
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "ddose_prev"),  
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "year", 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res3_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res3_ow.rds"))



# Trial 4 ##################

res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))
## incident #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))
### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = NULL,
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res4_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res4_ow.rds"))


# Trial 5 ##################
res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))
## incident #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs =  NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs =  NULL,
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs =  NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))
### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs =  NULL,
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs =  NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs =  NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs =  NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res5_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res5_ow.rds"))


# After 90 day ###############
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_gt90d/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/"

# Trial 1: 1-2 days #################
res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))

## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))


### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = "ddose_prev", 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")


### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev", 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))

### od ##############################


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year"),
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = NULL,
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")


### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"),  
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")




## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")



### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year"), 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


saveRDS(res_ipw, file = paste0(folder, "res1_12d_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res1_12d_ow.rds"))


# Trial 2: 3-4 days #################
## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

#res_ipw$od$incid = NA
res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                              covs = c("age_st", "sud",
                                       "phy_att_gr2", "social_assist"),
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

#res_ow$od$incid = NA
res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                             covs = c("age_st", "year",
                                      "ddose_prev"),
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights,
                                time_var="time_discont",
                                event_var="sens_discont",
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                               covs = c("age_st", "sud",
                                        "phy_att_gr2", "social_assist"),
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                              covs = c("age_st", "year",
                                       "ddose_prev"),
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")


# res_ipw$dis$incid = NA
# res_ow$dis$incid = NA
# 

## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "po_preOUD",
                                      "ddose_prev"),  
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = c("age_st", "year", "ddose_prev"), 
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("age_st", "year", "po_preOUD",
                                       "ddose_prev"),  
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")


baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs =  "ddose_prev", 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "ddose_prev", 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res2_34d_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res2_34d_ow.rds"))



# Trial 3 ########################

res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))

## incident #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")

baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_ow.rds"))


### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights,
                                time_var="time_od_acute",
                                event_var="sens_od_acute",
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                              covs = c("age_st", "year","po_preOUD",
                                       "sed_psy_med", "ddose_prev", "aud"),
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                             covs = c("age_st", "year","po_preOUD",
                                      "sed_psy_med", "ddose_prev",
                                      "udt", "aud"),
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights,
                                time_var="time_discont",
                                event_var="sens_discont",
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                               covs = c("age_st", "year","po_preOUD",
                                        "sed_psy_med", "ddose_prev", "aud"),
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                              covs = c("age_st", "year","po_preOUD",
                                       "sed_psy_med", "ddose_prev",
                                       "udt", "aud"),
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

# res_ow$od$incid = NA
# res_ow$dis$incid = NA
# res_ipw$od$incid = NA
# res_ipw$dis$incid = NA

## prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = c("age_st", "year", "ddose_prev"), 
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("age_st", "year", "ddose_prev"), 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "episode", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 #"ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year"), 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = c("year", "ddose_prev"),  
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "year", 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("year", "ddose_prev"), 
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res3_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res3_ow.rds"))



# Trial 4 ##################

res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))
## incident #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))
### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs = NULL,
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL, 
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res4_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res4_ow.rds"))


# Trial 5 ##################
res_ipw=list(od=list(incid=list(), prev=list(), rzero=list()), 
             dis=list(incid=list(), prev=list(), rzero=list()))

res_ow=list(od=list(incid=list(), prev=list(), rzero=list()), 
            dis=list(incid=list(), prev=list(), rzero=list()))
## incident #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_ow.rds"))

### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs =  NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs =  NULL,
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs =  NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

## prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))
### od ##############################

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs = NULL,
                             cluster_var = "studyid_rev",
                             wt = "wt_trunc")


res_ow$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                            covs =  NULL,
                            cluster_var = "studyid_rev",
                            wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")

res_ow$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs =  NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

## All prevalent #########################
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))

### od ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

res_ipw$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev",
                              wt = "wt_trunc")


res_ow$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                             covs =  NULL,
                             cluster_var = "studyid_rev",
                             wt = "ow_wt")

### dis ##############################
long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)

res_ipw$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev",
                               wt = "wt_trunc")

res_ow$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs =  NULL,
                              cluster_var = "studyid_rev",
                              wt = "ow_wt")

saveRDS(res_ipw, file = paste0(folder, "res5_iptw.rds"))
saveRDS(res_ow, file = paste0(folder, "res5_ow.rds"))




