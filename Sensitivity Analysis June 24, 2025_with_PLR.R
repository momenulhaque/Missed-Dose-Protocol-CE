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

# R functions #############
## Table function ###########
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
  
  dat1 <- dat %>% select(all_of(out_vars), all_of(clust_vars), all_of(exposure), all_of(baselinevars))
  
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


create_table = function(data, exposure = "trt", trial=1, 
                        incid=TRUE, baselinevars,
                        continuous_vars=c("age_st", "ddose_prev", "year", "ret_length")){
  
  # baselinevars0 <- c("sexM", "age_st", "age_lt30", "age_30_39", "age_40_49", "age_ge50", "rural", 
  #                    "year", "po_preOUD", "phy_att_gr2", "social_assist", "ret_length", "home_grp_b", "incar_grp_bl", 
  #                    "drug_acute", "cci_gt0", "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
  #                    "udt", "ddose_prev", "clin_size", 'virtual_care', "pss")
  #baselinevars = c(outcome, baselinevars) 
  
  # dat1 = select(data, c(studyid_rev, time0, exposure, baselinevars, outcome))
  # dat1 <- dat1 %>% mutate(across(-c(time0, age_st, ddose_prev, year, ret_length, studyid_rev), as.factor))
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
  
  # if(outcome=="sens_od_acute"){rename_vars[2] <- "Overdose-related acute care visits"}
  
  #  row.names(tr1_table1_before) = rename_vars
  
  if (incid) {
    design <- svydesign(id = ~1, weights = ~wt_trunc, data = dat1)
  } else {
    design <- svydesign(id = ~studyid_rev, weights = ~wt_trunc, data = dat1)
  }
  
  tr1_table1_after = svyCreateTableOne(vars = baselinevars, factorVars = setdiff(baselinevars, cont_vars), 
                                       strata = "trt", addOverall = F, data = design, smd = T, test = T)
  tr1_table1_after <- print(tr1_table1_after, smd = TRUE, showAllLevels = F, exact = "stage", 
                            quote = F, noSpaces = T, printToggle = FALSE, nonnormal = cont_vars, test = FALSE)
  
  #row.names(tr1_table1_after) = rename_vars
  
  smd_column <- tr1_table1_after[, "SMD"]
  table <- cbind(tr1_table1_before, SMD_after = smd_column)
  
  # if(trial %in% c(1, 2)){
  #   colnames(table) = c("Total", "No reduction, n(%)", "Reduction, n(%)", "Increase, n(%)", "SMD", "SMD after weighting")
  # }else if(trial == 4){
  #   colnames(table) = c("Total", "No reduction, n(%)", "Reduction, n(%)", "SMD", "SMD after weighting")
  # }else if (trial == 2 | trial == 3){
  #   colnames(table) = c("Total", "No reduction, n(%)", ">=50% reduction or <=40mg, n(%)", "1-<50% reduction, n(%)", "SMD", "SMD after weighting")
  # }else if (trial == 5){
  #   colnames(table) = c("Total", "No reduction, n(%)", ">=50% reduction or <=8mg, n(%)", "1-<50% reduction, n(%)", "SMD", "SMD after weighting")
  # } else if (trial == 6){
  #   colnames(table) = c("Total", "No reduction, n(%)", ">=50% reduction or <=4mg, n(%)", "1-<50% reduction, n(%)", "SMD", "SMD after weighting")
  # }
  # 
  # table[rename_vars[2], c("SMD", "SMD after weighting")] <- ""
  
  return(table)
  
}

check_imbalanced <- function(data, incid=TRUE, baselinevars,
                             continuous_vars=c("age_st", "ddose_prev", "year", "ret_length")){
  
  #Dat$home_grp_b = ifelse(Dat$home_grp_1 == 1 | Dat$home_grp_2 == 1, 1, 0)
  
  # baselinevars <- c("trt", "sexM", "age_st", "age_lt30", "age_30_39", "age_40_49", "age_ge50", "rural", 
  #                   "year", "po_preOUD", "phy_att_gr2", "social_assist", "ret_length", "home_grp_b", "incar_grp_bl", 
  #                   "drug_acute", "cci_gt0", "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
  #                   "udt", "ddose_prev", "clin_size", 'virtual_care', "virtual_oud", "pss")
  
  dat1 = data
  cont_vars <- continuous_vars[continuous_vars %in% baselinevars]
  
  if (incid) {
    design <- svydesign(id = ~1, weights = ~wt_trunc, data = dat1)
  } else {
    design <- svydesign(id = ~studyid_rev, weights = ~wt_trunc, data = dat1)
  }
  
  tab1 = svyCreateTableOne(vars = baselinevars, factorVars = setdiff(baselinevars, cont_vars), 
                           strata = "trt", addOverall = F, data = design, smd = T, test = T)
  
  smd_values <- print(tab1, smd = TRUE)[-c(1,2), "SMD"]
  smd_clean <- gsub("<", "", smd_values)  # Remove "<" symbols
  smd_df <- data.frame(
    Variable = names(smd_values), 
    SMD = as.numeric(smd_clean)
  )
  imbalanced_vars <- smd_df %>% filter(SMD > 0.1)
  imbalanced_var_names <- gsub(" \\(.*\\)| =.*", "", imbalanced_vars$Variable)
  return(imbalanced_var_names)
}

HR_estimation <- function(Dat, incid = TRUE, outcomes = c("time_od_acute", "sens_od_acute"), 
                          covs = baselinevars, exposure = c("no reversion", "increase", "reversion")) {
  library(survival)
  library(dplyr)
  
  # Keep only needed columns
  if (incid) {
    vars <- c("trt", covs, outcomes, "wt_trunc")
  } else {
    vars <- c("trt", covs, outcomes, "wt_trunc", "studyid_rev")
  }
  
  Dat <- Dat[, vars]
  
  # Build model formula
  covs1 <- paste0(c("trt", covs), collapse = "+")
  cox_formula <- as.formula(paste("Surv(", outcomes[1], ",", outcomes[2], ") ~", covs1))
  
  # Fit weighted Cox model
  fit <- if (incid) {
    coxph(cox_formula, data = Dat, weights = wt_trunc, robust = TRUE, model = TRUE)
  } else {
    coxph(cox_formula, data = Dat, weights = wt_trunc, cluster = studyid_rev, robust = TRUE, model = TRUE)
  }
  
  # Compute summary statistics
  summary_stats <- Dat %>%
    group_by(trt) %>%
    summarise(
      total_individuals = n(),  # unweighted
      events = sum(!!sym(outcomes[2]), na.rm = TRUE),  # unweighted
      weighted_person_years = round(sum(wt_trunc * !!sym(outcomes[1]), na.rm = TRUE) / 365.25, 0),  # rounded 0 digits
      rate_per_py = round((events / (sum(wt_trunc * !!sym(outcomes[1]), na.rm = TRUE) / 365.25)), 2),
      rate_per_10py = round((events / (sum(wt_trunc * !!sym(outcomes[1]), na.rm = TRUE) / 365.25)) * 10, 2),
      rate_per_100py = round((events / (sum(wt_trunc * !!sym(outcomes[1]), na.rm = TRUE) / 365.25)) * 100, 2)
    ) %>%
    ungroup()
  
  return(list(
    model = fit,
    stats = summary_stats
  ))
}



### hdPS related ###########
get_dat_all <- function(dat1){
  
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
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_diagx.sas7bdat", "diagx_trunc", "diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_proc.sas7bdat", "icode_trunc", "proc"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_diag.sas7bdat", "diagcd_trunc", "msp_diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_fee.sas7bdat", "fitm", "msp_fee"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Pnet_drug_rev.sas7bdat", "TC5_CD_trunc", "pnet")
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


get_dat_all_b <- function(dat1){
  
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
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_diagx_nw_b.sas7bdat", "diagx_trunc", "diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Dad_proc_nw_b.sas7bdat", "icode_trunc", "proc"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_diag_nw_b.sas7bdat", "diagcd_trunc", "msp_diag"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Msp_fee_nw_b.sas7bdat", "fitm", "msp_fee"),
    read_and_filter("E:/SAS Data/Missed dose protocol_v2/Pnet_drug_rev_nw_b.sas7bdat", "TC5_CD_trunc", "pnet")
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
  
  # Set up parallel LASSO fitting
  n_cores <- parallel::detectCores()
  cl <- makeCluster(n_cores - 1)
  registerDoParallel(cl)
  
  set.seed(123)
  fit.lasso <- cv.glmnet(x = X, y = Y, nfolds = nfolds, parallel = parallel, alpha = 1, family = "cox")
  stopCluster(cl)
  
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



# Restricted who completed induction ##################
## Trial 1 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")


resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr1.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"),
                                cluster_var = "studyid_rev")



#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"), 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"),
                               cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_induction/resTrial1_iptw.rds")



## Trial 2 ######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr2.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_induction/resTrial2_iptw.rds")


## Trial 3 ######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr3.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_induction/resTrial3_iptw.rds")




## Trial 4 ######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr4.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_induction/resTrial4_iptw.rds")


## Trial 5 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr5.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_induction/resTrial5_iptw.rds")





# High fentanyl area ##################

## Trial 1 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr1.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_high_fent_area/resTrial1_iptw.rds")



## Trial 2 ######################

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr2.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_high_fent_area/resTrial2_iptw.rds")

## Trial 3 ######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr3.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_high_fent_area/resTrial3_iptw.rds")

## Trial 4 ######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr4.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_high_fent_area/resTrial4_iptw.rds")


## Trial 5 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr5.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_high_fent_area/resTrial5_iptw.rds")



# Post fentanyl era: > April 1, 2012 ##################

## Trial 1 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr1.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/post_fentanyl/resTrial1_iptw.rds")



## Trial 2 ######################

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr2.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/post_fentanyl/resTrial2_iptw.rds")


## Trial 3 ######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr3.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/post_fentanyl/resTrial3_iptw.rds")

## Trial 4 ######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr4.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/post_fentanyl/resTrial4_iptw.rds")




## Trial 5 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr5.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/post_fentanyl/resTrial5_iptw.rds")



# Restricted on missed dose amount ##################
## Trial 1 #######################
### low dose < 80 mg ###########
#### incident users 
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr1.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_lt80_16/resTrial1_iptw.rds")


###  high dose >= 80 mg ###########
#### incident users 
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr1.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_ge80_16/resTrial1_iptw.rds")


## Trial 2 ######################
### low dose < 80 mg ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr2.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_lt80_16/resTrial2_iptw.rds")


### high dose >= 80 mg ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr2.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_ge80_16/resTrial2_iptw.rds")


## Trial 3 ######################
### low dose < 80 mg ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr3.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_lt80_16/resTrial3_iptw.rds")

### high dose >= 80 mg ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr3.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_ge80_16/resTrial3_iptw.rds")


## Trial 4 ######################
### low dose < 16 mg ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr4.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_lt80_16/resTrial4_iptw.rds")



### low dose >= 16 mg ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr4.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_ge80_16/resTrial4_iptw.rds")




## Trial 5 ######################
### low dose < 16 mg ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr5.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_lt80_16/resTrial5_iptw.rds")


### low dose >= 16 mg ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr5.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_result_dose_ge80_16/resTrial5_iptw.rds")




# 2017 Guidelines ##################
## Trial 1: 1-2 days #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust_17/ps_weight_tr1.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust_17/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust_17/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Result_robust_17/resTrial1_iptw.rds")



## Trial 2: 3-4 days ##############################################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")



ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust_17/ps_weight_tr2.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust_17/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust_17/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Result_robust_17/resTrial2_iptw.rds")



# Retention length  ###########################

## Trial 1 ###################################
###  retention length <=30 day ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr1.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_lt30d/resTrial1_iptw.rds")



###  retention length 31-90 day ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr1.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt30_lt90d/resTrial1_iptw.rds")



###  retention length >90 day ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr1.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt90d/resTrial1_iptw.rds")



## Trial 2 ######################
###  retention length <=30 day ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr2.rds")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)

# # resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
# #                                covs = c("age_st", "age_ge50", "ret_length",
# #                                         "copd_asthma", "sed_psy_med", 
# #                                         "po", "sud"),
# #                                cluster_var = "studyid_rev")
# # 
# # 
# # long_dat = expand_person_period(Dat=ps_fit$data_weights, 
# #                                 time_var="time_discont", 
# #                                 event_var="sens_discont", 
# #                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("age_st", "age_ge50", "ret_length",
#                                          "copd_asthma", "sed_psy_med", 
#                                          "po", "sud"), 
#                                 cluster_var = "studyid_rev")
# 

resTrial$od$incid <- NULL
resTrial$dis$incid <- NULL

### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_lt30d/resTrial2_iptw.rds")


###  retention length 31-90 day ###########

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr2.rds")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("age_st", "age_ge50", "ret_length",
#                                         "copd_asthma", "sed_psy_med", 
#                                         "po", "sud"),
#                                cluster_var = "studyid_rev")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_discont", 
#                                 event_var="sens_discont", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("age_st", "age_ge50", "ret_length",
#                                          "copd_asthma", "sed_psy_med", 
#                                          "po", "sud"), 
#                                 cluster_var = "studyid_rev")

resTrial$od$incid <- NULL
resTrial$dis$incid <- NULL

### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt30_lt90d/resTrial2_iptw.rds")


###  retention length >90 day ###########

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr2.rds")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("age_st", "age_ge50", "ret_length",
#                                         "copd_asthma", "sed_psy_med", 
#                                         "po", "sud"),
#                                cluster_var = "studyid_rev")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_discont", 
#                                 event_var="sens_discont", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("age_st", "age_ge50", "ret_length",
#                                          "copd_asthma", "sed_psy_med", 
#                                          "po", "sud"), 
#                                 cluster_var = "studyid_rev")
# 

resTrial$od$incid <- NULL
resTrial$dis$incid <- NULL


### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt90d/resTrial2_iptw.rds")


## Trial 3 ######################
### retention length <=90 day ###########

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr3.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_lt30d/resTrial3_iptw.rds")

### retention length 31-90 day ###########

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr3.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt30_lt90d/resTrial3_iptw.rds")


### retention length >90 day ###########

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr3.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt90d/resTrial3_iptw.rds")


## Trial 4 ######################
### retention length <=30 day ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr4.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_lt30d/resTrial4_iptw.rds")


### retention length 31-90 day ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr4.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt30_lt90d/resTrial4_iptw.rds")


### retention length >90 day ###########

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/resTrial4_iptw.rds")
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr4.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt90d/resTrial4_iptw.rds")





## Trial 5 ######################
### retention length <=30 day ###########


source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr5.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_lt30d/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_lt30d/resTrial5_iptw.rds")



### retention length 31-90 day ###########
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr5.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt30_lt90d/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt30_lt90d/resTrial5_iptw.rds")


### retention length >90 day ###########

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr5.rds")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_ret_length_gt90d/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_ret_length_gt90d/resTrial5_iptw.rds")

# Have > 1 y cum experience ##################

## Trial 1 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr1.rds")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = NULL, 
#                                cluster_var = "studyid_rev")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_discont", 
#                                 event_var="sens_discont", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = NULL, 
#                                 cluster_var = "studyid_rev")

resTrial$od$incid <- NULL
resTrial$dis$incid <- NULL

### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr1_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_retention_365/resTrial1_iptw.rds")


## Trial 2 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

# 
# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr2.rds")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("age_st", "age_ge50", "ret_length",
#                                         "copd_asthma", "sed_psy_med", 
#                                         "po", "sud"),
#                                cluster_var = "studyid_rev")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_discont", 
#                                 event_var="sens_discont", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("age_st", "age_ge50", "ret_length",
#                                          "copd_asthma", "sed_psy_med", 
#                                          "po", "sud"), 
#                                 cluster_var = "studyid_rev")
# 

resTrial$od$incid <- NULL
resTrial$dis$incid <- NULL


### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_retention_365/resTrial2_iptw.rds")



## Trial 3 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr3.rds")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("ret_length", "year", "ddose_prev"),  
#                                cluster_var = "studyid_rev")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_discont", 
#                                 event_var="sens_discont", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs =  c("ret_length", "year", "ddose_prev"),  
#                                 cluster_var = "studyid_rev")
# 

resTrial$od$incid <- NULL
resTrial$dis$incid <- NULL


### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_retention_365/resTrial3_iptw.rds")



## Trial 4 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr4.rds")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = NULL, 
#                                cluster_var = "studyid_rev")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_discont", 
#                                 event_var="sens_discont", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = NULL, 
#                                 cluster_var = "studyid_rev")

resTrial$od$incid <- NULL
resTrial$dis$incid <- NULL

### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_retention_365/resTrial4_iptw.rds")



## Trial 5 #######################

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr5.rds")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("year", "ret_length", "ddose_prev"),  
#                                cluster_var = "studyid_rev")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_discont", 
#                                 event_var="sens_discont", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("year", "ret_length", "ddose_prev"), 
#                                 cluster_var = "studyid_rev")

resTrial$od$incid <- NULL
resTrial$dis$incid <- NULL

### Prevalent new-user design 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_retention_365/resTrial5_iptw.rds")




# Outcomes: 30 days gap #######################

ps_folder = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_discont_30/"

## Trial 1 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont30", "time_discont30", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1.rds"))


# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = NULL, 
#                                cluster_var = "studyid_rev")


resTrial$od$incid <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_prev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                               covs = NULL, 
#                               cluster_var = "studyid_rev")

resTrial$od$prev <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_rprev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 30)
# 
# resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = NULL, 
#                                cluster_var = "studyid_rev")


resTrial$od$rzero <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial1_iptw.rds"))

## Trial 2 #######################

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont30", "time_discont30", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2.rds"))


# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("age_st", "age_ge50", "ret_length",
#                                         "copd_asthma", "sed_psy_med", 
#                                         "po", "sud"),
#                                cluster_var = "studyid_rev")


resTrial$od$incid <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_prev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                               covs = c("ret_length", "ddose_prev"),  
#                               cluster_var = "studyid_rev")

resTrial$od$prev <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_rprev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("ret_length", "ddose_prev"), 
#                                cluster_var = "studyid_rev")


resTrial$od$rzero <- NULL


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial2_iptw.rds"))



## Trial 3 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont30", "time_discont30", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3.rds"))


# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("ret_length", "year", "ddose_prev"),  
#                                cluster_var = "studyid_rev")

resTrial$od$incid <- NULL


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                               covs = c("ddose_prev"),
#                               cluster_var = "studyid_rev")


resTrial$od$prev <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("ddose_prev"),
#                                cluster_var = "studyid_rev")
# 


resTrial$od$rzero <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial3_iptw.rds"))



## Trial 4 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont30", "time_discont30", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4.rds"))
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = NULL, 
#                                cluster_var = "studyid_rev")


resTrial$od$incid <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                               covs = NULL, 
#                               cluster_var = "studyid_rev")
# 

resTrial$od$prev <- NULL


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("pss"), 
#                                cluster_var = "studyid_rev")
# 

resTrial$od$rzero <- NULL


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial4_iptw.rds"))



## Trial 5 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont30", "time_discont30", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5.rds"))


# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("year", "ret_length", "ddose_prev"),  
#                                cluster_var = "studyid_rev")
# 

resTrial$od$incid <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                               covs = NULL, 
#                               cluster_var = "studyid_rev")
# 


resTrial$od$prev <- NULL

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev.rds"))

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = NULL, 
#                                cluster_var = "studyid_rev")
# 

resTrial$od$rzero <- NULL


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont30", 
                                event_var="sens_discont30", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial5_iptw.rds"))



# Outcomes: All-cause acute care visits #######################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_outcome/"

## Trial 1 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial1_iptw.rds"))

## Trial 2 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial2_iptw.rds"))


## Trial 3 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial3_iptw.rds"))


## Trial 4 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial4_iptw.rds"))


## Trial 5 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_all_acute", 
                                event_var="sens_all_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial5_iptw.rds"))






# No truncation at 12m #######################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_no_truncation/"

## Trial 1 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "s_od_acute", "t_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial1_iptw.rds"))

#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



primary_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
sa_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/"

full11 <- readRDS(file=paste0(primary_folder, "resTrial1_iptw.rds"))

ps_folder = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_no_truncation/"
resTrial <- readRDS(file=paste0(folder, "resTrial1_iptw.rds"))

resTrial$od$rzero <- full11$od$rzero
resTrial$dis$rzero <- full11$dis$rzero
saveRDS(resTrial, file=paste0(folder, "resTrial1_iptw.rds"))

## Trial 2 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "s_od_acute", "t_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial2_iptw.rds"))


## Trial 3 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "s_od_acute", "t_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial3_iptw.rds"))


## Trial 4 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "s_od_acute", "t_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial4_iptw.rds"))


## Trial 5 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("s_discont", "t_discont", 
             "s_od_acute", "t_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_od_acute", 
                                event_var="s_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="t_discont", 
                                event_var="s_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial5_iptw.rds"))



# Restricted to highest fentanyl area ##################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_high_fent_area/"

## Trial 1 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial1_iptw.rds"))

## Trial 2 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial2_iptw.rds"))


## Trial 3 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial3_iptw.rds"))


## Trial 4 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial4_iptw.rds"))


## Trial 5 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial5_iptw.rds"))


# 30 day washout for prevalent users ##################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_30d_washout/"


## Trial 1 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial1_iptw.rds"))

## Trial 2 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial2_iptw.rds"))


## Trial 3 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial3_iptw.rds"))


## Trial 4 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial4_iptw.rds"))


## Trial 5 #######################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial5_iptw.rds"))


# Composite outcome #######################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

## Trial 1 #######################

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "1-3") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit = ps_model_ed(dat=dat1, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr1.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")


### Prevalent new-user design 
dat1_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "1-3") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

ps_fit = ps_model_ed(dat=dat1_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr1_prev.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr1_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "1-3") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

ps_fit = ps_model_ed(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", "year", 
                                       "ret_length", "prior_gap"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr1_rprev.rds"))


ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr1_rprev.rds"))




long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial1_iptw.rds"))


## Trial 2 #######################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"
dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "4") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "4") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
    
  ) %>%
  filter(trial_cat == "4") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()



resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit = ps_model_ed(dat=dat2, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr2.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr2.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit = ps_model_ed(dat=dat2_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr2_prev.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr2_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)



resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit = ps_model_ed(dat=dat2_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr2_rprev.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr2_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)



resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial2_iptw.rds"))



## Trial 3 #######################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"
dat3 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "5-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "5-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "5-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()



resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")


clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit = ps_model_ed(dat=dat3, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr3.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr3.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ret_length", "year", "ddose_prev"),  
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit = ps_model_ed(dat=dat3_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr3_prev.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr3_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)



resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit = ps_model_ed(dat=dat3_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr3_rprev.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr3_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)

resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ddose_prev"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial3_iptw.rds"))



## Trial 4 #######################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p4s"

dat4 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "1-5") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Prevalent new users 
dat4_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "1-5") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Multiple zero  
dat4_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "1-5") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()



resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit = ps_model_ed(dat=dat4, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr4.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr4.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 

ps_fit = ps_model_ed(dat=dat4_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr4_prev.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr4_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit = ps_model_ed(dat=dat4_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr4_rprev.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr4_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("pss"), 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial4_iptw.rds"))



## Trial 5 #######################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

trt_strategy = c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
trt_name = "trt_s5"

dat5 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "6-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "6-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(
    time_comp_outcome = pmin(time_od_acute, time_discont, na.rm = TRUE),
    sens_comp_outcome = pmax(sens_od_acute, sens_discont, na.rm = TRUE),
    time_od_comp_outcome = pmin(time_od_acute, time_od_death, na.rm = TRUE),
    sens_od_comp_outcome = pmax(sens_od_acute, sens_od_death, na.rm = TRUE),
    time_od_acute_30d = pmin(t_od_acute, 30),
    sens_od_acute_30d = ifelse(t_od_acute >= 30, 0, s_od_acute)
  ) %>%
  filter(trial_cat == "6-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()




resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit = ps_model_ed(dat=dat5, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr5.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr5.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")



### Prevalent new-user design 
ps_fit = ps_model_ed(dat=dat5_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr5_prev.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr5_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


#### Multiple zero 
ps_fit = ps_model_ed(dat=dat5_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr5_rprev.rds"))

ps_fit <- readRDS(file = paste0(folder, "ps_weight_tr5_rprev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_comp_outcome", 
                                event_var="sens_od_comp_outcome", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_comp_outcome", 
                                event_var="sens_comp_outcome", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial5_iptw.rds"))


# od acute in 30 days #######################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_outcome_od_acute_30d/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

## Trial 1 #######################

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 

out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


# # long_dat = expand_person_period(Dat=ps_fit$data_weights, 
# #                                 time_var="time_comp_outcome", 
# #                                 event_var="sens_comp_outcome", 
# #                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = NULL, 
#                                 cluster_var = "studyid_rev")


### Prevalent new-user design 

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = NULL, 
#                                cluster_var = "studyid_rev")
# 

#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")

# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 30)
# 
# 
# resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = NULL, 
#                                 cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial1_iptw.rds"))


## Trial 2 #######################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_outcome_od_acute_30d/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")


resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),
                               cluster_var = "studyid_rev")

# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("age_st", "age_ge50", "ret_length",
#                                          "copd_asthma", "sed_psy_med", 
#                                          "po", "sud"), 
#                                 cluster_var = "studyid_rev")
# 
# 

### Prevalent new-user design 

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# 
# resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("ret_length", "ddose_prev"), 
#                                cluster_var = "studyid_rev")


#### Multiple zero 

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# 
# resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("ret_length", "ddose_prev"), 
#                                 cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial2_iptw.rds"))



## Trial 3 #######################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_outcome_od_acute_30d/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")


resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")


clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "year", "ddose_prev"),  
                               cluster_var = "studyid_rev")


# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs =  c("ret_length", "year", "ddose_prev"),  
#                                 cluster_var = "studyid_rev")



### Prevalent new-user design 

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"),
                              cluster_var = "studyid_rev")


# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# 
# resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = c("ddose_prev"),
#                                cluster_var = "studyid_rev")
# 

#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ddose_prev"),
                               cluster_var = "studyid_rev")

# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("ddose_prev"), 
#                                 cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial3_iptw.rds"))



## Trial 4 #######################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_outcome_od_acute_30d/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")

# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = NULL, 
#                                 cluster_var = "studyid_rev")
# 


### Prevalent new-user design 

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")

# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = NULL, 
#                                cluster_var = "studyid_rev")


#### Multiple zero 

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("pss"), 
                               cluster_var = "studyid_rev")
# 
# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("pss"), 
#                                 cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial4_iptw.rds"))



## Trial 5 #######################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_comp_outcome/"

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_outcome_od_acute_30d/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()), 
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user design 
out_vars = c("sens_comp_outcome", "time_comp_outcome", 
             "sens_od_comp_outcome", "time_od_comp_outcome",
             "time_od_acute_30d", "sens_od_acute_30d",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute",
             "sens_discont", "time_discont")

clust_vars = c("studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")

# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = c("year", "ret_length", "ddose_prev"), 
#                                 cluster_var = "studyid_rev")



### Prevalent new-user design 

ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL, 
                              cluster_var = "studyid_rev")


# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                covs = NULL, 
#                                cluster_var = "studyid_rev")
# 

#### Multiple zero 
ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute_30d", 
                                event_var="sens_od_acute_30d", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")

# 
# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_comp_outcome", 
#                                 event_var="sens_comp_outcome", 
#                                 interval_length = 14)
# 
# 
# resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
#                                 covs = NULL, 
#                                 cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial5_iptw.rds"))




# Spline analysis #####################
library(survival)
library(splines)
library(marginaleffects)
library(ggplot2)
library(dplyr)


comparison_colors = c(
  "Dose reduction\nvs\nNo dose change"              = "#000000",  # black
  "Dose increase\nvs\nNo dose change"               = "#D62728",  # red
  "1-<50% reduction\nvs\nNo dose change"            = "#2CA02C",  # Green
  ">=50% reduction/<=40mg\nvs\nNo dose change"      = "#FF7F0E",  # orange
  ">=50% reduction/<=4mg\nvs\nNo dose change"       = "#CC79A7"   # Purple
)


# Retention time ###########
## Trial 1 #######################
### overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_rprev.rds")
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr1_rprev.rds")


### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr1_rprev.rds")



## Trial 2 ###############################
### Overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=40mg - no reversion" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=40mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr2_rprev.rds")


### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=40mg - no reversion" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=40mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr2_rprev.rds")



## Trial 3 ###############################
### Overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=40mg - no reversion" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=40mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr3_rprev.rds")

### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=40mg - no reversion" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=40mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr3_rprev.rds")


## Trial 4 #######################
### overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev.rds")
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr4_rprev.rds")


### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)


datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

mfx$ret_length <- datagrid_values$ret_length


mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr4_rprev.rds")



## Trial 5 ###############################
### Overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)

# Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)


mfx$ret_length <- datagrid_values$ret_length

# Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=4mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=4mg - no reversion" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=4mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr5_rprev.rds")


### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev.rds")


cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)


datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)


mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)


mfx$ret_length <- datagrid_values$ret_length

mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=4mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=4mg - no reversion" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=4mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr5_rprev.rds")



# hdPS analysis  ###########################

trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>%
  filter(trial_cat == "4") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>%
  filter(trial_cat == "4") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number(),
    unique_id = paste0(studyid_rev, "_", time0)  # Row number = order of gap_days within group
  ) %>%
  ungroup()


out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")



## Trial 1 #################################
resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ########
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")



ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr1.rds"))

# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars,
                                                  clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))
ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(ps_fit$data_weights,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=ps_fit$data_weights, 
                        exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr1.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  NULL,  
                                cluster_var = "studyid_rev")


### Prevalent new users ##########
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr1_prev.rds"))

dim(ps_fit$data_weights)

# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars,
                                                  clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))
ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(ps_fit$data_weights,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=ps_fit$data_weights, 
                        exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr2_prev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")

### All prevalent users #######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit = readRDS(file = "H:/SAS-Work/Missed doses Paper v2/hdPS/ps_weight_tr1_rprev.rds")

# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, "prior_gap",
                                                  clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))


ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# vars.empirical <- select_top200_empirical_covariates(ps_fit$data_weights,
#                                                      time_var = "time_discont",
#                                                      event_var = "sens_discont",
#                                                      max_vars = 200,
#                                                      nfolds = 5,
#                                                      parallel = FALSE,
#                                                      vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=ps_fit$data_weights, 
                        exposure = "trt", 
                        trunc_low = 0, trunc_high = 0.99,
                        out_vars=out_vars, clust_vars=clust_vars,
                        baselinevars=c(baselinevars, "prior_gap",
                                       vars.emp),
                        continuous_vars=c("age_st", "ddose_prev",
                                          "year", "ret_length",
                                          "prior_gap"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr1_rprev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights,
                                time_var="time_od_acute",
                                event_var="sens_od_acute",
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                               covs = c("ret_length", "ddose_prev"),
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights,
                                time_var="time_discont",
                                event_var="sens_discont",
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt",
                                covs = c("ret_length", "ddose_prev"),
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial1_iptw.rds"))


#resTrial <- readRDS(file=paste0(folder, "resTrial1_iptw.rds"))

Publish::publish(resTrial$od$rzero)


## Trial 2 ######################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")


resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user #######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr2.rds"))


# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))
ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(ps_fit$data_weights,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=ps_fit$data_weights, 
                        exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr2.rds"))



resTrial$od$incid = NA


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


### Prevalent new user #############
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr2_prev.rds"))
# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(ps_fit$data_weights,
                                   time_var = "time_discont",
                                   event_var = "sens_discont",
                                   max_vars = 200,
                                   nfolds = 5,
                                   vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=ps_fit$data_weights, 
                        exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr2_prev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")

### All prevalent users ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"

out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr2_rprev.rds"))

# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, "prior_gap",
                                                  clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(dat.all=ps_fit$data_weights,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_rprev.rds")

dim(ps_fit$data_weights)

ps_fit = ps_model_lasso(dat=ps_fit$data_weights, 
                        exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr2_rprev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"), 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")


saveRDS(resTrial, file=paste0(folder, "resTrial2_iptw.rds"))



## Trial 3 #######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat3 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>%
  filter(trial_cat == "5-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))

dat3_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>%
  filter(trial_cat == "5-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))

dat3_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>%
  filter(trial_cat == "5-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number(),
    unique_id = paste0(studyid_rev, "_", time0)   # Row number = order of gap_days within group
  ) %>%
  ungroup()



resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ########

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

dat.all <- get_dat_all(dat1=dat3)

# empirical covariates
vars.emp <- setdiff(names(dat.all), c(names(dat3)))


dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr3.rds"))

ps_fit = readRDS(file = paste0(folder, "ps_weight_tr3.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"),  
                                cluster_var = "studyid_rev")


### Prevalent new users ##########
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

dat.all <- get_dat_all(dat1=dat3_prev)

# empirical covariates
vars.emp <- setdiff(names(dat.all), c(names(dat3_prev)))


dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr3_prev.rds"))

ps_fit = readRDS(file = paste0(folder, "ps_weight_tr3_prev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ddose_prev"), 
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs =  c("ddose_prev"),  
                               cluster_var = "studyid_rev")

### All prevalent users #######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

dat.all <- get_dat_all(dat1=dat3_rprev)

# empirical covariates
vars.emp <- setdiff(names(dat.all), c(names(dat3_rprev)))


dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars=out_vars, clust_vars=clust_vars,
                        baselinevars=c(baselinevars, "prior_gap", vars.empirical),
                        continuous_vars=c("age_st", "ddose_prev",
                                          "year", "ret_length",
                                          "prior_gap"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr3_rprev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs =  c("ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs =  c("ddose_prev"),  
                                cluster_var = "studyid_rev")



saveRDS(resTrial, file=paste0(folder, "resTrial3_iptw.rds"))



## Trial 4 ######################
# trt_strategy = c("no reversion", "reversion", "increase")
# trt_name = "trt_p4s"
# dat4 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_nw_b.sas7bdat") %>%
#   filter(trial_cat == "1-5") %>%
#   rename(trt = !!trt_name) %>%
#   filter(trt %in% trt_strategy) %>%
#   mutate(trt = factor(trt, levels = trt_strategy),
#          unique_id = paste0(studyid_rev, "_", time0))
# 
# #### Prevalent new users
# dat4_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_nw_b.sas7bdat") %>%
#   filter(trial_cat == "1-5") %>%
#   rename(trt = !!trt_name) %>%
#   filter(trt %in% trt_strategy) %>%
#   mutate(trt = factor(trt, levels = trt_strategy),
#          unique_id = paste0(studyid_rev, "_", time0))
# 
# 
# #### Multiple zero
# dat4_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_nw_b.sas7bdat") %>%
#   filter(trial_cat == "1-5") %>%
#   rename(trt = !!trt_name) %>%
#   filter(trt %in% trt_strategy) %>%
#   mutate(trt = factor(trt, levels = trt_strategy)) %>%
#   group_by(studyid_rev, episode) %>%
#   mutate(
#     prior_gap = cumsum(gap_days),  # Cumulative sum including current row
#     prior_gap_num = row_number(),
#     unique_id = paste0(studyid_rev, "_", time0)   # Row number = order of gap_days within group
#   ) %>%
#   ungroup()



resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ########

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

# dat.all <- get_dat_all(dat1=dat4)
# dim(dat.all) # 460

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr4.rds"))



# dat.all <- dat.all %>%
#   mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
#          time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))
# 
# vars.empirical <- select_top200_empirical_covariates(dat.all, dat4,
#                                                      time_var = "time_discont",
#                                                      event_var = "sens_discont",
#                                                      max_vars = 200,
#                                                      nfolds = 5)

# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(names(dat2), 
                                                  "exposure", "ps", "wt_stab",
                                                  "wt_trunc"))
ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(ps_fit$data_weights,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr4.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL,  
                                cluster_var = "studyid_rev")


### Prevalent new users ##########
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

# dat.all <- get_dat_all(dat1=dat4_prev)
# 
# dat.all <- dat.all %>%
#   mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
#          time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))
# 
# vars.empirical <- select_top200_empirical_covariates(dat.all, dat4_prev,
#                                                      time_var = "time_discont",
#                                                      event_var = "sens_discont",
#                                                      max_vars = 200,
#                                                      nfolds = 5)

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev.rds"))


# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(names(dat2), 
                                                  "exposure", "ps", "wt_stab",
                                                  "wt_trunc"))
ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(ps_fit$data_weights,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)



ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr4_prev.rds"))



long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")

### All prevalent users #######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

# dat.all <- get_dat_all(dat1=dat4_rprev)
# 
# dat.all <- dat.all %>%
#   mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
#          time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))
# 
# vars.empirical <- select_top200_empirical_covariates(dat.all, dat4_rprev,
#                                                      time_var = "time_discont",
#                                                      event_var = "sens_discont",
#                                                      max_vars = 200,
#                                                      nfolds = 5)

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev.rds"))


# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(names(dat2_rprev), 
                                                  "exposure", "ps", "wt_stab",
                                                  "wt_trunc"))
ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(ps_fit$data_weights,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)



ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars=out_vars, clust_vars=clust_vars,
                        baselinevars=c(baselinevars, "prior_gap", vars.empirical),
                        continuous_vars=c("age_st", "ddose_prev",
                                          "year", "ret_length",
                                          "prior_gap"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr4_rprev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial4_iptw.rds"))



## Trial 5 #######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
trt_name = "trt_s5"

dat5 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5),
    ret_length_cat = factor(case_when(
      ret_length <= 30 ~ "0-30",
      ret_length > 30 & ret_length <= 90 ~ "31-90",
      ret_length > 90 ~ ">90"
    ), levels = c("0-30", "31-90", ">90"), ordered = TRUE),
    high_ddose = factor(case_when(
      ddose_prev >= 16 ~ ">=16",
      ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))

dat5_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5),
    ret_length_cat = factor(case_when(
      ret_length <= 30 ~ "0-30",
      ret_length > 30 & ret_length <= 90 ~ "31-90",
      ret_length > 90 ~ ">90"
    ), levels = c("0-30", "31-90", ">90"), ordered = TRUE),
    high_ddose = factor(case_when(
      ddose_prev >= 16 ~ ">=16",
      ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))

dat5_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5),
    ret_length_cat = factor(case_when(
      ret_length <= 30 ~ "0-30",
      ret_length > 30 & ret_length <= 90 ~ "31-90",
      ret_length > 90 ~ ">90"
    ), levels = c("0-30", "31-90", ">90"), ordered = TRUE),
    high_ddose = factor(case_when(
      ddose_prev >= 16 ~ ">=16",
      ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number(),
    unique_id = paste0(studyid_rev, "_", time0)   # Row number = order of gap_days within group
  ) %>%
  ungroup()




resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ########

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

dat.all <- get_dat_all_b(dat1=dat5)

# empirical covariates
vars.emp <- setdiff(names(dat.all), c(names(dat5)))
                    
                    
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr5.rds"))

ps_fit = readRDS(file = paste0(folder, "ps_weight_tr5.rds"))

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"),  
                                cluster_var = "studyid_rev")


### Prevalent new users ##########
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

dat.all <- get_dat_all_b(dat1=dat5_prev)

# empirical covariates
vars.emp <- setdiff(names(dat.all), c(names(dat5_prev)))


dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr5_prev.rds"))

ps_fit = readRDS(file = paste0(folder, "ps_weight_tr5_prev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")

### All prevalent users #######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

dat.all <- get_dat_all_b(dat1=dat5_rprev)

# empirical covariates
vars.emp <- setdiff(names(dat.all), c(names(dat5_rprev)))


dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

vars.empirical <- select_top200_empirical_covariates(dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200,
                                                     nfolds = 5,
                                                     vars.empirical=vars.emp)

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars=out_vars, clust_vars=clust_vars,
                        baselinevars=c(baselinevars, "prior_gap", vars.empirical),
                        continuous_vars=c("age_st", "ddose_prev",
                                          "year", "ret_length",
                                          "prior_gap"))

saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr5_rprev.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL, 
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL, 
                                cluster_var = "studyid_rev")

saveRDS(resTrial, file=paste0(folder, "resTrial5_iptw.rds"))



# hdPS with ML analysis  ###########################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/sens_hdPS_ML/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"

out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")



## Trial 1 #################################
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr1.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr1.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL,  
                                cluster_var = "studyid_rev")

### Prevalent new users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr1_prev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr1_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = c("ret_length", "ddose_prev"),  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")

### All prevalent users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_rprev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, "prior_gap",
                                                  clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt",
                                                  "wt_stab", "prior_gap_num",
                                                  "wt_trunc"))




ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "prior_gap", 
                                              "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr1_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"),  
                                cluster_var = "studyid_rev")


saveRDS(resTrial, file = paste0(folder, "resTrial1_iptw.rds"))


## Trial 2 ##################
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"
resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

### Incident users ####################

out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr2.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr2.rds"))


resTrial$od$incid = NA


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("age_st", "age_ge50", "ret_length",
                                         "copd_asthma", "sed_psy_med", 
                                         "po", "sud"),  
                                cluster_var = "studyid_rev")



### Prevalent new users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr2_prev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr2_prev.rds"))


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev",  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "ddose_prev",  
                               cluster_var = "studyid_rev")



### All prevalent users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr2_rprev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, "prior_gap",
                                                  clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))


ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "prior_gap", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr2_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "udt",  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = "udt",  
                                cluster_var = "studyid_rev")


saveRDS(resTrial, file = paste0(folder, "resTrial2_iptw.rds"))





## Trial 3 ######################
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

### Incident users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr3.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr3.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"),  
                                cluster_var = "studyid_rev")

### Prevalent new users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr3_prev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr3_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev",  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "ddose_prev",  
                               cluster_var = "studyid_rev")

### All prevalent users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr3_rprev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, "prior_gap",
                                                  clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))


ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "prior_gap", 
                                              "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr3_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "ddose_prev",  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "ddose_prev",  
                               cluster_var = "studyid_rev")


saveRDS(resTrial, file = paste0(folder, "resTrial3_iptw.rds"))


## Trial 4 ######################
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr4.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr4.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = NULL,  
                                cluster_var = "studyid_rev")

### Prevalent new users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr4_prev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr4_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev")

### All prevalent users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr4_rprev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, "prior_gap",
                                                  clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))


ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "prior_gap", 
                                              "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr4_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = "pss",  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = "pss",  
                               cluster_var = "studyid_rev")


saveRDS(resTrial, file = paste0(folder, "resTrial4_iptw.rds"))





## Trial 5 ######################
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"

source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr5.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr5.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("year", "ret_length", "ddose_prev"),  
                               cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("year", "ret_length", "ddose_prev"),  
                                cluster_var = "studyid_rev")

### Prevalent new users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr5_prev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))

ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr5_prev.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$prev = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev")

### All prevalent users ####################
out_vars = c("sens_discont", "time_discont",
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")

clust_vars = c("studyid_rev", "time0")
baselinevars <- c("sexM", "age_st", "age_lt30", "age_30_39",
                  "age_40_49", "age_ge50", "rural", "year", "po_preOUD",
                  "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                  "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma",
                  "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
                  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr5_rprev.rds"))

vars.emp <- setdiff(names(ps_fit$data_weights), c(out_vars, "prior_gap",
                                                  clust_vars, baselinevars,  
                                                  "exposure", "trt", "ps", "wt_stab",
                                                  "wt_trunc"))


ps_fit_ml <- ps_model_xgb(dat=ps_fit$data_weights,
                          exposure = "trt", 
                          trunc_low = 0, 
                          trunc_high = 0.99, 
                          out_vars=out_vars, 
                          clust_vars=clust_vars,
                          baselinevars=baselinevars,
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "prior_gap", 
                                              "ret_length"),
                          emp_vars=vars.emp)

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr5_rprev.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                              covs = NULL,  
                              cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


resTrial$dis$rzero = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = NULL,  
                               cluster_var = "studyid_rev")


saveRDS(resTrial, file = paste0(folder, "resTrial5_iptw.rds"))


