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

# Trial 1: 1-3 days #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>% 
  filter(trial_cat == "1-3") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 



#### Prevalent new users 
dat1_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "1-3") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
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

ps_fit = ps_model_ow(dat=dat1, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                           out_vars=out_vars, clust_vars=clust_vars, 
                           baselinevars=baselinevars, 
                           continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_ow.rds"))

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

ps_fit = ps_model_ow(dat=dat1_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_prev_ow.rds"))

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

ps_fit = ps_model_ow(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_rprev_ow.rds"))

# Trial 2 #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
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

ps_fit = ps_model_ow(dat=dat2, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_ow.rds"))



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

ps_fit = ps_model_ow(dat=dat2_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_rprev_ow.rds"))

# Trial 3 #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat3 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
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

ps_fit = ps_model_ow(dat=dat3, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_ow.rds"))

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

ps_fit = ps_model_ow(dat=dat3_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))

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

ps_fit = ps_model_ow(dat=dat3_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))

# Trial 1: 1-2 days #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Prevalent new users 
dat1_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()


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

ps_fit = ps_model_ow(dat=dat1, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat1_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))


# Trial 2: 3-4 days #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"
dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat2, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))




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

ps_fit = ps_model_ow(dat=dat2_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))


# Trial 4 #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p4s"
dat4 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Prevalent new users 
dat4_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat4_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
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


ps_fit = ps_model_ow(dat=dat4, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_ow.rds"))



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


ps_fit = ps_model_ow(dat=dat4_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))



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

ps_fit = ps_model_ow(dat=dat4_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))



# Trial 5 #######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
trt_name = "trt_s5"

dat5 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
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

ps_fit = ps_model_ow(dat=dat5, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat5_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat5_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))

saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))


# Less than 30 days ##############
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/"

# Trial 1:  #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>% 
  filter(trial_cat == "1-3" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 



#### Prevalent new users 
dat1_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "1-3" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "1-3" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat1, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_ow.rds"))

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

ps_fit = ps_model_ow(dat=dat1_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_rprev_ow.rds"))

# Trial 2 #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()


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

ps_fit = ps_model_ow(dat=dat2, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_rprev_ow.rds"))


# Trial 3 #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat3 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat3, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_ow.rds"))

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


ps_fit = ps_model_ow(dat=dat3_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat3_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))


# Trial 1: 1-2 days #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2 & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Prevalent new users 
dat1_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2 & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2 & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat1, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat1_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))



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

ps_fit = ps_model_ow(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))

# Trial 2: 3-4 days #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"
dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4 & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4 & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4 & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()



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

ps_fit = ps_model_ow(dat=dat2, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))

# Trial 4 #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p4s"
dat4 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Prevalent new users 
dat4_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat4_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5" & ret_length_cat == "0-30") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat4, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat4_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat4_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))

# Trial 5 #######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
trt_name = "trt_s5"

dat5 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
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
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
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
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
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
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()


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

ps_fit = ps_model_ow(dat=dat5, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat5_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat5_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))



# 31 to 90 days ##############
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_30_lt90d/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/"
# Trial 1:  #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>% 
  filter(trial_cat == "1-3" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 



#### Prevalent new users 
dat1_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "1-3" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "1-3" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat1, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_ow.rds"))

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

ps_fit = ps_model_ow(dat=dat1_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_rprev_ow.rds"))

# Trial 2 #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()


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

# dat2 for ret 30 to 90 has no ">=50%/<=40mg". 
# ps_fit = ps_model_ow(dat=dat2, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
#                      out_vars=out_vars, clust_vars=clust_vars, 
#                      baselinevars=baselinevars, 
#                      continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))
# 
# 
# saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_rprev_ow.rds"))


# Trial 3 #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat3 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat3, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_ow.rds"))

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


ps_fit = ps_model_ow(dat=dat3_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat3_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))


# Trial 1: 1-2 days #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2 & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Prevalent new users 
dat1_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2 & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2 & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat1, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat1_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))



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

ps_fit = ps_model_ow(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))

# Trial 2: 3-4 days #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"
dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4 & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4 & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4 & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()



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

ps_fit = ps_model_ow(dat=dat2, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))

# Trial 4 #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p4s"
dat4 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Prevalent new users 
dat4_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat4_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat4, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat4_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat4_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))

# Trial 5 #######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
trt_name = "trt_s5"

dat5 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14" & ret_length_cat == "31-90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()


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

ps_fit = ps_model_ow(dat=dat5, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat5_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat5_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))




# More than 90 days ##############
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_gt90d/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/"

# Trial 1:  #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>% 
  filter(trial_cat == "1-3" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 



#### Prevalent new users 
dat1_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "1-3" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "1-3" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat1, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_ow.rds"))

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

ps_fit = ps_model_ow(dat=dat1_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_rprev_ow.rds"))

# Trial 2 #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "4" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()


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

ps_fit = ps_model_ow(dat=dat2, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_rprev_ow.rds"))


# Trial 3 #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat3 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat3_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_cat == "5-14" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat3, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_ow.rds"))

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


ps_fit = ps_model_ow(dat=dat3_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat3_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))


# Trial 1: 1-2 days #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2 & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Prevalent new users 
dat1_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2 & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:2 & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat1, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat1_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))



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

ps_fit = ps_model_ow(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))

# Trial 2: 3-4 days #######################
#### incident users 
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"
dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4 & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4 & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 3:4 & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()



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

ps_fit = ps_model_ow(dat=dat2, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat2_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))

# Trial 4 #######################
#### incident users 
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p4s"
dat4 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

#### Prevalent new users 
dat4_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 


#### Multiple zero  
dat4_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "1-5" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()

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

ps_fit = ps_model_ow(dat=dat4, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat4_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat4_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))

# Trial 5 #######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
trt_name = "trt_s5"

dat5 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) 

dat5_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5
  )) %>% 
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 16 ~ ">=16", 
    ddose_prev < 16 ~ "<16"), levels = c("<16", ">=16"), ordered = TRUE)) %>%
  filter(trial_cat == "6-14" & ret_length_cat == ">90") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()


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

ps_fit = ps_model_ow(dat=dat5, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat5_prev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=baselinevars, 
                     continuous_vars=c("age_st", "ddose_prev", "year", "ret_length"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))


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

ps_fit = ps_model_ow(dat=dat5_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99, 
                     out_vars=out_vars, clust_vars=clust_vars, 
                     baselinevars=c(baselinevars, "prior_gap"), 
                     continuous_vars=c("age_st", "ddose_prev", 
                                       "year", "ret_length", "prior_gap"))


saveRDS(ps_fit, file = paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))

