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

# Incident users ###########################
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
## Trial 1: 1-3 days #################################

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
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))



dat.all <- get_dat_all_incid_m(dat1=dat1)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat1))
length(vars.empirical)
# 167

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))
saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1.rds")


## Trial 2: 4 days ######################
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



dat.all <- get_dat_all_incid_m(dat1=dat2)
dim(dat.all)
# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat2))
# 3 empirical covariates were selected

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2.rds")


## Trial 3 ######################
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


dat.all <- get_dat_all_incid_m(dat1=dat3)
dim(dat.all)
# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat3))
# 13 empirical covariates were selected
ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr3.rds")

## Trial 4 ######################
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p4s"
dat4 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  filter(trial_cat == "1-5") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))

table(dat4$trt)

dat.all <- get_dat_all_incid_b(dat1=dat4)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat4))
length(vars.empirical)
# 239
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)



ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr4.rds")


## Trial 5 ######################
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


dat.all <- get_dat_all_incid_b(dat1=dat5)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat5))
length(vars.empirical)
# 239
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)



ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr5.rds")

## Trial 1: 1-2 days #################################

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
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))



dat.all <- get_dat_all_incid_m(dat1=dat1)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat1))
length(vars.empirical)
# 167

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))
saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_12d_ow.rds")


## Trial 2: 3-4 days ######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2 <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>%
  filter(trial_gap %in% 3:4) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))



dat.all <- get_dat_all_incid_m(dat1=dat2)
dim(dat.all)
# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat2))
# 3 empirical covariates were selected

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_34d_ow.rds")



# Prevalent new users ###########################
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

## Trial 1: 1-3 days #################################
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
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
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))


dat.all <- get_dat_all_prev_m(dat1=dat1_prev)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat1_prev))
length(vars.empirical)
# 772
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_prev.rds")


## Trial 2: 4 days ######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>%
  filter(trial_cat == "4") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))



dat.all <- get_dat_all_prev_m(dat1=dat2_prev)
dim(dat.all)
# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat2_prev))
length(vars.empirical)
# 239
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_prev.rds")


## Trial 3 ######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

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


dat.all <- get_dat_all_prev_m(dat1=dat3_prev)
dim(dat.all)
# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat3_prev))
length(vars.empirical)

dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr3_prev.rds")



## Trial 4 ######################
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p4s"
dat4_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  filter(trial_cat == "1-5") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))

dat.all <- get_dat_all_incid_b(dat1=dat4_prev)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat4_prev))
length(vars.empirical)
# 239
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr4_prev.rds")


## Trial 5 ######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
trt_name = "trt_s5"

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


dat.all <- get_dat_all_incid_b(dat1=dat5_prev)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat5_prev))
length(vars.empirical)
# 239
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)



ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr5_prev.rds")

## Trial 1: 1-2 days #################################
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
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
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))


dat.all <- get_dat_all_prev_m(dat1=dat1_prev)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat1_prev))
length(vars.empirical)
# 772
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_12d_prev_ow.rds")


## Trial 2: 3-4 days ######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2_prev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>%
  filter(trial_gap %in% 3:4) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy),
         unique_id = paste0(studyid_rev, "_", time0))



dat.all <- get_dat_all_prev_m(dat1=dat2_prev)
dim(dat.all)
# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat2_prev))
length(vars.empirical)
# 239
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_34d_prev_ow.rds")


# All prevalent new users ###########################
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

## Trial 1: 1-3 days #################################
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
dat1_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE),
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80",
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  filter(trial_gap %in% 1:3) %>%
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


# dat.all <- get_dat_all_rprev_m(dat1=dat1_rprev)
# 
# # empirical covariates
# vars.empirical <- setdiff(names(dat.all), names(dat1_rprev))
# length(vars.empirical)
# # 772
# dat.all <- dat.all %>%
#   mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
#          time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))
# 
# # Select top 200 empirical covariate
# vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
#                                                      time_var = "time_discont",
#                                                      event_var = "sens_discont",
#                                                      max_vars = 200, nfolds = 5,
#                                                      vars.empirical=vars.empirical)
# 
# 
# ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
#                         out_vars, clust_vars,
#                         baselinevars = c(vars.empirical, baselinevars),
#                         continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))



ps_fit = ps_model_lasso(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, "prior_gap", "ret_length_cat",
                                         "maint_ddose", "maint_dur"),
                        continuous_vars = c("age_st", "ddose_prev", 
                                            "year", "ret_length", "prior_gap", 
                                            "maint_ddose", "maint_dur"))


saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_rprev.rds")


## Trial 2: 4 days ######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

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



dat.all <- get_dat_all_prev_m(dat1=dat2_rprev)
dim(dat.all)
# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat2_rprev))
length(vars.empirical)
# 34
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_rprev.rds")


## Trial 3 ######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

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


dat.all <- get_dat_all_prev_m(dat1=dat3_rprev)
dim(dat.all)
# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat3_rprev))
length(vars.empirical)

dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr3_rprev.rds")



## Trial 4 ######################
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p4s"
dat4_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  filter(trial_cat == "1-5") %>%
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

table(dat4_rprev$trt)

dat.all <- get_dat_all_rprev_b(dat1=dat4_rprev)

table(dat.all$trt)
dim(dat.all)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat4_rprev))
length(vars.empirical)

vars.empirical <- vars.empirical[1:200]
# 
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# # Select top 200 empirical covariate
# vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
#                                                      time_var = "time_discont",
#                                                      event_var = "sens_discont",
#                                                      max_vars = 200, nfolds = 5,
#                                                      vars.empirical=vars.empirical)
# 

ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars, "prior_gap"),
                        continuous_vars = c("age_st", "ddose_prev", 
                                            "year", "ret_length", "prior_gap"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr4_rprev.rds")


## Trial 5 ######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
trt_name = "trt_s5"

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
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number(),
    unique_id = paste0(studyid_rev, "_", time0)   # Row number = order of gap_days within group
  ) %>%
  ungroup()


dat.all <- get_dat_all_rprev_b(dat1=dat5_rprev)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat5_rprev))
length(vars.empirical)
# 239
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)



ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr5_rprev.rds")

## Trial 1: 1-2 days #################################
trt_strategy = c("no reversion", "reversion", "increase")
trt_name = "trt_p1"
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
    prior_gap_num = row_number(),
    unique_id = paste0(studyid_rev, "_", time0)   # Row number = order of gap_days within group
  ) %>%
  ungroup()


dat.all <- get_dat_all_rprev_m(dat1=dat1_rprev)

# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat1_rprev))
length(vars.empirical)
# 772
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(vars.empirical, "prior_gap", baselinevars),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length", "prior_gap"))



# ps_fit = ps_model_lasso(dat=dat1_rprev, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
#                         out_vars, clust_vars,
#                         baselinevars = c(baselinevars, "prior_gap", "ret_length_cat",
#                                          "maint_ddose", "maint_dur"),
#                         continuous_vars = c("age_st", "ddose_prev", 
#                                             "year", "ret_length", "prior_gap", 
#                                             "maint_ddose", "maint_dur"))


saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_12d_rprev_ow.rds")


## Trial 2: 3-4 days ######################
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2_rprev <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  )) %>%
  filter(trial_gap %in% 3:4) %>%
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



dat.all <- get_dat_all_prev_m(dat1=dat2_rprev)
dim(dat.all)
# empirical covariates
vars.empirical <- setdiff(names(dat.all), names(dat2_rprev))
length(vars.empirical)
# 34
dat.all <- dat.all %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# Select top 200 empirical covariate
vars.empirical <- select_top200_empirical_covariates(dat.all=dat.all,
                                                     time_var = "time_discont",
                                                     event_var = "sens_discont",
                                                     max_vars = 200, nfolds = 5,
                                                     vars.empirical=vars.empirical)


ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0.01, trunc_high = 0.99,
                        out_vars, clust_vars,
                        baselinevars = c(baselinevars, vars.empirical),
                        continuous_vars = c("age_st", "ddose_prev", "year", "ret_length"))

saveRDS(ps_fit, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_34d_rprev_ow.rds")


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



## Trial 1: 1-3 days #################################
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
ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr1_prev.rds"))

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
ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr1_rprev.rds"))

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


## Trial 2: 4 days ######################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")


resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user #######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")


ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr2.rds"))


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

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr2_prev.rds"))

ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

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

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr2_rprev.rds"))

ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))


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

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ########

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr3.rds"))

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
ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev.rds"))


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

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev.rds"))

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
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
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


resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ########

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr4.rds"))
dim(ps_fit$data_weights)

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

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev.rds"))

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

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev.rds"))

dim(ps_fit$data_weights)

table(ps_fit$data_weights$trt)

# empirical covariates
vars.emp <- setdiff(names(ps_fit$data_weights), c(names(dat4_rprev), 
                                                  "exposure", "ps", "wt_stab",
                                                  "wt_trunc"))
ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

# vars.empirical <- select_top200_empirical_covariates(ps_fit$data_weights,
#                                                      time_var = "time_discont",
#                                                      event_var = "sens_discont",
#                                                      max_vars = 200,
#                                                      nfolds = 5,
#                                                      vars.empirical=vars.emp)

# ps_fit = ps_model_lasso(dat=dat.all, exposure = "trt", trunc_low = 0, trunc_high = 0.99,
#                         out_vars=out_vars, clust_vars=clust_vars,
#                         baselinevars=c(baselinevars, "prior_gap", vars.empirical),
#                         continuous_vars=c("age_st", "ddose_prev",
#                                           "year", "ret_length",
#                                           "prior_gap"))
# 
# saveRDS(ps_fit, file = paste0(folder, "ps_weight_tr4_rprev.rds"))
# 

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

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ########

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr5.rds"))

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
ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev.rds"))


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

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev.rds"))

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


## Trial 1: 1-2 days #################################
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



ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))


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
ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))

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
ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))

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

saveRDS(resTrial, file=paste0(folder, "resTrial1_12d_iptw.rds"))


## Trial 2: 3-4 days ######################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")


resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident user #######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS/"
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")


ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                covs = c("ret_length", "ddose_prev"), 
                                cluster_var = "studyid_rev")

#resTrial$od$incid = NA


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

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))

ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))

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

ps_fit = readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))

ps_fit$data_weights <- ps_fit$data_weights %>%
  mutate(time_od_acute = ifelse(time_od_acute <= 0, 1e-5, time_od_acute),
         time_discont = ifelse(time_discont <= 0, 1e-5, time_discont))


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


saveRDS(resTrial, file=paste0(folder, "resTrial2_34d_iptw.rds"))




# hdPS with ML analysis  ###########################
## Trial 1: 1-3 days #################################
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"

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
ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr1_rprev.rds"))

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
                          baselinevars=c(baselinevars, "prior_gap"),
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "prior_gap",
                                              "ret_length", "maint_ddose",  
                                              "maint_dur"),
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


## Trial 2: 4 days ##################
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"
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

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ####################
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
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"
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


### Incident users ####################
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
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"
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
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"

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



resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ####################
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



## Trial 1: 1-2 days #################################
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"

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


ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr1_12d_ow.rds"))

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

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr1_12d_ow.rds"))

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
ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr1_12d_prev_ow.rds"))

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

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr1_12d_prev_ow.rds"))

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
ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr1_12d_rprev_ow.rds"))

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
                          baselinevars=c(baselinevars, "prior_gap"),
                          continuous_vars = c("age_st", "ddose_prev",
                                              "year", "prior_gap",
                                              "ret_length", "maint_ddose",  
                                              "maint_dur"),
                          emp_vars=vars.emp)


saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr1_rprev_ow.rds"))

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


saveRDS(resTrial, file = paste0(folder, "resTrial1_12d_iptw.rds"))


## Trial 2: 3-4 days ##################
ps_hdPS = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/"
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/"
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

resTrial=list(od=list(incid=list(), prev=list(), rzero=list()),
              dis=list(incid=list(), prev=list(), rzero=list()))

### Incident users ####################
ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr2_34d_ow.rds"))

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

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr2_34d_ow.rds"))

long_dat = expand_person_period(Dat=ps_fit_ml$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

resTrial$od$incid = fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                               covs = c("age_st", "age_ge50", "ret_length",
                                        "copd_asthma", "sed_psy_med", 
                                        "po", "sud"),  
                               cluster_var = "studyid_rev")



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
ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr2_34d_prev_ow.rds"))

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

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr2_34d_prev_ow.rds"))


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
ps_fit <- readRDS(file = paste0(ps_hdPS, "ps_weight_tr2_34d_rprev_ow.rds"))


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

saveRDS(ps_fit_ml, file = paste0(folder, "ps_weight_tr2_34d_rprev_ow.rds"))

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


saveRDS(resTrial, file = paste0(folder, "resTrial2_34d_iptw.rds"))



