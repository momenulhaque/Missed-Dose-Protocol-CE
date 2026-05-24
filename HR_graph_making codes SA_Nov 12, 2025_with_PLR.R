# Reading results ################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
primary_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
sa_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/"
## Primary analysis ##############
full11 <- readRDS(file=paste0(primary_folder, "res1_12d_iptw.rds"))
full12 <- readRDS(file=paste0(primary_folder, "res2_34d_iptw.rds"))
full13 <- readRDS(file=paste0(primary_folder, "resTrial3_iptw.rds"))
full14 <- readRDS(file=paste0(primary_folder, "resTrial4_iptw.rds"))
full15 <- readRDS(file=paste0(primary_folder, "resTrial5_iptw.rds"))

## Completed induction ##################
ind11 <- readRDS(file=paste0(sa_folder, "Sens_result_induction/resTrial1_12d_iptw.rds"))
ind12 <- readRDS(file=paste0(sa_folder, "Sens_result_induction/resTrial2_34d_iptw.rds"))
ind13 <- readRDS(file=paste0(sa_folder, "Sens_result_induction/resTrial3_iptw.rds"))
ind14 <- readRDS(file=paste0(sa_folder, "Sens_result_induction/resTrial4_iptw.rds"))
ind15 <- readRDS(file=paste0(sa_folder, "Sens_result_induction/resTrial5_iptw.rds"))

## Missed low dose ##################
lowdose11 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_lt80_16/resTrial1_12d_iptw.rds"))
lowdose12 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_lt80_16/resTrial2_34d_iptw.rds"))
lowdose13 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_lt80_16/resTrial3_iptw.rds"))
lowdose14 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_lt80_16/resTrial4_iptw.rds"))
lowdose15 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_lt80_16/resTrial5_iptw.rds"))

## Missed high dose ##################
highdose11 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_ge80_16/resTrial1_12d_iptw.rds"))
highdose12 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_ge80_16/resTrial2_34d_iptw.rds"))
highdose13 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_ge80_16/resTrial3_iptw.rds"))
highdose14 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_ge80_16/resTrial4_iptw.rds"))
highdose15 <- readRDS(file=paste0(sa_folder, "Sens_result_dose_ge80_16/resTrial5_iptw.rds"))

## 2023 Guideline for Trial 1, and 2 ##################
guide23_11 <- readRDS(file=paste0(primary_folder, "resTrial1_iptw.rds"))
guide23_12 <- readRDS(file=paste0(primary_folder, "resTrial2_iptw.rds"))

## Retention length <=30 days ##################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"
ret_le30_11 <- readRDS(file=paste0(folder, "res1_12d_iptw.rds"))
ret_le30_12 <- readRDS(file=paste0(folder, "res2_34d_iptw.rds"))
ret_le30_13 <- readRDS(file=paste0(folder, "res3_iptw.rds"))
ret_le30_14 <- readRDS(file=paste0(folder, "res4_iptw.rds"))
ret_le30_15 <- readRDS(file=paste0(folder, "res5_iptw.rds"))

## Retention length 31-90 days ##################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_30_lt90d/"
ret_gt30_le90_11 <- readRDS(file=paste0(folder, "res1_12d_iptw.rds"))
ret_gt30_le90_12 <- readRDS(file=paste0(folder, "res2_34d_iptw.rds"))
ret_gt30_le90_13 <- readRDS(file=paste0(folder, "res3_iptw.rds"))
ret_gt30_le90_14 <- readRDS(file=paste0(folder, "res4_iptw.rds"))
ret_gt30_le90_15 <- readRDS(file=paste0(folder, "res5_iptw.rds"))


## Retention length >90 days ##################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_gt90d/"
ret_gt90_11 <- readRDS(file=paste0(folder, "res1_12d_iptw.rds"))
ret_gt90_12 <- readRDS(file=paste0(folder, "res2_34d_iptw.rds"))
ret_gt90_13 <- readRDS(file=paste0(folder, "res3_iptw.rds"))
ret_gt90_14 <- readRDS(file=paste0(folder, "res4_iptw.rds"))
ret_gt90_15 <- readRDS(file=paste0(folder, "res5_iptw.rds"))


## > 1 year experience #########
year1_exp_11 <- readRDS(file=paste0(sa_folder, "Sens_retention_365/resTrial1_12d_iptw.rds"))
year1_exp_12 <- readRDS(file=paste0(sa_folder, "Sens_retention_365/resTrial2_34d_iptw.rds"))
year1_exp_13 <- readRDS(file=paste0(sa_folder, "Sens_retention_365/resTrial3_iptw.rds"))
year1_exp_14 <- readRDS(file=paste0(sa_folder, "Sens_retention_365/resTrial4_iptw.rds"))
year1_exp_15 <- readRDS(file=paste0(sa_folder, "Sens_retention_365/resTrial5_iptw.rds"))

## 30 days discontinuation #########
discont_30d_11 <- readRDS(file=paste0(sa_folder, "Sens_discont_30/resTrial1_12d_iptw.rds"))
discont_30d_12 <- readRDS(file=paste0(sa_folder, "Sens_discont_30/resTrial2_34d_iptw.rds"))
discont_30d_13 <- readRDS(file=paste0(sa_folder, "Sens_discont_30/resTrial3_iptw.rds"))
discont_30d_14 <- readRDS(file=paste0(sa_folder, "Sens_discont_30/resTrial4_iptw.rds"))
discont_30d_15 <- readRDS(file=paste0(sa_folder, "Sens_discont_30/resTrial5_iptw.rds"))

## All-acute care visits ###########
all_acute_11 <- readRDS(file=paste0(sa_folder, "Sens_outcome/resTrial1_12d_iptw.rds"))
all_acute_12 <- readRDS(file=paste0(sa_folder, "Sens_outcome/resTrial2_34d_iptw.rds"))
all_acute_13 <- readRDS(file=paste0(sa_folder, "Sens_outcome/resTrial3_iptw.rds"))
all_acute_14 <- readRDS(file=paste0(sa_folder, "Sens_outcome/resTrial4_iptw.rds"))
all_acute_15 <- readRDS(file=paste0(sa_folder, "Sens_outcome/resTrial5_iptw.rds"))

## od-acute care visits in 30d ###########
# folder_sd = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/Sens_outcome_od_acute_30d/"
# od_acute30d_11 <- readRDS(file=paste0(folder_sd, "resTrial1_iptw.rds"))
# od_acute30d_12 <- readRDS(file=paste0(folder_sd, "resTrial2_iptw.rds"))
# od_acute30d_13 <- readRDS(file=paste0(folder_sd, "resTrial3_iptw.rds"))
# od_acute30d_14 <- readRDS(file=paste0(folder_sd, "Sens_outcome_od_acute_30d/resTrial4_iptw.rds"))
# od_acute30d_15 <- readRDS(file=paste0(folder_sd, "Sens_outcome_od_acute_30d/resTrial5_iptw.rds"))



## No truncation at 12m ###########

no_trunc_11 <- readRDS(file=paste0(sa_folder, "Sens_no_truncation/resTrial1_12d_iptw.rds"))
no_trunc_12 <- readRDS(file=paste0(sa_folder, "Sens_no_truncation/resTrial2_34d_iptw.rds"))
no_trunc_13 <- readRDS(file=paste0(sa_folder, "Sens_no_truncation/resTrial3_iptw.rds"))
no_trunc_14 <- readRDS(file=paste0(sa_folder, "Sens_no_truncation/resTrial4_iptw.rds"))
no_trunc_15 <- readRDS(file=paste0(sa_folder, "Sens_no_truncation/resTrial5_iptw.rds"))


## High fentanyl area ################
high_fentanyl_11 <- readRDS(file=paste0(sa_folder, "Sens_high_fent_area/resTrial1_12d_iptw.rds"))
high_fentanyl_12 <- readRDS(file=paste0(sa_folder, "Sens_high_fent_area/resTrial2_34d_iptw.rds"))
high_fentanyl_13 <- readRDS(file=paste0(sa_folder, "Sens_high_fent_area/resTrial3_iptw.rds"))
high_fentanyl_14 <- readRDS(file=paste0(sa_folder, "Sens_high_fent_area/resTrial4_iptw.rds"))
high_fentanyl_15 <- readRDS(file=paste0(sa_folder, "Sens_high_fent_area/resTrial5_iptw.rds"))

## Post-fentanyl era: after April 1, 2012 ###############

post_fentanyl_11 <- readRDS(file=paste0(sa_folder, "post_fentanyl/resTrial1_12d_iptw.rds"))
post_fentanyl_12 <- readRDS(file=paste0(sa_folder, "post_fentanyl/resTrial2_34d_iptw.rds"))
post_fentanyl_13 <- readRDS(file=paste0(sa_folder, "post_fentanyl/resTrial3_iptw.rds"))
post_fentanyl_14 <- readRDS(file=paste0(sa_folder, "post_fentanyl/resTrial4_iptw.rds"))
post_fentanyl_15 <- readRDS(file=paste0(sa_folder, "post_fentanyl/resTrial5_iptw.rds"))


## 30d washout ################
d30_Washout_11 <- readRDS(file=paste0(sa_folder, "Sens_30d_washout/resTrial1_12d_iptw.rds"))
d30_Washout_12 <- readRDS(file=paste0(sa_folder, "Sens_30d_washout/resTrial2_34d_iptw.rds"))
d30_Washout_13 <- readRDS(file=paste0(sa_folder, "Sens_30d_washout/resTrial3_iptw.rds"))
d30_Washout_14 <- readRDS(file=paste0(sa_folder, "Sens_30d_washout/resTrial4_iptw.rds"))
d30_Washout_15 <- readRDS(file=paste0(sa_folder, "Sens_30d_washout/resTrial5_iptw.rds"))

## Comp. outcome: od = od acute + od death ################
## and dis  = discont + od acute 
comp_outcome_11 <- readRDS(file=paste0(sa_folder, "Sens_comp_outcome/resTrial1_12d_iptw.rds"))
comp_outcome_12 <- readRDS(file=paste0(sa_folder, "Sens_comp_outcome/resTrial2_34d_iptw.rds"))
comp_outcome_13 <- readRDS(file=paste0(sa_folder, "Sens_comp_outcome/resTrial3_iptw.rds"))
comp_outcome_14 <- readRDS(file=paste0(sa_folder, "Sens_comp_outcome/resTrial4_iptw.rds"))
comp_outcome_15 <- readRDS(file=paste0(sa_folder, "Sens_comp_outcome/resTrial5_iptw.rds"))


## hdPS  ##################
hdPS_11 <- readRDS(file=paste0(sa_folder, "Sens_hdPS/resTrial1_12d_iptw.rds"))
hdPS_12 <- readRDS(file=paste0(sa_folder, "Sens_hdPS/resTrial2_34d_iptw.rds"))
hdPS_13 <- readRDS(file=paste0(sa_folder, "Sens_hdPS/resTrial3_iptw.rds"))
hdPS_14 <- readRDS(file=paste0(sa_folder, "Sens_hdPS/resTrial4_iptw.rds"))
hdPS_15 <- readRDS(file=paste0(sa_folder, "Sens_hdPS/resTrial5_iptw.rds"))


## hdPS with ML  ##################
hdPS_ML_11 <- readRDS(file=paste0(sa_folder, "Sens_hdPS_ML/resTrial1_12d_iptw.rds"))
hdPS_ML_12 <- readRDS(file=paste0(sa_folder, "Sens_hdPS_ML/resTrial2_34d_iptw.rds"))
hdPS_ML_13 <- readRDS(file=paste0(sa_folder, "Sens_hdPS_ML/resTrial3_iptw.rds"))
hdPS_ML_14 <- readRDS(file=paste0(sa_folder, "Sens_hdPS_ML/resTrial4_iptw.rds"))
hdPS_ML_15 <- readRDS(file=paste0(sa_folder, "Sens_hdPS_ML/resTrial5_iptw.rds"))

## Overlap weights ###################
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
res1_full_ow <- readRDS(file=paste0(folder, "res1_12d_ow.rds"))
res2_full_ow <- readRDS(file=paste0(folder, "res2_34d_ow.rds"))
res3_full_ow <- readRDS(file=paste0(folder, "res3_ow.rds"))
res4_full_ow <- readRDS(file=paste0(folder, "res4_ow.rds"))
res5_full_ow <- readRDS(file=paste0(folder, "res5_ow.rds"))

## IV ############################
IV_1 <- readRDS(file=paste0(sa_folder, "sens_IV/resTrial1_12d_iv.rds"))
IV_2 <- readRDS(file=paste0(sa_folder, "sens_IV/resTrial2_34d_iv.rds"))
IV_3 <- readRDS(file=paste0(sa_folder, "sens_IV/resTrial3_iv.rds"))
IV_4 <- readRDS(file=paste0(sa_folder, "sens_IV/resTrial4_iv.rds"))
IV_5 <- readRDS(file=paste0(sa_folder, "sens_IV/resTrial5_iv.rds"))




extract_all_results <- function(res_obj, trial_label = "Trial 1", model_type = "IPTW") {
  designs <- c("incid" = "Incident users", 
               "prev" = "Prevalent new users", 
               "rzero" = "Prevalent new* users")
  outcomes <- c("od" = "Overdose-related acute care visits", 
                "dis" = "OAT discontinuation")
  
  bind_rows(
    lapply(names(outcomes), function(outcome) {
      lapply(names(designs), function(design) {
        result <- res_obj[[outcome]][[design]]
        
        # Skip if result is NULL or NA
        if (is.null(result) || all(is.na(result))) return(NULL)
        
        data.frame(
          Comparison = result$names,
          HR = result$OR,
          lower = result$lower,
          upper = result$upper,
          Design = designs[[design]],
          Outcome = outcomes[[outcome]],
          Model = model_type,
          Trial = trial_label,
          stringsAsFactors = FALSE
        )
      })
    }) |> unlist(recursive = FALSE)
  )
}


extract_all_results_v1 <- function(res_crude = NULL, res_iptw = NULL, trial_label = "Trial 1") {
  designs <- c("incid" = "Incident users", "prev" = "Prevalent new users", "rzero" = "Prevalent new* users")
  outcomes <- c("od" = "Overdose-related acute care visits", "dis" = "OAT discontinuation")
  
  extract <- function(res_list, model_type) {
    if (is.null(res_list)) return(NULL)
    
    bind_rows(
      lapply(names(outcomes), function(outcome) {
        lapply(names(designs), function(design) {
          result <- res_list[[outcome]][[design]]
          if (is.null(result)) return(NULL)
          
          # Check if result is a list with "model" OR is a coxph directly
          model <- if (inherits(result, "coxph")) {
            result
          } else if (is.list(result) && "model" %in% names(result) && inherits(result$model, "coxph")) {
            result$model
          } else {
            NULL
          }
          
          if (is.null(model)) return(NULL)
          
          coefs <- tryCatch(summary(model)$coefficients, error = function(e) NULL)
          if (is.null(coefs)) return(NULL)
          
          data.frame(
            Comparison = rownames(coefs),
            HR = exp(coefs[, "coef"]),
            lower = exp(coefs[, "coef"] - 1.96 * coefs[, "robust se"]),
            upper = exp(coefs[, "coef"] + 1.96 * coefs[, "robust se"]),
            Design = designs[[design]],
            Outcome = outcomes[[outcome]],
            Model = model_type,
            Trial = trial_label,
            stringsAsFactors = FALSE
          )
        })
      }) |> unlist(recursive = FALSE)
    )
  }
  
  # Only call extract if list is not NULL
  df1 <- if (!is.null(res_crude)) extract(res_crude, "Unadjusted") else NULL
  df2 <- if (!is.null(res_iptw))  extract(res_iptw, "IPTW") else NULL
  
  bind_rows(df1, df2)
}

extract_all_results_IV <- function(res_obj, trial_label = "Trial 1", model_type = "IPTW") {
  designs <- c(
    "incid" = "Incident users", 
    "prev"  = "Prevalent new users", 
    "rzero" = "Prevalent new* users", 
    "rzerp" = "Prevalent new* users"
  )
  
  outcomes <- c(
    "od"  = "Overdose-related acute care visits", 
    "dis" = "OAT discontinuation"
  )
  
  out_list <- list()
  
  for (group_name in names(res_obj)) {
    group <- res_obj[[group_name]]
    
    for (outcome in names(outcomes)) {
      if (!outcome %in% names(group)) next
      outcome_block <- group[[outcome]]
      
      for (design in names(designs)) {
        if (!design %in% names(outcome_block)) next
        result <- outcome_block[[design]]
        
        if (is.null(result) || nrow(result) == 0) next
        
        out_list[[length(out_list) + 1]] <- data.frame(
          Group      = group_name,
          Comparison = result$names,
          HR         = result$OR,
          lower      = result$lower,
          upper      = result$upper,
          Design     = designs[[design]],
          Outcome    = outcomes[[outcome]],
          Model      = model_type,
          Trial      = trial_label,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  do.call(rbind, out_list)
}

## Extracting Results #####################
### Trial 1 ###########
full <- extract_all_results(res_obj=full11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "1. Primary analysis"
  )

ind <- extract_all_results(res_obj=ind11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "2. Completed induction"
  )


low_dose <- extract_all_results(res_obj=lowdose11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "3. Missed <80mg of methadone"
  )


high_dose <- extract_all_results(res_obj=highdose11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "4. Missed >=80mg of methadone"
  )


ret_le30 <- extract_all_results(res_obj=ret_le30_11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "5. Days since OAT initiation: <=30 days"
  )

ret_gt30_le90 <- extract_all_results(res_obj=ret_gt30_le90_11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "6. Days since OAT initiation: 31-90 days"
  )

ret_gt90 <- extract_all_results(res_obj=ret_gt90_11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "7. Days since OAT initiation: >90 days"
  )

guide23 <- extract_all_results(res_obj=guide23_11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "8. 2023 Guidelines: 1-3 days gap"
  )

high_fentanyl <- extract_all_results(res_obj=high_fentanyl_11, 
                                     trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "9. Restricted to high fentanyl area"
  )


post_fentanyl <- extract_all_results(res_obj=post_fentanyl_11, 
                                     trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)"
  )



year1_exp <- extract_all_results(res_obj=year1_exp_11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "11. 1-year or longer cumulative OAT experience"
  )

no_trunc <- extract_all_results(res_obj=no_trunc_11, 
                                trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "12. No truncation at 12 months of follow-up"
  )

d30_washout <- extract_all_results(res_obj=d30_Washout_11, 
                                   trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "13. 30 days washout for prevalent new users"
  )



discont_30d <- extract_all_results(res_obj=discont_30d_11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "14. Outcome: No OAT dispensations for 30 days"
  )


all_acute <- extract_all_results(res_obj=all_acute_11, 
                                 trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease"),
         Outcome == "Overdose-related acute care visits") %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "15. Outcome: All-cause acute care visits"
  )



comp_outcome <- extract_all_results(res_obj=comp_outcome_11, 
                                    trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = ifelse(Outcome == "Overdose-related acute care visits", 
                      "16. Outcome: Overdose-related acute care \nvisits or overdose-related death",
                      "17. Outcome: Overdose-related acute \ncare visits or OAT discontinuation")
  )

hdPS <- extract_all_results(res_obj = hdPS_11, trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "18. hdPS analysis"
  )

hdPS_ML <- extract_all_results(res_obj=hdPS_ML_11, 
                               trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "19. hdPS with machine learning algorithm"
  )

OW <- extract_all_results(res_obj=res1_full_ow, 
                               trial_label = "Trial 1") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "20. Overlap weights based analysis"
  )


IV <- extract_all_results_IV(res_obj=IV_1, 
                             trial_label = "Trial 1") %>%
  select(-Comparison) %>% 
  filter(Group %in% c("grp_red", "grp_incr")) %>%
  mutate(
    Comparison = recode(Group,
                        "grp_red" = "Dose reduction\nvs\nNo dose change",
                        "grp_incr"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "21. Prescriber's preference based IV estimation"
  ) %>% select(Comparison, HR, lower, upper, Design, Outcome, Model, Trial, SA_group)


res_com <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, guide23, high_fentanyl,
  post_fentanyl, year1_exp, no_trunc, d30_washout, all_acute, discont_30d, comp_outcome, hdPS,
  hdPS_ML, OW, IV) %>%
  mutate(
    SA_group = factor(SA_group), 
    Trial = recode(Trial, "Trial 1" = "Trial 1: Methadone interruption of 1-2 days")
  ) %>% 
  mutate(
    SA_group = factor(
      SA_group,
      levels = unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  ) %>%
  arrange(Design)

saveRDS(res_com, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral1_12d_PLR.rds")


### Trial 2 ###########
full <- extract_all_results(res_obj=full12, trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "1. Primary analysis"
  )

ind <- extract_all_results(res_obj=ind12, trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "2. Completed induction"
  )



low_dose <- extract_all_results(res_obj=lowdose12, trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "3. Missed <80mg of methadone"
  )

high_dose <- extract_all_results(res_obj=highdose12, trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "4. Missed >=80mg of methadone"
  )




ret_le30 <- extract_all_results(res_obj=ret_le30_12, 
                                trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "5. Days since OAT initiation: <=30 days"
  )

ret_gt30_le90 <- extract_all_results(res_obj=ret_gt30_le90_12, 
                                     trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "6. Days since OAT initiation: 31-90 days"
  )

ret_gt90 <- extract_all_results(res_obj=ret_gt90_12, 
                                trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "7. Days since OAT initiation: >90 days"
  )

guide23 <- extract_all_results(res_obj=guide23_12, trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "8. 2023 Guidelines: 4 days gap"
  )


high_fentanyl <- extract_all_results(res_obj=high_fentanyl_12, 
                                     trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "9. Restricted to high fentanyl area"
  )

post_fentanyl <- extract_all_results(res_obj=post_fentanyl_12, 
                                     trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)"
  )

year1_exp <- extract_all_results(res_obj=year1_exp_12, 
                                 trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "11. 1-year or longer cumulative OAT experience"
  )

no_trunc <- extract_all_results(res_obj=no_trunc_12, 
                                trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "12. No truncation at 12 months of follow-up"
  )


d30_washout <- extract_all_results(res_obj=d30_Washout_12, 
                                   trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "13. 30 days washout for prevalent new users"
  )


discont_30d <- extract_all_results(res_obj=discont_30d_12, 
                                   trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "14. Outcome: No OAT dispensations for 30 days"
  )

all_acute <- extract_all_results(res_obj=all_acute_12, 
                                 trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease"),
         Outcome == "Overdose-related acute care visits") %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "15. Outcome: All-cause acute care visits"
  )


comp_outcome <- extract_all_results(res_obj=comp_outcome_12, 
                                    trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = ifelse(Outcome == "Overdose-related acute care visits", 
                      "16. Outcome: Overdose-related acute care \nvisits or overdose-related death",
                      "17. Outcome: Overdose-related acute \ncare visits or OAT discontinuation")
  )


hdPS <- extract_all_results(res_obj=hdPS_12, trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "18. hdPS analysis"
  )

hdPS_ML <- extract_all_results(res_obj=hdPS_ML_12, 
                               trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "19. hdPS with machine learning algorithm"
  )

OW <- extract_all_results(res_obj=res2_full_ow, 
                          trial_label = "Trial 2") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "20. Overlap weights based analysis"
  )


IV <- extract_all_results_IV(res_obj=IV_2, 
                             trial_label = "Trial 2") %>%
  select(-Comparison) %>% 
  filter(Group %in% c("grp_lt50", "grp_ge50", "grp_incr")) %>%
  mutate(
    Comparison = recode(Group,
                        "grp_lt50" = "1-<50% reduction\nvs\nNo dose change",
                        "grp_ge50" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "grp_incr" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "21. Prescriber's preference based IV estimation"
  ) %>% select(Comparison, HR, lower, upper, Design, Outcome, Model, Trial, SA_group)






res_com <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, guide23, high_fentanyl,
  post_fentanyl, year1_exp, no_trunc, d30_washout, all_acute, discont_30d, comp_outcome, hdPS, hdPS_ML, OW, IV
) %>%
  mutate(
    SA_group = factor(SA_group), 
    Trial = recode(Trial, "Trial 2" = "Trial 2: Methadone interruption of 3-4 days")
  ) %>% 
  mutate(
    SA_group = factor(
      SA_group,
      levels = unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  ) %>%
  arrange(Design)

saveRDS(res_com, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral2_34d_PLR.rds")



### Trial 3 ###########
full <- extract_all_results(res_obj=full13, trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "1. Primary analysis"
  )

ind <- extract_all_results(res_obj=ind13, trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "2. Completed induction"
  )



low_dose <- extract_all_results(res_obj=lowdose13, trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "3. Missed <80mg of methadone"
  )

high_dose <- extract_all_results(res_obj=highdose13, trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "4. Missed >=80mg of methadone"
  )


ret_le30 <- extract_all_results(res_obj=ret_le30_13, trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "5. Days since OAT initiation: <=30 days"
  )

ret_gt30_le90 <- extract_all_results(res_obj= 
                                       ret_gt30_le90_13, trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "6. Days since OAT initiation: 31-90 days"
  )

ret_gt90 <- extract_all_results(res_obj=ret_gt90_13, trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "7. Days since OAT initiation: >90 days"
  )

high_fentanyl <- extract_all_results(res_obj=high_fentanyl_13, 
                                     trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "9. Restricted to high fentanyl area"
  )

post_fentanyl <- extract_all_results(res_obj=post_fentanyl_13, 
                                     trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)"
  )


year1_exp <- extract_all_results(res_obj=year1_exp_13, 
                                 trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "11. 1-year or longer cumulative OAT experience"
  )

no_trunc <- extract_all_results(res_obj=no_trunc_13, 
                                trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "12. No truncation at 12 months of follow-up"
  )

d30_washout <- extract_all_results(res_obj=d30_Washout_13, 
                                   trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "13. 30 days washout for prevalent new users"
  )


discont_30d <- extract_all_results(res_obj=discont_30d_13, 
                                   trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "14. Outcome: No OAT dispensations for 30 days"
  )

all_acute <- extract_all_results(res_obj=all_acute_13, 
                                 trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease"),
         Outcome == "Overdose-related acute care visits") %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "15. Outcome: All-cause acute care visits"
  )



comp_outcome <- extract_all_results(res_obj=comp_outcome_13, 
                                    trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = ifelse(Outcome == "Overdose-related acute care visits", 
                      "16. Outcome: Overdose-related acute care \nvisits and overdose-related death",
                      "17. Outcome: Overdose-related acute \ncare visits and OAT discontinuation")
  )





hdPS <- extract_all_results(res_obj=hdPS_13, trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "18. hdPS analysis"
  )



hdPS_ML <- extract_all_results(res_obj=hdPS_ML_13, trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "19. hdPS with machine learning algorithm"
  )

OW <- extract_all_results(res_obj=res3_full_ow, 
                          trial_label = "Trial 3") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "20. Overlap weights based analysis"
  )


IV <- extract_all_results_IV(res_obj=IV_3, 
                             trial_label = "Trial 3") %>%
  select(-Comparison) %>% 
  filter(Group %in% c("grp_lt50", "grp_ge50", "grp_incr")) %>%
  mutate(
    Comparison = recode(Group,
                        "grp_lt50" = "1-<50% reduction\nvs\nNo dose change",
                        "grp_ge50" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "grp_incr" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "21. Prescriber's preference based IV estimation"
  ) %>% select(Comparison, HR, lower, upper, Design, Outcome, Model, Trial, SA_group)





res_com <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, high_fentanyl,
  post_fentanyl, year1_exp, no_trunc, d30_washout, all_acute, discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>%
  mutate(
    SA_group = factor(SA_group), 
    Trial = recode(Trial, "Trial 3" = "Trial 3: Methadone interruption of 5-14 days")
  ) %>% 
  mutate(
    SA_group = factor(
      SA_group,
      levels = unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  ) %>%
  arrange(Design)

#View(res_com)
saveRDS(res_com, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral3_PLR.rds")


### Trial 4 ###########

full <- extract_all_results(res_obj=full14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "1. Primary analysis"
  )

ind <- extract_all_results(res_obj=ind14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "2. Completed induction"
  )


low_dose <- extract_all_results(res_obj=lowdose14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "3. Missed <16/4mg of buprenorphine/naloxone"
  )


high_dose <- extract_all_results(res_obj=highdose14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "4. Missed >=16/4mg of buprenorphine/naloxone"
  )

ret_le30 <- extract_all_results(res_obj=ret_le30_14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "5. Days since OAT initiation: <=30 days"
  )

ret_gt30_le90 <- extract_all_results(res_obj=ret_gt30_le90_14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "6. Days since OAT initiation: 31-90 days"
  )

ret_gt90 <- extract_all_results(res_obj=ret_gt90_14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "7. Days since OAT initiation: >90 days"
  )



high_fentanyl <- extract_all_results(res_obj=high_fentanyl_14, 
                                     trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "9. Restricted to high fentanyl area"
  )

post_fentanyl <- extract_all_results(res_obj=post_fentanyl_14, 
                                     trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "10. Restricted to post fentanyl era \n(from April 1, 2012 - December 31, 2022)"
  )



year1_exp <- extract_all_results(res_obj=year1_exp_14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "11. 1-year or longer cumulative OAT experience"
  )

no_trunc <- extract_all_results(res_obj=no_trunc_14, 
                                trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "12. No truncation at 12 months of follow-up"
  )



d30_washout <- extract_all_results(res_obj=d30_Washout_14, 
                                   trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "13. 30 days washout for prevalent new users"
  )



discont_30d <- extract_all_results(res_obj=discont_30d_14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "14. Outcome: No OAT dispensations for 30 days"
  )


all_acute <- extract_all_results(res_obj=all_acute_14, 
                                 trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease"),
         Outcome == "Overdose-related acute care visits") %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "15. Outcome: All-cause acute care visits"
  )

comp_outcome <- extract_all_results(res_obj=comp_outcome_14, 
                                    trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = ifelse(Outcome == "Overdose-related acute care visits", 
                      "16. Outcome: Overdose-related acute care \nvisits and overdose-related death",
                      "17. Outcome: Overdose-related acute \ncare visits and OAT discontinuation")
  )

hdPS <- extract_all_results(res_obj=hdPS_14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "18. hdPS analysis"
  )

hdPS_ML <- extract_all_results(res_obj=hdPS_ML_14, trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "19. hdPS with machine learning algorithm"
  )


OW <- extract_all_results(res_obj=res4_full_ow, 
                          trial_label = "Trial 4") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "20. Overlap weights based analysis"
  )


IV <- extract_all_results_IV(res_obj=IV_4, 
                             trial_label = "Trial 4") %>%
  select(-Comparison) %>% 
  filter(Group %in% c("grp_red", "grp_incr")) %>%
  mutate(
    Comparison = recode(Group,
                        "grp_red" = "Dose reduction\nvs\nNo dose change",
                        "grp_incr"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )), 
    SA_group = "21. Prescriber's preference based IV estimation"
  ) %>% select(Comparison, HR, lower, upper, Design, Outcome, Model, Trial, SA_group)




res_com <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, high_fentanyl,
  post_fentanyl, year1_exp, no_trunc, d30_washout, all_acute, discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>%
  mutate(
    SA_group = factor(SA_group), 
    Trial = recode(Trial, "Trial 4" = "Trial 4: Buprenorphine/naloxone interruption of 1-5 days")
  ) %>%  
  mutate(
    SA_group = factor(
      SA_group,
      levels = unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  ) %>%
  arrange(Design)


saveRDS(res_com, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral4_PLR.rds")


### Trial 5 ###########
full <- extract_all_results(res_obj=full15, trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "1. Primary analysis"
  )

ind <- extract_all_results(res_obj=ind15, trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "2. Completed induction"
  )

low_dose <- extract_all_results(res_obj=lowdose15, trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "3. Missed <16/4mg of buprenorphine/naloxone"
  )

high_dose <- extract_all_results(res_obj=highdose15, trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "4. Missed >=16/4mg of buprenorphine/naloxone"
  )


ret_le30 <- extract_all_results(res_obj=ret_le30_15, trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "5. Days since OAT initiation: <=30 days"
  )

ret_gt30_le90 <- extract_all_results(res_obj= 
                                       ret_gt30_le90_15, 
                                     trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "6. Days since OAT initiation: 31-90 days"
  )

ret_gt90 <- extract_all_results(res_obj=ret_gt90_15, trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "7. Days since OAT initiation: >90 days"
  )


high_fentanyl <- extract_all_results(res_obj=high_fentanyl_15, 
                                     trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "9. Restricted to high fentanyl area"
  )

post_fentanyl <- extract_all_results(res_obj=post_fentanyl_15, 
                                     trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)"
  )

year1_exp <- extract_all_results(res_obj=year1_exp_15, 
                                 trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "11. 1-year or longer cumulative OAT experience"
  )

no_trunc <- extract_all_results(res_obj=no_trunc_15, 
                                trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "12. No truncation at 12 months of follow-up"
  )

d30_washout <- extract_all_results(res_obj=d30_Washout_15, 
                                   trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "13. 30 days washout for prevalent new users"
  )




discont_30d <- extract_all_results(res_obj=discont_30d_15, 
                                   trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "14. Outcome: No OAT dispensations for 30 days"
  )

all_acute <- extract_all_results(res_obj=all_acute_15, 
                                 trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease"),
         Outcome == "Overdose-related acute care visits") %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "15. Outcome: All-cause acute care visits"
  )




comp_outcome <- extract_all_results(res_obj=comp_outcome_15, 
                                    trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = ifelse(Outcome == "Overdose-related acute care visits", 
                      "16. Outcome: Overdose-related acute care \nvisits and overdose-related death",
                      "17. Outcome: Overdose-related acute \ncare visits and OAT discontinuation")
    
  )



hdPS <- extract_all_results(res_obj=hdPS_15, 
                            trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "18. hdPS analysis"
  )

hdPS_ML <- extract_all_results(res_obj=hdPS_ML_15, 
                               trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "19. hdPS with machine learning algorithm"
  )

OW <- extract_all_results(res_obj=res5_full_ow, 
                          trial_label = "Trial 5") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "20. Overlap weights based analysis"
  )


IV <- extract_all_results_IV(res_obj=IV_5, 
                             trial_label = "Trial 5") %>%
  select(-Comparison) %>% 
  filter(Group %in% c("grp_lt50", "grp_ge50", "grp_incr")) %>%
  mutate(
    Comparison = recode(Group,
                        "grp_lt50" = "1-<50% reduction\nvs\nNo dose change",
                        "grp_ge50" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "grp_incr" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    SA_group = "21. Prescriber's preference based IV estimation"
  ) %>% select(Comparison, HR, lower, upper, Design, Outcome, Model, Trial, SA_group)



res_com <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, high_fentanyl,
  post_fentanyl, year1_exp, no_trunc, d30_washout, all_acute, discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>%
  mutate(
    SA_group = factor(SA_group), 
    Trial = recode(Trial, "Trial 5" = "Trial 5: Buprenorphine/naloxone interruption of 6-14 days")
  ) %>%  
  mutate(
    SA_group = factor(
      SA_group,
      levels = unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  ) %>%
  arrange(Design)




saveRDS(res_com, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral5_PLR.rds")


# Create graphs: Nov 12, 2025 ##############
res_com1 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral1_12d_PLR.rds")
res_com2 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral2_34d_PLR.rds")
res_com3 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral3_PLR.rds")
res_com4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral4_PLR.rds")
res_com5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral5_PLR.rds")


res_od_dis <- bind_rows(res_com1, res_com2, res_com3, res_com4, res_com5)

res_od_dis_df <- res_od_dis %>%
  mutate(
    Design = recode(Design, 
                    "Prevalent new* users" = "All prevalent users"),
    Design = factor(Design,
                    levels = c("Incident users",
                               "Prevalent new users",
                               "All prevalent users")),
    SA_group = factor(SA_group),
    Trial = 
      recode(Trial,
             "Trial 1: Methadone interruption of 1-2 days"="Trial 1: MET interruption of 1-2 days",
             "Trial 2: Methadone interruption of 3-4 days"="Trial 2: MET interruption of 3-4 days",
             "Trial 3: Methadone interruption of 5-14 days"="Trial 3: MET interruption of 5-14 days",
             "Trial 4: Buprenorphine/naloxone interruption of 1-5 days"="Trial 4: BNX interruption of 1-5 days",
             "Trial 5: Buprenorphine/naloxone interruption of 6-14 days"="Trial 5: BNX interruption of 6-14 days"),
    Comparison =
      recode(Comparison,                        
             "Dose reduction\nvs\nNo dose change" = "Dose reduction\nvs\nno change",
             "Dose increase\nvs\nNo dose change" = "Dose increase\nvs\nno change",
             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50%\nvs\nno change",
             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose*\nvs\nno change",
             ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose*\nvs\nno change"),
    Comparison = factor(Comparison, levels = rev(c("Dose reduction\nvs\nno change",
                                               "Dose reduction to <50%\nvs\nno change",
                                               "Dose reduction to >=50%/initial dose*\nvs\nno change",
                                               "Dose increase\nvs\nno change")))
    
  )


res_od_dis_df <- res_od_dis_df %>% 
  mutate(
    SA_group = factor(
      SA_group,
      levels = unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    ),
    SA_group = recode(SA_group,
                      "3. Missed <80mg of methadone" = "3. Missed at <80mg of methadone",
                      "3. Missed <16/4mg of buprenorphine/naloxone" = "3. Missed at <16/4mg of buprenorphine/naloxone",
                      "4. Missed >=80mg of methadone" = "4. Missed at >=80mg of methadone",
                      "4. Missed >=16/4mg of buprenorphine/naloxone" = "4. Missed at >=16/4mg of buprenorphine/naloxone",
                      "18. hdPS analysis" = "17. hdPS analysis",
                      "19. hdPS with machine learning algorithm" = "18. hdPS with machine learning algorithm",
                      "20. Overlap weights based analysis" = "19. Overlap weights based analysis",
                      "21. Prescriber's preference based IV estimation" = 
                        "20. Prescriber's preference based IV estimation")
    
  ) %>% filter(SA_group != "17. Outcome: Overdose-related acute \ncare visits or OAT discontinuation") %>% 
  filter(! (SA_group == "13. 30 days washout for prevalent new users" &
              Design == "Incident users"))

levels(res_od_dis_df$Comparison)

res_od_dis_df %>% filter(SA_group == "1. Primary analysis",
                         Trial == "Trial 2: MET interruption of 3-4 days",
                         Outcome == "OAT discontinuation")


comparison_colors_v1 = c(
  "Dose reduction\nvs\nno change"              = "#000000",  # Blue
  "Dose reduction to <50%\nvs\nno change"            = "#2CA02C",  # Green
  "Dose reduction to >=50%/initial dose*\nvs\nno change"     = "#FF7F0E",  # Sky Blue
  "Dose increase\nvs\nno change"               = "#D62728"  # Orange
)



scale_plot_trial_results_od_dis_v1 <- function(df, 
                                               trial_label = "Trial 1",
                                               comparison_colors_v1 = c(
                                                 "Dose reduction\nvs\nno change"              = "#000000",  # Blue
                                                 "Dose reduction to <50%\nvs\nno change"            = "#2CA02C",  # Green
                                                 "Dose reduction to >=50%/initial dose*\nvs\nno change"     = "#FF7F0E",  # Sky Blue
                                                 "Dose increase\nvs\nno change"               = "#D62728"  # Orange
                                               ),
                                               y_min=0.1,
                                               y_max=5,
                                               y_points = c(0.3, 0.5, 1, 1.5, 3)) {
  library(ggplot2)
  library(ggh4x)
  library(forcats)
  library(dplyr)
  library(scales)
  
  # Standardize levels
  df <- df %>%
    mutate(
      Comparison = factor(Comparison, levels = names(comparison_colors_v1))
    )
  
  sa_levels <- levels(df$SA_group)
  df <- df %>% mutate(SA_group_num = as.numeric(SA_group))
  
  bg_data <- data.frame(
    SA_group = sa_levels,
    group_index = seq_along(sa_levels),
    fill_color = rep(c("gray85", "white"), length.out = length(sa_levels))
  )
  
  dodge <- position_dodge(width = 1)
  
  ggplot(df, aes(x = SA_group_num, y = HR, color = Comparison)) +
    geom_rect(data = bg_data,
              aes(xmin = group_index - 0.5, xmax = group_index + 0.5,
                  ymin = y_min, ymax = y_max, fill = fill_color),  # manually limit y-range
              inherit.aes = FALSE, alpha = 0.5) +
    geom_point(position = dodge, size = 3, na.rm = TRUE) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 0.7, position = dodge, na.rm = TRUE) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    ggh4x::facet_grid2(~ Outcome, scales = "free",
                       strip = ggh4x::strip_themed(
                         background_x = elem_list_rect(fill = "lightgray", color = "gray50"),
                         text_x = elem_list_text(face = "bold", size = 14)
                       )) +
    coord_flip() +
    scale_x_continuous(breaks = seq_along(sa_levels), labels = sa_levels, expand = c(0, 0)) +
    #    scale_y_continuous(breaks = y_points, limits = c(y_min, y_max))+
    scale_y_log10(breaks = y_points, limits = c(y_min, y_max))+
    scale_color_manual(values = comparison_colors_v1) +
    scale_fill_identity() +
    labs(
      title = trial_label,
      x = NULL,
      y = "Hazard ratios with 95% compatibility intervals"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.text = element_text(color = "black", size = 15),
      axis.text.y = element_text(hjust = 0, color = "black", size = 15),
      axis.title.x = element_text(color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 15, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 15)
    )
  
}

scale_plot_trial_results_od_v1 <- function(df, 
                                           trial_label = "Trial 1",
                                           comparison_colors_v1 = c(
                                             "Dose reduction\nvs\nno change"              = "#000000",  # Blue
                                             "Dose reduction to <50%\nvs\nno change"            = "#2CA02C",  # Green
                                             "Dose reduction to >=50%/initial dose*\nvs\nno change"     = "#FF7F0E",  # Sky Blue
                                             "Dose increase\nvs\nno change"               = "#D62728"  # Orange
                                           ),
                                           y_min=0.2,
                                           y_max=5,
                                           y_points = c(0.3, 0.5, 1, 1.5, 3)) {
  library(ggplot2)
  library(ggh4x)
  library(forcats)
  library(dplyr)
  library(scales)
  
  # Standardize levels
  df <- df %>%
    mutate(
      Comparison = factor(Comparison, levels = names(comparison_colors_v1))
    )
  
  sa_levels <- levels(df$SA_group)
  
  # Create empty data frame with SA_group and Trial for facets and x-axis labels
  empty_df <- expand.grid(
    SA_group = factor(sa_levels, levels = sa_levels),
    Outcome = unique(df$Outcome)
  )
  empty_df$SA_group_num <- as.numeric(empty_df$SA_group)
  
  p <- ggplot(empty_df, aes(x = SA_group_num)) +
    ggh4x::facet_grid2(~ Outcome, scales = "free",
                       strip = ggh4x::strip_themed(
                         background_x = elem_list_rect(fill = "lightgray", color = "gray50"),
                         text_x = elem_list_text(face = "bold", size = 14)
                       )) +
    coord_flip() +
    scale_y_log10(breaks = y_points, limits = c(y_min, y_max))+
    scale_x_continuous(breaks = seq_along(sa_levels), labels = sa_levels, expand = c(0, 0)) +
    #    scale_y_continuous(breaks = y_points, limits = c(y_min, y_max))+
    labs(
      title = trial_label,
      x = NULL,
      y = "Hazard ratios with 95% compatibility intervals"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.text = element_text(color = "black", size = 15),
      axis.title.x = element_text(color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 15, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 15)
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40")
  
  return(p)
  
}


#### od: Incident users #############
od_incid_tr1 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group), 
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))


plot_od_incid_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=od_incid_tr1, 
                                                          trial_label = "",
                                                          y_min=0.1)


od_incid_tr2 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))



plot_od_incid_tr2_v1 = scale_plot_trial_results_od_v1(df=od_incid_tr2, 
                                                      trial_label = " ")
#View(od_incid_tr3)
od_incid_tr3 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))



plot_od_incid_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=od_incid_tr3, 
                                                          trial_label = " ")


od_incid_tr4 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_od_incid_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=od_incid_tr4, 
                                                          trial_label = " ")


od_incid_tr5 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  )  %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))



plot_od_incid_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=od_incid_tr5, 
                                                          trial_label = " ")



#### od: prevalent new users #############
od_prev_tr1 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                       Design == "Prevalent new users",
                                       Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_od_prev_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=od_prev_tr1, 
                                                         trial_label = " ")


od_prev_tr2 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                       Design == "Prevalent new users",
                                       Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_od_prev_tr2_v1 = scale_plot_trial_results_od_dis_v1(df=od_prev_tr2, 
                                                         trial_label = " ")



od_prev_tr3 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                       Design == "Prevalent new users",
                                       Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == ">=50% reduction/initial dose*\nvs\nNo dose change"))%>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))




plot_od_prev_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=od_prev_tr3, 
                                                         trial_label = " ")




od_prev_tr4 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                       Design == "Prevalent new users",
                                       Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_od_prev_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=od_prev_tr4, 
                                                         trial_label = " ")



od_prev_tr5 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                       Design == "Prevalent new users",
                                       Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_od_prev_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=od_prev_tr5, 
                                                         trial_label = " ")


#### od: All prevalent users #############
od_all_prev_tr1 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_od_all_prev_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr1, 
                                                             trial_label = " ")



od_all_prev_tr2 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))



plot_od_all_prev_tr2_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr2, 
                                                             trial_label = " ")


od_all_prev_tr3 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))



plot_od_all_prev_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr3, 
                                                             trial_label = " ")

od_all_prev_tr4 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))


plot_od_all_prev_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr4, 
                                                             trial_label = " ")


od_all_prev_tr5 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change"))  %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))



plot_od_all_prev_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr5, 
                                                             trial_label = " ")

#### dis: Incident users #############
dis_incid_tr1 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_dis_incid_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr1, 
                                                           trial_label = " ")



dis_incid_tr2 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))



plot_dis_incid_tr2_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr2, 
                                                           trial_label = " ")


dis_incid_tr3 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))



plot_dis_incid_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr3, 
                                                           trial_label = " ")



dis_incid_tr4 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_dis_incid_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr4, 
                                                           trial_label = "")


dis_incid_tr5 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))


plot_dis_incid_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr5, 
                                                           trial_label = "") 


#### dis: prevalent new users #############
dis_prev_tr1 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                        Design == "Prevalent new users",
                                        Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))
plot_dis_prev_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=dis_prev_tr1, 
                                                          trial_label = "")


dis_prev_tr2 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                        Design == "Prevalent new users",
                                        Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    ) 
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))


plot_dis_prev_tr2_v1 = scale_plot_trial_results_od_dis_v1(df=dis_prev_tr2, 
                                                          trial_label = "")


dis_prev_tr3 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                        Design == "Prevalent new users",
                                        Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))



plot_dis_prev_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=dis_prev_tr3, 
                                                          trial_label = "")


dis_prev_tr4 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                        Design == "Prevalent new users",
                                        Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    ) 
  )  %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))
plot_dis_prev_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=dis_prev_tr4, 
                                                          trial_label = "")



dis_prev_tr5 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                        Design == "Prevalent new users",
                                        Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))


plot_dis_prev_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=dis_prev_tr5, 
                                                          trial_label = "")



#### dis: All prevalent users #############
dis_all_prev_tr1 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    ) 
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))
plot_dis_all_prev_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr1, 
                                                              trial_label = "")


dis_all_prev_tr2 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  filter(HR < 3,
         upper < 5) %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))


plot_dis_all_prev_tr2_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr2, 
                                                              trial_label = "")


dis_all_prev_tr3 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(HR < 3,
         upper < 5,
         lower > 0.02) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))


plot_dis_all_prev_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr3, 
                                                              trial_label = "")



dis_all_prev_tr4 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  filter(HR < 3,
         upper < 5) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_dis_all_prev_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr4, 
                                                              trial_label = "")




dis_all_prev_tr5 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  filter(HR < 3,
         upper < 5) %>%  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reduction to >=50%/initial dose*\nvs\nno change")) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))


plot_dis_all_prev_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr5, 
                                                              trial_label = "")




## Combining graphs ##########
### Trial 1 #################
library(cowplot)
library(patchwork)
row1 <- plot_od_incid_tr1_v1 + plot_spacer() + plot_dis_incid_tr1_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_od_prev_tr1_v1 + plot_spacer() + plot_dis_prev_tr1_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_od_all_prev_tr1_v1 + plot_spacer() + plot_dis_all_prev_tr1_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr1_v1 <- plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() / # spacer between row1 and row2
  row2 / 
  plot_spacer() / # spacer between row2 and row3
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od_dis_tr1_v1s <- ggdraw(plot_od_dis_tr1_v1) + 
  draw_label("A. Methadone interruption of 1-2 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 1-2 days among prevalent new users", 
             x = 0.02, y = 0.66, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Methadone interruption of 1-2 days among all prevalent users", 
             x = 0.02, y = 0.34, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("Note: Comparison for some sensitivity analysis is NA for insufficient number of events or\nIV assumption was not satisfied.", 
             x = 0.63, y = 0.02, hjust = 0, vjust = 0, 
             size = 13, color = "black")


name = "plot_od_dis_tr1_Nov_12"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7800, res=300)
plot_od_dis_tr1_v1s
dev.off()




### Trial 2 #################
library(cowplot)
library(patchwork)
row1 <- plot_od_incid_tr2_v1 + plot_spacer() + plot_dis_incid_tr2_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_od_prev_tr2_v1 + plot_spacer() + plot_dis_prev_tr2_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_od_all_prev_tr2_v1 + plot_spacer() + plot_dis_all_prev_tr2_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr2_v1 <-plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() / # spacer between row1 and row2
  row2 / 
  plot_spacer() / # spacer between row2 and row3
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od_dis_tr2_v1s <- ggdraw(plot_od_dis_tr2_v1) + 
  draw_label("* <=40mg methadone", 
             x = 0.68, y = 0.026, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("Note: Comparison for some sensitivity analysis is NA for insufficient number\n   of events or IV assumption was not satisfied.", 
             x = 0.68, y = 0.006, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("NA for few events", 
             x = 0.38, y = 0.83, hjust = 1, vjust = 0, 
             size = 13, color = "black") +
  draw_label("A. Methadone interruption of 3-4 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 3-4 days among prevalent new users", 
             x = 0.02, y = 0.66, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Methadone interruption of 3-4 days among all prevalent users", 
             x = 0.02, y = 0.34, hjust = 0, vjust = 0, 
             size = 18, color = "black")



name = "plot_od_dis_tr2_Nov_12"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=8500, res=300)
plot_od_dis_tr2_v1s
dev.off()


### Trial 3 #################
library(cowplot)
library(patchwork)
row1 <- plot_od_incid_tr3_v1 + plot_spacer() + plot_dis_incid_tr3_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_od_prev_tr3_v1 + plot_spacer() + plot_dis_prev_tr3_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_od_all_prev_tr3_v1 + plot_spacer() + plot_dis_all_prev_tr3_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr3_v1 <- plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() / # spacer between row1 and row2
  row2 / 
  plot_spacer() / # spacer between row2 and row3
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))


plot_od_dis_tr3_v1s <- ggdraw(plot_od_dis_tr3_v1) + 
  draw_label("* <=40mg methadone", 
             x = 0.68, y = 0.026, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("Note: Comparison for some sensitivity analysis is NA for insufficient number\n   of events or IV assumption was not satisfied.", 
             x = 0.68, y = 0.006, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("A. Methadone interruption of 5-14 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 5-14 days among prevalent new users", 
             x = 0.02, y = 0.66, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Methadone interruption of 5-14 days among all prevalent users", 
             x = 0.02, y = 0.34, hjust = 0, vjust = 0, 
             size = 18, color = "black")


name = "plot_od_dis_tr3_Nov_12"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=8500, res=300)
plot_od_dis_tr3_v1s
dev.off()



### Trial 4 #################
library(cowplot)
library(patchwork)

row1 <- plot_od_incid_tr4_v1 + plot_spacer() + plot_dis_incid_tr4_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_od_prev_tr4_v1 + plot_spacer() + plot_dis_prev_tr4_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_od_all_prev_tr4_v1 + plot_spacer() + plot_dis_all_prev_tr4_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr4_v1 <- plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() /# spacer between row1 and row2
  row2 / 
  plot_spacer() /# spacer between row2 and row3
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od_dis_tr4_v1s <- ggdraw(plot_od_dis_tr4_v1) + 
  draw_label("A. Buprenorphine/naloxone interruption of 1-5 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone interruption of 1-5 days among prevalent new users", 
             x = 0.02, y = 0.66, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Buprenorphine/naloxone interruption of 1-5 days among all prevalent users", 
             x = 0.02, y = 0.34, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("Note: Comparison for some sensitivity analysis is NA for insufficient number of events or\nIV assumption was not satisfied.", 
             x = 0.63, y = 0.006, hjust = 0, vjust = 0, 
             size = 13, color = "black")


name = "plot_od_dis_tr4_Nov_12"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7800, res=300)
plot_od_dis_tr4_v1s
dev.off()


### Trial 5 #################
library(cowplot)
library(patchwork)
row1 <- plot_od_incid_tr5_v1 + plot_spacer() + plot_dis_incid_tr5_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_od_prev_tr5_v1 + plot_spacer() + plot_dis_prev_tr5_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_od_all_prev_tr5_v1 + plot_spacer() + plot_dis_all_prev_tr5_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr5_v1 <- plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() / # spacer between row1 and row2
  row2 / 
  plot_spacer() / # spacer between row2 and row3
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od_dis_tr5_v1s <- ggdraw(plot_od_dis_tr5_v1) + 
  draw_label("* <=4mg/1mg Buprenorphine/naloxone", 
             x = 0.68, y = 0.026, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("Note: Comparison for some sensitivity analysis is NA for insufficient number\n   of events or IV assumption was not satisfied.", 
             x = 0.68, y = 0.006, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("A. Buprenorphine/naloxone interruption of 6-14 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone interruption of 6-14 days among prevalent new users", 
             x = 0.02, y = 0.66, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Buprenorphine/naloxone interruption of 6-14 days among all prevalent users", 
             x = 0.02, y = 0.34, hjust = 0, vjust = 0, 
             size = 18, color = "black")


name = "plot_od_dis_tr5_Nov_12"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=8500, res=300)
plot_od_dis_tr5_v1s
dev.off()





# Cumulative Incidence: Nov 12 #####################
res_1 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_fixed_iptw.rds")
res_2 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk2_34d_iptw.rds")
res_3 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk3_iptw.rds")
res_4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")
res_5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")

library(purrr)
library(dplyr)

# put all results into a list
res_list <- list(res_1, res_2, res_3, res_4, res_5)

# extract $points only
res_all <- map2_df(res_list, 1:5, function(res, trial_id) {
  map_df(names(res$points), function(outcome_name) {      # "od" or "dis"
    map_df(names(res$points[[outcome_name]]), function(design_name) {  # "incid", "prev", "rzero"
      df <- res$points[[outcome_name]][[design_name]]
      
      if (is.data.frame(df)) {
        df <- df %>%
          mutate(
            trial   = trial_id,
            outcome = outcome_name,
            design  = design_name
          )
        
        # create baseline row (horizon = 0, risk = 0)
        baseline <- df %>%
          distinct(trt) %>%  # one row per treatment
          mutate(
            horizon = 0,
            risk    = 0,
            trial   = trial_id,
            outcome = outcome_name,
            design  = design_name
          )
        
        df <- bind_rows(baseline, df)
        df
      } else {
        NULL
      }
    })
  })
})


res_all <- res_all %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "reversion" = "Dose reduction",
                  "1-<50%" = "Dose reduction to <50%",
                  ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose*",
                  ">=50%/<=4mg" = "Dose reduction to >=50%/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "Dose reduction", "Dose reduction to <50%",
                                          "Dose reduction to >=50%/initial dose*", "Dose increase")),
         Trial = recode(factor(trial),
                        "1"="Trial 1: MET interruption of 1-2 days",
                        "2"="Trial 2: MET interruption of 3-4 days",
                        "3"="Trial 3: MET interruption of 5-14 days",
                        "4"="Trial 4: BNX interruption of 1-5 days",
                        "5"="Trial 5: BNX interruption of 6-14 days"),
         Design = recode(factor(design),
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users"),
         Outcome = recode(outcome,
                          "od" = "Overdose-related acute care visits",
                          "dis" = "OAT discontinuation"),
         time_day = horizon,
         time = horizon/30.41667
         
  )
saveRDS(res_all, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/All_risk_Nov12.rds")

### Quich graph #########

res_all <- readRS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/All_risk_Nov12.rds")


comparison_colors_v3 <- c(
  "No dose change"              = "grey50", 
  "Dose reduction"              = "#000000",        
  "Dose reduction to <50%"            = "#2CA02C",  
  "Dose reduction to >=50%/initial dose*"  = "#FF7F0E",
  "Dose increase"               = "#D62728"  
)

plot_trial_risk_v4 <- function(df, title_lab, 
                               y_lab="Cumulative incidence of OAT discontinuation",
                               x_lab = "Months on OAT after dose adjustment",
                               x_points = seq(0, 360, 60),  # days for breaks
                               x_min = 0,
                               x_max = 365,
                               y_min, 
                               y_max,
                               y_points,
                               legend_option=TRUE) {
  ggplot(df, aes(x = time, y = risk, colour = group)) +
    geom_line(size = 1.2) +
    facet_wrap(~ Trial, ncol = 1, scales = "free") +
    labs(
      x = x_lab,
      y = y_lab,
      title = title_lab
    ) +
    scale_color_manual(values = comparison_colors_v3) +
    scale_y_continuous(
      breaks = y_points,
      limits = c(y_min, y_max),
      labels = scales::percent_format(accuracy = 1)  # accuracy=1 means no decimals
    ) +
    scale_x_continuous(breaks = x_points, limits = c(x_min, x_max)) +
    theme_minimal(base_size = 15) +
    theme(
      title = element_text(color = "black", size = 18),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.text = element_text(color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      legend.position = if(legend_option) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 15, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 15)
    )
}




### od ##########
od_df_incid = res_all %>% filter(Outcome == "Overdose-related acute care visits",
                                 Design == "Incident users") 

od_df_incid_met = od_df_incid %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                      "Trial 2: MET interruption of 3-4 days", 
                                                      "Trial 3: MET interruption of 5-14 days"))

od_df_incid_bnx = od_df_incid %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                      "Trial 5: BNX interruption of 6-14 days"))


plot_risk_incid_met = plot_trial_risk_v4(df=od_df_incid_met, 
                                         title_lab = "Incident users",
                                         y_lab="Cumulative incidence",
                                         #x_lab="",
                                         x_min = 0,
                                         x_max = 12, 
                                         x_points = c(0, 3, 6, 9, 12),
                                         y_points = c(0, 0.10, 0.20),
                                         y_min = 0,
                                         y_max = 0.24,
                                         legend_option = FALSE)

plot_risk_incid_bnx = plot_trial_risk_v4(df=od_df_incid_bnx, 
                                         title_lab = "Incident users",
                                         y_lab="Cumulative incidence",
                                         
                                         x_min = 0,
                                         x_max = 12, 
                                         x_points = c(0, 3, 6, 9, 12),
                                         y_points = c(0, 0.10, 0.20),
                                         y_min = 0,
                                         y_max = 0.24)



od_df_prev = res_all %>% filter(Outcome == "Overdose-related acute care visits",
                                Design == "Prevalent new users")

od_df_prev_met = od_df_prev %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                    "Trial 2: MET interruption of 3-4 days", 
                                                    "Trial 3: MET interruption of 5-14 days"))

od_df_prev_bnx = od_df_prev %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                    "Trial 5: BNX interruption of 6-14 days"))




plot_risk_prev_met = plot_trial_risk_v4(df=od_df_prev_met, 
                                        title_lab = "Prevalent new users",
                                        y_lab="",
                                        #x_lab="",
                                        x_min = 0,
                                        x_max = 12, 
                                        x_points = c(0, 3, 6, 9, 12),
                                        y_points = c(0, 0.10, 0.20),
                                        y_min = 0,
                                        y_max = 0.24,
                                        legend_option = FALSE)

plot_risk_prev_bnx = plot_trial_risk_v4(df=od_df_prev_bnx, 
                                        title_lab = "Prevalent new users",
                                        y_lab="",
                                        x_min = 0,
                                        x_max = 12, 
                                        x_points = c(0, 3, 6, 9, 12),
                                        y_points = c(0, 0.10, 0.20),
                                        y_min = 0,
                                        y_max = 0.24)




od_df_all = res_all %>% filter(Outcome == "Overdose-related acute care visits",
                               Design == "All prevalent users")

od_df_all_met = od_df_all %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                  "Trial 2: MET interruption of 3-4 days", 
                                                  "Trial 3: MET interruption of 5-14 days"))

od_df_all_bnx = od_df_all %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                  "Trial 5: BNX interruption of 6-14 days"))


plot_risk_all_met = plot_trial_risk_v4(df=od_df_all_met, 
                                       title_lab = "All prevalent users",
                                       y_lab="",
                                       #x_lab="",
                                       x_min = 0,
                                       x_max = 12, 
                                       x_points = c(0, 3, 6, 9, 12),
                                       y_points = c(0, 0.10, 0.20),
                                       y_min = 0,
                                       y_max = 0.24,
                                       legend_option = FALSE)

plot_risk_all_bnx = plot_trial_risk_v4(df=od_df_all_bnx, 
                                       title_lab = "All prevalent users",
                                       y_lab="",
                                       x_min = 0,
                                       x_max = 12, 
                                       x_points = c(0, 3, 6, 9, 12),
                                       y_points = c(0, 0.10, 0.20),
                                       y_min = 0,
                                       y_max = 0.375)



library(patchwork)
library(cowplot)
plot_od_met <- plot_risk_incid_met + plot_spacer() + plot_risk_prev_met + plot_spacer() + plot_risk_all_met + # use `+` to arrange side by side
  plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1)) 


plot_od_bnx <- (plot_risk_incid_bnx + plot_spacer() + plot_risk_prev_bnx +
                  plot_spacer() + plot_risk_all_bnx +
                  plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1), 
                              guides = "collect")) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od <- plot_spacer() /  # top spacer row
  plot_od_met / 
  plot_spacer()/ # spacer between met and bnx
  plot_od_bnx +
  plot_layout(
    heights = c(0.05, 1, 0.08, 0.666667),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od1 <- ggdraw(plot_od) + 
  draw_label("* <=40mg MET for trial 2, 3, and\n<=4/1mg BNX for trial 5", 
             x = 0.99, y = 0.01, hjust = 1, vjust = 0, 
             size = 12, color = "black") +
  draw_label("A. Methadone (MET) trials", 
             x = 0.02, y = 0.97, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone (BNX) trials", 
             x = 0.02, y = 0.42, hjust = 0, vjust = 0, 
             size = 18, color = "black")


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
name = "od_risk_Nov_12"
figure=paste0(folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=5000, height=6000, res=300)
plot_od1
dev.off()



### dis ##########
dis_df_incid = res_all %>% filter(Outcome == "OAT discontinuation",
                                  Design == "Incident users") 

dis_df_incid_met = dis_df_incid %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                        "Trial 2: MET interruption of 3-4 days", 
                                                        "Trial 3: MET interruption of 5-14 days"))

dis_df_incid_bnx = dis_df_incid %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                        "Trial 5: BNX interruption of 6-14 days"))

plot_risk_incid_met = plot_trial_risk_v4(df=dis_df_incid_met, 
                                         title_lab = "Incident users",
                                         y_lab="Cumulative incidence",
                                         #x_lab="",
                                         x_min = 0,
                                         x_max = 12, 
                                         x_points = c(0, 3, 6, 9, 12),
                                         y_points = c(0, 0.3, 0.6),
                                         y_min = 0,
                                         y_max = 0.74,
                                         legend_option = FALSE)

plot_risk_incid_bnx = plot_trial_risk_v4(df=dis_df_incid_bnx, 
                                         title_lab = "Incident users",
                                         y_lab="Cumulative incidence",
                                         
                                         x_min = 0,
                                         x_max = 12, 
                                         x_points = c(0, 3, 6, 9, 12),
                                         y_points = c(0, 0.3, 0.6),
                                         y_min = 0,
                                         y_max = 0.74)


dis_df_prev = res_all %>% filter(Outcome == "OAT discontinuation",
                                 Design == "Prevalent new users")

dis_df_prev_met = dis_df_prev %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                      "Trial 2: MET interruption of 3-4 days", 
                                                      "Trial 3: MET interruption of 5-14 days"))

dis_df_prev_bnx = dis_df_prev %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                      "Trial 5: BNX interruption of 6-14 days"))




plot_risk_prev_met = plot_trial_risk_v4(df=dis_df_prev_met, 
                                        title_lab = "Prevalent new users",
                                        y_lab="",
                                        #x_lab="",
                                        x_min = 0,
                                        x_max = 12, 
                                        x_points = c(0, 3, 6, 9, 12),
                                        y_points = c(0, 0.3, 0.6),
                                        y_min = 0,
                                        y_max = 0.74,
                                        legend_option = FALSE)

plot_risk_prev_bnx = plot_trial_risk_v4(df=dis_df_prev_bnx, 
                                        title_lab = "Prevalent new users",
                                        y_lab="",
                                        x_min = 0,
                                        x_max = 12, 
                                        x_points = c(0, 3, 6, 9, 12),
                                        y_points = c(0, 0.3, 0.6),
                                        y_min = 0,
                                        y_max = 0.74)




dis_df_all = res_all %>% filter(Outcome == "OAT discontinuation",
                                Design == "All prevalent users")

dis_df_all_met = dis_df_all %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                    "Trial 2: MET interruption of 3-4 days", 
                                                    "Trial 3: MET interruption of 5-14 days"))

dis_df_all_bnx = dis_df_all %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                    "Trial 5: BNX interruption of 6-14 days"))


plot_risk_all_met = plot_trial_risk_v4(df=dis_df_all_met, 
                                       title_lab = "All prevalent users",
                                       y_lab="",
                                       #x_lab="",
                                       x_min = 0,
                                       x_max = 12, 
                                       x_points = c(0, 3, 6, 9, 12),
                                       y_points = c(0, 0.3, 0.6),
                                       y_min = 0,
                                       y_max = 0.74,
                                       legend_option = FALSE)


plot_risk_all_bnx = plot_trial_risk_v4(df=dis_df_all_bnx, 
                                       title_lab = "All prevalent users",
                                       y_lab="",
                                       x_min = 0,
                                       x_max = 12, 
                                       x_points = c(0, 3, 6, 9, 12),
                                       y_points = c(0, 0.3, 0.6),
                                       y_min = 0,
                                       y_max = 0.74)



library(patchwork)
plot_dis_met <- plot_risk_incid_met + plot_spacer() + plot_risk_prev_met + plot_spacer() + plot_risk_all_met + # use `+` to arrange side by side
  plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1)) 


plot_dis_bnx <- (plot_risk_incid_bnx + plot_spacer() + plot_risk_prev_bnx +
                   plot_spacer() + plot_risk_all_bnx +
                   plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1), 
                               guides = "collect")) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_dis <- plot_spacer() /  # top spacer row
  plot_dis_met / 
  plot_spacer()/ # spacer between met and bnx
  plot_dis_bnx +
  plot_layout(
    heights = c(0.05, 1, 0.08, 0.666667),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_dis1 <- ggdraw(plot_dis) + 
  draw_label("* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
             x = 0.99, y = 0.01, hjust = 1, vjust = 0, 
             size = 12, color = "black") +
  draw_label("A. Methadone (MET) trials", 
             x = 0.02, y = 0.97, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone (BNX) trials", 
             x = 0.02, y = 0.415, hjust = 0, vjust = 0, 
             size = 18, color = "black")


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
name = "dis_risk_Nov_12"
figure=paste0(folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=5000, height=6000, res=300)
plot_dis1
dev.off()



# Risk Table ###########################
res_1 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_fixed_iptw.rds")
res_2 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk2_34d_iptw.rds")
res_3 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk3_iptw.rds")
res_4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")
res_5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")

library(purrr)
library(dplyr)

# put all results into a list
res_list <- list(res_1, res_2, res_3, res_4, res_5)

library(dplyr)
library(tidyr)
library(glue)

smart_round <- function(x){
  sapply(x, function(v){
    if(is.na(v)) return(NA)
    if(abs(v) < 10) {
      return(round(v, 1))
    } else {
      return(round(v, 1))
    }
  })
}

smart_round_rr <- function(x){
  sapply(x, function(v){
    if(is.na(v)) return(NA)
    if(abs(v) < 10) {
      return(round(v, 2))
    } else {
      return(round(v, 1))
    }
  })
}


make_table <- function(res, trial_id){
  
  map_df(names(res$points_ci), function(outcome_name){       # "od" or "dis"
    map_df(names(res$points_ci[[outcome_name]]), function(design_name){  # "incid", "prev", "rzero"
      
      entry <- res$points_ci[[outcome_name]][[design_name]]
      if(!is.list(entry)) return(NULL)
      
      # --- extract and process risks ---
      risks <- entry$risks %>%
        mutate(
          risk_pct       = smart_round(risk * 100),
          risk_lower_pct = smart_round(risk_lower * 100),
          risk_upper_pct = smart_round(risk_upper * 100),
          `risk%` = glue("{risk_pct}% ({risk_lower_pct}%, {risk_upper_pct}%)")
        ) %>%
        select(trt, horizon, `risk%`) %>%
        pivot_wider(names_from = trt, values_from = `risk%`)
      
      # --- extract and process contrasts ---
      contrasts <- entry$contrasts %>%
        mutate(
          # RD in %
          RD_pct       = smart_round(RD * 100),
          RD_lower_pct = smart_round(RD_lower * 100),
          RD_upper_pct = smart_round(RD_upper * 100),
          `RD%` = glue("{RD_pct}% ({RD_lower_pct}%, {RD_upper_pct}%)"),
          
          # RR (numeric, rounded, no % symbol since it's ratio)
          RR_rnd       = smart_round_rr(RR),
          RR_lower_rnd = smart_round_rr(RR_lower),
          RR_upper_rnd = smart_round_rr(RR_upper),
          `RR` = glue("{RR_rnd} ({RR_lower_rnd}, {RR_upper_rnd})")
        ) %>%
        select(horizon, trt, `RD%`, `RR`)
      
      # --- merge risks + contrasts ---
      df <- left_join(risks, contrasts, by = "horizon")
      
      df %>%
        mutate(
          trial   = trial_id,
          outcome = outcome_name,
          design  = design_name
        )
    })
  })
}


# apply to all trials
all_tables <- map2_df(res_list, 1:5, ~ make_table(.x, .y))



library(dplyr)
library(openxlsx)

tables_by_trial <- split(all_tables, all_tables$trial)



trt_map <- c(
  "no reversion" ,
  "reversion"  ,
  "1-<50%"   ,
  ">=50%/<=40mg" ,
  ">=50%/<=4mg" ,
  "increase" 
)


wb <- createWorkbook()
for(trial_id in names(tables_by_trial)){
  
  df <- tables_by_trial[[trial_id]]
  
  # Remove columns that are entirely NA
  df <- df %>% select(where(~ !all(is.na(.))))
  
  
  df <- df %>%
    select(horizon, any_of(trt_map), everything()) 
  
  if(trial_id %in% c("3", "5")) df <- df %>% filter(horizon %in% c(28, 30, 60, 63, 90, 91, 180, 182, 365))
  
  df <- df %>% mutate(outcome = factor(outcome, levels = c("od", "dis")),
                      trt = factor(trt, levels = trt_map)) %>% 
    arrange(outcome, design, trt, horizon)
  
  # Add a worksheet named by trial
  addWorksheet(wb, paste0("Trial_", trial_id))
  
  # Write the data
  writeData(wb, sheet = paste0("Trial_", trial_id), x = df)
}


saveWorkbook(wb, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/All_risk_tables_v6.xlsx", overwrite = TRUE)




## For writing ##################
res_1 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_iptw.rds")
res_2 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk2_iptw.rds")
res_3 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk3_iptw.rds")
res_4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")
res_5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")

library(purrr)
library(dplyr)

# put all results into a list
res_list <- list(res_1, res_2, res_3, res_4, res_5)

library(dplyr)
library(tidyr)
library(glue)

make_table <- function(res, trial_id){
  
  map_df(names(res$points_ci), function(outcome_name){       # "od" or "dis"
    map_df(names(res$points_ci[[outcome_name]]), function(design_name){  # "incid", "prev", "rzero"
      
      entry <- res$points_ci[[outcome_name]][[design_name]]
      if(!is.list(entry)) return(NULL)
      
      contrasts <- entry$contrasts %>%
        mutate(
          # RD in %
          RD_pct       = smart_round (RD * 100),
          RD_lower_pct = smart_round (RD_lower * 100),
          RD_upper_pct = smart_round (RD_upper * 100),
          # RR (numeric, rounded)
          RR_rnd       = smart_round (RR),
          RR_lower_rnd = smart_round (RR_lower),
          RR_upper_rnd = smart_round (RR_upper)
        ) 
      
      contrasts %>%
        mutate(
          trial   = trial_id,
          outcome = outcome_name,
          design  = design_name
        )
    })
  })
}

# 1st sentence

all_tables <- map2_df(res_list, 1:5, ~ make_table(.x, .y))
all_tables <- all_tables %>% filter(horizon == 365, 
                                    trial %in% 1:3, 
                                    trt != "increase",
                                    RD_upper < 0, 
                                    outcome == "od") %>% 
  select(trt, horizon, RD, RD_lower, RD_upper, design, outcome, trial)

all_tables

res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs.rds")
res <- res %>% filter( 
  Trial %in% c("Trial 1: MET interruption of 1-3 days",
               "Trial 2: MET interruption of 4 days",
               "Trial 3: MET interruption of 5-14 days"), 
  Comparison != "Dose increase\nvs\nNo dose change",
  #lower > 1,
  Outcome == "Overdose-related acute care visits") %>% 
  select(Comparison, HR, lower, upper, Design, Trial)

res

# 2nd sentence 
all_tables <- map2_df(res_list, 1:5, ~ make_table(.x, .y))


table1 <- all_tables %>% filter(horizon == 365, 
                                trial %in% 1:3, 
                                trt == "increase",
                                RD_upper < 0, 
                                outcome == "od") %>% 
  select(trt, horizon, RD, RD_lower, RD_upper, design, outcome, trial)

table1

res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs.rds")
res <- res %>% filter( 
  Trial %in% c("Trial 1: MET interruption of 1-3 days",
               "Trial 2: MET interruption of 4 days",
               "Trial 3: MET interruption of 5-14 days"), 
  Comparison == "Dose increase\nvs\nNo dose change",
  upper < 1,
  Outcome == "Overdose-related acute care visits") %>% 
  select(Comparison, HR, lower, upper, Design, Trial)

res






# 3rd sentence 
all_tables <- map2_df(res_list, 1:5, ~ make_table(.x, .y))


table1 <- all_tables %>% filter(horizon == 365, 
                                trial %in% 1:3, 
                                trt == "increase",
                                RD_lower > 0, 
                                outcome == "od") %>% 
  select(trt, horizon, RD, RD_lower, RD_upper, design, outcome, trial)

table1

res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs.rds")
res <- res %>% filter( 
  Trial %in% c("Trial 1: MET interruption of 1-3 days",
               "Trial 2: MET interruption of 4 days",
               "Trial 3: MET interruption of 5-14 days"), 
  Comparison == "Dose increase\nvs\nNo dose change",
  lower > 1,
  Outcome == "Overdose-related acute care visits") %>% 
  select(Comparison, HR, lower, upper, Design, Trial)

res








# 4th sentence 
all_tables <- map2_df(res_list, 1:5, ~ make_table(.x, .y))


table1 <- all_tables %>% filter(horizon == 365, 
                                trial %in% 4:5, 
                                trt != "increase",
                                RD_lower > 0, 
                                outcome == "od") %>% 
  select(trt, horizon, RD, RD_lower, RD_upper, design, outcome, trial)

table1

res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs.rds")
res <- res %>% filter( 
  Trial %in% c("Trial 4: BNX interruption of 1-5 days",
               "Trial 5: BNX interruption of 6-14 days"), 
  Comparison != "Dose increase\nvs\nNo dose change",
  lower > 1,
  Outcome == "Overdose-related acute care visits") %>% 
  select(Comparison, HR, lower, upper, Design, Trial)

res




# 5th sentence 
table1 <- all_tables %>% filter(horizon == 365, 
                                trial %in% 4:5, 
                                trt != "increase",
                                #RD_lower > 0, 
                                outcome == "od") %>% 
  select(trt, horizon, RD, RD_lower, RD_upper, design, outcome, trial)

table1

res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs.rds")
res <- res %>% filter( 
  Trial %in% c("Trial 4: BNX interruption of 1-5 days",
               "Trial 5: BNX interruption of 6-14 days"), 
  Comparison != "Dose increase\nvs\nNo dose change",
  #lower > 1, 
  Outcome == "Overdose-related acute care visits") %>% 
  select(Comparison, HR, lower, upper, Design, Trial)

res





# 6th sentence 
table1 <- all_tables %>% filter(horizon == 365, 
                                trial %in% 1:3, 
                                trt != "increase",
                                RD_lower > 0, 
                                outcome == "dis") %>% 
  select(trt, horizon, RD, RD_lower, RD_upper, design, outcome, trial)

table1

res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs.rds")
res <- res %>% filter( 
  Trial %in% c("Trial 4: BNX interruption of 1-5 days",
               "Trial 5: BNX interruption of 6-14 days"), 
  Comparison == "Dose increase\nvs\nNo dose change",
  lower > 1, 
  Outcome == "Overdose-related acute care visits") %>% 
  select(Comparison, HR, lower, upper, Design, Trial)

res







# 6th sentence 
res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs.rds")

res1 <- res %>% filter( 
  Trial %in% c("Trial 1: MET interruption of 1-3 days",
               "Trial 2: MET interruption of 4 days",
               "Trial 3: MET interruption of 5-14 days"), 
  Comparison == "Dose increase\nvs\nNo dose change",
  #upper > 1,
  Outcome == "OAT discontinuation") %>% 
  select(Comparison, HR, lower, upper, Design, Trial)

res1



res1 <- res %>% filter( 
  Trial %in% c("Trial 1: MET interruption of 1-3 days",
               "Trial 2: MET interruption of 4 days",
               "Trial 3: MET interruption of 5-14 days"), 
  Comparison != "Dose increase\nvs\nNo dose change",
  upper < 1,
  Outcome == "OAT discontinuation") %>% 
  select(Comparison, HR, lower, upper, Design, Trial)

res1





res1 <- res %>% filter( 
  Trial %in% c("Trial 4: BNX interruption of 1-5 days",
               "Trial 5: BNX interruption of 6-14 days"), 
  Comparison == "Dose increase\nvs\nNo dose change",
  #upper < 1,
  Outcome == "OAT discontinuation") %>% 
  select(Comparison, HR, lower, upper, Design, Trial)

res1




