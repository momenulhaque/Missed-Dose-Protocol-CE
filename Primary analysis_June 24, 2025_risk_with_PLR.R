# Cum Incid from PLR ###################################
## Trial 1 ################################
rm(list = ls())
gc()
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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

### incident users ########################
risk_Trial=list(od=list(incid=list(), prev=list(), rzero=list()), 
                dis=list(incid=list(), prev=list(), rzero=list()))

risk_points=list(od=list(incid=list(), prev=list(), rzero=list()), 
                 dis=list(incid=list(), prev=list(), rzero=list()))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

risk_points$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("ret_length", "ddose_prev"), 
                                        horizons=c(seq(14, 357, by=14), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 50, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("ret_length", "ddose_prev"), 
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 100, conf_int = TRUE,
                                       ref_group = "no reversion")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


risk_points$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("ret_length", "ddose_prev"),  
                                        horizons=c(seq(14, 357, by=14), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 50, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("ret_length", "ddose_prev"),  
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = TRUE,
                                       ref_group = "no reversion")

saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_iptw.rds")

### Prevalent new users ########################

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

risk_points$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("ret_length", "ddose_prev"),
                                      horizons=c(seq(14, 357, by=14), 365),
                                      cluster_var="studyid_rev", 
                                      samples = 50, conf_int = FALSE,
                                      ref_group = "no reversion")


risk_Trial$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                     trt_var="trt", 
                                     covs =  c("ret_length", "ddose_prev"),
                                     horizons=c(30, 60, 90, 180, 365),
                                     cluster_var="studyid_rev", 
                                     samples = 50, conf_int = TRUE,
                                     ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


risk_points$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs =  c("ret_length", "ddose_prev"),
                                       horizons=c(seq(14, 357, by=14), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = FALSE,
                                       ref_group = "no reversion")


risk_Trial$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt",
                                      covs =  c("ret_length", "ddose_prev"),
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 50, conf_int = TRUE,
                                      ref_group = "no reversion")

saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_iptw.rds")

res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_iptw.rds")
risk_points <- res$points
risk_Trial <- res$points_ci

### All prevalent users ################
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_rprev.rds")

od_results <- predict_cuminc_coxPH(Dat = ps_fit$data_weights,
                                time_var = "time_od_acute",
                                event_var = "sens_od_acute",
                                covs = c("ret_length", "ddose_prev"),
                                weight_var = "wt_trunc", 
                                cluster_var = "studyid_rev",
                                horizons_risk = c(seq(14, 357, by=14), 365),
                                horizons_ci = c(30, 60, 90, 180, 365))



risk_points$od$rzero <- tibble(od_results$points %>% select(trt, horizon, risk))

risk_Trial$od$rzero$risks <- od_results$risk
risk_Trial$od$rzero$contrasts <- od_results$contrasts
risk_Trial$od$rzero$ref_group <- od_results$ref_group


dis_results <- predict_cuminc_coxPH(Dat = ps_fit$data_weights,
                                time_var = "time_discont",
                                event_var = "sens_discont",
                                covs = c("ret_length", "ddose_prev"),
                                weight_var = "wt_trunc", 
                                cluster_var = "studyid_rev",
                                horizons_risk = c(seq(14, 357, by=14), 365),
                                horizons_ci = c(30, 60, 90, 180, 365))




risk_points$dis$rzero <- tibble(dis_results$points %>% select(trt, horizon, risk))

risk_Trial$dis$rzero$risks <- dis_results$risk
risk_Trial$dis$rzero$contrasts <- dis_results$contrasts
risk_Trial$dis$rzero$ref_group <- dis_results$ref_group

saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_iptw.rds")
  


## Trial 2 ################################
rm(list = ls())
gc()
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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

### incident users ########################
risk_Trial=list(od=list(incid=list(), prev=list(), rzero=list()), 
                dis=list(incid=list(), prev=list(), rzero=list()))

risk_points=list(od=list(incid=list(), prev=list(), rzero=list()), 
                 dis=list(incid=list(), prev=list(), rzero=list()))



ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$incid <- NA


risk_Trial$od$incid <- NA

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$incid <- predict_cuminc( 
                                        Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("age_st", "age_ge50", "ret_length",
                                                 "copd_asthma", "sed_psy_med", 
                                                 "po", "sud"),
                                        horizons=c(seq(7, 357, by=7), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 100, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$incid <- predict_cuminc(
                                       Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("age_st", "age_ge50", "ret_length",
                                                "copd_asthma", "sed_psy_med", 
                                                "po", "sud"),
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 100, conf_int = TRUE,
                                       ref_group = "no reversion")



### Prevalent new users ########################

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$prev <- predict_cuminc(
                                      Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("ddose_prev"),
                                      horizons=c(seq(7, 357, by=7), 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = FALSE,
                                      ref_group = "no reversion")


risk_Trial$od$prev <- predict_cuminc( 
                                     Dat_expanded=long_dat, 
                                     trt_var="trt", 
                                     covs = c("ddose_prev"),
                                     horizons=c(30, 60, 90, 180, 365),
                                     cluster_var="studyid_rev", 
                                     samples = 100, conf_int = TRUE,
                                     ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$prev <- predict_cuminc(
                                       Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("ddose_prev"),
                                       horizons=c(seq(7, 357, by=7), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 100, conf_int = FALSE,
                                       ref_group = "no reversion")


risk_Trial$dis$prev <- predict_cuminc(
                                      Dat_expanded=long_dat, 
                                      trt_var="trt",
                                      covs = c("ddose_prev"),
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = TRUE,
                                      ref_group = "no reversion")


### All prevalent users ################
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$rzero <- predict_cuminc( 
                                       Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("udt"),
                                       horizons=c(seq(7, 357, by=7), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 100, conf_int = FALSE,
                                       ref_group = "no reversion")


risk_Trial$od$rzero <- predict_cuminc(
                                      Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("udt"),
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = TRUE,
                                      ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$rzero <- predict_cuminc(
                                        Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("udt"),
                                        horizons=c(seq(7, 357, by=7), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 100, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$rzero <- predict_cuminc(
                                       Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("udt"),
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 100, conf_int = TRUE,
                                       ref_group = "no reversion")

saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk2_iptw.rds")



## Trial 3 ################################
rm(list = ls())
gc()
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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

### incident users ########################
risk_Trial=list(od=list(incid=list(), prev=list(), rzero=list()), 
                dis=list(incid=list(), prev=list(), rzero=list()))

risk_points=list(od=list(incid=list(), prev=list(), rzero=list()), 
                 dis=list(incid=list(), prev=list(), rzero=list()))



ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("ret_length", "year", "ddose_prev"), 
                                      horizons=c(seq(7, 357, by=7), 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = FALSE,
                                      ref_group = "no reversion")


risk_Trial$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("ret_length", "year", "ddose_prev"), 
                                      horizons=c(seq(7, 357, by=7), 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = TRUE,
                                      ref_group = "no reversion")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("ret_length", "year", "ddose_prev"), 
                                        horizons=c(seq(7, 357, by=7), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 100, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("ret_length", "year", "ddose_prev"), 
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = TRUE,
                                      ref_group = "no reversion")


### Prevalent new users ########################

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("ddose_prev"),
                                      horizons=c(seq(7, 357, by=7), 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = FALSE,
                                      ref_group = "no reversion")


risk_Trial$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("ddose_prev"),
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = TRUE,
                                      ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("ddose_prev"),
                                        horizons=c(seq(7, 357, by=7), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 100, conf_int = FALSE,
                                        ref_group = "no reversion")
  

risk_Trial$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt",
                                      covs = c("ddose_prev"),
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = TRUE,
                                      ref_group = "no reversion")


### All prevalent users ################
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("ddose_prev"),
                                      horizons=c(seq(7, 357, by=7), 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = FALSE,
                                      ref_group = "no reversion")


risk_Trial$od$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("ddose_prev"),
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = TRUE,
                                      ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("ddose_prev"),
                                        horizons=c(seq(7, 357, by=7), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 100, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                       covs = c("ddose_prev"),
                                        horizons=c(30, 60, 90, 180, 365),
                                        cluster_var="studyid_rev", 
                                        samples = 100, conf_int = TRUE,
                                        ref_group = "no reversion")


saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk3_iptw.rds")




## Trial 4 ################################
rm(list = ls())
gc()
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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

### incident users ########################
risk_Trial=list(od=list(incid=list(), prev=list(), rzero=list()), 
                dis=list(incid=list(), prev=list(), rzero=list()))

risk_points=list(od=list(incid=list(), prev=list(), rzero=list()), 
                 dis=list(incid=list(), prev=list(), rzero=list()))



ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

risk_points$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = NULL, 
                                        horizons=c(seq(14, 357, by=14), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 50, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = NULL, 
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = TRUE,
                                       ref_group = "no reversion")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


risk_points$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = NULL, 
                                        horizons=c(seq(14, 357, by=14), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 50, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = NULL, 
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = TRUE,
                                       ref_group = "no reversion")


### Prevalent new users ########################

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

risk_points$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = NULL,
                                      horizons=c(seq(14, 357, by=14), 365),
                                      cluster_var="studyid_rev", 
                                      samples = 50, conf_int = FALSE,
                                      ref_group = "no reversion")


risk_Trial$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                     trt_var="trt", 
                                     covs = NULL,
                                     horizons=c(30, 60, 90, 180, 365),
                                     cluster_var="studyid_rev", 
                                     samples = 50, conf_int = TRUE,
                                     ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


risk_points$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = NULL,
                                       horizons=c(seq(14, 357, by=14), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = FALSE,
                                       ref_group = "no reversion")


risk_Trial$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt",
                                      covs = NULL,
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 50, conf_int = TRUE,
                                      ref_group = "no reversion")


### All prevalent users ################
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

risk_points$od$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = "pss",
                                       horizons=c(seq(14, 357, by=14), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = FALSE,
                                       ref_group = "no reversion")


risk_Trial$od$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = "pss",
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 50, conf_int = TRUE,
                                      ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


risk_points$dis$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = "pss",
                                        horizons=c(seq(14, 357, by=14), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 50, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = "pss",
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = TRUE,
                                       ref_group = "no reversion")


saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")


## Trial 5 ################################
rm(list = ls())
gc()
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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

### incident users ########################
risk_Trial=list(od=list(incid=list(), prev=list(), rzero=list()), 
                dis=list(incid=list(), prev=list(), rzero=list()))

risk_points=list(od=list(incid=list(), prev=list(), rzero=list()), 
                 dis=list(incid=list(), prev=list(), rzero=list()))



ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("year", "ret_length", "ddose_prev"),
                                        horizons=c(seq(7, 357, by=7), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 50, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("year", "ret_length", "ddose_prev"), 
                                       horizons=c(seq(7, 357, by=7), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = TRUE,
                                       ref_group = "no reversion")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("year", "ret_length", "ddose_prev"), 
                                        horizons=c(seq(7, 357, by=7), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 50, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("year", "ret_length", "ddose_prev"), 
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = TRUE,
                                       ref_group = "no reversion")


### Prevalent new users ########################

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_prev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = NULL,
                                      horizons=c(seq(7, 357, by=7), 365),
                                      cluster_var="studyid_rev", 
                                      samples = 50, conf_int = FALSE,
                                      ref_group = "no reversion")


risk_Trial$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                     trt_var="trt", 
                                     covs = NULL,
                                     horizons=c(30, 60, 90, 180, 365),
                                     cluster_var="studyid_rev", 
                                     samples = 50, conf_int = TRUE,
                                     ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = NULL,
                                       horizons=c(seq(7, 357, by=7), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = FALSE,
                                       ref_group = "no reversion")


risk_Trial$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt",
                                      covs = NULL,
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 50, conf_int = TRUE,
                                      ref_group = "no reversion")


### All prevalent users ################
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = NULL,
                                       horizons=c(seq(7, 357, by=7), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = FALSE,
                                       ref_group = "no reversion")


risk_Trial$od$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = NULL,
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 50, conf_int = TRUE,
                                      ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = NULL,
                                        horizons=c(seq(7, 357, by=7), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 50, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = NULL,
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = TRUE,
                                       ref_group = "no reversion")


saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")




## Trial 1: 1-2 days ################################
rm(list = ls())
gc()
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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

### incident users ########################
risk_Trial=list(od=list(incid=list(), prev=list(), rzero=list()), 
                dis=list(incid=list(), prev=list(), rzero=list()))

risk_points=list(od=list(incid=list(), prev=list(), rzero=list()), 
                 dis=list(incid=list(), prev=list(), rzero=list()))



ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_ow.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

risk_points$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("year", "ret_length"),  
                                        horizons=c(seq(14, 357, by=14), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 100, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("ret_length", "ddose_prev"), 
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 100, conf_int = TRUE,
                                       ref_group = "no reversion")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


risk_points$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs = c("ret_length", "ddose_prev"),  
                                        horizons=c(seq(14, 357, by=14), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 100, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$incid <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("ret_length", "ddose_prev"),  
                                       horizons=c(30, 60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 100, conf_int = TRUE,
                                       ref_group = "no reversion")

saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_iptw.rds")

### Prevalent new users ########################

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_prev_ow.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

risk_points$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs = c("year", "ret_length", "ddose_prev"),
                                      horizons=c(seq(14, 357, by=14), 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = FALSE,
                                      ref_group = "no reversion")


risk_Trial$od$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                     trt_var="trt", 
                                     covs =  c("year", "ret_length", "ddose_prev"),
                                     horizons=c(30, 60, 90, 180, 365),
                                     cluster_var="studyid_rev", 
                                     samples = 100, conf_int = TRUE,
                                     ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


risk_points$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs =  c("year", "ret_length", "ddose_prev"),
                                       horizons=c(seq(14, 357, by=14), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 100, conf_int = FALSE,
                                       ref_group = "no reversion")


risk_Trial$dis$prev <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt",
                                      covs =  c("year", "ret_length", "ddose_prev"),
                                      horizons=c(30, 60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 100, conf_int = TRUE,
                                      ref_group = "no reversion")

saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_iptw.rds")

res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_iptw.rds")
risk_points <- res$points
risk_Trial <- res$points_ci

### All prevalent users ################
library(stdReg)
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")

od_results <- predict_cuminc_coxPH(Dat = ps_fit$data_weights,
                                   time_var = "time_od_acute",
                                   event_var = "sens_od_acute",
                                   covs = c("ret_length", "ddose_prev"), 
                                   weight_var = "wt_trunc", 
                                   cluster_var = "studyid_rev",
                                   horizons_risk = c(seq(14, 357, by=14), 365),
                                   horizons_ci = c(30, 60, 90, 180, 365))



risk_points$od$rzero <- tibble(od_results$points %>% select(trt, horizon, risk))

risk_Trial$od$rzero$risks <- od_results$risk
risk_Trial$od$rzero$contrasts <- od_results$contrasts
risk_Trial$od$rzero$ref_group <- od_results$ref_group


dis_results <- predict_cuminc_coxPH(Dat = ps_fit$data_weights,
                                    time_var = "time_discont",
                                    event_var = "sens_discont",
                                    covs = c("ret_length", "ddose_prev"), 
                                    weight_var = "wt_trunc", 
                                    cluster_var = "studyid_rev",
                                    horizons_risk = c(seq(14, 357, by=14), 365),
                                    horizons_ci = c(30, 60, 90, 180, 365))




risk_points$dis$rzero <- tibble(dis_results$points %>% select(trt, horizon, risk))

risk_Trial$dis$rzero$risks <- dis_results$risk
risk_Trial$dis$rzero$contrasts <- dis_results$contrasts
risk_Trial$dis$rzero$ref_group <- dis_results$ref_group

saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_iptw.rds")

#### Fixed part of All prevalent users  ##################
res <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_iptw.rds")
risk_points <- res$points
risk_Trial <- res$points_ci

library(stdReg)
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")



long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 60)

risk_points$od$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt", 
                                       covs = c("ret_length"),
                                       horizons=c(seq(60, 357, by=30), 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = FALSE,
                                       ref_group = "no reversion")


risk_Trial$od$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                      trt_var="trt", 
                                      covs =  c("ret_length"),
                                      horizons=c(60, 90, 180, 365),
                                      cluster_var="studyid_rev", 
                                      samples = 50, conf_int = TRUE,
                                      ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 60)


risk_points$dis$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                        trt_var="trt", 
                                        covs =  c("ret_length"),
                                        horizons=c(seq(60, 357, by=30), 365),
                                        cluster_var="studyid_rev", 
                                        samples = 50, conf_int = FALSE,
                                        ref_group = "no reversion")


risk_Trial$dis$rzero <- predict_cuminc(Dat_expanded=long_dat, 
                                       trt_var="trt",
                                       covs =  c("ret_length"),
                                       horizons=c(60, 90, 180, 365),
                                       cluster_var="studyid_rev", 
                                       samples = 50, conf_int = TRUE,
                                       ref_group = "no reversion")


saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_fixed_iptw.rds")



## Trial 2: 3-4 days ################################
rm(list = ls())
gc()
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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

### incident users ########################
risk_Trial=list(od=list(incid=list(), prev=list(), rzero=list()), 
                dis=list(incid=list(), prev=list(), rzero=list()))

risk_points=list(od=list(incid=list(), prev=list(), rzero=list()), 
                 dis=list(incid=list(), prev=list(), rzero=list()))



ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_ow.rds")

risk_points$od$incid <- NA


risk_Trial$od$incid <- NA

# long_dat = expand_person_period(Dat=ps_fit$data_weights, 
#                                 time_var="time_od_acute", 
#                                 event_var="sens_od_acute", 
#                                 interval_length = 14)
# 
# risk_points$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
#                                         trt_var="trt", 
#                                         covs = c("age_st", "sud",
#                                                  "phy_att_gr2", "social_assist"),
#                                         horizons=c(seq(14, 357, by=14), 365),
#                                         cluster_var="studyid_rev", 
#                                         samples = 100, conf_int = FALSE,
#                                         ref_group = "no reversion")
# 
# 
# risk_Trial$od$incid <-  predict_cuminc(Dat_expanded=long_dat, 
#                                        trt_var="trt", 
#                                        covs = c("age_st", "sud",
#                                                 "phy_att_gr2", "social_assist"),
#                                        horizons=c(30, 60, 90, 180, 365),
#                                        cluster_var="studyid_rev", 
#                                        samples = 100, conf_int = TRUE,
#                                        ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$incid <- predict_cuminc( 
  Dat_expanded=long_dat, 
  trt_var="trt", 
  covs = c("age_st", "sud",
           "phy_att_gr2", "social_assist"),
  horizons=c(seq(7, 357, by=7), 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = FALSE,
  ref_group = "no reversion")


risk_Trial$dis$incid <- predict_cuminc(
  Dat_expanded=long_dat, 
  trt_var="trt", 
  covs = c("age_st", "sud",
           "phy_att_gr2", "social_assist"),
  horizons=c(30, 60, 90, 180, 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = TRUE,
  ref_group = "no reversion")



### Prevalent new users ########################

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_prev_ow.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)

risk_points$od$prev <- predict_cuminc(
  Dat_expanded=long_dat, 
  trt_var="trt", 
  covs = c("ddose_prev"),
  horizons=c(seq(7, 357, by=7), 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = FALSE,
  ref_group = "no reversion")


risk_Trial$od$prev <- predict_cuminc( 
  Dat_expanded=long_dat, 
  trt_var="trt", 
  covs = c("ddose_prev"),
  horizons=c(30, 60, 90, 180, 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = TRUE,
  ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$prev <- predict_cuminc(
  Dat_expanded=long_dat, 
  trt_var="trt", 
  covs = c("ddose_prev"),
  horizons=c(seq(7, 357, by=7), 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = FALSE,
  ref_group = "no reversion")


risk_Trial$dis$prev <- predict_cuminc(
  Dat_expanded=long_dat, 
  trt_var="trt",
  covs = c("ddose_prev"),
  horizons=c(30, 60, 90, 180, 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = TRUE,
  ref_group = "no reversion")


### All prevalent users ################
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_rprev_ow.rds")

long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)

risk_points$od$rzero <- predict_cuminc( 
  Dat_expanded=long_dat, 
  trt_var="trt", 
  covs = "ddose_prev", 
  horizons=c(seq(14, 357, by=14), 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = FALSE,
  ref_group = "no reversion")


risk_Trial$od$rzero <- predict_cuminc(
  Dat_expanded=long_dat, 
  trt_var="trt", 
  covs = "ddose_prev", 
  horizons=c(30, 60, 90, 180, 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = TRUE,
  ref_group = "no reversion")


long_dat = expand_person_period(Dat=ps_fit$data_weights, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


risk_points$dis$rzero <- predict_cuminc(
  Dat_expanded=long_dat, 
  trt_var="trt", 
  covs = "ddose_prev", 
  horizons=c(seq(7, 357, by=7), 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = FALSE,
  ref_group = "no reversion")


risk_Trial$dis$rzero <- predict_cuminc(
  Dat_expanded=long_dat, 
  trt_var="trt", 
  covs = "ddose_prev", 
  horizons=c(30, 60, 90, 180, 365),
  cluster_var="studyid_rev", 
  samples = 100, conf_int = TRUE,
  ref_group = "no reversion")

saveRDS(list(points=risk_points, points_ci=risk_Trial), file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk2_34d_iptw.rds")

