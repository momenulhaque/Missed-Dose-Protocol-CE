# Trial 1 ###################
## inequality: pres ###############
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction  #############
res_valid = list(grp_red=list(od = list(incid=list(), prev=list(), rzero=list()),
                              dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()



res_valid$grp_red$od$incid <- inequality_test_binaryIV(dat1, 
                                                    iv_cont = "iv_var",
                                                    exposure = "trt",
                                                    outcome = "sens_od_acute")

res_valid$grp_red$dis$incid <- inequality_test_binaryIV(dat1, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_discont")

dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$prev <- inequality_test_binaryIV(dat1_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$prev <- inequality_test_binaryIV(dat1_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                      iv_cont = "iv_var",
                                                      exposure = "trt",
                                                      outcome = "sens_od_acute")

res_valid$grp_red$dis$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_discont")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat1, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat1, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat1_prev, 
                                                      iv_cont = "iv_var",
                                                      exposure = "trt",
                                                      outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat1_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_discont")

dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr1.rds")




## inequality: faci ###############
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction  #############
res_valid = list(grp_red=list(od = list(incid=list(), prev=list(), rzero=list()),
                              dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_red$od$incid <- inequality_test_binaryIV(dat1, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$incid <- inequality_test_binaryIV(dat1, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$prev <- inequality_test_binaryIV(dat1_prev, 
                                                      iv_cont = "iv_var",
                                                      exposure = "trt",
                                                      outcome = "sens_od_acute")

res_valid$grp_red$dis$prev <- inequality_test_binaryIV(dat1_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_discont")

dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat1, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat1, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat1_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat1_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr1.rds")



## bias comp: pres ###################
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction  #############

trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_red$incid <- compute_bias_components_robust(data=dat1, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_red$prev <- compute_bias_components_robust(data=dat1_prev, iv_bin = "IV_bin", 
                                                         exposure = "trt", covars=covariates)



dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_red$rzero <- compute_bias_components_robust(data=dat1_rprev, 
                                                          iv_bin = "IV_bin", 
                                                          exposure = "trt", 
                                                          covars=c(covariates, "prior_gap_b"))





#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat1, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat1_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat1_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr1.rds")


## bias comp: faci ###################
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction  #############

trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_red$incid <- compute_bias_components_robust(data=dat1, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_red$prev <- compute_bias_components_robust(data=dat1_prev, iv_bin = "IV_bin", 
                                                         exposure = "trt", covars=covariates)



dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_red$rzero <- compute_bias_components_robust(data=dat1_rprev, 
                                                          iv_bin = "IV_bin", 
                                                          exposure = "trt", 
                                                          covars=c(covariates, "prior_gap_b"))





#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat1, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat1_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat1_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_bias_comp_tr1.rds")




# Trial 2 ####################################
## inequality: pres #################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = "reversion"
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to 1-<50% ##############
res_valid = list(grp_lt50=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))



trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_lt50$od$incid <- inequality_test_binaryIV(dat2, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$incid <- inequality_test_binaryIV(dat2, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")




dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$prev <- inequality_test_binaryIV(dat2_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_lt50$dis$prev <- inequality_test_binaryIV(dat2_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat2, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat2, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat2_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat2_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr2.rds")



## inequality: faci #################

out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = "reversion"
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to 1-<50% ##############
res_valid = list(grp_lt50=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))



trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_lt50$od$incid <- inequality_test_binaryIV(dat2, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$incid <- inequality_test_binaryIV(dat2, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")




dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$prev <- inequality_test_binaryIV(dat2_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_lt50$dis$prev <- inequality_test_binaryIV(dat2_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat2, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat2, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat2_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat2_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr2.rds")


## bias comp: pres ###################

res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to <50%  #############

trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_lt50$incid <- compute_bias_components_robust(data=dat2, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_lt50$prev <- compute_bias_components_robust(data=dat2_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_lt50$rzero <- compute_bias_components_robust(data=dat2_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat2, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat2_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat2_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr2.rds")


## bias comp: faci ###################
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to <50%  #############

trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_lt50$incid <- compute_bias_components_robust(data=dat2, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_lt50$prev <- compute_bias_components_robust(data=dat2_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_lt50$rzero <- compute_bias_components_robust(data=dat2_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat2, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat2_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat2_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_bias_comp_tr2.rds")




# Trial 3 ###################################
## inequality: pres #################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = "reversion"
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to 1-<50% ##############
res_valid = list(grp_lt50=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))



trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_lt50$od$incid <- inequality_test_binaryIV(dat3, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$incid <- inequality_test_binaryIV(dat3, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")




dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$prev <- inequality_test_binaryIV(dat3_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_lt50$dis$prev <- inequality_test_binaryIV(dat3_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat3_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr3.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$rzero <- inequality_test_binaryIV(dat3_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$rzero <- inequality_test_binaryIV(dat3_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat3, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat3, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat3_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat3_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


dat3_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr3.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat3_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat3_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr3.rds")


## inequality: faci #################


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = "reversion"
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to 1-<50% ##############
res_valid = list(grp_lt50=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))



trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_lt50$od$incid <- inequality_test_binaryIV(dat3, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$incid <- inequality_test_binaryIV(dat3, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")




dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$prev <- inequality_test_binaryIV(dat3_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_lt50$dis$prev <- inequality_test_binaryIV(dat3_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat3_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr3.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$rzero <- inequality_test_binaryIV(dat3_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$rzero <- inequality_test_binaryIV(dat3_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat3, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat3, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat3_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat3_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


dat3_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr3.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat3_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat3_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr3.rds")

## bias comp: pres ###################
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to <50%  #############

trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_lt50$incid <- compute_bias_components_robust(data=dat3, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_lt50$prev <- compute_bias_components_robust(data=dat3_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat3_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr3.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_lt50$rzero <- compute_bias_components_robust(data=dat3_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat3, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat3_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat3_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr3.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat3_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr3.rds")


## bias comp: faci ###################
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to <50%  #############

trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_lt50$incid <- compute_bias_components_robust(data=dat3, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_lt50$prev <- compute_bias_components_robust(data=dat3_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat3_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr3.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_lt50$rzero <- compute_bias_components_robust(data=dat3_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat3, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat3_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat3_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr3.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat3_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_bias_comp_tr3.rds")





# Trial 4 ###################
## inequality: pres #################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=4mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=4mg"
# grp5 = "increase"

#### Reduction  #############
res_valid = list(grp_red=list(od = list(incid=list(), prev=list(), rzero=list()),
                              dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))


trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_red$od$incid <- inequality_test_binaryIV(dat4, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$incid <- inequality_test_binaryIV(dat4, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$prev <- inequality_test_binaryIV(dat4_prev, 
                                                      iv_cont = "iv_var",
                                                      exposure = "trt",
                                                      outcome = "sens_od_acute")

res_valid$grp_red$dis$prev <- inequality_test_binaryIV(dat4_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_discont")

dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$rzero <- inequality_test_binaryIV(dat4_rprev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$rzero <- inequality_test_binaryIV(dat4_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat4, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat4, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat4_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat4_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat4_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat4_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr4.rds")



## inequality: faci #################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=4mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=4mg"
# grp5 = "increase"

#### Reduction  #############
res_valid = list(grp_red=list(od = list(incid=list(), prev=list(), rzero=list()),
                              dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))


trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_red$od$incid <- inequality_test_binaryIV(dat4, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$incid <- inequality_test_binaryIV(dat4, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$prev <- inequality_test_binaryIV(dat4_prev, 
                                                      iv_cont = "iv_var",
                                                      exposure = "trt",
                                                      outcome = "sens_od_acute")

res_valid$grp_red$dis$prev <- inequality_test_binaryIV(dat4_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_discont")

dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$rzero <- inequality_test_binaryIV(dat4_rprev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$rzero <- inequality_test_binaryIV(dat4_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat4, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat4, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat4_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat4_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat4_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat4_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr4.rds")



## bias comp: pres ###################
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=4mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=4mg"
# grp5 = "increase"

#### Reduction  #############

trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_red$incid <- compute_bias_components_robust(data=dat4, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_red$prev <- compute_bias_components_robust(data=dat4_prev, iv_bin = "IV_bin", 
                                                         exposure = "trt", covars=covariates)



dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_red$rzero <- compute_bias_components_robust(data=dat4_rprev, 
                                                          iv_bin = "IV_bin", 
                                                          exposure = "trt", 
                                                          covars=c(covariates, "prior_gap_b"))





#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat4, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat4_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat4_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr4.rds")


## bias comp: faci ###################
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=4mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=4mg"
# grp5 = "increase"

#### Reduction  #############

trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_red$incid <- compute_bias_components_robust(data=dat4, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_red$prev <- compute_bias_components_robust(data=dat4_prev, iv_bin = "IV_bin", 
                                                         exposure = "trt", covars=covariates)



dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_red$rzero <- compute_bias_components_robust(data=dat4_rprev, 
                                                          iv_bin = "IV_bin", 
                                                          exposure = "trt", 
                                                          covars=c(covariates, "prior_gap_b"))





#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat4, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat4_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 16, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat4_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_bias_comp_tr4.rds")




# Trial 5 ###################################
## inequality: pres #################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = "reversion"
# grp3 = "1-<50%"
# grp4 = ">=50%/<=4mg"
# grp5 = "increase"

#### Reduction to 1-<50% ##############
res_valid = list(grp_lt50=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))



trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_lt50$od$incid <- inequality_test_binaryIV(dat5, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$incid <- inequality_test_binaryIV(dat5, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")




dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$prev <- inequality_test_binaryIV(dat5_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_lt50$dis$prev <- inequality_test_binaryIV(dat5_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat5_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr5.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$rzero <- inequality_test_binaryIV(dat5_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$rzero <- inequality_test_binaryIV(dat5_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


#### Reduction to >=50%/<=4mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat5, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat5, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat5_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat5_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


dat5_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr5.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat5_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat5_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr5.rds")



## inequality: faci #################

out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = "reversion"
# grp3 = "1-<50%"
# grp4 = ">=50%/<=4mg"
# grp5 = "increase"

#### Reduction to 1-<50% ##############
res_valid = list(grp_lt50=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))



trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_lt50$od$incid <- inequality_test_binaryIV(dat5, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$incid <- inequality_test_binaryIV(dat5, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")




dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$prev <- inequality_test_binaryIV(dat5_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_lt50$dis$prev <- inequality_test_binaryIV(dat5_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat5_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr5.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$rzero <- inequality_test_binaryIV(dat5_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$rzero <- inequality_test_binaryIV(dat5_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


#### Reduction to >=50%/<=4mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat5, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat5, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat5_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat5_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


dat5_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr5.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat5_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat5_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr5.rds")




## bias comp: pres ###################
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=4mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=4mg"
# grp5 = "increase"

#### Reduction to <50%  #############

trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_lt50$incid <- compute_bias_components_robust(data=dat5, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_lt50$prev <- compute_bias_components_robust(data=dat5_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat5_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr5.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_lt50$rzero <- compute_bias_components_robust(data=dat5_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



#### Reduction to >=50%/<=4mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat5, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat5_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat5_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr5.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat5_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr5.rds")


## bias comp: faci ###################
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=4mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=4mg"
# grp5 = "increase"

#### Reduction to <50%  #############

trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_lt50$incid <- compute_bias_components_robust(data=dat5, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_lt50$prev <- compute_bias_components_robust(data=dat5_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat5_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr5.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_lt50$rzero <- compute_bias_components_robust(data=dat5_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



#### Reduction to >=50%/<=4mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat5, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat5_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat5_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr5.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat5_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_bias_comp_tr5.rds")


# Trial 1: 1-2 days ###################
## inequality: pres ###############
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction  #############
res_valid = list(grp_red=list(od = list(incid=list(), prev=list(), rzero=list()),
                              dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()



res_valid$grp_red$od$incid <- inequality_test_binaryIV(dat1, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$incid <- inequality_test_binaryIV(dat1, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$prev <- inequality_test_binaryIV(dat1_prev, 
                                                      iv_cont = "iv_var",
                                                      exposure = "trt",
                                                      outcome = "sens_od_acute")

res_valid$grp_red$dis$prev <- inequality_test_binaryIV(dat1_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_discont")

dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat1, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat1, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat1_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat1_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr1_12d_ow.rds")




## inequality: faci ###############
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction  #############
res_valid = list(grp_red=list(od = list(incid=list(), prev=list(), rzero=list()),
                              dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_red$od$incid <- inequality_test_binaryIV(dat1, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$incid <- inequality_test_binaryIV(dat1, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$prev <- inequality_test_binaryIV(dat1_prev, 
                                                      iv_cont = "iv_var",
                                                      exposure = "trt",
                                                      outcome = "sens_od_acute")

res_valid$grp_red$dis$prev <- inequality_test_binaryIV(dat1_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_discont")

dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_red$od$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_red$dis$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat1, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat1, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat1_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat1_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat1_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr1_12d_ow.rds")



## bias comp: pres ###################
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction  #############

trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_red$incid <- compute_bias_components_robust(data=dat1, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_red$prev <- compute_bias_components_robust(data=dat1_prev, iv_bin = "IV_bin", 
                                                         exposure = "trt", covars=covariates)



dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_red$rzero <- compute_bias_components_robust(data=dat1_rprev, 
                                                          iv_bin = "IV_bin", 
                                                          exposure = "trt", 
                                                          covars=c(covariates, "prior_gap_b"))





#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat1, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat1_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat1_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr1_12d_ow.rds")


## bias comp: faci ###################
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction  #############

trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_red$incid <- compute_bias_components_robust(data=dat1, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_red$prev <- compute_bias_components_robust(data=dat1_prev, iv_bin = "IV_bin", 
                                                         exposure = "trt", covars=covariates)



dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_red$rzero <- compute_bias_components_robust(data=dat1_rprev, 
                                                          iv_bin = "IV_bin", 
                                                          exposure = "trt", 
                                                          covars=c(covariates, "prior_gap_b"))





#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat1, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat1_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat1_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_bias_comp_tr1_12d_ow.rds")






# Trial 2: 3-4 days ####################################
## inequality: pres #################
out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = "reversion"
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to 1-<50% ##############
res_valid = list(grp_lt50=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))



trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_lt50$od$incid <- inequality_test_binaryIV(dat2, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$incid <- inequality_test_binaryIV(dat2, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")




dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$prev <- inequality_test_binaryIV(dat2_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_lt50$dis$prev <- inequality_test_binaryIV(dat2_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2_34d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat2, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat2, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat2_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat2_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2_34d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr2_34d_ow.rds")



## inequality: faci #################

out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = "reversion"
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to 1-<50% ##############
res_valid = list(grp_lt50=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od = list(incid=list(), prev=list(), rzero=list()),
                               dis = list(incid=list(), prev=list(), rzero=list())))



trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_lt50$od$incid <- inequality_test_binaryIV(dat2, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$incid <- inequality_test_binaryIV(dat2, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")




dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$prev <- inequality_test_binaryIV(dat2_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_lt50$dis$prev <- inequality_test_binaryIV(dat2_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")

dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2_34d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_lt50$od$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_lt50$dis$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$incid <- inequality_test_binaryIV(dat2, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$incid <- inequality_test_binaryIV(dat2, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()

res_valid$grp_incr$od$prev <- inequality_test_binaryIV(dat2_prev, 
                                                       iv_cont = "iv_var",
                                                       exposure = "trt",
                                                       outcome = "sens_od_acute")

res_valid$grp_incr$dis$prev <- inequality_test_binaryIV(dat2_prev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_discont")


dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2_34d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev)) %>% as.data.frame()


res_valid$grp_incr$od$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                        iv_cont = "iv_var",
                                                        exposure = "trt",
                                                        outcome = "sens_od_acute")

res_valid$grp_incr$dis$rzero <- inequality_test_binaryIV(dat2_rprev, 
                                                         iv_cont = "iv_var",
                                                         exposure = "trt",
                                                         outcome = "sens_discont")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr2_34d_ow.rds")


## bias comp: pres ###################

res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to <50%  #############

trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_lt50$incid <- compute_bias_components_robust(data=dat2, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_lt50$prev <- compute_bias_components_robust(data=dat2_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2_34d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_lt50$rzero <- compute_bias_components_robust(data=dat2_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat2, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat2_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2_34d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat2_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr2_34d_ow.rds")


## bias comp: faci ###################
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


out_vars = c("sens_discont", "time_discont", 
             "sens_od_acute", "time_od_acute",
             "sens_all_acute", "time_all_acute")
clust_vars = c("prac_studyid_rev", "studyid_rev", "time0")
baselinevars = c("sexM", "age_st", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")
iv_vars = c("grp1_eps_pct", "grp2_eps_pct", "grp3_eps_pct", "grp4_eps_pct", "grp5_eps_pct")
grp_vars = c("grp1", "grp2", "grp3", "grp4", "grp5")

# grp2 = c("1-<50%", ">=50%/<=40mg")
# grp3 = "1-<50%"
# grp4 = ">=50%/<=40mg"
# grp5 = "increase"

#### Reduction to <50%  #############

trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_lt50$incid <- compute_bias_components_robust(data=dat2, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_lt50$prev <- compute_bias_components_robust(data=dat2_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2_34d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_lt50$rzero <- compute_bias_components_robust(data=dat2_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



#### Reduction to >=50%/<=40mg ###########

# First condition doesn't hold, therefore won't check the 2nd condition

#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()

covariates = c("sexM", "age_lt30", "age_30_39", 
               "age_40_49", "age_ge50", "rural", "year_1",  "year_2", 
               "po_preOUD", "phy_att_gr2", "social_assist", "ret_1","ret_2",
               "home_grp_b", "incar_grp_b2", "drug_acute", "cci_gt0", "copd_asthma", 
               "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
               "udt", "ddose_prev_b", "clin_size_b", "virtual_care", "pss")


res_valid$grp_incr$incid <- compute_bias_components_robust(data=dat2, iv_bin = "IV_bin", 
                                                           exposure = "trt", covars=covariates)



dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()


res_valid$grp_incr$prev <- compute_bias_components_robust(data=dat2_prev, iv_bin = "IV_bin", 
                                                          exposure = "trt", covars=covariates)



dat2_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2_34d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", 
                  iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev),
         year_1 = ifelse(year >= 2016, 1, 0),
         year_2 = ifelse(year >= 2020, 1, 0),
         ret_1 = ifelse(ret_length <= 30, 1, 0),
         ret_2 = ifelse(ret_length <= 90, 1, 0),
         incar_grp_b2 = ifelse(incar_grp_bl=="1_pastyear", 1, 0), 
         ddose_prev_b = ifelse(ddose_prev >= 80, 1, 0),
         prior_gap_b = ifelse(prior_gap <= median(prior_gap, na.rm=TRUE), 1, 0),
         IV_bin = ifelse(iv_var <= median(iv_var, na.rm=TRUE), 0, 1)
  ) %>% as.data.frame()



res_valid$grp_incr$rzero <- compute_bias_components_robust(data=dat2_rprev, 
                                                           iv_bin = "IV_bin", 
                                                           exposure = "trt", 
                                                           covars=c(covariates, "prior_gap_b"))



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_bias_comp_tr2_34d_ow.rds")




## Summary: ineq values ########

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



ineq1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr1.rds")
ineq2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr2.rds")
ineq3 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr3.rds")
ineq4 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr4.rds")
ineq5 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr5.rds")
ineq1_12d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr1_12d_ow.rds")
ineq2_34d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond2_3_tr2_34d_ow.rds")



max_val_1 <- extract_max_values(ineq1, trial_label = "Trial 1") %>% 
  mutate(IV_type = "Prescriber's preference above 50% vs below 50%")
max_val_2 <- extract_max_values(ineq1, trial_label = "Trial 2") %>% 
  mutate(IV_type = "Prescriber's preference above 50% vs below 50%")
max_val_3 <- extract_max_values(ineq1, trial_label = "Trial 3") %>% 
  mutate(IV_type = "Prescriber's preference above 50% vs below 50%")
max_val_4 <- extract_max_values(ineq1, trial_label = "Trial 4") %>% 
  mutate(IV_type = "Prescriber's preference above 50% vs below 50%")
max_val_5 <- extract_max_values(ineq1, trial_label = "Trial 5") %>% 
  mutate(IV_type = "Prescriber's preference above 50% vs below 50%")

max_val_1_12d <- extract_max_values(ineq1_12d, trial_label = "Trial 1_12d") %>% 
  mutate(IV_type = "Prescriber's preference above 50% vs below 50%")
max_val_2_34d <- extract_max_values(ineq2_34d, trial_label = "Trial 2_34d") %>% 
  mutate(IV_type = "Prescriber's preference above 50% vs below 50%")




faci_ineq1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr1.rds")
faci_ineq2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr2.rds")
faci_ineq3 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr3.rds")
faci_ineq4 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr4.rds")
faci_ineq5 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr5.rds")

faci_ineq1_12d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr1_12d_ow.rds")
faci_ineq2_34d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond2_3_tr2_34d_ow.rds")



faci_max_val_1 <- extract_max_values(faci_ineq1, trial_label = "Trial 1") %>% 
  mutate(IV_type = "Facility-level's preference above 50% vs below 50%")
faci_max_val_2 <- extract_max_values(faci_ineq1, trial_label = "Trial 2") %>% 
  mutate(IV_type = "Facility-level's preference above 50% vs below 50%")
faci_max_val_3 <- extract_max_values(faci_ineq1, trial_label = "Trial 3") %>% 
  mutate(IV_type = "Facility-level's preference above 50% vs below 50%")
faci_max_val_4 <- extract_max_values(faci_ineq1, trial_label = "Trial 4") %>% 
  mutate(IV_type = "Facility-level's preference above 50% vs below 50%")
faci_max_val_5 <- extract_max_values(faci_ineq1, trial_label = "Trial 5") %>% 
  mutate(IV_type = "Facility-level's preference above 50% vs below 50%")

faci_max_val_1_12d <- extract_max_values(faci_ineq1_12d, trial_label = "Trial 1_12d") %>% 
  mutate(IV_type = "Facility-level's preference above 50% vs below 50%")
faci_max_val_2_34d <- extract_max_values(faci_ineq2_34d, trial_label = "Trial 2_34d") %>% 
  mutate(IV_type = "Facility-level's preference above 50% vs below 50%")

ineq_values = bind_rows(max_val_1, max_val_2, max_val_1_12d, max_val_2_34d, max_val_3, max_val_4, max_val_5, 
                        faci_max_val_1, faci_max_val_2, faci_max_val_1_12d, faci_max_val_2_34d, faci_max_val_3, faci_max_val_4, faci_max_val_5)


saveRDS(ineq_values, file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/ineq_values.rds")






## Summary: bias plot ########
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



bias_1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr1.rds")
bias_2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr2.rds")
bias_3 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr3.rds")
bias_4 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr4.rds")
bias_5 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr5.rds")

bias_1_12d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr1_12d_ow.rds")
bias_2_34d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_comp_tr2_34d_ow.rds")


res1 <- extract_bias_results(res_obj=bias_1) %>% mutate(Trial = "Trial 1")
res2 <- extract_bias_results(res_obj=bias_2) %>% mutate(Trial = "Trial 2")
res3 <- extract_bias_results(res_obj=bias_3) %>% mutate(Trial = "Trial 3")
res4 <- extract_bias_results(res_obj=bias_4) %>% mutate(Trial = "Trial 4")
res5 <- extract_bias_results(res_obj=bias_5) %>% mutate(Trial = "Trial 5")
res1_12d <- extract_bias_results(res_obj=bias_1_12d) %>% mutate(Trial = "Trial 1_12d")
res2_34d <- extract_bias_results(res_obj=bias_2_34d) %>% mutate(Trial = "Trial 2_34d")


res_bias <- bind_rows(res1, res2, res1_12d, res2_34d, res3, res4, res5) 


bias_plot <- function(df, trial_label = "Trial 1",
                      variable_labels = new.names,
                      comparison_colors_v1 = c(
                        "By treatment"     = "#000000",  # black
                        "By IV (unscaled)" = "#D62728",  # red
                        "By IV (scaled)"   = "#2CA02C"   # green
                      ), legend = FALSE, 
                      y_low = -6,
                      y_high = 6,
                      y_breaks = c(-4, -2, 0, 2, 4),
                      y_text = TRUE,
                      x_lab="Prevalence difference with 95% compatibility intervals") { 
  
  library(ggplot2)
  library(ggh4x)
  library(forcats)
  library(dplyr)
  library(scales)
  
  # Label variables
  df <- df %>%
    mutate(variable_label = recode(variable, !!!variable_labels)) %>%
    filter(!is.na(variable_label)) %>%
    # Reverse the factor levels
    mutate(variable_label = factor(variable_label, levels = rev(variable_labels)),
           Comparison = factor(biascomponent, levels = names(comparison_colors_v1)))
  
  df <- df[!is.na(df$variable_label), ]
  
  # Background stripes
  bg_data <- data.frame(
    variable_label = levels(df$variable_label),
    group_index = seq_along(levels(df$variable_label)),
    fill_color = rep(c("gray85", "white"), length.out = length(levels(df$variable_label)))
  )
  
  dodge <- position_dodge(width = 1)
  
  ggplot(df, aes(x = variable_label, y = Mean, color = Comparison)) +
    geom_rect(data = bg_data, 
              aes(xmin = group_index - 0.5, xmax = group_index + 0.5,
                  ymin = -Inf, ymax = Inf, fill = fill_color),
              inherit.aes = FALSE, alpha = 0.5) +
    geom_point(position = dodge, size = 3, na.rm = TRUE) +
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                  width = 0.1, linewidth = 0.7, position = dodge, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_y_continuous(
      limits = c(y_low, y_high),   # set the range of the x-axis
      breaks = y_breaks,  # customize tick marks
      expand = c(0, 0)         # avoid extra padding
    ) +
    coord_flip() +
    scale_color_manual(values = comparison_colors_v1) +
    scale_fill_identity() +
    labs(
      title = trial_label,
      x = NULL,
      y = x_lab
    ) + 
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.text = element_text(color = "black", size = 14),
      axis.text.y = if(y_text) element_text(hjust = 0, color = "black", size = 15) else element_blank(),
      axis.title.x = element_text(color = "black", size = 15),
      legend.position = if(legend) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 15, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 15),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
}


new.names1 <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year_1        = "Year of reinitiation date: >=2016",
  year_2        = "Year of reinitiation date: >=2020",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_1    = "Duration of OAT episode before time zero: <=30 days",
  ret_2    = "Duration of OAT episode before time zero: <=90 days",
  home_grp_b    = "Unstable housing in past 5 years",
  incar_grp_bl  = "Incarcerated in the past year",
  drug_acute    = "Drug-related acute care visits in the past month",
  cci_gt0       = "Charlson Comorbidity Index (>0)",
  copd_asthma   = "Asthma or COPD in the past year",
  sed_psy_med   = "Dispensations of psychiatric or sedative medication in the past month",
  po            = "Dispensations opioid other than OAT in the past month",
  sud           = "Other substance use disorder",
  aud           = "Alcohol use disorder",
  smd           = "Serious mental disorder",
  hcv           = "Hepatitis C",
  tob           = "Tobacco use disorder (past year)",
  udt           = "Urine drug test in the past month",
  ddose_prev_b    = "Last daily dose: >=80mg",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_gap_b     = "Total OAT interruptions since initiation: above 50%"
)




new.names <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year_1        = "Year of reinitiation date: >=2016",
  year_2        = "Year of reinitiation date: >=2020",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_1    = "Duration of OAT episode before time zero: <=30 days",
  ret_2    = "Duration of OAT episode before time zero: <=90 days",
  home_grp_b    = "Unstable housing in past 5 years",
  incar_grp_bl  = "Incarcerated in the past year",
  drug_acute    = "Drug-related acute care visits in the past month",
  cci_gt0       = "Charlson Comorbidity Index (>0)",
  copd_asthma   = "Asthma or COPD in the past year",
  sed_psy_med   = "Dispensations of psychiatric or sedative medication in the past month",
  po            = "Dispensations opioid other than OAT in the past month",
  sud           = "Other substance use disorder",
  aud           = "Alcohol use disorder",
  smd           = "Serious mental disorder",
  hcv           = "Hepatitis C",
  tob           = "Tobacco use disorder (past year)",
  udt           = "Urine drug test in the past month",
  ddose_prev_b    = "Last daily dose: >=80mg",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply"
)

## Trial 1 #############

### Incident users

res_dat <- res_bias %>% filter(Trial == "Trial 1",
                               Design == "Incident users",
                               Group == "grp_red")

plot_incid_1 <- bias_plot(df=res_dat, trial_label = "Dose reduction\nvs\nNo dose change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 1",
                               Design == "Incident users",
                               Group == "grp_incr")

plot_incid_2 <- bias_plot(df=res_dat, trial_label = "Dose increase\nvs\nNo dose change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL,
                          y_text = FALSE)

### Prevalent new users

res_dat <- res_bias %>% filter(Trial == "Trial 1",
                               Design == "Prevalent new users",
                               Group == "grp_red")

plot_prev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 1",
                               Design == "Prevalent new users",
                               Group == "grp_incr")

plot_prev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL,
                         y_text = FALSE)


### All prevalent users

res_dat <- res_bias %>% filter(Trial == "Trial 1",
                               Design == "All prevalent users",
                               Group == "grp_red")

plot_rprev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab="Prevalence difference with 95% compatibility intervals")

res_dat <- res_bias %>% filter(Trial == "Trial 1",
                               Design == "All prevalent users",
                               Group == "grp_incr")

plot_rprev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          y_text = FALSE,
                          x_lab="Prevalence difference with 95% compatibility intervals")




library(patchwork)
library(cowplot)

row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_prev_1 + plot_spacer() + plot_prev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_1 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer()/
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

library(cowplot)


plot_1s <- ggdraw(plot_1) + 
  draw_label("A. Methadone interruption of 1-3 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 1-3 days among prevalent new users",  
             x = 0.02, y = 0.64, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Methadone interruption of 1-3 days among all prevalent users",  
             x = 0.02, y = 0.33, hjust = 0, vjust = 0, 
             size = 18, color = "black") 

ggsave("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_plot_tr1_Sep22.pdf",
       plot = plot_1s, width = 22, height = 26, units = "in", dpi = 600)


name = "bias_plot_tr1_Sep22"
figure=paste0("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7000, res=300)
plot_1s
dev.off()




## Trial 2 #############

### Incident users

res_dat <- res_bias %>% filter(Trial == "Trial 2",
                               Design == "Incident users",
                               Group == "grp_lt50")

plot_incid_1 <- bias_plot(df=res_dat, trial_label = "1-<50% reduction\nvs\nNo dose change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 2",
                               Design == "Incident users",
                               Group == "grp_incr")

plot_incid_2 <- bias_plot(df=res_dat, trial_label = "Dose increase\nvs\nNo dose change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL,
                          y_text = FALSE)

### Prevalent new users

res_dat <- res_bias %>% filter(Trial == "Trial 2",
                               Design == "Prevalent new users",
                               Group == "grp_lt50")

plot_prev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 2",
                               Design == "Prevalent new users",
                               Group == "grp_incr")

plot_prev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL,
                         y_text = FALSE)


### All prevalent users

res_dat <- res_bias %>% filter(Trial == "Trial 2",
                               Design == "All prevalent users",
                               Group == "grp_lt50")

plot_rprev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab="Prevalence difference with 95% compatibility intervals")

res_dat <- res_bias %>% filter(Trial == "Trial 2",
                               Design == "All prevalent users",
                               Group == "grp_incr")

plot_rprev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          y_text = FALSE,
                          x_lab="Prevalence difference with 95% compatibility intervals")




library(patchwork)
library(cowplot)

row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_prev_1 + plot_spacer() + plot_prev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_2 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer()/
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

library(cowplot)


plot_2s <- ggdraw(plot_2) + 
  draw_label("A. Methadone interruption of 4 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 4 days among prevalent new users",  
             x = 0.02, y = 0.64, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Methadone interruption of 4 days among all prevalent users",  
             x = 0.02, y = 0.33, hjust = 0, vjust = 0, 
             size = 18, color = "black") 

ggsave("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_plot_tr2_Sep22.pdf",
       plot = plot_2s, width = 22, height = 26, units = "in", dpi = 600)


name = "bias_plot_tr2_Sep22"
figure=paste0("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7000, res=300)
plot_2s
dev.off()




## Trial 3 #############

### Incident users

res_dat <- res_bias %>% filter(Trial == "Trial 3",
                               Design == "Incident users",
                               Group == "grp_lt50")

plot_incid_1 <- bias_plot(df=res_dat, trial_label = "1-<50% reduction\nvs\nNo dose change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 3",
                               Design == "Incident users",
                               Group == "grp_incr")

plot_incid_2 <- bias_plot(df=res_dat, trial_label = "Dose increase\nvs\nNo dose change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL,
                          y_text = FALSE)

### Prevalent new users

res_dat <- res_bias %>% filter(Trial == "Trial 3",
                               Design == "Prevalent new users",
                               Group == "grp_lt50")

plot_prev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 3",
                               Design == "Prevalent new users",
                               Group == "grp_incr")

plot_prev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL,
                         y_text = FALSE)


### All prevalent users

res_dat <- res_bias %>% filter(Trial == "Trial 3",
                               Design == "All prevalent users",
                               Group == "grp_lt50")

plot_rprev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab="Prevalence difference with 95% compatibility intervals")

res_dat <- res_bias %>% filter(Trial == "Trial 3",
                               Design == "All prevalent users",
                               Group == "grp_incr")

plot_rprev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          y_text = FALSE,
                          x_lab="Prevalence difference with 95% compatibility intervals")




library(patchwork)
library(cowplot)

row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_prev_1 + plot_spacer() + plot_prev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_3 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer()/
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

library(cowplot)


plot_3s <- ggdraw(plot_3) + 
  draw_label("A. Methadone interruption of 5-14 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 5-14 days among prevalent new users",  
             x = 0.02, y = 0.64, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Methadone interruption of 5-14 days among all prevalent users",  
             x = 0.02, y = 0.33, hjust = 0, vjust = 0, 
             size = 18, color = "black") 

ggsave("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_plot_tr3_Sep22.pdf",
       plot = plot_2s, width = 22, height = 26, units = "in", dpi = 600)


name = "bias_plot_tr3_Sep22"
figure=paste0("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7000, res=300)
plot_3s
dev.off()







## Trial 4 #############

new.names1_bnx <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year_1        = "Year of reinitiation date: >=2016",
  year_2        = "Year of reinitiation date: >=2020",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_1    = "Duration of OAT episode before time zero: <=30 days",
  ret_2    = "Duration of OAT episode before time zero: <=90 days",
  home_grp_b    = "Unstable housing in past 5 years",
  incar_grp_bl  = "Incarcerated in the past year",
  drug_acute    = "Drug-related acute care visits in the past month",
  cci_gt0       = "Charlson Comorbidity Index (>0)",
  copd_asthma   = "Asthma or COPD in the past year",
  sed_psy_med   = "Dispensations of psychiatric or sedative medication in the past month",
  po            = "Dispensations opioid other than OAT in the past month",
  sud           = "Other substance use disorder",
  aud           = "Alcohol use disorder",
  smd           = "Serious mental disorder",
  hcv           = "Hepatitis C",
  tob           = "Tobacco use disorder (past year)",
  udt           = "Urine drug test in the past month",
  ddose_prev_b    = "Last daily dose: >=16mg/4mg",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_gap_b     = "Total OAT interruptions since initiation: above 50%"
)




new.names_bnx <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year_1        = "Year of reinitiation date: >=2016",
  year_2        = "Year of reinitiation date: >=2020",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_1    = "Duration of OAT episode before time zero: <=30 days",
  ret_2    = "Duration of OAT episode before time zero: <=90 days",
  home_grp_b    = "Unstable housing in past 5 years",
  incar_grp_bl  = "Incarcerated in the past year",
  drug_acute    = "Drug-related acute care visits in the past month",
  cci_gt0       = "Charlson Comorbidity Index (>0)",
  copd_asthma   = "Asthma or COPD in the past year",
  sed_psy_med   = "Dispensations of psychiatric or sedative medication in the past month",
  po            = "Dispensations opioid other than OAT in the past month",
  sud           = "Other substance use disorder",
  aud           = "Alcohol use disorder",
  smd           = "Serious mental disorder",
  hcv           = "Hepatitis C",
  tob           = "Tobacco use disorder (past year)",
  udt           = "Urine drug test in the past month",
  ddose_prev_b    = "Last daily dose: >=16mg/4mg",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply"
)

### Incident users

res_dat <- res_bias %>% filter(Trial == "Trial 4",
                               Design == "Incident users",
                               Group == "grp_red")

plot_incid_1 <- bias_plot(df=res_dat, trial_label = "Dose reduction\nvs\nNo dose change",
                          variable_labels=new.names_bnx, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 4",
                               Design == "Incident users",
                               Group == "grp_incr")

plot_incid_2 <- bias_plot(df=res_dat, trial_label = "Dose increase\nvs\nNo dose change",
                          variable_labels=new.names_bnx, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL,
                          y_text = FALSE)

### Prevalent new users

res_dat <- res_bias %>% filter(Trial == "Trial 4",
                               Design == "Prevalent new users",
                               Group == "grp_red")

plot_prev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names_bnx, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 4",
                               Design == "Prevalent new users",
                               Group == "grp_incr")

plot_prev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names_bnx, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL,
                         y_text = FALSE)


### All prevalent users

res_dat <- res_bias %>% filter(Trial == "Trial 4",
                               Design == "All prevalent users",
                               Group == "grp_red")

plot_rprev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1_bnx, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab="Prevalence difference with 95% compatibility intervals")

res_dat <- res_bias %>% filter(Trial == "Trial 4",
                               Design == "All prevalent users",
                               Group == "grp_incr")

plot_rprev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1_bnx, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          y_text = FALSE,
                          x_lab="Prevalence difference with 95% compatibility intervals")




library(patchwork)
library(cowplot)

row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_prev_1 + plot_spacer() + plot_prev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_4 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer()/
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

library(cowplot)


plot_4s <- ggdraw(plot_4) + 
  draw_label("A. Buprenorphine/naloxone interruption of 1-5 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone interruption of 1-5 days among prevalent new users",  
             x = 0.02, y = 0.64, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Buprenorphine/naloxone interruption of 1-5 days among all prevalent users",  
             x = 0.02, y = 0.33, hjust = 0, vjust = 0, 
             size = 18, color = "black") 

ggsave("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_plot_tr4_Sep22.pdf",
       plot = plot_1s, width = 22, height = 26, units = "in", dpi = 600)


name = "bias_plot_tr4_Sep22"
figure=paste0("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7000, res=300)
plot_4s
dev.off()






## Trial 5 #############

### Incident users

res_dat <- res_bias %>% filter(Trial == "Trial 5",
                               Design == "Incident users",
                               Group == "grp_lt50")

plot_incid_1 <- bias_plot(df=res_dat, trial_label = "1-<50% reduction\nvs\nNo dose change",
                          variable_labels=new.names_bnx, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 5",
                               Design == "Incident users",
                               Group == "grp_incr")

plot_incid_2 <- bias_plot(df=res_dat, trial_label = "Dose increase\nvs\nNo dose change",
                          variable_labels=new.names_bnx, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL,
                          y_text = FALSE)

### Prevalent new users

res_dat <- res_bias %>% filter(Trial == "Trial 5",
                               Design == "Prevalent new users",
                               Group == "grp_lt50")

plot_prev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names_bnx, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 5",
                               Design == "Prevalent new users",
                               Group == "grp_incr")

plot_prev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names_bnx, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL,
                         y_text = FALSE)


### All prevalent users

res_dat <- res_bias %>% filter(Trial == "Trial 5",
                               Design == "All prevalent users",
                               Group == "grp_lt50")

plot_rprev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1_bnx, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab="Prevalence difference with 95% compatibility intervals")

res_dat <- res_bias %>% filter(Trial == "Trial 5",
                               Design == "All prevalent users",
                               Group == "grp_incr")

plot_rprev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1_bnx, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          y_text = FALSE,
                          x_lab="Prevalence difference with 95% compatibility intervals")




library(patchwork)
library(cowplot)

row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_prev_1 + plot_spacer() + plot_prev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_5 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer()/
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

library(cowplot)


plot_5s <- ggdraw(plot_5) + 
  draw_label("A. Buprenophine/naloxone interruption of 5-14 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenophine/naloxone interruption of 5-14 days among prevalent new users",  
             x = 0.02, y = 0.64, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Buprenophine/naloxone interruption of 5-14 days among all prevalent users",  
             x = 0.02, y = 0.33, hjust = 0, vjust = 0, 
             size = 18, color = "black") 

ggsave("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/bias_plot_tr5_Sep22.pdf",
       plot = plot_5s, width = 22, height = 26, units = "in", dpi = 600)


name = "bias_plot_tr5_Sep22"
figure=paste0("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7000, res=300)
plot_5s
dev.off()



## Trial 1: 1-2 days #############

### Incident users

res_dat <- res_bias %>% filter(Trial == "Trial 1_12d",
                               Design == "Incident users",
                               Group == "grp_red")

plot_incid_1 <- bias_plot(df=res_dat, trial_label = "Dose reduction\nvs\nno change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 1_12d",
                               Design == "Incident users",
                               Group == "grp_incr")

plot_incid_2 <- bias_plot(df=res_dat, trial_label = "Dose increase\nvs\nno change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL,
                          y_text = FALSE)

### Prevalent new users

res_dat <- res_bias %>% filter(Trial == "Trial 1_12d",
                               Design == "Prevalent new users",
                               Group == "grp_red")

plot_prev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 1_12d",
                               Design == "Prevalent new users",
                               Group == "grp_incr")

plot_prev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL,
                         y_text = FALSE)


### All prevalent users

res_dat <- res_bias %>% filter(Trial == "Trial 1_12d",
                               Design == "All prevalent users",
                               Group == "grp_red")

plot_rprev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab="Prevalence difference with 95% compatibility intervals")

res_dat <- res_bias %>% filter(Trial == "Trial 1_12d",
                               Design == "All prevalent users",
                               Group == "grp_incr")

plot_rprev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          y_text = FALSE,
                          x_lab="Prevalence difference with 95% compatibility intervals")




library(patchwork)
library(cowplot)

row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_prev_1 + plot_spacer() + plot_prev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_1 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer()/
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

library(cowplot)


plot_1s <- ggdraw(plot_1) + 
  draw_label("A. Methadone interruption of 1-2 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 1-2 days among prevalent new users",  
             x = 0.02, y = 0.64, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Methadone interruption of 1-2 days among all prevalent users",  
             x = 0.02, y = 0.33, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


name = "bias_plot_tr1_12d_Sep22"
figure=paste0("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7800, res=300)
plot_1s
dev.off()




## Trial 2: 1-2 days #############

### Incident users

res_dat <- res_bias %>% filter(Trial == "Trial 2_34d",
                               Design == "Incident users",
                               Group == "grp_lt50")

plot_incid_1 <- bias_plot(df=res_dat, trial_label = "Dose reduction to <50%\nvs\nno change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 2_34d",
                               Design == "Incident users",
                               Group == "grp_incr")

plot_incid_2 <- bias_plot(df=res_dat, trial_label = "Dose increase\nvs\nno change",
                          variable_labels=new.names, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab=NULL,
                          y_text = FALSE)

### Prevalent new users

res_dat <- res_bias %>% filter(Trial == "Trial 2_34d",
                               Design == "Prevalent new users",
                               Group == "grp_lt50")

plot_prev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL)

res_dat <- res_bias %>% filter(Trial == "Trial 2_34d",
                               Design == "Prevalent new users",
                               Group == "grp_incr")

plot_prev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                         variable_labels=new.names, legend = FALSE,
                         y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                         x_lab=NULL,
                         y_text = FALSE)


### All prevalent users

res_dat <- res_bias %>% filter(Trial == "Trial 2_34d",
                               Design == "All prevalent users",
                               Group == "grp_lt50")

plot_rprev_1 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          x_lab="Prevalence difference with 95% compatibility intervals")

res_dat <- res_bias %>% filter(Trial == "Trial 2_34d",
                               Design == "All prevalent users",
                               Group == "grp_incr")

plot_rprev_2 <- bias_plot(df=res_dat, trial_label = NULL,
                          variable_labels=new.names1, legend = FALSE,
                          y_low = -6, y_high = 6, y_breaks = c(-4, -2, 0, 2, 4),
                          y_text = FALSE,
                          x_lab="Prevalence difference with 95% compatibility intervals")




library(patchwork)
library(cowplot)

row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + 
  plot_layout(widths = c(1, 0.07, 1))

row2 <- plot_prev_1 + plot_spacer() + plot_prev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_2 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer()/
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

library(cowplot)


plot_2s <- ggdraw(plot_2) + 
  draw_label("A. Methadone interruption of 3-4 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 3-4 days among prevalent new users",  
             x = 0.02, y = 0.64, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("C. Methadone interruption of 3-4 days among all prevalent users",  
             x = 0.02, y = 0.33, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


name = "bias_plot_tr2_34d_Sep22"
figure=paste0("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7800, res=300)
plot_2s
dev.off()



