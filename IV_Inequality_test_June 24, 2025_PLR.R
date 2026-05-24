## Trial 1 ###################
### Condition 2 & 3: pres #################
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




### Condition 2 & 3: faci #################
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



### bias comp: pres ###################
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


### bias comp: faci ###################
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




## Trial 2 ####################################
### ### Condition 2 & 3: pres #################
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




### ### Condition 2 & 3: faci #################
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

### bias comp: pres ###################
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


### bias comp: faci ###################
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




## Trial 3 ###################################
### ### Condition 2 & 3: pres #################
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




### ### Condition 2 & 3: faci #################
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





## Trial 4 ###################

### Condition 2 & 3: pres #################
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




### Condition 2 & 3: faci #################
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



### bias comp: pres ###################
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


### bias comp: faci ###################
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




## Trial 5 ###################################
### ### Condition 2 & 3: pres #################
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




### ### Condition 2 & 3: faci #################
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




### bias comp: pres ###################
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


### bias comp: faci ###################
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

