# Condition 1 #############
## Trial 1 #################
### predicted prob : pres #####
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"
dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_logistic(data = dat1, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_logistic(data = dat1_prev, 
                                          trt_var = "trt", 
                                          iv_var = "iv_var", 
                                          spline_df = 3)


dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_logistic(data = dat1_rprev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <-fit_iv_logistic(data = dat1, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat1_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat1_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr1.rds")


### predicted prob : faci #####
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"
dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_logistic(data = dat1, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_logistic(data = dat1_prev, 
                                          trt_var = "trt", 
                                          iv_var = "iv_var", 
                                          spline_df = 3)


dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_logistic(data = dat1_rprev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <-fit_iv_logistic(data = dat1, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat1_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat1_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr1.rds")



### Condition 1: pres #################
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"


dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_F_lmm(data = dat1,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_F_lmm(data = dat1_prev,
                                       trt_var = "trt",
                                       iv_var = "iv_var",
                                       baselinevars = baselinevars,
                                       prescriber_id = "prac_studyid_rev",
                                       spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_F_lmm(data = dat1_rprev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = c(baselinevars, "prior_gap"),
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat1,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat1_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat1_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)
saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr1.rds")

### Condition 1: faci #################
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"


dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_F_lmm(data = dat1,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_F_lmm(data = dat1_prev,
                                       trt_var = "trt",
                                       iv_var = "iv_var",
                                       baselinevars = baselinevars,
                                       prescriber_id = "prac_studyid_rev",
                                       spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_F_lmm(data = dat1_rprev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = c(baselinevars, "prior_gap"),
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat1,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat1_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat1_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr1.rds")



## Trial 2 ################

### predicted prob: pres #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_logistic(data = dat2, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_logistic(data = dat2_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_logistic(data = dat2_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)

#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <- fit_iv_logistic(data = dat2, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_logistic(data = dat2_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_logistic(data = dat2_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_logistic(data = dat2, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <-  fit_iv_logistic(data = dat2_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <-  fit_iv_logistic(data = dat2_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr2.rds")

### predicted prob: faci #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <-  fit_iv_logistic(data = dat2, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <-  fit_iv_logistic(data = dat2_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <-  fit_iv_logistic(data = dat2_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)

#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <-  fit_iv_logistic(data = dat2, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <-  fit_iv_logistic(data = dat2_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <-  fit_iv_logistic(data = dat2_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <-  fit_iv_logistic(data = dat2, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat2_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat2_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr2.rds")




### Condition 1 #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

#### Reduction to >=50%/<=4mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr2.rds")

### Condition 1: faci #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

#### Reduction to >=50%/<=4mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr2.rds")

## Trial 3 ###################


### predicted prob: pres #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_logistic(data = dat3, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_logistic(data = dat3_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_logistic(data = dat3_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)

#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <- fit_iv_logistic(data = dat3, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_logistic(data = dat3_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_logistic(data = dat3_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_logistic(data = dat3, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)

dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <-  fit_iv_logistic(data = dat3_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <-  fit_iv_logistic(data = dat3_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr3.rds")

### predicted prob: faci #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <-  fit_iv_logistic(data = dat3, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <-  fit_iv_logistic(data = dat3_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <-  fit_iv_logistic(data = dat3_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)

#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <-  fit_iv_logistic(data = dat3, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <-  fit_iv_logistic(data = dat3_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <-  fit_iv_logistic(data = dat3_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <-  fit_iv_logistic(data = dat3, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)

dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat3_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat3_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr3.rds")



### Condition 1: pres #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_F_lmm(data = dat3,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_F_lmm(data = dat3_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_F_lmm(data = dat3_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <- fit_iv_F_lmm(data = dat3,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_F_lmm(data = dat3_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_F_lmm(data = dat3_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat3,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat3_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat3_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr3.rds")


### Condition 1: faci #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_F_lmm(data = dat3,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_F_lmm(data = dat3_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_F_lmm(data = dat3_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <- fit_iv_F_lmm(data = dat3,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_F_lmm(data = dat3_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_F_lmm(data = dat3_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat3,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat3_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat3_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr3.rds")

## Trial 4 #######################


### predicted prob : pres #####
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"
dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_logistic(data = dat4, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_logistic(data = dat4_prev, 
                                          trt_var = "trt", 
                                          iv_var = "iv_var", 
                                          spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_logistic(data = dat4_rprev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <-fit_iv_logistic(data = dat4, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat4_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat4_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr4.rds")


### predicted prob : faci #####
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"
dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_logistic(data = dat4, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_logistic(data = dat4_prev, 
                                          trt_var = "trt", 
                                          iv_var = "iv_var", 
                                          spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_logistic(data = dat4_rprev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <-fit_iv_logistic(data = dat4, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat4_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat4_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr4.rds")




### Condition 1: pres #####
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"
dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_F_lmm(data = dat4,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_F_lmm(data = dat4_prev,
                                       trt_var = "trt",
                                       iv_var = "iv_var",
                                       baselinevars = baselinevars,
                                       prescriber_id = "prac_studyid_rev",
                                       spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_F_lmm(data = dat4_rprev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = c(baselinevars, "prior_gap"),
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat4,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat4_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat4_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr4.rds")


### Condition 1: faci #####
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"
dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_F_lmm(data = dat4,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_F_lmm(data = dat4_prev,
                                       trt_var = "trt",
                                       iv_var = "iv_var",
                                       baselinevars = baselinevars,
                                       prescriber_id = "prac_studyid_rev",
                                       spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_F_lmm(data = dat4_rprev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = c(baselinevars, "prior_gap"),
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat4,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat4_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat4_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr4.rds")

## Trial 5 ##########################

### predicted prob: pres #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_logistic(data = dat5, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_logistic(data = dat5_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_logistic(data = dat5_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)

#### Reduction to >=50%/<=4mg ###########
trt_strategy = ">=50%/<=4mg"
iv_name = "grp3_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

table(dat5$trt)

res_valid$grp_ge50$incid <- fit_iv_logistic(data = dat5, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_logistic(data = dat5_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_logistic(data = dat5_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_logistic(data = dat5, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)

dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <-  fit_iv_logistic(data = dat5_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <-  fit_iv_logistic(data = dat5_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr5.rds")



### predicted prob: faci #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <-  fit_iv_logistic(data = dat5, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <-  fit_iv_logistic(data = dat5_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <-  fit_iv_logistic(data = dat5_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)

#### Reduction to >=50%/<=4mg ###########
trt_strategy = ">=50%/<=4mg"
iv_name = "grp3_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <-  fit_iv_logistic(data = dat5, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <-  fit_iv_logistic(data = dat5_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <-  fit_iv_logistic(data = dat5_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <-  fit_iv_logistic(data = dat5, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)

dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat5_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat5_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr5.rds")






### Condition 1: pres #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_F_lmm(data = dat5,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_F_lmm(data = dat5_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_F_lmm(data = dat5_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

#### Reduction to >=50%/<=4mg ###########
trt_strategy = ">=50%/<=4mg"
iv_name = "grp3_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_ge50$incid <- fit_iv_F_lmm(data = dat5,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_F_lmm(data = dat5_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_F_lmm(data = dat5_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat5,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat5_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat5_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr5.rds")


### Condition 1: faci #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_F_lmm(data = dat5,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_F_lmm(data = dat5_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_F_lmm(data = dat5_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

#### Reduction to >=50%/<=4mg ###########
trt_strategy = ">=50%/<=4mg"
iv_name = "grp3_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_ge50$incid <- fit_iv_F_lmm(data = dat5,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_F_lmm(data = dat5_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_F_lmm(data = dat5_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat5,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat5_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat5_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr5.rds")


## Trial 1: 1-2 days #################
### predicted prob : pres #####
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"
dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_logistic(data = dat1, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_logistic(data = dat1_prev, 
                                          trt_var = "trt", 
                                          iv_var = "iv_var", 
                                          spline_df = 3)


dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_logistic(data = dat1_rprev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <- fit_iv_logistic(data = dat1, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat1_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat1_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr1_12d_ow.rds")


### predicted prob : faci #####
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"
dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_logistic(data = dat1, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_logistic(data = dat1_prev, 
                                          trt_var = "trt", 
                                          iv_var = "iv_var", 
                                          spline_df = 3)


dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_logistic(data = dat1_rprev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <-fit_iv_logistic(data = dat1, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat1_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


dat1_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1_12d_ow.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat1_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr1_12d_ow.rds")



### Condition 1: pres #################
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"


dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_F_lmm(data = dat1,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_F_lmm(data = dat1_prev,
                                       trt_var = "trt",
                                       iv_var = "iv_var",
                                       baselinevars = baselinevars,
                                       prescriber_id = "prac_studyid_rev",
                                       spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_F_lmm(data = dat1_rprev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = c(baselinevars, "prior_gap"),
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat1,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat1_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat1_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)
saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr1_12d_ow.rds")

### Condition 1: faci #################
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
res_valid = list(grp_red=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"


dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_red$incid <- fit_iv_F_lmm(data = dat1,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$prev <- fit_iv_F_lmm(data = dat1_prev,
                                       trt_var = "trt",
                                       iv_var = "iv_var",
                                       baselinevars = baselinevars,
                                       prescriber_id = "prac_studyid_rev",
                                       spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_red$rzero <- fit_iv_F_lmm(data = dat1_rprev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = c(baselinevars, "prior_gap"),
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat1,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat1_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat1_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr1_12d_ow.rds")



## Trial 2: 3-4 days ################

### predicted prob: pres #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_logistic(data = dat2, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_logistic(data = dat2_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_logistic(data = dat2_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)

#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <- fit_iv_logistic(data = dat2, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_logistic(data = dat2_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_logistic(data = dat2_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_logistic(data = dat2, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <-  fit_iv_logistic(data = dat2_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <-  fit_iv_logistic(data = dat2_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr2_34d_ow.rds")

### predicted prob: faci #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <-  fit_iv_logistic(data = dat2, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <-  fit_iv_logistic(data = dat2_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <-  fit_iv_logistic(data = dat2_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)

#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <-  fit_iv_logistic(data = dat2, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <-  fit_iv_logistic(data = dat2_prev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <-  fit_iv_logistic(data = dat2_rprev, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <-  fit_iv_logistic(data = dat2, 
                                             trt_var = "trt", 
                                             iv_var = "iv_var", 
                                             spline_df = 3)

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_logistic(data = dat2_prev, 
                                           trt_var = "trt", 
                                           iv_var = "iv_var", 
                                           spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_logistic(data = dat2_rprev, 
                                            trt_var = "trt", 
                                            iv_var = "iv_var", 
                                            spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr2_34d_ow.rds")




### Condition 1 #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr2_34d_ow.rds")

### Condition 1: faci #####
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
res_valid = list(grp_lt50=list(incid=list(), prev=list(), rzero=list()),
                 grp_ge50=list(incid=list(), prev=list(), rzero=list()),
                 grp_incr=list(incid=list(), prev=list(), rzero=list()))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


res_valid$grp_lt50$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_lt50$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

#### Reduction to >=50%/<=4mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_ge50$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)



#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$incid <- fit_iv_F_lmm(data = dat2,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = baselinevars,
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)

dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$prev <- fit_iv_F_lmm(data = dat2_prev,
                                        trt_var = "trt",
                                        iv_var = "iv_var",
                                        baselinevars = baselinevars,
                                        prescriber_id = "prac_studyid_rev",
                                        spline_df = 3)


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
         prac_studyid_rev = as.factor(prac_studyid_rev))

res_valid$grp_incr$rzero <- fit_iv_F_lmm(data = dat2_rprev,
                                         trt_var = "trt",
                                         iv_var = "iv_var",
                                         baselinevars = c(baselinevars, "prior_gap"),
                                         prescriber_id = "prac_studyid_rev",
                                         spline_df = 3)


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr2_34d_ow.rds")




## histogram ###############
### Trial 1 ############
pres_1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr1.rds")
#### pres #############

pres_incid_1 = plot_iv_logistic(fit = pres_1$grp_red$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_prev_1 = plot_iv_logistic(fit = pres_1$grp_red$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1, 
                               plot_title = "Prevalent new users")

pres_rprev_1 = plot_iv_logistic(fit = pres_1$grp_red$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "All prevalent users")



pres_incid_2 = plot_iv_logistic(fit = pres_1$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_2 = plot_iv_logistic(fit = pres_1$grp_incr$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_1$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr1.rds")

faci_incid_1 = plot_iv_logistic(fit = faci_1$grp_red$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_1 = plot_iv_logistic(fit = faci_1$grp_red$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_1$grp_red$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_1$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_2 = plot_iv_logistic(fit = faci_1$grp_incr$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_1$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_1 + plot_spacer() + pres_prev_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row2 <- pres_incid_2 + plot_spacer() + pres_prev_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))


row3 <- faci_incid_1 + plot_spacer() + faci_prev_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))


row4 <- faci_incid_2 + plot_spacer() + faci_prev_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))


plot_1 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer() + 
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.16, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_1s <- ggdraw(plot_1) + 
  draw_label("A. Prescriber's facicribing preference for dose reduction versus no dose change after methadone interruption of 1-3 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose increase versus no dose change after methadone interruption of 1-3 days",  
             x = 0.02, y = 0.73, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Facility-level prescribing preference for dose reduction versus no dose change after methadone interruption of 1-3 days", 
             x = 0.02, y = 0.44, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("D. Facility-level prescribing preference for dose increase versus no dose change after methadone interruption of 1-3 days",  
             x = 0.02, y = 0.20, hjust = 0, vjust = 0, 
             size = 18, color = "black")  


name = "histogram_1"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=6000, res=300)
plot_1s
dev.off()


### Trial 2 ############
pres_2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr2.rds")
#### pres #############
pres_incid_0 = plot_iv_logistic(fit = pres_2$grp_lt50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_prev_0 = plot_iv_logistic(fit = pres_2$grp_lt50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1, 
                               plot_title = "Prevalent new users")

pres_rprev_0 = plot_iv_logistic(fit = pres_2$grp_lt50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2), 
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "All prevalent users")

pres_incid_1 = plot_iv_logistic(fit = pres_2$grp_ge50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_1 = plot_iv_logistic(fit = pres_2$grp_ge50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_1 = plot_iv_logistic(fit = pres_2$grp_ge50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



pres_incid_2 = plot_iv_logistic(fit = pres_2$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_2 = plot_iv_logistic(fit = pres_2$grp_incr$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_2$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr2.rds")


faci_incid_0 = plot_iv_logistic(fit = faci_2$grp_lt50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_0 = plot_iv_logistic(fit = faci_2$grp_lt50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_0 = plot_iv_logistic(fit = faci_2$grp_lt50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


faci_incid_1 = plot_iv_logistic(fit = faci_2$grp_ge50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_1 = plot_iv_logistic(fit = faci_2$grp_ge50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_2$grp_ge50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_2$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_2 = plot_iv_logistic(fit = faci_2$grp_incr$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_2$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_0 + plot_spacer() + pres_prev_0 + plot_spacer() + pres_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row2 <- pres_incid_1 + plot_spacer() + pres_prev_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row3 <- pres_incid_2 + plot_spacer() + pres_prev_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row4 <- faci_incid_0 + plot_spacer() + faci_prev_0 + plot_spacer() + faci_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row5 <- faci_incid_1 + plot_spacer() + faci_prev_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row6 <- faci_incid_2 + plot_spacer() + faci_prev_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

plot_2 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer()/
  row5/
  plot_spacer()/
  row6/
  plot_spacer() +
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_2s <- ggdraw(plot_2) + 
  draw_label("A. Prescriber's facicribing preference for dose reduction to <50% versus no dose change after methadone interruption of 4 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose reduction to >=50%/initial dose* versus no dose change after methadone interruption of 4 days",  
             x = 0.02, y = 0.815, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Prescriber's facicribing preference for dose increase versus no dose change after methadone interruption of 4 days",  
             x = 0.02, y = 0.656, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("D. Facility-level prescribing preference for dose reduction to <50% versus no dose change after methadone interruption of 4 days", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("E. Facility-level prescribing preference for dose reduction to >=50%/initial dose* versus no dose change after methadone interruption of 4 days",  
             x = 0.02, y = 0.325, hjust = 0, vjust = 0, 
             size = 18, color = "black") + 
  draw_label("F. Facility-level prescribing preference for dose increase versus no dose change after methadone interruption of 4 days",  
             x = 0.02, y = 0.168, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


name = "histogram_2"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7000, res=300)
plot_2s
dev.off()


### Trial 3 ############
pres_3 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr3.rds")
#### pres #############
pres_incid_0 = plot_iv_logistic(fit = pres_3$grp_lt50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_prev_0 = plot_iv_logistic(fit = pres_3$grp_lt50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1, 
                               plot_title = "Prevalent new users")

pres_rprev_0 = plot_iv_logistic(fit = pres_3$grp_lt50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2), 
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "All prevalent users")

pres_incid_1 = plot_iv_logistic(fit = pres_3$grp_ge50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_1 = plot_iv_logistic(fit = pres_3$grp_ge50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_1 = plot_iv_logistic(fit = pres_3$grp_ge50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



pres_incid_2 = plot_iv_logistic(fit = pres_3$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_2 = plot_iv_logistic(fit = pres_3$grp_incr$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_3$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_3 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr3.rds")


faci_incid_0 = plot_iv_logistic(fit = faci_3$grp_lt50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_0 = plot_iv_logistic(fit = faci_3$grp_lt50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_0 = plot_iv_logistic(fit = faci_3$grp_lt50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


faci_incid_1 = plot_iv_logistic(fit = faci_3$grp_ge50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_1 = plot_iv_logistic(fit = faci_3$grp_ge50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_3$grp_ge50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_3$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_2 = plot_iv_logistic(fit = faci_3$grp_incr$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_3$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_0 + plot_spacer() + pres_prev_0 + plot_spacer() + pres_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row2 <- pres_incid_1 + plot_spacer() + pres_prev_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row3 <- pres_incid_2 + plot_spacer() + pres_prev_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row4 <- faci_incid_0 + plot_spacer() + faci_prev_0 + plot_spacer() + faci_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row5 <- faci_incid_1 + plot_spacer() + faci_prev_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row6 <- faci_incid_2 + plot_spacer() + faci_prev_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

plot_3 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer()/
  row5/
  plot_spacer()/
  row6/
  plot_spacer() +
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_3s <- ggdraw(plot_3) + 
  draw_label("A. Prescriber's facicribing preference for dose reduction to <50% versus no dose change after methadone interruption of 5-14 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose reduction to >=50%/initial dose* versus no dose change after methadone interruption of 5-14 days",  
             x = 0.02, y = 0.815, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Prescriber's facicribing preference for dose increase versus no dose change after methadone interruption of 5-14 days",  
             x = 0.02, y = 0.652, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("D. Facility-level prescribing preference for dose reduction to <50% versus no dose change after methadone interruption of 5-14 days", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("E. Facility-level prescribing preference for dose reduction to >=50%/initial dose* versus no dose change after methadone interruption of 5-14 days",  
             x = 0.02, y = 0.325, hjust = 0, vjust = 0, 
             size = 18, color = "black") + 
  draw_label("F. Facility-level prescribing preference for dose increase versus no dose change after methadone interruption of 5-14 days",  
             x = 0.02, y = 0.165, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


name = "histogram_3"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7000, res=300)
plot_3s
dev.off()



### Trial 4 ############
pres_4 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr4.rds")
#### pres #############

pres_incid_1 = plot_iv_logistic(fit = pres_4$grp_red$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_prev_1 = plot_iv_logistic(fit = pres_4$grp_red$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1, 
                               plot_title = "Prevalent new users")

pres_rprev_1 = plot_iv_logistic(fit = pres_4$grp_red$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "All prevalent users")



pres_incid_2 = plot_iv_logistic(fit = pres_4$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_2 = plot_iv_logistic(fit = pres_4$grp_incr$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_4$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_4 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr4.rds")

faci_incid_1 = plot_iv_logistic(fit = faci_4$grp_red$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_1 = plot_iv_logistic(fit = faci_4$grp_red$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_4$grp_red$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_4$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_2 = plot_iv_logistic(fit = faci_4$grp_incr$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_4$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_1 + plot_spacer() + pres_prev_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row2 <- pres_incid_2 + plot_spacer() + pres_prev_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))


row3 <- faci_incid_1 + plot_spacer() + faci_prev_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))


row4 <- faci_incid_2 + plot_spacer() + faci_prev_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))


plot_4 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer() + 
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_4s <- ggdraw(plot_4) + 
  draw_label("A. Prescriber's facicribing preference for dose reduction versus no dose change after buprenorphine/naloxone interruption of 1-5 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose increase versus no dose change after buprenorphine/naloxone interruption of 1-5 days",  
             x = 0.02, y = 0.73, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Facility-level prescribing preference for dose reduction versus no dose change after buprenorphine/naloxone interruption of 1-5 days", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("D. Facility-level prescribing preference for dose increase versus no dose change after buprenorphine/naloxone interruption of 1-5 days",  
             x = 0.02, y = 0.24, hjust = 0, vjust = 0, 
             size = 18, color = "black")  


name = "histogram_4"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=6000, res=300)
plot_4s
dev.off()



### Trial 5 ############
pres_5 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr5.rds")
#### pres #############
pres_incid_0 = plot_iv_logistic(fit = pres_5$grp_lt50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_prev_0 = plot_iv_logistic(fit = pres_5$grp_lt50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1, 
                               plot_title = "Prevalent new users")

pres_rprev_0 = plot_iv_logistic(fit = pres_5$grp_lt50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "All prevalent users")

pres_incid_1 = plot_iv_logistic(fit = pres_5$grp_ge50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_1 = plot_iv_logistic(fit = pres_5$grp_ge50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_1 = plot_iv_logistic(fit = pres_5$grp_ge50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



pres_incid_2 = plot_iv_logistic(fit = pres_5$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_2 = plot_iv_logistic(fit = pres_5$grp_incr$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_5$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_5 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr5.rds")


faci_incid_0 = plot_iv_logistic(fit = faci_5$grp_lt50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_0 = plot_iv_logistic(fit = faci_5$grp_lt50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_0 = plot_iv_logistic(fit = faci_5$grp_lt50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


faci_incid_1 = plot_iv_logistic(fit = faci_5$grp_ge50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_1 = plot_iv_logistic(fit = faci_5$grp_ge50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_5$grp_ge50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_5$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_2 = plot_iv_logistic(fit = faci_5$grp_incr$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_5$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_0 + plot_spacer() + pres_prev_0 + plot_spacer() + pres_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row2 <- pres_incid_1 + plot_spacer() + pres_prev_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row3 <- pres_incid_2 + plot_spacer() + pres_prev_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row4 <- faci_incid_0 + plot_spacer() + faci_prev_0 + plot_spacer() + faci_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row5 <- faci_incid_1 + plot_spacer() + faci_prev_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row6 <- faci_incid_2 + plot_spacer() + faci_prev_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

plot_5 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer()/
  row5/
  plot_spacer()/
  row6/
  plot_spacer() +
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_5s <- ggdraw(plot_5) + 
  draw_label("A. Prescriber's facicribing preference for dose reduction to <50% versus no dose change after buprenorphine/naloxone interruption of 6-14 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose reduction to >=50%/initial dose* versus no dose change after buprenorphine/naloxone interruption of 6-14 days",  
             x = 0.02, y = 0.815, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Prescriber's facicribing preference for dose increase versus no dose change after buprenorphine/naloxone interruption of 6-14 days",  
             x = 0.02, y = 0.656, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("D. Facility-level prescribing preference for dose reduction to <50% versus no dose change after buprenorphine/naloxone interruption of 6-14 days", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("E. Facility-level prescribing preference for dose reduction to >=50%/initial dose* versus no dose change after buprenorphine/naloxone interruption of 6-14 days",  
             x = 0.02, y = 0.325, hjust = 0, vjust = 0, 
             size = 18, color = "black") + 
  draw_label("F. Facility-level prescribing preference for dose increase versus no dose change after buprenorphine/naloxone interruption of 6-14 days",  
             x = 0.02, y = 0.168, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


name = "histogram_5"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7000, res=300)
plot_5s
dev.off()

### Trial 1: 1-2 days ############
pres_1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr1_12d_ow.rds")
#### pres #############

pres_incid_1 = plot_iv_logistic(fit = pres_1$grp_red$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_prev_1 = plot_iv_logistic(fit = pres_1$grp_red$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1, 
                               plot_title = "Prevalent new users")

pres_rprev_1 = plot_iv_logistic(fit = pres_1$grp_red$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "All prevalent users")



pres_incid_2 = plot_iv_logistic(fit = pres_1$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_2 = plot_iv_logistic(fit = pres_1$grp_incr$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_1$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr1_12d_ow.rds")

faci_incid_1 = plot_iv_logistic(fit = faci_1$grp_red$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_1 = plot_iv_logistic(fit = faci_1$grp_red$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_1$grp_red$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_1$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_2 = plot_iv_logistic(fit = faci_1$grp_incr$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_1$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_1 + plot_spacer() + pres_prev_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row2 <- pres_incid_2 + plot_spacer() + pres_prev_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))


row3 <- faci_incid_1 + plot_spacer() + faci_prev_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))


row4 <- faci_incid_2 + plot_spacer() + faci_prev_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))


plot_1 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer() + 
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.16, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_1s <- ggdraw(plot_1) + 
  draw_label("A. Prescriber's facicribing preference for dose reduction versus no dose change after methadone interruption of 1-2 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose increase versus no dose change after methadone interruption of 1-2 days",  
             x = 0.02, y = 0.73, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Facility-level prescribing preference for dose reduction versus no dose change after methadone interruption of 1-2 days", 
             x = 0.02, y = 0.48, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("D. Facility-level prescribing preference for dose increase versus no dose change after methadone interruption of 1-2 days",  
             x = 0.02, y = 0.24, hjust = 0, vjust = 0, 
             size = 18, color = "black")  




name = "histogram_1_12d"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=6000, res=300)
plot_1s
dev.off()


### Trial 2: 3-4 days ############
pres_2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr2_34d_ow.rds")
#### pres #############
pres_incid_0 = plot_iv_logistic(fit = pres_2$grp_lt50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_prev_0 = plot_iv_logistic(fit = pres_2$grp_lt50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1, 
                               plot_title = "Prevalent new users")

pres_rprev_0 = plot_iv_logistic(fit = pres_2$grp_lt50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2), 
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "All prevalent users")

pres_incid_1 = plot_iv_logistic(fit = pres_2$grp_ge50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_1 = plot_iv_logistic(fit = pres_2$grp_ge50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_1 = plot_iv_logistic(fit = pres_2$grp_ge50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



pres_incid_2 = plot_iv_logistic(fit = pres_2$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_2 = plot_iv_logistic(fit = pres_2$grp_incr$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_2$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr2_34d_ow.rds")


faci_incid_0 = plot_iv_logistic(fit = faci_2$grp_lt50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_0 = plot_iv_logistic(fit = faci_2$grp_lt50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_0 = plot_iv_logistic(fit = faci_2$grp_lt50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


faci_incid_1 = plot_iv_logistic(fit = faci_2$grp_ge50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_1 = plot_iv_logistic(fit = faci_2$grp_ge50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_2$grp_ge50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_2$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_2 = plot_iv_logistic(fit = faci_2$grp_incr$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_2$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_0 + plot_spacer() + pres_prev_0 + plot_spacer() + pres_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row2 <- pres_incid_1 + plot_spacer() + pres_prev_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row3 <- pres_incid_2 + plot_spacer() + pres_prev_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row4 <- faci_incid_0 + plot_spacer() + faci_prev_0 + plot_spacer() + faci_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row5 <- faci_incid_1 + plot_spacer() + faci_prev_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row6 <- faci_incid_2 + plot_spacer() + faci_prev_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

plot_2 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer()/
  row5/
  plot_spacer()/
  row6/
  plot_spacer() +
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_2s <- ggdraw(plot_2) + 
  draw_label("A. Prescriber's facicribing preference for dose reduction to <50% versus no dose change after methadone interruption of 3-4 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose reduction to >=50%/initial dose* versus no dose change after methadone interruption of 3-4 days",  
             x = 0.02, y = 0.815, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Prescriber's facicribing preference for dose increase versus no dose change after methadone interruption of 3-4 days",  
             x = 0.02, y = 0.656, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("D. Facility-level prescribing preference for dose reduction to <50% versus no dose change after methadone interruption of 3-4 days", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("E. Facility-level prescribing preference for dose reduction to >=50%/initial dose* versus no dose change after methadone interruption of 3-4 days",  
             x = 0.02, y = 0.325, hjust = 0, vjust = 0, 
             size = 18, color = "black") + 
  draw_label("F. Facility-level prescribing preference for dose increase versus no dose change after methadone interruption of 3-4 days",  
             x = 0.02, y = 0.168, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


name = "histogram_2_34d"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6000, height=7000, res=300)
plot_2s
dev.off()



## summary of F-value ###############
pres_1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr1.rds")
faci_1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr1.rds")

pres_2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr2.rds")
faci_2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr2.rds")

pres_3 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr3.rds")
faci_3 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr3.rds")

pres_4 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr4.rds")
faci_4 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr4.rds")

pres_5 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr5.rds")
faci_5 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr5.rds")


pres_1_12d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr1_12d_ow.rds")
faci_1_12d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr1_12d_ow.rds")

pres_2_34d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/cond1_tr2_34d_ow.rds")
faci_2_34d <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_cond1_tr2_34d_ow.rds")


extract_F_values <- function(res_obj, trial_label = "Trial") {
  designs <- c("incid" = "Incident users", 
               "prev"  = "Prevalent new users", 
               "rzero" = "All prevalent users")
  
  # loop over whatever groups exist in the object
  bind_rows(
    lapply(names(res_obj), function(g) {
      lapply(names(designs), function(d) {
        result <- res_obj[[g]][[d]]
        
        # skip if missing or invalid
        if (is.null(result) || is.null(result$F_value) || 
            is.na(result$F_value) || is.nan(result$F_value)) return(NULL)
        
        data.frame(
          Group   = g,
          Design  = designs[[d]],
          F_value = result$F_value,
          Trial   = trial_label,
          stringsAsFactors = FALSE
        )
      })
    }) |> unlist(recursive = FALSE)
  )
}




F_value_pres_1 <- extract_F_values(res_obj=pres_1, trial_label = "Trial 1")
F_value_pres_2 <- extract_F_values(res_obj=pres_2, trial_label = "Trial 2")
F_value_pres_3 <- extract_F_values(res_obj=pres_3, trial_label = "Trial 3")
F_value_pres_4 <- extract_F_values(res_obj=pres_4, trial_label = "Trial 4")
F_value_pres_5 <- extract_F_values(res_obj=pres_5, trial_label = "Trial 5")
F_value_pres_1_12d <- extract_F_values(res_obj=pres_1_12d, trial_label = "Trial 1_12d")
F_value_pres_2_34d <- extract_F_values(res_obj=pres_2_34d, trial_label = "Trial 2_34d")


F_value_faci_1 <- extract_F_values(res_obj=faci_1, trial_label = "Trial 1")
F_value_faci_2 <- extract_F_values(res_obj=faci_2, trial_label = "Trial 2")
F_value_faci_3 <- extract_F_values(res_obj=faci_3, trial_label = "Trial 3")
F_value_faci_4 <- extract_F_values(res_obj=faci_4, trial_label = "Trial 4")
F_value_faci_5 <- extract_F_values(res_obj=faci_5, trial_label = "Trial 5")
F_value_faci_1_12d <- extract_F_values(res_obj=faci_1_12d, trial_label = "Trial 1_12d")
F_value_faci_2_34d <- extract_F_values(res_obj=faci_2_34d, trial_label = "Trial 2_34d")


pres = bind_rows(F_value_pres_1, F_value_pres_1_12d, F_value_pres_2, F_value_pres_2_34d, F_value_pres_3, F_value_pres_4, F_value_pres_5) %>% 
  mutate(IV_level = "pres")
faci = bind_rows(F_value_faci_1, F_value_faci_2, F_value_faci_1_12d, F_value_faci_2_34d, F_value_faci_3, F_value_faci_4, F_value_faci_5) %>% 
  mutate(IV_level = "faci")

F_value <- bind_rows(pres, faci)
saveRDS(F_value, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/F_value.rds")

F_value <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/F_value.rds")
write.csv(F_value, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/F_value.csv")


