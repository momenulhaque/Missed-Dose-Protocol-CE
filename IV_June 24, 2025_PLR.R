
# HR from PLR ############################
## Trial 1 ################################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")
### IV calculation: pres ############
dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                              oat=iv_cohort_alls_m, 
                                              cond="trial_cat == '1-3'",
                                              trt_name="trt_p23")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_cat == '1-3'",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds")


dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_gap %in% 2:3",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1.rds")


### IV calculation: faci ############
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")

prov_clinics_v2 <- read_sas("W:/SUD/Data/2022/derived/prov_clinics_v2.sas7bdat")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)



faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_cat == '1-3'",
                                        trt_name="trt_p23")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1.rds")



dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)


faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_cat == '1-3'",
                                        trt_name='trt_p23')


saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1.rds")


dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)


faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_gap %in% 2:3",
                                        trt_name='trt_p23')

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1.rds")


### 2nd stage: pres #####
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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
res_valid = list(grp_red=list(od=list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(dis=list(incid=list(), prev=list(), rzero=list())))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"
dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat1, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_red$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_red$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat1_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_red$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                          covs = c(baselinevars, "pearson") , 
                                          cluster_var = "studyid_rev")



long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_red$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")



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

resid_dat <- fit1st_iv_logistic(data = dat1_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)



res_valid$grp_red$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")



long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


res_valid$grp_red$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "prior_gap", "pearson") , 
                                            cluster_var = "studyid_rev")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat1, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat1_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_incr$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")



long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_incr$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")



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

resid_dat <- fit1st_iv_logistic(data = dat1_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)


res_valid$grp_incr$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")




long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


res_valid$grp_incr$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "prior_gap", "pearson") , 
                                             cluster_var = "studyid_rev")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/resTrial1_iv.rds")



## Trial 2 ################################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")
### IV calculation: pres ############
dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_cat == 4",
                                          trt_name="trt_p23")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_cat == 4",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds")



dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_cat == 4",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2.rds")

### IV calculation: faci ############
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")
prov_clinics_v2 <- read_sas("W:/SUD/Data/2022/derived/prov_clinics_v2.sas7bdat")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)



faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_cat == 4",
                                        trt_name="trt_p23")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2.rds")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)

faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_cat == 4",
                                        trt_name='trt_p23')

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2.rds")



dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)

faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_cat == 4",
                                        trt_name='trt_p23')

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2.rds")



### 2nd stage: pres #####
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
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
res_valid = list(grp_lt50=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())),
                 grp_ge50=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


resid_dat <- fit1st_iv_logistic(data = dat2, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_lt50$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_lt50$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat2_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_lt50$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_lt50$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat2_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)


res_valid$grp_lt50$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_lt50$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


resid_dat <- fit1st_iv_logistic(data = dat2, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_ge50$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_ge50$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat2_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_ge50$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_ge50$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat2_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_ge50$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_ge50$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat2, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat2_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat2_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)


res_valid$grp_incr$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


res_valid$grp_incr$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/resTrial2_iv.rds")


## Trial 3 ################################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")

### IV calculation: pres ############
dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_cat == '5-14'",
                                          trt_name="trt_p23")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_cat == '5-14'",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds")



dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_cat == '5-14'",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr3.rds")

### IV calculation: faci ############
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")

prov_clinics_v2 <- read_sas("W:/SUD/Data/2022/derived/prov_clinics_v2.sas7bdat")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)


faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_cat == '5-14'",
                                        trt_name="trt_p23")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr3.rds")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)

faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_cat == '5-14'",
                                        trt_name='trt_p23')

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr3.rds")



dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)

faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_cat == '5-14'",
                                        trt_name='trt_p23')

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr3.rds")


### 2nd stage: pres #####
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
res_valid = list(grp_lt50=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())),
                 grp_ge50=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


resid_dat <- fit1st_iv_logistic(data = dat3, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_lt50$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_lt50$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat3_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_lt50$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_lt50$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat3_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_lt50$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_lt50$dis$rzerp = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


resid_dat <- fit1st_iv_logistic(data = dat3, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_ge50$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_ge50$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat3_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_ge50$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_ge50$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat3_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_ge50$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_ge50$dis$rzerp = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat3 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat3, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat3_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr3.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat3_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat3_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_incr$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_incr$dis$rzerp = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/resTrial3_iv.rds")



## Trial 4 ################################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
iv_cohort_alls_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_b.sas7bdat")
### IV calculation: pres ############
dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt_b(dat=dat_b, 
                                          oat=iv_cohort_alls_b, 
                                          cond="trial_cat == '1-5'",
                                          trt_name="trt_p45")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds")

dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt_b(dat=dat_b, 
                                          oat=iv_cohort_alls_b, 
                                          cond="trial_cat == '1-5'",
                                          trt_name="trt_p45")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds")


dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt_b(dat=dat_b, 
                                          oat=iv_cohort_alls_b, 
                                          cond="trial_cat == '1-5'",
                                          trt_name="trt_p45")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr4.rds")

### IV calculation: faci ############
iv_cohort_alls_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_b.sas7bdat")

prov_clinics_v2 <- read_sas("W:/SUD/Data/2022/derived/prov_clinics_v2.sas7bdat")


dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)


faciiv_iv_dt <- compute_facility_pct_dt_b(dat=dat_b, 
                                        oat=iv_cohort_alls_b, 
                                        cond="trial_cat == '1-5'",
                                        trt_name="trt_p45")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr4.rds")

dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)


faciiv_iv_dt <- compute_facility_pct_dt_b(dat=dat_b, 
                                        oat=iv_cohort_alls_b, 
                                        cond="trial_cat == '1-5'",
                                        trt_name="trt_p45")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr4.rds")


dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)


faciiv_iv_dt <- compute_facility_pct_dt_b(dat=dat_b, 
                                        oat=iv_cohort_alls_b, 
                                        cond="trial_cat == '1-5'",
                                        trt_name="trt_p45")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr4.rds")

### 2nd stage: pres #####
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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
res_valid = list(grp_red=list(od=list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(dis=list(incid=list(), prev=list(), rzero=list())))


trt_strategy = c("1-<50%", ">=50%/<=4mg")
iv_name = "grp2_eps_pct"
dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat4, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_red$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_red$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")



dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat4_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_red$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                          covs = c(baselinevars, "pearson") , 
                                          cluster_var = "studyid_rev")



long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_red$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")



dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat4_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)



res_valid$grp_red$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")



long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


res_valid$grp_red$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "prior_gap", "pearson") , 
                                            cluster_var = "studyid_rev")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat4 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat4, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")



dat4_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr4.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat4_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_incr$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")



long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_incr$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")



dat4_rprev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr4.rds") %>%
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days)
  ) %>%
  ungroup() %>% 
  select(all_of(c(clust_vars, out_vars, baselinevars, "prior_gap", iv_vars, grp_vars,  "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat4_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)


res_valid$grp_incr$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")




long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


res_valid$grp_incr$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "prior_gap", "pearson") , 
                                             cluster_var = "studyid_rev")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/resTrial4_iv.rds")


## Trial 5 ################################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
iv_cohort_alls_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_b.sas7bdat")
### IV calculation: pres ############
dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt_b(dat=dat_b, 
                                          oat=iv_cohort_alls_b, 
                                          cond="trial_cat == '6-14'",
                                          trt_name="trt_p45")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds")

dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt_b(dat=dat_b, 
                                          oat=iv_cohort_alls_b, 
                                          cond="trial_cat == '6-14'",
                                          trt_name="trt_p45")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds")


dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt_b(dat=dat_b, 
                                          oat=iv_cohort_alls_b, 
                                          cond="trial_cat == '6-14'",
                                          trt_name="trt_p45")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_b_tr5.rds")


### IV calculation: faci ############
iv_cohort_alls_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_b.sas7bdat")

prov_clinics_v2 <- read_sas("W:/SUD/Data/2022/derived/prov_clinics_v2.sas7bdat")

dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)



faciiv_iv_dt <- compute_facility_pct_dt_b(dat=dat_b, 
                                        oat=iv_cohort_alls_b, 
                                        cond="trial_cat == '6-14'",
                                        trt_name="trt_p45")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_b_tr5.rds")

dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)



faciiv_iv_dt <- compute_facility_pct_dt_b(dat=dat_b, 
                                        oat=iv_cohort_alls_b, 
                                        cond="trial_cat == '6-14'",
                                        trt_name="trt_p45")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_b_tr5.rds")


dat_b <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  rename(trt_p45 = trt_s5) %>%
  mutate(trt_p45 = factor(trt_p45),
         time0_12m_prior = time0 %m-% months(12),
         trt_p45 = recode(trt_p45, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p45 == "no reversion"),
         grp2 = (trt_p45 %in% c("1-<50%", ">=50%/<=4mg")),
         grp3 = (trt_p45 == "1-<50%"),
         grp4 = (trt_p45 == ">=50%/<=4mg"),
         grp5 = (trt_p45 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)



faciiv_iv_dt <- compute_facility_pct_dt_b(dat=dat_b, 
                                        oat=iv_cohort_alls_b, 
                                        cond="trial_cat == '6-14'",
                                        trt_name="trt_p45")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_b_tr5.rds")



### 2nd stage: pres #####
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
res_valid = list(grp_lt50=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())),
                 grp_ge50=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


resid_dat <- fit1st_iv_logistic(data = dat5, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_lt50$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_lt50$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat5_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_lt50$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_lt50$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat5_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_lt50$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_lt50$dis$rzerp = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


#### Reduction to >=50%/<=4mg ###########
trt_strategy = ">=50%/<=4mg"
iv_name = "grp3_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


resid_dat <- fit1st_iv_logistic(data = dat5, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_ge50$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_ge50$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat5_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_ge50$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_ge50$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat5_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_ge50$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_ge50$dis$rzerp = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat5 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat5, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat5_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_b_tr5.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p45"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p45 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p45 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat5_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat5_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_incr$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_incr$dis$rzerp = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/resTrial5_iv.rds")






## Trial 1: 1-2 days ################################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
#iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")
### IV calculation: pres ############
dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_gap %in% 1:2",
                                          trt_name="trt_p23")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_gap %in% 1:2",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds")


dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_gap %in% 1:2",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr1_12d_ow.rds")


### IV calculation: faci ############
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
#iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")

#prov_clinics_v2 <- read_sas("W:/SUD/Data/2022/derived/prov_clinics_v2.sas7bdat")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)



faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_gap %in% 1:2",
                                        trt_name="trt_p23")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr1_12d_ow.rds")



dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)


faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_gap %in% 1:2",
                                        trt_name='trt_p23')


saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr1_12d_ow.rds")


dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)


faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_gap %in% 1:2",
                                        trt_name='trt_p23')

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr1_12d_ow.rds")


### 2nd stage: pres #####
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")

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
res_valid = list(grp_red=list(od=list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(dis=list(incid=list(), prev=list(), rzero=list())))


trt_strategy = c("1-<50%", ">=50%/<=40mg")
iv_name = "grp2_eps_pct"
dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat1, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_red$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_red$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat1_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_red$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                          covs = c(baselinevars, "pearson") , 
                                          cluster_var = "studyid_rev")



long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_red$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")



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

resid_dat <- fit1st_iv_logistic(data = dat1_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)



res_valid$grp_red$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")



long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


res_valid$grp_red$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "prior_gap", "pearson") , 
                                            cluster_var = "studyid_rev")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat1 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat1, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")



dat1_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr1_12d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars,  "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter( trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse( trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat1_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_incr$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")



long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_incr$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")



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

resid_dat <- fit1st_iv_logistic(data = dat1_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")

long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 30)


res_valid$grp_incr$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")




long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 30)


res_valid$grp_incr$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "prior_gap", "pearson") , 
                                             cluster_var = "studyid_rev")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/resTrial1_12d_iv.rds")





## Trial 2: 3-4 days ################################
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
#iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")
### IV calculation: pres ############
dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 


praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_gap %in% 3:4",
                                          trt_name="trt_p23")

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_gap %in% 3:4",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds")



dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) 

praciv_iv_dt <- compute_prescriber_pct_dt(dat=dat_m, 
                                          oat=iv_cohort_alls_m, 
                                          cond="trial_gap %in% 3:4",
                                          trt_name='trt_p23')

saveRDS(praciv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/praciv_rprev_m_tr2_34d_ow.rds")

### IV calculation: faci ############
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
#iv_cohort_alls_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/iv_cohort_alls_m.sas7bdat")
#prov_clinics_v2 <- read_sas("W:/SUD/Data/2022/derived/prov_clinics_v2.sas7bdat")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_incid_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)



faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_gap %in% 3:4",
                                        trt_name="trt_p23")

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_incid_m_tr2_34d_ow.rds")

dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_prev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)

faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_gap %in% 3:4",
                                        trt_name='trt_p23')

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_prev_m_tr2_34d_ow.rds")



dat_m <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  rename(trt_p23 = trt_s456) %>%
  mutate(trt_p23 = factor(trt_p23),
         time0_12m_prior = time0 %m-% months(12),
         trt_p23 = recode(trt_p23, 
                          "1-<25%" = "1-<50%",
                          "25-<50%" = "1-<50%"),
         incar_grp_bls = (incar_grp_bl == "1_pastyear"),
         grp1 = (trt_p23 == "no reversion"),
         grp2 = (trt_p23 %in% c("1-<50%", ">=50%/<=40mg")),
         grp3 = (trt_p23 == "1-<50%"),
         grp4 = (trt_p23 == ">=50%/<=40mg"),
         grp5 = (trt_p23 == "increase")) %>% 
  mutate(
    year = year(time0),
    month = month(time0)
  ) %>%
  left_join(
    prov_clinics_v2,
    by = c("year" = "year",
           "month" = "month",
           "prac_studyid_rev" = "prac_studyid_rev")
  ) %>%
  arrange(studyid_rev)

faciiv_iv_dt <- compute_facility_pct_dt(dat=dat_m, 
                                        oat=iv_cohort_alls_m, 
                                        cond="trial_gap %in% 3:4",
                                        trt_name='trt_p23')

saveRDS(faciiv_iv_dt, file = "E:/SAS Data/Missed dose protocol_v2/faciiv_rprev_m_tr2_34d_ow.rds")



### 2nd stage: pres #####
source(file = "H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/func_PLR.R")
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
res_valid = list(grp_lt50=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())),
                 grp_ge50=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())),
                 grp_incr=list(od=list(incid=list(), prev=list(), rzero=list()),
                               dis=list(incid=list(), prev=list(), rzero=list())))


trt_strategy = "1-<50%"
iv_name = "grp3_eps_pct"
dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


resid_dat <- fit1st_iv_logistic(data = dat2, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_lt50$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_lt50$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat2_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_lt50$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_lt50$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat2_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_lt50$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_lt50$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


#### Reduction to >=50%/<=40mg ###########
trt_strategy = ">=50%/<=40mg"
iv_name = "grp3_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))


resid_dat <- fit1st_iv_logistic(data = dat2, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_ge50$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_ge50$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat2_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_ge50$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_ge50$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat2_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_ge50$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_ge50$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


#### Dose increase ###########
trt_strategy = "increase"
iv_name = "grp5_eps_pct"

dat2 <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_incid_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat2, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$incid = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


dat2_prev <- readRDS(file = "E:/SAS Data/Missed dose protocol_v2/praciv_prev_m_tr2_34d_ow.rds") %>%
  select(all_of(c(clust_vars, out_vars, baselinevars, iv_vars, grp_vars, "trt_p23"))) %>%
  rename(iv_var = all_of(iv_name)) %>%  
  filter(trt_p23 %in% c("no reversion", trt_strategy)) %>%
  mutate(trt = ifelse(trt_p23 %in% trt_strategy, 1, 0),
         prac_studyid_rev = as.factor(prac_studyid_rev))

resid_dat <- fit1st_iv_logistic(data = dat2_prev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 7)


res_valid$grp_incr$od$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                           covs = c(baselinevars, "pearson") , 
                                           cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 7)


res_valid$grp_incr$dis$prev = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


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

resid_dat <- fit1st_iv_logistic(data = dat2_rprev, 
                                trt_var = "trt", 
                                covs = baselinevars,
                                iv_var = "iv_var")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_od_acute", 
                                event_var="sens_od_acute", 
                                interval_length = 14)


res_valid$grp_incr$od$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                            covs = c(baselinevars, "pearson") , 
                                            cluster_var = "studyid_rev")


long_dat = expand_person_period(Dat=resid_dat$residuals, 
                                time_var="time_discont", 
                                event_var="sens_discont", 
                                interval_length = 14)


res_valid$grp_incr$dis$rzero = iv_fit_PLR_HR(Dat_expanded=long_dat, trt_var = "trt", 
                                             covs = c(baselinevars, "pearson") , 
                                             cluster_var = "studyid_rev")


saveRDS(res_valid, "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/resTrial2_34d_iv.rds")


