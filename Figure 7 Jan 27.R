# Figure 8 #################

res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction of <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction of >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction of >=50%/initial dose* vs no change")
  ) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design)


res_df2 %>% count(res_type)

### od: Less than 30 days, iptw ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"


#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_od_red <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 1: MET interruption of 1-2 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr)


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 2: MET interruption of 3-4 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) 


#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) 


#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_od_red <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 4: BNX interruption of 1-5 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr)

#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) 


#### met forestplot ##########################
lt30d_met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                                    tr2_full_iptw_od,
                                    tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))


lt30d_met_full_iptw_od$`OD-related acute care visits` <- "" 
lt30d_met_full_iptw_od$`HR(95% CI)` <- lt30d_met_full_iptw_od$HR_display
lt30d_met_full_iptw_od$`Dose adjustment after interruption` <- lt30d_met_full_iptw_od$Comparison
lt30d_met_full_iptw_od$`User design` <- lt30d_met_full_iptw_od$Design

length(lt30d_met_full_iptw_od$HR)

#### bnx forestplot ##########################
lt30d_bnx_full_iptw_od <- bind_rows(tr4_full_iptw_od,
                                    tr5_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))




lt30d_bnx_full_iptw_od$`OD-related acute care visits` <- "" 
lt30d_bnx_full_iptw_od$`HR(95% CI)` <- lt30d_bnx_full_iptw_od$HR_display
lt30d_bnx_full_iptw_od$`Dose adjustment after interruption` <- lt30d_bnx_full_iptw_od$Comparison
lt30d_bnx_full_iptw_od$`User design` <- lt30d_bnx_full_iptw_od$Design


### dis: Less than 30 days, iptw ######

#### Trial 1: dis: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr1_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 1: MET interruption of 1-2 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr)

#### Trial 2: dis: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) 

#### Trial 3: dis: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) 

#### Trial 4: dis: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 4: BNX interruption of 1-5 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr)

#### Trial 5: dis: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) 



#### met forestplot ##########################
lt30d_met_full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                                     tr2_full_iptw_dis,
                                     tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))


lt30d_met_full_iptw_dis$`OAT discontinuation` <- "" 
lt30d_met_full_iptw_dis$`HR(95% CI)` <- lt30d_met_full_iptw_dis$HR_display
lt30d_met_full_iptw_dis$`Dose adjustment after interruption` <- lt30d_met_full_iptw_dis$Comparison
lt30d_met_full_iptw_dis$`User design` <- lt30d_met_full_iptw_dis$Design


#### bnx forestplot ##########################
lt30d_bnx_full_iptw_dis <- bind_rows(tr4_full_iptw_dis,
                                     tr5_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))


lt30d_bnx_full_iptw_dis$`OAT discontinuation` <- "" 
lt30d_bnx_full_iptw_dis$`HR(95% CI)` <- lt30d_bnx_full_iptw_dis$HR_display
lt30d_bnx_full_iptw_dis$`Dose adjustment after interruption` <- lt30d_bnx_full_iptw_dis$Comparison
lt30d_bnx_full_iptw_dis$`User design` <- lt30d_bnx_full_iptw_dis$Design


### od: 31 to 90 days, iptw ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_30_lt90d/"

res_df2 %>% count(res_type)

#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_od_red <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 1: MET interruption of 1-2 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) 


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 2: MET interruption of 3-4 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) 

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) 


#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_od_red <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 4: BNX interruption of 1-5 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) 

#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) 

#### met forestplot ##########################
f31to90_met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                                      tr2_full_iptw_od,
                                      tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))


f31to90_met_full_iptw_od$`OD-related acute care visits` <- "" 
f31to90_met_full_iptw_od$`HR(95% CI)` <- f31to90_met_full_iptw_od$HR_display
f31to90_met_full_iptw_od$`Dose adjustment after interruption` <- f31to90_met_full_iptw_od$Comparison
f31to90_met_full_iptw_od$`User design` <- f31to90_met_full_iptw_od$Design


#### bnx forestplot ##########################
f31to90_bnx_full_iptw_od <- bind_rows(tr4_full_iptw_od,
                                      tr5_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))


f31to90_bnx_full_iptw_od$`OD-related acute care visits` <- "" 
f31to90_bnx_full_iptw_od$`HR(95% CI)` <- f31to90_bnx_full_iptw_od$HR_display
f31to90_bnx_full_iptw_od$`Dose adjustment after interruption` <- f31to90_bnx_full_iptw_od$Comparison
f31to90_bnx_full_iptw_od$`User design` <- f31to90_bnx_full_iptw_od$Design


### dis: 31 to 90 days, iptw ######

#### Trial 1: dis: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr1_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 1: MET interruption of 1-2 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) 

#### Trial 2: dis: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr)


#### Trial 3: dis: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr)

#### Trial 4: dis: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 4: BNX interruption of 1-5 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr)

#### Trial 5: dis: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) 




#### met forestplot ##########################
f31to90_met_full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                                       tr2_full_iptw_dis,
                                       tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))


f31to90_met_full_iptw_dis$`OAT discontinuation` <- "" 
f31to90_met_full_iptw_dis$`HR(95% CI)` <- f31to90_met_full_iptw_dis$HR_display
f31to90_met_full_iptw_dis$`Dose adjustment after interruption` <- f31to90_met_full_iptw_dis$Comparison
f31to90_met_full_iptw_dis$`User design` <- f31to90_met_full_iptw_dis$Design


#### bnx forestplot ##########################
f31to90_bnx_full_iptw_dis <- bind_rows(tr4_full_iptw_dis,
                                       tr5_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )%>% filter(Design %in% c("", "All prevalent users"))


f31to90_bnx_full_iptw_dis$`OAT discontinuation` <- "" 
f31to90_bnx_full_iptw_dis$`HR(95% CI)` <- f31to90_bnx_full_iptw_dis$HR_display
f31to90_bnx_full_iptw_dis$`Dose adjustment after interruption` <- f31to90_bnx_full_iptw_dis$Comparison
f31to90_bnx_full_iptw_dis$`User design` <- f31to90_bnx_full_iptw_dis$Design



### od: after 90 days, iptw ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_gt90d/"

res_df2 %>% count(res_type)

#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_od_red <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 1: MET interruption of 1-2 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_iptw_od_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) 

#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 2: MET interruption of 3-4 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) 

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) 

#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_od_red <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 4: BNX interruption of 1-5 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_od_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr)

#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) 




#### met forestplot ##########################
af90d_met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                                    tr2_full_iptw_od,
                                    tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))


af90d_met_full_iptw_od$`OD-related acute care visits` <- "" 
af90d_met_full_iptw_od$`HR(95% CI)` <- af90d_met_full_iptw_od$HR_display
af90d_met_full_iptw_od$`Dose adjustment after interruption` <- af90d_met_full_iptw_od$Comparison
af90d_met_full_iptw_od$`User design` <- af90d_met_full_iptw_od$Design


#### bnx forestplot ##########################
af90d_bnx_full_iptw_od <- bind_rows(tr4_full_iptw_od,
                                    tr5_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))


af90d_bnx_full_iptw_od$`OD-related acute care visits` <- "" 
af90d_bnx_full_iptw_od$`HR(95% CI)` <- af90d_bnx_full_iptw_od$HR_display
af90d_bnx_full_iptw_od$`Dose adjustment after interruption` <- af90d_bnx_full_iptw_od$Comparison
af90d_bnx_full_iptw_od$`User design` <- af90d_bnx_full_iptw_od$Design


### dis: after 90 days, iptw ######

#### Trial 1: dis: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_dis_red <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr1_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 1: MET interruption of 1-2 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) 

#### Trial 2: dis: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) 

#### Trial 3: dis: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) 


#### Trial 4: dis: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_dis_red <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 4: BNX interruption of 1-5 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) 

#### Trial 5: dis: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction of <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction of >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) 



#### met forestplot ##########################
af90d_met_full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                                     tr2_full_iptw_dis,
                                     tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  ) %>% filter(Design %in% c("", "All prevalent users"))


af90d_met_full_iptw_dis$`OAT discontinuation` <- "" 
af90d_met_full_iptw_dis$`HR(95% CI)` <- af90d_met_full_iptw_dis$HR_display
af90d_met_full_iptw_dis$`Dose adjustment after interruption` <- af90d_met_full_iptw_dis$Comparison
af90d_met_full_iptw_dis$`User design` <- af90d_met_full_iptw_dis$Design

#### bnx forestplot ##########################
af90d_bnx_full_iptw_dis <- bind_rows(tr4_full_iptw_dis,
                                     tr5_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )%>% filter(Design %in% c("", "All prevalent users"))


af90d_bnx_full_iptw_dis$`OAT discontinuation` <- "" 
af90d_bnx_full_iptw_dis$`HR(95% CI)` <- af90d_bnx_full_iptw_dis$HR_display
af90d_bnx_full_iptw_dis$`Dose adjustment after interruption` <- af90d_bnx_full_iptw_dis$Comparison
af90d_bnx_full_iptw_dis$`User design` <- af90d_bnx_full_iptw_dis$Design





#### Creating forestplot ###############

library(forestploter)
tm_met <- forest_theme(ci_lty = "solid",
                       ci_pch = c(17, 15, 19),
                       base_size = 14,
                       #ci_pch = 15,
                       ci_lwd = 1.8,
                       ci_col = c("springgreen3", "black", "black"),
                       legend_name = NULL,
                       legend_value = c(" Within 30 day of OAT initiation ", 
                                        " Within 31 to 90 day of OAT initiation ", 
                                        " After 90 day of OAT initiation "),
                       legend_position = "none",
                       legend_ncol = 1,
                       refline_gp = gpar(lwd = 1.2, lty = "dashed", col = "red"),
                       footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                                          fontface="plain", col = "black"),
                       arrow_length = 0.06,
                       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
                       xlab_adjust = c("refline", "center"),
                       core = list(padding = unit(c(4, 3), "mm")))

tm <- forest_theme(ci_lty = "solid",
                   ci_pch = c(17, 15, 19),
                   base_size = 14,
                   #ci_pch = 15,
                   ci_lwd = 1.8,
                   ci_col = c("springgreen3", "black", "black"),
                   
                   legend_name = NULL,
                   legend_value = c(" Within 30 day of OAT initiation ", 
                                    " Within 31 to 90 day of OAT initiation ", 
                                    " After 90 day of OAT initiation "),
                   legend_position = "bottom",
                   legend_ncol = 1,
                   refline_gp = gpar(lwd = 1.2, lty = "dashed", col = "red"),
                   footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                                      fontface="plain", col = "black"),
                   arrow_length = 0.06,
                   arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
                   xlab_adjust = c("refline", "center"),
                   core = list(padding = unit(c(4, 3), "mm")))



met_full_iptw <- lt30d_met_full_iptw_od %>% select(`Dose adjustment after interruption`,
                                                   `User design`) %>% 
  mutate(`Dose adjustment after interruption`=recode(`Dose adjustment after interruption`,
                                                     "   Dose reduction of >=50%/initial dose* vs no change"=
                                                       "   Dose reduction of >=50%/initial dose* vs\n   no change"))

met_full_iptw <- met_full_iptw %>%  mutate(aHR_od=ifelse(`User design` %in% c(""), "", 
                                                         paste0(lt30d_met_full_iptw_od$`HR(95% CI)`, "\n",
                                                                f31to90_met_full_iptw_od$`HR(95% CI)`, "\n",
                                                                af90d_met_full_iptw_od$`HR(95% CI)`)),
                                           aHR_dis=ifelse(`User design` %in% c(""), "", 
                                                          paste0(lt30d_met_full_iptw_dis$`HR(95% CI)`, "\n",
                                                                 f31to90_met_full_iptw_dis$`HR(95% CI)`, "\n",
                                                                 af90d_met_full_iptw_dis$`HR(95% CI)`)),
                                           `User design` = paste0(" ", `User design`)
                                           
)

met_full_iptw$`aHR(95% CI)**` <-  met_full_iptw$aHR_od
met_full_iptw$`  aHR(95% CI)` <-  met_full_iptw$aHR_dis
met_full_iptw$`Overdose-related acute care visits` <- ""
met_full_iptw$`OAT discontinuation` <- ""


p_od_dis_met <- forest(
  met_full_iptw[, c(1, 5, 7, 6, 8)],
  est   = list(lt30d_met_full_iptw_od$HR,
               lt30d_met_full_iptw_dis$HR,
               f31to90_met_full_iptw_od$HR,
               f31to90_met_full_iptw_dis$HR,
               af90d_met_full_iptw_od$HR,
               af90d_met_full_iptw_dis$HR),
  lower = list(lt30d_met_full_iptw_od$lower,
               lt30d_met_full_iptw_dis$lower,
               f31to90_met_full_iptw_od$lower,
               f31to90_met_full_iptw_dis$lower, 
               af90d_met_full_iptw_od$lower,
               af90d_met_full_iptw_dis$lower),
  upper = list(lt30d_met_full_iptw_od$upper,
               lt30d_met_full_iptw_dis$upper,
               f31to90_met_full_iptw_od$upper,
               f31to90_met_full_iptw_dis$upper, 
               af90d_met_full_iptw_od$upper,
               af90d_met_full_iptw_dis$upper),
  ci_column = c(3, 5),
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 3),
  ticks_at = c(0.5, 1, 1.5),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "A. Methadone (MET) trials among prevalent new users",
  footnote = NULL,
  nudge_y = 0.3,
  theme = tm_met)
# p_od_dis_met <- add_grob(p_od_dis_met, 
#                      row = 12, 
#                      col = 1:2,
#                      order = "background",
#                      gb_fn = gridtext::richtext_grob,
#                      text = txt,
#                      gp = gpar(fontsize = 11),
#                      hjust = 0, vjust = 1, halign = 0, valign = 1,
#                      x = unit(0, "npc"), y = unit(1, "npc"))


# 
# folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"
# name = "All_prev_full_od_dis_iptw_Nov12"
# figure=paste0(folder, name, ".png")
# ppi = 300
# png(figure, width=3800, height=2800, res=300)
# p_od_dis_met
# dev.off()
# 



bnx_full_iptw <- lt30d_bnx_full_iptw_od %>% select(`Dose adjustment after interruption`,
                                                   `User design`) %>% 
  mutate(`Dose adjustment after interruption`=recode(`Dose adjustment after interruption`,
                                                     "   Dose reduction of >=50%/initial dose* vs no change"=
                                                       "   Dose reduction of >=50%/initial dose* vs\n   no change"))

bnx_full_iptw <- bnx_full_iptw %>%  mutate(aHR_od=ifelse(`User design` %in% c(""), "", 
                                                         paste0(lt30d_bnx_full_iptw_od$`HR(95% CI)`, "\n",
                                                                f31to90_bnx_full_iptw_od$`HR(95% CI)`, "\n",
                                                                af90d_bnx_full_iptw_od$`HR(95% CI)`)),
                                           aHR_dis=ifelse(`User design` %in% c(""), "", 
                                                          paste0(lt30d_bnx_full_iptw_dis$`HR(95% CI)`, "\n",
                                                                 f31to90_bnx_full_iptw_dis$`HR(95% CI)`, "\n",
                                                                 af90d_bnx_full_iptw_dis$`HR(95% CI)`)),
                                           `User design` = paste0(" ", `User design`)
                                           
)

bnx_full_iptw$`aHR(95% CI)` <-  bnx_full_iptw$aHR_od
bnx_full_iptw$`  aHR(95% CI)` <-  bnx_full_iptw$aHR_dis
bnx_full_iptw$`Overdose-related acute care visits` <- ""
bnx_full_iptw$`OAT discontinuation` <- ""

names(bnx_full_iptw)
p_od_dis_bnx <- forest(
  bnx_full_iptw[, c(1, 5, 7, 6, 8)],
  est   = list(lt30d_bnx_full_iptw_od$HR,
               lt30d_bnx_full_iptw_dis$HR,
               f31to90_bnx_full_iptw_od$HR,
               f31to90_bnx_full_iptw_dis$HR,
               af90d_bnx_full_iptw_od$HR,
               af90d_bnx_full_iptw_dis$HR),
  lower = list(lt30d_bnx_full_iptw_od$lower,
               lt30d_bnx_full_iptw_dis$lower,
               f31to90_bnx_full_iptw_od$lower,
               f31to90_bnx_full_iptw_dis$lower, 
               af90d_bnx_full_iptw_od$lower,
               af90d_bnx_full_iptw_dis$lower),
  upper = list(lt30d_bnx_full_iptw_od$upper,
               lt30d_bnx_full_iptw_dis$upper,
               f31to90_bnx_full_iptw_od$upper,
               f31to90_bnx_full_iptw_dis$upper, 
               af90d_bnx_full_iptw_od$upper,
               af90d_bnx_full_iptw_dis$upper),
  ci_column = c(3, 5),
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 3),
  ticks_at = c(0.5, 1, 1.5),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials among prevalent new users",
  footnote = NULL,
  nudge_y = 0.3,
  theme = tm)


txt <- "<sup>*</sup> <=40mg MET for Trial 2, 3 and <=4/1mg BNX for Trial 5<br><sup>**</sup> adjusted hazard ratio with 95% compatibility intervals"

p_od_dis_bnx <- add_grob(p_od_dis_bnx, 
                         row = 8, 
                         col = 1:2,
                         order = "background",
                         gb_fn = gridtext::richtext_grob,
                         text = txt,
                         gp = gpar(fontsize = 11),
                         hjust = 0, vjust = 1, halign = 0, valign = 1,
                         x = unit(0, "npc"), y = unit(1, "npc"))




library(patchwork)
plot_od_dis <- plot_spacer() /  # top spacer row
  p_od_dis_met / 
  plot_spacer() /
  p_od_dis_bnx/ plot_spacer() +
  plot_layout(
    heights = c(0.18, 1, 0.12, 1, 0.01),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"
name = "Figure8_iptw_Jan27_MK"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3700, height=4700, res=300)
plot_od_dis
dev.off()

