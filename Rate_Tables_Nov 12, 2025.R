library(dplyr)
library(tidyr)

compute_summary_stats <- function(Dat, 
                                  time_var = "time_od_acute", 
                                  event_var = "sens_od_acute",
                                  weight_var = "wt_trunc",
                                  trt_var = "trt",
                                  weighted = TRUE) {
  
  Dat %>%
    group_by(!!sym(trt_var)) %>%
    summarise(
      Event = sum(!!sym(event_var), na.rm = TRUE),
      Episode = n(),
      Total_PY = if (weighted) {
        round(sum(!!sym(weight_var) * !!sym(time_var), na.rm = TRUE) / 365.25, 0)
      } else {
        round(sum(!!sym(time_var), na.rm = TRUE) / 365.25, 0)
      },
      Rate_PY = if (weighted) {
        round(Event / (sum(!!sym(weight_var) * !!sym(time_var), na.rm = TRUE) / 365.25), 2)
      } else {
        round(Event / (sum(!!sym(time_var), na.rm = TRUE) / 365.25), 2)
      },
      Rate_10PY = Rate_PY * 10,
      Rate_100PY = Rate_PY * 100
    ) %>%
    ungroup()
}

#res_HRs %>% count(SA_group)

# Trial 1 ##################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral1_12d_PLR.rds")
## overdose ################
### Primary analysis ####################
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())


### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <80mg of methadone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=80mg of methadone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())

res_HRs %>% count(SA_group)

### Guideline 2023 ####################
hr_wide <- res_HRs %>% filter(SA_group == "8. 2023 Guidelines: 1-3 days gap", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr1.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr1_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr1_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

guide23 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "8. 2023 Guidelines: 1-3 days gap") %>%
  select(SA_group, Design = Design, everything())



### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_12d_rprev_ow.rds")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")

ps_fit$data_weights <- ps_fit$data_weights %>% mutate(
  wt_trunc = wt_trunc + rnorm(n(), mean = 0, sd = 1e-3))

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())



### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))
#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr1_rprev_ow.rds")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")
ps_fit$data_weights <- ps_fit$data_weights %>% mutate(
  wt_trunc = wt_trunc + rnorm(n(), mean = 0, sd = 1e-2))

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())




### Minimum 1-year of cumulative OAT experience ####################
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr1_12d_ow.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
#                                         time_var = "time_od_acute", 
#                                         event_var = "sens_od_acute",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>% 
#   mutate(Design = rep("Incident users", 3))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. 1-year or longer cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())




### All-acute care visits ####################
hr_wide <- res_HRs %>% filter(SA_group == "15. Outcome: All-cause acute care visits", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

all_acute <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "15. Outcome: All-cause acute care visits") %>%
  select(SA_group, Design = Design, everything())



### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_od_acute",
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_od_acute", 
                                       event_var = "s_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_od_acute", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())



### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())


### Restricted to post fentanyl era ####################
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())


### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())



### Comp outcome: od related ####################
hr_wide <- res_HRs %>% filter(SA_group == "16. Outcome: Overdose-related acute care \nvisits or overdose-related death", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_comp_outcome",
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_comp_outcome", 
                                       event_var = "sens_od_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_comp_outcome", 
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "16. Outcome: Overdose-related acute care visits or overdose-related related death") %>%
  select(SA_group, Design = Design, everything())

### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())


### IV estimation #########

hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())




SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, guide23, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  all_acute, 
  #discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design) 


write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral1_12d_od_PLR_ow.csv")


## OAT discontinuation ################
### Primary analysis ####################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral1_12d_PLR.rds")
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)



# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%  
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())



### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr1_12d_ow.rds")

dim(ps_fit$data_weights)

StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <80mg of methadone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=80mg of methadone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())


### Guideline 2023 ####################
hr_wide <- res_HRs %>% filter(SA_group == "8. 2023 Guidelines: 1-3 days gap", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr1.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr1_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr1_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

guide23 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "8. 2023 Guidelines: 1-3 days gap") %>%
  select(SA_group, Design = Design, everything())



### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr1_12d_rprev_ow.rds")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")
ps_fit$data_weights <- ps_fit$data_weights %>% mutate(
  wt_trunc = wt_trunc + rnorm(n(), mean = 0, sd = 1e-2))

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())


### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr1_rprev_ow.rds")
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")
ps_fit$data_weights <- ps_fit$data_weights %>% mutate(
  wt_trunc = wt_trunc + rnorm(n(), mean = 0, sd = 1e-2))


StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())



### Minimum 1-year of cumulative OAT experience ####################
#res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder/HRs_Tiral1_PLR.rds")
# table(res_HRs$SA_group)
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr1_12d_ow.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
#                                         time_var = "time_discont", 
#                                         event_var = "sens_discont",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>% 
#   mutate(Design = rep("Incident users", 3))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. 1-year or longer cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())


### 30d discontinuation ####################
hr_wide <- res_HRs %>% filter(SA_group == "14. Outcome: No OAT dispensations for 30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont30",
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont30", 
                                       event_var = "sens_discont30",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont30", 
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

discont_30d <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "14. Outcome: No OAT dispensations for 30 days") %>%
  select(SA_group, Design = Design, everything())


### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_discont",
                                        event_var = "s_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_discont", 
                                       event_var = "s_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_discont", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())



### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())


### Restricted to post fentanyl era ####################
#res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder/HRs_Tiral1_PLR.rds")
# table(res_HRs$SA_group)
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())



### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())



### Comp outcome: dis related ####################
hr_wide <- res_HRs %>% filter(SA_group == "17. Outcome: Overdose-related acute \ncare visits or OAT discontinuation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_comp_outcome",
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_comp_outcome", 
                                       event_var = "sens_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_comp_outcome", 
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "17. Outcome: Overdose-related acute care visits or OAT discontinuation") %>%
  select(SA_group, Design = Design, everything())


### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())




### IV ####################
hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont",
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_12d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())





SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, guide23, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  #all_acute, 
  discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design) 


write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral1_12d_dis_PLR_ow.csv")

# Trial 2 ##################
## overdose ################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral2_34d_PLR.rds")
res_HRs %>% count(SA_group)
### Primary analysis ####################
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())

### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <80mg of methadone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=80mg of methadone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())


### Guideline 2023 ####################
hr_wide <- res_HRs %>% filter(SA_group == "8. 2023 Guidelines: 4 days gap", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr2.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr2_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr2_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

guide23 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "8. 2023 Guidelines: 4 days gap") %>%
  select(SA_group, Design = Design, everything())



### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_34d_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_34d_rprev_ow.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())

### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr2_34d_rprev_ow.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())



### Minimum 1-year of cumulative OAT experience ####################
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr2_34d.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
#                                         time_var = "time_od_acute", 
#                                         event_var = "sens_od_acute",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>% 
#   mutate(Design = rep("Incident users", 4))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. 1-year or longer cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())




### All-acute care visits ####################
hr_wide <- res_HRs %>% filter(SA_group == "15. Outcome: All-cause acute care visits", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_all_acute",
                                        event_var = "sens_all_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_all_acute", 
                                       event_var = "sens_all_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_all_acute", 
                                        event_var = "sens_all_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

all_acute <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "15. Outcome: All-cause acute care visits") %>%
  select(SA_group, Design = Design, everything())



### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_od_acute",
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_od_acute", 
                                       event_var = "s_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_od_acute", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())



### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())

### Restricted to post fentanyl era ####################
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())


### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())



### Comp outcome: od related ####################
hr_wide <- res_HRs %>% filter(SA_group == "16. Outcome: Overdose-related acute care \nvisits or overdose-related death", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_comp_outcome",
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_comp_outcome", 
                                       event_var = "sens_od_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_comp_outcome", 
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "16. Outcome: Overdose-related acute care visits or overdose-related related death") %>%
  select(SA_group, Design = Design, everything())


### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())




### IV estimation #########

hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())




SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, guide23, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  all_acute, 
  #discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design) 

write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral2_34d_od_PLR_ow.csv")


## discontinuation ################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral2_34d_PLR.rds")

### Primary analysis ####################
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())

### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <80mg of methadone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=80mg of methadone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())


### Guideline 2023 ####################
hr_wide <- res_HRs %>% filter(SA_group == "8. 2023 Guidelines: 4 days gap", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr2.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr2_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_robust/ps_weight_tr2_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

guide23 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "8. 2023 Guidelines: 4 days gap") %>%
  select(SA_group, Design = Design, everything())



### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_34d_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr2_34d_rprev_ow.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())

### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr2_34d_rprev_ow.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())



### Minimum 1-year of cumulative OAT experience ####################
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr2_34d.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
#                                         time_var = "time_discont", 
#                                         event_var = "sens_discont",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>% 
#   mutate(Design = rep("Incident users", 4))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. 1-year or longer cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())




### 30d discontinuation ####################
hr_wide <- res_HRs %>% filter(SA_group == "14. Outcome: No OAT dispensations for 30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont30",
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont30", 
                                       event_var = "sens_discont30",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont30", 
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

discont_30d <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "14. Outcome: No OAT dispensations for 30 days") %>%
  select(SA_group, Design = Design, everything())




### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_od_acute",
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_od_acute", 
                                       event_var = "s_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_od_acute", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())



### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())

### Restricted to post fentanyl era ####################
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())


### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())



### Comp outcome: od related ####################
hr_wide <- res_HRs %>% filter(SA_group == "17. Outcome: Overdose-related acute \ncare visits or OAT discontinuation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_comp_outcome",
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_comp_outcome", 
                                       event_var = "sens_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_comp_outcome", 
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "17. Outcome: Overdose-related acute care visits or OAT discontinuation") %>%
  select(SA_group, Design = Design, everything())


### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())




### IV estimation #########

hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_34d_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())




SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, guide23, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  #all_acute, 
  discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design) 

write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral2_34d_dis_PLR_ow.csv")




# Trial 3 ##################
## overdose ################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral3_PLR.rds")
res_HRs %>% count(SA_group)
### Primary analysis ####################
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())



### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <80mg of methadone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=80mg of methadone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())


### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())

### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr3_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())



### Minimum 1-year of cumulative OAT experience ####################
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr3.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
#                                         time_var = "time_od_acute", 
#                                         event_var = "sens_od_acute",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>% 
#   mutate(Design = rep("Incident users", 4))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. 1-year or longer cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())




### All-acute care visits ####################
hr_wide <- res_HRs %>% filter(SA_group == "15. Outcome: All-cause acute care visits", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

all_acute <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "15. Outcome: All-cause acute care visits") %>%
  select(SA_group, Design = Design, everything())



### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_od_acute",
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_od_acute", 
                                       event_var = "s_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_od_acute", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())



### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())

### Restricted to post fentanyl era ####################
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())



### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())



### Comp outcome: od related ####################
hr_wide <- res_HRs %>% filter(SA_group == "16. Outcome: Overdose-related acute care \nvisits and overdose-related death", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_comp_outcome",
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_comp_outcome", 
                                       event_var = "sens_od_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_comp_outcome", 
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "16. Outcome: Overdose-related acute care visits or overdose-related related death") %>%
  select(SA_group, Design = Design, everything())

### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())



### IV estimation #########
hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())





SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  all_acute, 
  #discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design)  

write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral3_od_PLR_ow.csv")


res_HRs %>% count(SA_group)
## OAT discontinuation ################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral3_PLR.rds")

### Primary analysis ####################
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%  
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())



### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <80mg of methadone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=80mg of methadone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=80mg of methadone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())


### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr3_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())


### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr3_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())



### Minimum 1-year of cumulative OAT experience ####################
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr3.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
#                                         time_var = "time_discont",
#                                         event_var = "sens_discont",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>%
#   mutate(Design = rep("Incident users", 4))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. 1-year or longer cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())


### 30d discontinuation ####################
hr_wide <- res_HRs %>% filter(SA_group == "14. Outcome: No OAT dispensations for 30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont30",
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont30", 
                                       event_var = "sens_discont30",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont30", 
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

discont_30d <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "14. Outcome: No OAT dispensations for 30 days") %>%
  select(SA_group, Design = Design, everything())


### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_discont",
                                        event_var = "s_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_discont", 
                                       event_var = "s_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_discont", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())



### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())


### Restricted to post fentanyl era ####################
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())


### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())



### Comp outcome: dis related ####################
hr_wide <- res_HRs %>% filter(SA_group == "17. Outcome: Overdose-related acute \ncare visits and OAT discontinuation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_comp_outcome",
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr3_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_comp_outcome", 
                                       event_var = "sens_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_comp_outcome", 
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "17. Outcome: Overdose-related acute care visits or OAT discontinuation") %>%
  select(SA_group, Design = Design, everything())


### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())




### IV ####################
hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_prev.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont",
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())




SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  #all_acute, 
  discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design) 


write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral3_dis_PLR_ow.csv")


res_HRs %>% count(SA_group)
# Trial 4 ##################
## overdose ################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral4_PLR.rds")
### Primary analysis ####################
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)



# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())



### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <16/4mg of buprenorphine/naloxone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <16/4mg of buprenorphine/naloxone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=16/4mg of buprenorphine/naloxone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=16/4mg of buprenorphine/naloxone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())




### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr4_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())


### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr4_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())




### Minimum 1-year of cumulative OAT experience ####################
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr4.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
#                                         time_var = "time_od_acute", 
#                                         event_var = "sens_od_acute",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>% 
#   mutate(Design = rep("Incident users", 3))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. Minimum 1-year of cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())




### All-acute care visits ####################
hr_wide <- res_HRs %>% filter(SA_group == "15. Outcome: All-cause acute care visits", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

all_acute <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "15. Outcome: All-cause acute care visits") %>%
  select(SA_group, Design = Design, everything())



### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_od_acute",
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_od_acute", 
                                       event_var = "s_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_od_acute", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())



### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())


### Restricted to post fentanyl era ####################
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(from April 1, 2012 - December 31, 2022)", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())


### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())



### Comp outcome: od related ####################
hr_wide <- res_HRs %>% filter(SA_group == "16. Outcome: Overdose-related acute care \nvisits and overdose-related death", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_comp_outcome",
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_comp_outcome", 
                                       event_var = "sens_od_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_comp_outcome", 
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "16. Outcome: Overdose-related acute care visits or overdose-related related death") %>%
  select(SA_group, Design = Design, everything())


### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())



### IV estimation #########
hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())



SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  all_acute, 
  #discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design) 

write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral4_od_PLR_ow.csv")



## OAT discontinuation ################
### Primary analysis ####################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral4_PLR.rds")
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%  
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())



### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <16/4mg of buprenorphine/naloxone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <16/4mg of buprenorphine/naloxone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=16/4mg of buprenorphine/naloxone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=16/4mg of buprenorphine/naloxone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())


### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())



### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr4_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())


### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr4_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())



### Minimum 1-year of cumulative OAT experience ####################
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr4.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
#                                         time_var = "time_discont", 
#                                         event_var = "sens_discont",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>% 
#   mutate(Design = rep("Incident users", 3))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. 1-year or longer cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())

### 30d discontinuation ####################
hr_wide <- res_HRs %>% filter(SA_group == "14. Outcome: No OAT dispensations for 30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont30",
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont30", 
                                       event_var = "sens_discont30",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont30", 
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_discont <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "14. Outcome: No OAT dispensations for 30 days") %>%
  select(SA_group, Design = Design, everything())



### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_discont",
                                        event_var = "s_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_discont", 
                                       event_var = "s_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_discont", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())




### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())


### Restricted to post fentanyl era ####################
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(from April 1, 2012 - December 31, 2022)", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())


### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())


### Comp outcome: dis related ####################
hr_wide <- res_HRs %>% filter(SA_group == "17. Outcome: Overdose-related acute \ncare visits and OAT discontinuation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr4.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_comp_outcome",
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr4_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_comp_outcome", 
                                       event_var = "sens_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr4_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_comp_outcome", 
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "17. Outcome: Overdose-related acute care visits or OAT discontinuation") %>%
  select(SA_group, Design = Design, everything())

### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())





### IV ####################
hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont",
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 3))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 3))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())




SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  #all_acute, 
  d30_discont,
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design) 


write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral4_dis_PLR_ow.csv")




# Trial 5 ##################
## overdose ################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral5_PLR.rds")

### Primary analysis ####################
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev_ow.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())



### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <16/4mg of buprenorphine/naloxone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <16/4mg of buprenorphine/naloxone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=16/4mg of buprenorphine/naloxone", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=16/4mg of buprenorphine/naloxone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr5_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr5_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr5_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())


### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr5_prev.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())

### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr5_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())



### Minimum 1-year of cumulative OAT experience ####################
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
#ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr5.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
#                                         time_var = "time_od_acute", 
#                                         event_var = "sens_od_acute",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>% 
#   mutate(Design = rep("Incident users", 4))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. 1-year or longer cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())




### All-acute care visits ####################
hr_wide <- res_HRs %>% filter(SA_group == "15. Outcome: All-cause acute care visits", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_outcome/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

all_acute <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "15. Outcome: All-cause acute care visits") %>%
  select(SA_group, Design = Design, everything())



### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_od_acute",
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_od_acute", 
                                       event_var = "s_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_od_acute", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())





### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())

### Restricted to post fentanyl era ####################
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())


### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())



### Comp outcome: od related ####################
hr_wide <- res_HRs %>% filter(SA_group == "16. Outcome: Overdose-related acute care \nvisits and overdose-related death", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_comp_outcome",
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_comp_outcome", 
                                       event_var = "sens_od_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_comp_outcome", 
                                        event_var = "sens_od_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "16. Outcome: Overdose-related acute care visits or overdose-related related death") %>%
  select(SA_group, Design = Design, everything())


### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 3))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 3))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())



### IV estimation #########
hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "Overdose-related acute care visits") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_od_acute",
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_od_acute", 
                                       event_var = "sens_od_acute",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_od_acute", 
                                        event_var = "sens_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())



SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  all_acute, 
  #discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design) 

#write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder2/Rates_Tiral5_od_PLR.csv")



write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral5_od_PLR_ow.csv")



## OAT discontinuation ################
res_HRs <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral5_PLR.rds")
### Primary analysis ####################
hr_wide <- res_HRs %>% filter(SA_group == "1. Primary analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)



# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%  
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

full <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "1. Primary analysis") %>%
  select(SA_group, Design = Design, everything())



### Completion induction ####################
hr_wide <- res_HRs %>% filter(SA_group == "2. Completed induction", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_induction/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ind <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "2. Completed induction") %>%
  select(SA_group, Design = Design, everything())


### Missed low dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "3. Missed <16/4mg of buprenorphine/naloxone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_lt80_16/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

low_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "3. Missed <16/4mg of buprenorphine/naloxone") %>%
  select(SA_group, Design = Design, everything())


### Missed high dose ####################
hr_wide <- res_HRs %>% filter(SA_group == "4. Missed >=16/4mg of buprenorphine/naloxone", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_result_dose_ge80_16/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_dose <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "4. Missed >=16/4mg of buprenorphine/naloxone") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: <=30 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "5. Days since OAT initiation: <=30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr5_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_lt30d/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_le30 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "5. Days since OAT initiation: <=30 days") %>%
  select(SA_group, Design = Design, everything())



### Days since OAT initiation: 31-90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr5_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_30_lt90d/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt30_le90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "6. Days since OAT initiation: 31-90 days") %>%
  select(SA_group, Design = Design, everything())


### Days since OAT initiation: >90 days ####################
hr_wide <- res_HRs %>% filter(SA_group == "7. Days since OAT initiation: >90 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

#StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr5_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ret_gt90d/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

ret_gt90 <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "7. Days since OAT initiation: >90 days") %>%
  select(SA_group, Design = Design, everything())


### hdPS ####################
hr_wide <- res_HRs %>% filter(SA_group == "18. hdPS analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr5_prev.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_wt/ps_weight_tr5_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "18. hdPS analysis") %>%
  select(SA_group, Design = Design, everything())

### hdPS with ML ####################
hr_wide <- res_HRs %>% filter(SA_group == "19. hdPS with machine learning algorithm", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_hdPS_ML/ps_weight_tr5_rprev.rds")

StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

hdPS_ML <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "19. hdPS with machine learning algorithm") %>%
  select(SA_group, Design = Design, everything())


### Minimum 1-year of cumulative OAT experience ####################
hr_wide <- res_HRs %>% filter(SA_group == "11. 1-year or longer cumulative OAT experience", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
# ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr5.rds")
# StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
#                                         time_var = "time_discont",
#                                         event_var = "sens_discont",
#                                         weight_var = "wt_trunc",
#                                         trt_var = "trt")%>%
#   mutate(Design = rep("Incident users", 4))

StatTrial$incid <- NULL

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_retention_365/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

year1_exp <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "11. 1-year or longer cumulative OAT experience") %>%
  select(SA_group, Design = Design, everything())


### 30d discontinuation ####################
hr_wide <- res_HRs %>% filter(SA_group == "14. Outcome: No OAT dispensations for 30 days", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont30",
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont30", 
                                       event_var = "sens_discont30",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_discont_30/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont30", 
                                        event_var = "sens_discont30",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

discont_30d <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "14. Outcome: No OAT dispensations for 30 days") %>%
  select(SA_group, Design = Design, everything())


### No truncation at 12 months of follow-up ####################
hr_wide <- res_HRs %>% filter(SA_group == "12. No truncation at 12 months of follow-up", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "t_discont",
                                        event_var = "s_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "t_discont", 
                                       event_var = "s_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_no_truncation/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "t_discont", 
                                        event_var = "s_od_acute",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

no_trunc <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "12. No truncation at 12 months of follow-up") %>%
  select(SA_group, Design = Design, everything())



### Restricted to high fentanyl area ####################
hr_wide <- res_HRs %>% filter(SA_group == "9. Restricted to high fentanyl area", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_high_fent_area/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

high_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "9. Restricted to high fentanyl area") %>%
  select(SA_group, Design = Design, everything())

### Restricted to post fentanyl era ####################
hr_wide <- res_HRs %>% filter(SA_group == "10. Restricted to post fentanyl era \n(04/01/2012-12/31/2022)", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/post_fentanyl/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)


post_fentanyl <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "10. Restricted to post fentanyl era (04/01/2012-12/31/2022)") %>%
  select(SA_group, Design = Design, everything())



### 30 days washout for prevalent new users ####################
hr_wide <- res_HRs %>% filter(SA_group == "13. 30 days washout for prevalent new users", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_30d_washout/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

d30_washout <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "13. 30 days washout for prevalent new users") %>%
  select(SA_group, Design = Design, everything())



### Comp outcome: dis related ####################
hr_wide <- res_HRs %>% filter(SA_group == "17. Outcome: Overdose-related acute \ncare visits and OAT discontinuation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr5.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_comp_outcome",
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt")%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr5_prev.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_comp_outcome", 
                                       event_var = "sens_comp_outcome",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_comp_outcome/ps_weight_tr5_rprev.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_comp_outcome", 
                                        event_var = "sens_comp_outcome",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

comp_outcome <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "17. Outcome: Overdose-related acute care visits or OAT discontinuation") %>%
  select(SA_group, Design = Design, everything())

### Overlap weights ###################
hr_wide <- res_HRs %>% filter(SA_group == "20. Overlap weights based analysis", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt")%>% 
  mutate(Design = rep("Incident users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_prev_ow.rds")
StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont", 
                                       event_var = "sens_discont",
                                       weight_var = "ow_wt",
                                       trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new users", 4))

ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont", 
                                        event_var = "sens_discont",
                                        weight_var = "ow_wt",
                                        trt_var = "trt") %>% 
  mutate(Design = rep("Prevalent new* users", 4))



rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

OW <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "20. Overlap weights based analysis") %>%
  select(SA_group, Design = Design, everything())




### IV ####################
hr_wide <- res_HRs %>% filter(SA_group == "21. Prescriber's preference based IV estimation", 
                              Model == "IPTW",
                              Outcome == "OAT discontinuation") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  ) %>%
  select(Design, Comparison, HR_CI) %>%
  pivot_wider(names_from = Comparison, values_from = HR_CI)


# Read rates
StatTrial=list(incid=list(), prev=list(), rzero=list())
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_ow.rds")
StatTrial$incid = compute_summary_stats(Dat=ps_fit$data_weights,
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE)%>%
  mutate(Design = rep("Incident users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_prev_ow.rds")

StatTrial$prev = compute_summary_stats(Dat=ps_fit$data_weights, 
                                       time_var = "time_discont",
                                       event_var = "sens_discont",
                                       weight_var = "wt_trunc",
                                       trt_var = "trt",
                                       weighted = FALSE) %>% 
  mutate(Design = rep("Prevalent new users", 4))


ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev_ow.rds")
StatTrial$rzero = compute_summary_stats(Dat=ps_fit$data_weights, 
                                        time_var = "time_discont",
                                        event_var = "sens_discont",
                                        weight_var = "wt_trunc",
                                        trt_var = "trt",
                                        weighted = FALSE) %>%
  mutate(Design = rep("Prevalent new* users", 4))

rate_wide <- bind_rows(StatTrial) %>%
  select(Design, trt, Event, Episode, Total_PY, Rate_PY) %>%
  pivot_longer(cols = c(Event, Episode, Total_PY, Rate_PY),
               names_to = "Metric", values_to = "Value") %>%
  unite("trt_metric", trt, Metric, sep = "_") %>%
  pivot_wider(names_from = trt_metric, values_from = Value)

IV <- rate_wide %>%
  left_join(hr_wide, by = "Design") %>%
  mutate(SA_group = "21. Prescriber's preference based IV estimation") %>%
  select(SA_group, Design = Design, everything())







SA_rate <- bind_rows(
  full, ind, low_dose, high_dose, ret_le30, ret_gt30_le90, ret_gt90, 
  high_fentanyl, post_fentanyl,
  year1_exp, no_trunc, d30_washout, 
  #all_acute, 
  discont_30d, 
  comp_outcome, hdPS, hdPS_ML, OW, IV
) %>% mutate(
  SA_group = factor(SA_group),
  SA_group = factor(
    SA_group,
    levels = rev(  # reverse the order here
      unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    )
  )
) %>% arrange(Design) 



write.csv(SA_rate, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral5_dis_PLR_ow.csv")








## Hiding low frequency #######
tr1 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral1_12d_dis_PLR_ow.csv")
tr2 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral2_34d_dis_PLR_ow.csv")
tr3 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral3_dis_PLR_ow.csv")
tr4 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral4_dis_PLR_ow.csv")
tr5 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral5_dis_PLR_ow.csv")

mask_small_counts <- function(df) {
  # Columns with event counts
  event_cols <- grep("_Event$", names(df), value = TRUE)
  episode_cols <- grep("_Episode$", names(df), value = TRUE)
  rate_cols <- grep("_Rate_PY$", names(df), value = TRUE)
  
  for (i in seq_along(event_cols)) {
    e_col <- event_cols[i]
    ep_col <- episode_cols[i]
    rate_col <- rate_cols[i]
    
    small_rows <- df[[e_col]] < 10 
    df[[e_col]][small_rows] <- "<10"
    df[[rate_col]][small_rows] <- "-"
    
    small_rows_epi <- df[[ep_col]] < 10
    df[[ep_col]][small_rows_epi] <- "<10"
    df[[rate_col]][small_rows_epi] <- "-"
    
  }
  
  return(df)
}

# Apply to all trial datasets
tr1 <- mask_small_counts(tr1)
tr2 <- mask_small_counts(tr2)
tr3 <- mask_small_counts(tr3)
tr4 <- mask_small_counts(tr4)
tr5 <- mask_small_counts(tr5)

write.csv(tr1, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral1_12d_dis_PLR_ow_v1.csv")
write.csv(tr2, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral2_34d_dis_PLR_ow_v1.csv")
write.csv(tr3, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral3_dis_PLR_ow_v1.csv")
write.csv(tr4, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral4_dis_PLR_ow_v1.csv")
write.csv(tr5, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral5_dis_PLR_ow_v1.csv")






## Hiding low frequency #######
tr1 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral1_12d_od_PLR_ow.csv")
tr2 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral2_34d_od_PLR_ow.csv")
tr3 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral3_od_PLR_ow.csv")
tr4 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral4_od_PLR_ow.csv")
tr5 <- read.csv(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral5_od_PLR_ow.csv")


# Apply to all trial datasets
tr1 <- mask_small_counts(tr1)
tr2 <- mask_small_counts(tr2)
tr3 <- mask_small_counts(tr3)
tr4 <- mask_small_counts(tr4)
tr5 <- mask_small_counts(tr5)

write.csv(tr1, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral1_12d_od_PLR_ow_v1.csv")
write.csv(tr2, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral2_34d_od_PLR_ow_v1.csv")
write.csv(tr3, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral3_od_PLR_ow_v1.csv")
write.csv(tr4, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral4_od_PLR_ow_v1.csv")
write.csv(tr5, file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/Rates_Tiral5_od_PLR_ow_v1.csv")

