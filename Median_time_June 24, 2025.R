compute_summary_stats <- function(Dat, 
                                  time_var = "time_od_acute", 
                                  event_var = "sens_od_acute",
                                  trt_var = "trt") {
  library(dplyr)
  
  Dat %>%
    group_by(!!sym(trt_var)) %>%
    summarise(
      Median_time = round(median(!!sym(time_var), na.rm = TRUE), 0)
    ) %>%
    ungroup()
}

library(survminer)



# Overdose ################
compute_summary_stats <- function(Dat, 
                                  time_var = "time_od_acute", 
                                  event_var = "sens_od_acute",
                                  trt_var = "trt", 
                                  Trial = "Trial 1",
                                  Design = "Incident users") {
  library(dplyr)
  
  Dat %>%
    group_by(!!sym(trt_var)) %>%
    summarise(
      Median_time = round(median(!!sym(time_var), na.rm = TRUE), 0),
      Events_pct  = round(mean(!!sym(event_var) == 1, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    ) %>% 
    mutate(Trial = Trial,
           Design = Design)
}


## Trial 1 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_rprev.rds")$data_weights

med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 1",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 1",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 1",
                                        Design = "All prevalent users")


(tr1 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 2 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_rprev.rds")$data_weights

med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 2",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 2",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 2",
                                        Design = "All prevalent users")

(tr2 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 3 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev.rds")$data_weights
med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 3",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 3",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 3",
                                        Design = "All prevalent users")

(tr3 = bind_rows(med_time, med_time_prev, med_time_rprev))

met = bind_rows(tr1, tr2, tr3) %>% 
  mutate(trt = recode(trt, 
                      "1-<50%" = "reversion",
                      ">=50%/<=40mg" = "reversion")) %>% 
  group_by(trt) %>%
  summarise(
    min = min(Events_pct, na.rm = TRUE),
    max = max(Events_pct, na.rm = TRUE)
  ) %>%
  ungroup()


## Trial 4 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev.rds")$data_weights

med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 4",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 4",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 4",
                                        Design = "All prevalent users")

(tr4 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 5 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev.rds")$data_weights
med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 5",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 5",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 5",
                                        Design = "All prevalent users")

(tr5 = bind_rows(med_time, med_time_prev, med_time_rprev))


bnx = bind_rows(tr4, tr5) %>% 
  mutate(trt = recode(trt, 
                      "1-<50%" = "reversion",
                      ">=50%/<=4mg" = "reversion")) %>% 
  group_by(trt) %>%
  summarise(
    min = min(Events_pct, na.rm = TRUE),
    max = max(Events_pct, na.rm = TRUE)
  ) %>%
  ungroup()



# OAT discontinuation ################
form_discont <- as.formula(paste("Surv(time_discont, sens_discont) ~ trt"))
## Trial 1 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_rprev.rds")$data_weights
fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=reversion"="Dose reduction",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 1") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=reversion"="Dose reduction",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 1") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=reversion"="Dose reduction",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 1") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr1 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 2 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_rprev.rds")$data_weights
fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=1-<50%"="Dose reduction to <50%",
                                                     "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 2") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=1-<50%"="Dose reduction to <50%",
                                                               "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 2") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=1-<50%"="Dose reduction to <50%",
                                                                 "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 2") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr2 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 3 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev.rds")$data_weights
fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=1-<50%"="Dose reduction to <50%",
                                                     "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 3") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=1-<50%"="Dose reduction to <50%",
                                                               "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 3") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=1-<50%"="Dose reduction to <50%",
                                                                 "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 3") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr3 = bind_rows(med_time, med_time_prev, med_time_rprev))


## Trial 4 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev.rds")$data_weights
fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=reversion"="Dose reduction",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 4") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=reversion"="Dose reduction",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 4") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=reversion"="Dose reduction",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 4") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr4 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 5 ####
dat <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5.rds")$data_weights
dat_prev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_prev.rds")$data_weights
dat_rprev <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev.rds")$data_weights
fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=1-<50%"="Dose reduction to <50%",
                                                     "trt=>=50%/<=4mg"="Dose reduction to >=50%/4mg",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 5") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=1-<50%"="Dose reduction to <50%",
                                                               "trt=>=50%/<=4mg"="Dose reduction to >=50%/4mg",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 5") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=1-<50%"="Dose reduction to <50%",
                                                                 "trt=>=50%/<=4mg"="Dose reduction to >=50%/4mg",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 5") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr5 = bind_rows(med_time, med_time_prev, med_time_rprev))


dis_median = bind_rows(tr1, tr2, tr3, tr4, tr5)


dis_median %>% filter(Trial %in% c("Trial 1", "Trial 2", "Trial 3"),
                      trt == "No dose change")

dis_median %>% filter(Trial %in% c("Trial 1", "Trial 2", "Trial 3"),
                      trt != "No dose change", 
                      trt != "Dose increase")

dis_median %>% filter(Trial %in% c("Trial 1", "Trial 2", "Trial 3"),
                      trt == "Dose increase")


dis_median %>% filter(Trial %in% c("Trial 4", "Trial 5"),
                      trt == "No dose change")

dis_median %>% filter(Trial %in% c("Trial 4", "Trial 5"),
                      trt != "No dose change", 
                      trt != "Dose increase")

dis_median %>% filter(Trial %in% c("Trial 4", "Trial 5"),
                      trt == "Dose increase")




dis_median %>% filter(Trial %in% c("Trial 1"))






# Updates on Nov 12, 2025 ##############
## Trial 1: od #############################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))$data_weights


med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 1",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 1",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 1",
                                        Design = "All prevalent users")

(tr1 = bind_rows(med_time, med_time_prev, med_time_rprev))

## Trial 2: od ###############
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))$data_weights


med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 2",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 2",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 2",
                                        Design = "All prevalent users")

(tr2 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 3 ####
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr3_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))$data_weights

med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 3",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 3",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 3",
                                        Design = "All prevalent users")

(tr3 = bind_rows(med_time, med_time_prev, med_time_rprev))

met = bind_rows(tr1, tr2, tr3) %>% 
  mutate(trt = recode(trt, 
                      "1-<50%" = "reversion",
                      ">=50%/<=40mg" = "reversion")) %>% 
  group_by(trt) %>%
  summarise(
    min = min(Events_pct, na.rm = TRUE),
    max = max(Events_pct, na.rm = TRUE)
  ) %>%
  ungroup()


## Trial 4 ####
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr4_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))$data_weights

med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 4",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 4",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 4",
                                        Design = "All prevalent users")

(tr4 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 5 ####
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr5_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))$data_weights


med_time <- compute_summary_stats(dat, 
                                  Trial = "Trial 5",
                                  Design = "Incident users")
med_time_prev <- compute_summary_stats(dat_prev, 
                                       Trial = "Trial 5",
                                       Design = "Prevalent new users")
med_time_rprev <- compute_summary_stats(dat_rprev, 
                                        Trial = "Trial 5",
                                        Design = "All prevalent users")

(tr5 = bind_rows(med_time, med_time_prev, med_time_rprev))


bnx = bind_rows(tr4, tr5) %>% 
  mutate(trt = recode(trt, 
                      "1-<50%" = "reversion",
                      ">=50%/<=4mg" = "reversion")) %>% 
  group_by(trt) %>%
  summarise(
    min = min(Events_pct, na.rm = TRUE),
    max = max(Events_pct, na.rm = TRUE)
  ) %>%
  ungroup()
met


## OAT discontinuation ################
form_discont <- as.formula(paste("Surv(time_discont, sens_discont) ~ trt"))
## Trial 1: dis ####
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))$data_weights

fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=reversion"="Dose reduction",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 1") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=reversion"="Dose reduction",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 1") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=reversion"="Dose reduction",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 1") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr1 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 2: dis ####
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))$data_weights

fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=1-<50%"="Dose reduction to <50%",
                                                     "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 2") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=1-<50%"="Dose reduction to <50%",
                                                               "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 2") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=1-<50%"="Dose reduction to <50%",
                                                                 "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 2") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr2 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 3 ####
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr3_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))$data_weights
fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=1-<50%"="Dose reduction to <50%",
                                                     "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 3") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=1-<50%"="Dose reduction to <50%",
                                                               "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 3") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=1-<50%"="Dose reduction to <50%",
                                                                 "trt=>=50%/<=40mg"="Dose reduction to >=50%/40mg",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 3") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr3 = bind_rows(med_time, med_time_prev, med_time_rprev))


## Trial 4 ####
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr4_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))$data_weights
fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=reversion"="Dose reduction",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 4") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=reversion"="Dose reduction",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 4") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=reversion"="Dose reduction",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 4") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr4 = bind_rows(med_time, med_time_prev, med_time_rprev))



## Trial 5 ####
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr5_ow.rds"))$data_weights
dat_prev <- readRDS(file=paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))$data_weights
dat_rprev <- readRDS(file=paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))$data_weights

fit <- surv_fit(form_discont, data = dat)
med_time <- surv_median(fit) %>% mutate(trt = recode(strata, 
                                                     "trt=no reversion"="No dose change",
                                                     "trt=1-<50%"="Dose reduction to <50%",
                                                     "trt=>=50%/<=4mg"="Dose reduction to >=50%/4mg",
                                                     "trt=increase"="Dose increase"),
                                        Design = "Incident users", 
                                        Trial = "Trial 5") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_prev <- surv_fit(form_discont, data = dat_prev)
med_time_prev <- surv_median(fit_prev) %>% mutate(trt = recode(strata, 
                                                               "trt=no reversion"="No dose change",
                                                               "trt=1-<50%"="Dose reduction to <50%",
                                                               "trt=>=50%/<=4mg"="Dose reduction to >=50%/4mg",
                                                               "trt=increase"="Dose increase"),
                                                  Design = "Prevalent new users", 
                                                  Trial = "Trial 5") %>% 
  select(median, lower, upper, trt, Design, Trial)


fit_rprev <- surv_fit(form_discont, data = dat_rprev)
med_time_rprev <- surv_median(fit_rprev) %>% mutate(trt = recode(strata, 
                                                                 "trt=no reversion"="No dose change",
                                                                 "trt=1-<50%"="Dose reduction to <50%",
                                                                 "trt=>=50%/<=4mg"="Dose reduction to >=50%/4mg",
                                                                 "trt=increase"="Dose increase"),
                                                    Design = "All prevalent users", 
                                                    Trial = "Trial 5") %>% 
  select(median, lower, upper, trt, Design, Trial)


(tr5 = bind_rows(med_time, med_time_prev, med_time_rprev))


dis_median = bind_rows(tr1, tr2, tr3, tr4, tr5)


dis_median %>% filter(Trial %in% c("Trial 1", "Trial 2", "Trial 3"),
                      trt == "No dose change")

dis_median %>% filter(Trial %in% c("Trial 1", "Trial 2", "Trial 3"),
                      trt != "No dose change", 
                      trt != "Dose increase")

dis_median %>% filter(Trial %in% c("Trial 1", "Trial 2", "Trial 3"),
                      trt == "Dose increase")


dis_median %>% filter(Trial %in% c("Trial 4", "Trial 5"),
                      trt == "No dose change")

dis_median %>% filter(Trial %in% c("Trial 4", "Trial 5"),
                      trt != "No dose change", 
                      trt != "Dose increase")

dis_median %>% filter(Trial %in% c("Trial 4", "Trial 5"),
                      trt == "Dose increase")




dis_median %>% filter(Trial %in% c("Trial 1"))

