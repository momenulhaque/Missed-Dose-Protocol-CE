install.packages("interactionRCS")
library(splines)
library(interactionRCS)
library(survival)
library(ggplot2)
library(scales)
library(cobalt)
library(glmtoolbox)
source("H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/R_functions.R")

# Retention time ###########
## Trial 1 #######################
### overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_rprev.rds")
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr1_rprev.rds")


### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr1_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr1_rprev.rds")



## Trial 2 ###############################
### Overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

hr_df <- data_or_hr(
  model = cox_spline_model,
  x = "ret_length",
  by = "trt",
  xvals = seq(0, ret_max, by = 5)
)

# Step 3: Recode contrasts for plotting
hr_df <- hr_df %>%
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    ))
  )



# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=40mg - no reversion" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=40mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr2_rprev.rds")


### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr2_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=40mg - no reversion" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=40mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr2_rprev.rds")



## Trial 3 ###############################
### Overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=40mg - no reversion" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=40mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr3_rprev.rds")

### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr3_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=40mg - no reversion" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=40mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr3_rprev.rds")


## Trial 4 #######################
### overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev.rds")
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr4_rprev.rds")


### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr4_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)


datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

mfx$ret_length <- datagrid_values$ret_length


mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr4_rprev.rds")



## Trial 5 ###############################
### Overdose ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev.rds")

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)

# Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)


mfx$ret_length <- datagrid_values$ret_length

# Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=4mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=4mg - no reversion" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=4mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr5_rprev.rds")


### Discont ##############
ps_fit <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/ps_weight_tr5_rprev.rds")


cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)


datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)


mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)


mfx$ret_length <- datagrid_values$ret_length

mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=4mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "1-<50% reduction\nvs\nNo dose change",
                      ">=50%/<=4mg - no reversion" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                      "increase - no reversion"  = "Dose increase\nvs\nNo dose change"
    ),
    contrast = factor(contrast, levels = c(
      "1-<50% reduction\nvs\nNo dose change",
      ">=50% reduction/<=4mg\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr5_rprev.rds")



# Updates: Nov 12, 2025 ##########
install.packages("interactionRCS")
library(splines)
library(interactionRCS)
library(survival)
library(ggplot2)
library(scales)
library(cobalt)
library(glmtoolbox)
source("H:/SAS-Work/Missed doses Paper v1/Analysis June 24, 2025/R_functions.R")

## Trial 1: 1-2 day #######################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"
ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))

### overdose ##############
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nno change",
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr1_12d_rprev_ow.rds")


### Discont ##############
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nno change",
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr1_12d_rprev_ow.rds")



## Trial 2: 3-4 days ###############################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"
ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))

### Overdose ##############
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)


# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "Dose reduction to <50%\nvs\nno change",
                      ">=50%/<=40mg - no reversion" = "Dose reduction to >=50%/initial dose*\nvs\nno change" ,
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction to <50%\nvs\nno change",
      "Dose reduction to >=50%/initial dose*\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr2_34d_rprev_ow.rds")


### Discont ##############

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "Dose reduction to <50%\nvs\nno change",
                      ">=50%/<=40mg - no reversion" = "Dose reduction to >=50%/initial dose*\nvs\nno change" ,
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction to <50%\nvs\nno change",
      "Dose reduction to >=50%/initial dose*\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr2_34d_rprev_ow.rds")


## Trial 3 ###############################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"
ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))


### Overdose ##############
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "Dose reduction to <50%\nvs\nno change",
                      ">=50%/<=40mg - no reversion" = "Dose reduction to >=50%/initial dose*\nvs\nno change" ,
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction to <50%\nvs\nno change",
      "Dose reduction to >=50%/initial dose*\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr3_rprev_ow.rds")

### Discont ##############

# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.95)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=40mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "Dose reduction to <50%\nvs\nno change",
                      ">=50%/<=40mg - no reversion" = "Dose reduction to >=50%/initial dose*\nvs\nno change" ,
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction to <50%\nvs\nno change",
      "Dose reduction to >=50%/initial dose*\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))


saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr3_rprev_ow.rds")


## Trial 4 #######################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"
ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))

### overdose ##############
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)

# Step 2: Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Step 3: Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

# Step 4: Add ret_length manually since slopes() drops it
mfx$ret_length <- datagrid_values$ret_length

# Step 5: Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nno change",
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))



saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr4_rprev_ow.rds")


### Discont ##############
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)


datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)

mfx$ret_length <- datagrid_values$ret_length


mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("reversion - no reversion", "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "reversion - no reversion" = "Dose reduction\nvs\nno change",
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))



saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr4_rprev_ow.rds")



## Trial 5 ###############################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"
ps_fit <- readRDS(file=paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))


### Overdose ##############
# Step 1: Fit Cox model with interaction and spline
cox_spline_model <- coxph(
  Surv(time_od_acute, sens_od_acute) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)

# Create data grid for prediction
datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)

# Estimate log-HR (linear predictor)
mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)


mfx$ret_length <- datagrid_values$ret_length

# Filter to desired comparisons and compute HR
mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=4mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "Dose reduction to <50%\nvs\nno change",
                      ">=50%/<=4mg - no reversion" = "Dose reduction to >=50%/initial dose*\nvs\nno change" ,
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction to <50%\nvs\nno change",
      "Dose reduction to >=50%/initial dose*\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))



saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr5_rprev_ow.rds")


### Discont ##############

cox_spline_model <- coxph(
  Surv(time_discont, sens_discont) ~ trt * ns(ret_length, df = 3),
  data = ps_fit$data_weights,
  weights = wt_trunc
)
max=quantile(ps_fit$data_weights$ret_length, 0.75)


datagrid_values <- datagrid(
  model = cox_spline_model,
  ret_length = seq(0, max, by = 5),
  trt = levels(ps_fit$data_weights$trt)
)


mfx <- slopes(
  model = cox_spline_model,
  newdata = datagrid_values,
  variables = "trt",
  type = "lp"
)


mfx$ret_length <- datagrid_values$ret_length

mfx_filtered <- as.data.frame(mfx) %>%
  filter(
    term == "trt",
    contrast %in% c("1-<50% - no reversion", 
                    ">=50%/<=4mg - no reversion", 
                    "increase - no reversion"),
    !is.na(ret_length)
  ) %>%
  mutate(
    HR = exp(estimate),
    HR_low = exp(conf.low),
    HR_high = exp(conf.high)
  ) %>% 
  mutate(
    contrast = recode(contrast,
                      "1-<50% - no reversion" = "Dose reduction to <50%\nvs\nno change",
                      ">=50%/<=4mg - no reversion" = "Dose reduction to >=50%/initial dose*\nvs\nno change" ,
                      "increase - no reversion"  = "Dose increase\nvs\nno change"
    ),
    contrast = factor(contrast, levels = c(
      "Dose reduction to <50%\nvs\nno change",
      "Dose reduction to >=50%/initial dose*\nvs\nno change",
      "Dose increase\nvs\nno change"
    )))



saveRDS(mfx_filtered, file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr5_rprev_ow.rds")

