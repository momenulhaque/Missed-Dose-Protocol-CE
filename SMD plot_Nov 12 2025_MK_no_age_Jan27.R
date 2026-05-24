# Primary: SMD: iptw: with notes #################
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"
library(cobalt)
library(tidyverse)

## Trial 1: 1-2d, iptw #############

x_low = -.6; x_high = 0.6
x_breaks = c(-.4, -.2, 0, 0.2, 0.4)

baselinevars = c("sexM", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 "ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

new.names1 <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_gap     = "Total OAT interruptions since initiation"
)




new.names <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply"
)


var_levels <- names(new.names1)

bg_data1 <- data.frame(
  var = var_levels,
  var_index = seq_along(var_levels),
  fill_color = rep(c("gray65", "white"), length.out = length(var_levels))
)

bg_data <- bg_data1[bg_data1$var != "prior_gap", ]


### incident users ##########
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))$data_weights
dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "reversion" = "Dose reduction",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)


fmla <- as.formula(
  paste("trt ~", paste(baselinevars, collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "Dose reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, .1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose reduction\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 

plot_incid_1 + 
  geom_vline(xintercept = c(-0.2, 0.2)) 


dat_12 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose increase\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 

un_smd_11 <- round((prep_11$Balance$Diff.Un), 3)
adj_smd_11 <- round((prep_11$Balance$Diff.Adj), 3)
un_smd_12 <- round((prep_12$Balance$Diff.Un), 3)
adj_smd_12 <- round((prep_12$Balance$Diff.Adj), 3)
row_names <- row.names(prep_11$Balance)

smd_incid_1 <- data.frame(un_smd_11, un_smd_12, adj_smd_11, adj_smd_12) %>%
  mutate(
    un_smd = apply(cbind(un_smd_11, un_smd_12), 1, 
                   function(x) x[which.max(abs(x))]),
    adj_smd = apply(cbind(adj_smd_11, adj_smd_12), 1, 
                    function(x) x[which.max(abs(x))]),    
    Design = "Incident users",
    Trial = "Trial 1",
    cov=row_names
  )

### prevalent new users #################

dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))$data_weights
dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "reversion" = "Dose reduction",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)



fmla <- as.formula(
  paste("trt ~", paste(c(baselinevars, "prior_gap"), collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "Dose reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_rprev_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "bottom",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_12 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)


plot_rprev_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 

library(patchwork)
library(cowplot)

row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_1 = plot_spacer() / 
  row1/
  plot_spacer() /
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))

library(cowplot)
library(patchwork)


plot_1s <- ggdraw(plot_1) + 
  draw_label("A. Methadone interruption of 1-2 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 1-2 days among prevalent new users",  
             x = 0.02, y = 0.485, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("Note: Confounders with SMD >= 0.1 were included in the \noutcome models as covariates.", 
             x = 0.65, y = 0.01, hjust = 0, vjust = 0, 
             size = 14, color = "black")

name = "smd_plot_tr1_12d_Jan27_iptw_MK"

figure=paste0(ps_folder, name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=5000, height=6000, res=300)
plot_1s
dev.off()


## Trial 2: 3-4d, iptw #############
baselinevars = c("sexM", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", 
                 "ret_length", 
                 "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

new.names1 <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_gap     = "Total OAT interruptions since initiation"
)




new.names <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply"
)


var_levels <- names(new.names1)

bg_data1 <- data.frame(
  var = var_levels,
  var_index = seq_along(var_levels),
  fill_color = rep(c("gray65", "white"), length.out = length(var_levels))
)

bg_data <- bg_data1[bg_data1$var != "prior_gap", ]

### incident users ##########
x_low = -.6; x_high = 0.6
x_breaks = c(-.4, -.2, 0, 0.2, 0.4)

dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))$data_weights

dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "1-<50%" = "1-<50% reduction",
                                   ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)



fmla <- as.formula(
  paste("trt ~", paste(baselinevars, collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose reductions of <50%\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 





dat_12 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose reductions of >=50%/initial dose*\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_13 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_13$trt <- droplevels(dat_13$trt)

prep_13 <- bal.tab(
  fmla,
  data = dat_13,
  weights = dat_13$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_3 <- love.plot(
  prep_13,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose increase\nvs\nno change"   # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  theme(
    title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


### All prevalent users #################

dat <- readRDS(paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))$data_weights
dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "1-<50%" = "1-<50% reduction",
                                   ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)

fmla <- as.formula(
  paste("trt ~", paste(c(baselinevars, "prior_gap"), collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_rprev_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "bottom",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_12 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_rprev_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_13 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_13$trt <- droplevels(dat_13$trt)

prep_13 <- bal.tab(
  fmla,
  data = dat_13,
  weights = dat_13$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)




plot_rprev_3 <- love.plot(
  prep_13,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "bottom",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 

library(patchwork)

row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + plot_spacer() + plot_incid_3 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 + plot_spacer() + plot_rprev_3 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

plot_2 = plot_spacer() / 
  row1/
  plot_spacer() /
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))

library(cowplot)


plot_2s <- ggdraw(plot_2) + 
  draw_label("A. Methadone interruption of 3-4 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 3-4 days among prevalent new users",  
             x = 0.02, y = 0.485, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("* <=40mg methadone\nNote: Confounders with SMD >= 0.1 were included in the \noutcome models as covariates.", 
             x = 0.73, y = 0.01, hjust = 0, vjust = 0, 
             size = 13, color = "black") 


name = "smd_plot_tr2_34d_Jan27_iptw_MK"
figure=paste0(ps_folder, name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6200, height=6000, res=300)
plot_2s
dev.off()



## Trial 3: 5-14d, iptw #############
baselinevars = c("sexM", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

new.names1 <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_gap     = "Total OAT interruptions since initiation"
)




new.names <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply"
)


var_levels <- names(new.names1)

bg_data1 <- data.frame(
  var = var_levels,
  var_index = seq_along(var_levels),
  fill_color = rep(c("gray65", "white"), length.out = length(var_levels))
)

bg_data <- bg_data1[bg_data1$var != "prior_gap", ]

### incident users ##########
x_low = -.6; x_high = 0.6
x_breaks = c(-.4, -.2, 0, 0.2, 0.4)
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr3_ow.rds"))$data_weights
dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "1-<50%" = "1-<50% reduction",
                                   ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)



fmla <- as.formula(
  paste("trt ~", paste(baselinevars, collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose reductions of <50%\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 





dat_12 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose reductions of >=50%/initial dose*\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_13 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_13$trt <- droplevels(dat_13$trt)

prep_13 <- bal.tab(
  fmla,
  data = dat_13,
  weights = dat_13$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_3 <- love.plot(
  prep_13,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose increase\nvs\nno change"   # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 

### All prevalent users #################
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))$data_weights
dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "1-<50%" = "1-<50% reduction",
                                   ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)

fmla <- as.formula(
  paste("trt ~", paste(c(baselinevars, "prior_gap"), collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_rprev_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "bottom",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_12 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_rprev_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_13 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_13$trt <- droplevels(dat_13$trt)

prep_13 <- bal.tab(
  fmla,
  data = dat_13,
  weights = dat_13$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)




plot_rprev_3 <- love.plot(
  prep_13,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "bottom",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 

library(patchwork)


row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + plot_spacer() + plot_incid_3 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 + plot_spacer() + plot_rprev_3 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

plot_3 = plot_spacer() / 
  row1/
  plot_spacer() /
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))

library(cowplot)


plot_3s <- ggdraw(plot_3) + 
  draw_label("A. Methadone interruption of 5-14 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 5-14 days among prevalent new users",  
             x = 0.02, y = 0.485, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("* <=40mg methadone\nNote: Confounders with SMD >= 0.1 were included in the \noutcome models as covariates.", 
             x = 0.73, y = 0.01, hjust = 0, vjust = 0, 
             size = 13, color = "black") 



name = "smd_plot_tr3_Jan27_iptw_MK"
figure=paste0(ps_folder, name,".png")
ppi = 300
png(figure, width=6200, height=6000, res=300)
plot_3s
dev.off()



## Trial 4: bnx 1-5d, iptw #############
x_low = -.6; x_high = 0.6
x_breaks = c(-.4, -.2, 0, 0.2, 0.4)

baselinevars = c("sexM", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

new.names1 <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_gap     = "Total OAT interruptions since initiation"
)




new.names <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply"
)


var_levels <- names(new.names1)

bg_data1 <- data.frame(
  var = var_levels,
  var_index = seq_along(var_levels),
  fill_color = rep(c("gray65", "white"), length.out = length(var_levels))
)

bg_data <- bg_data1[bg_data1$var != "prior_gap", ]


### incident users ##########
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr4_ow.rds"))$data_weights

dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "reversion" = "Dose reduction",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)


fmla <- as.formula(
  paste("trt ~", paste(baselinevars, collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "Dose reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose reduction\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_12 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose increase\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


### All prevalent users #################
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))$data_weights

dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "reversion" = "Dose reduction",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)



fmla <- as.formula(
  paste("trt ~", paste(c(baselinevars, "prior_gap"), collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "Dose reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_rprev_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "bottom",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_12 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)


plot_rprev_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_4 = plot_spacer() / 
  row1/
  plot_spacer() /
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))

library(cowplot)


plot_4s <- ggdraw(plot_4) + 
  draw_label("A. Buprenorphine/naloxone interruption of 1-5 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone interruption of 1-5 days among prevalent new users",  
             x = 0.02, y = 0.485, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("Note: Confounders with SMD >= 0.1 were included in the \noutcome models as covariates.", 
             x = 0.65, y = 0.01, hjust = 0, vjust = 0, 
             size = 14, color = "black")



name = "smd_plot_tr4_Jan27_iptw"
figure=paste0(ps_folder, name,".png")
ppi = 300
png(figure, width=5000, height=6000, res=300)
plot_4s
dev.off()



## Trial 5: bnx 6-14d, iptw #############
baselinevars = c("sexM", "age_lt30", "age_30_39", 
                 "age_40_49", "age_ge50", "rural", "year", "po_preOUD", 
                 "phy_att_gr2", "social_assist", "ret_length", "home_grp_b",
                 "incar_grp_bl", "drug_acute", "cci_gt0", "copd_asthma", 
                 "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob", 
                 "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss")

new.names1 <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply",
  prior_gap     = "Total OAT interruptions since initiation"
)




new.names <- c(
  sexM          = "Sex: Male",
  age_lt30      = "Age category: 18-29",
  age_30_39     = "Age category: 30-39",
  age_40_49     = "Age category: 40-49",
  age_ge50      = "Age category: 50+",
  rural         = "Region: Rural",
  year          = "Year of reinitiation date",
  po_preOUD     = "Opioid dispensations other than OAT prior to the first OUD indication",
  phy_att_gr2   = "Attachment to the prescriber who reinitiated OAT",
  social_assist = "Receipt of income assistance in the past year",
  ret_length    = "Duration of OAT episode before time zero",
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
  ddose_prev    = "Last daily dose",
  clin_size_b   = "Large practice size (>300 clients)",
  virtual_care  = "Virtual care in past month",
  pss           = "Prescriber safer supply"
)


var_levels <- names(new.names1)

bg_data1 <- data.frame(
  var = var_levels,
  var_index = seq_along(var_levels),
  fill_color = rep(c("gray65", "white"), length.out = length(var_levels))
)

bg_data <- bg_data1[bg_data1$var != "prior_gap", ]

### incident users ##########
x_low = -.6; x_high = 0.6
x_breaks = c(-.4, -.2, 0, 0.2, 0.4)

dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr5_ow.rds"))$data_weights

dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "1-<50%" = "1-<50% reduction",
                                   ">=50%/<=4mg" = ">=50% reduction/initial dose*",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)



fmla <- as.formula(
  paste("trt ~", paste(baselinevars, collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose reductions of <50%\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 





dat_12 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose reductions of >=50%/initial dose*\nvs\nno change"  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_13 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_13$trt <- droplevels(dat_13$trt)

prep_13 <- bal.tab(
  fmla,
  data = dat_13,
  weights = dat_13$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_incid_3 <- love.plot(
  prep_13,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = "Dose increase\nvs\nno change"   # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 

### All prevalent users #################
dat <- readRDS(file=paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))$data_weights

dat <- dat %>% mutate(trt = recode(trt,
                                   "no reversion" = "No dose change",
                                   "1-<50%" = "1-<50% reduction",
                                   ">=50%/<=4mg" = ">=50% reduction/initial dose*",
                                   "increase" = "Dose increase"),
                      incar_grp_bl=ifelse(incar_grp_bl=="1_pastyear", 1, 0) 
)

fmla <- as.formula(
  paste("trt ~", paste(c(baselinevars, "prior_gap"), collapse = " + "))
)

dat_11 <- subset(dat, trt %in% c("No dose change", "1-<50% reduction"))
dat_11$trt <- droplevels(dat_11$trt)

prep_11 <- bal.tab(
  fmla,
  data = dat_11,
  weights = dat_11$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_rprev_1 <- love.plot(
  prep_11,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "bottom",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_12 <- subset(dat, trt %in% c("No dose change", ">=50% reduction/initial dose*"))
dat_12$trt <- droplevels(dat_12$trt)

prep_12 <- bal.tab(
  fmla,
  data = dat_12,
  weights = dat_12$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)

plot_rprev_2 <- love.plot(
  prep_12,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names1,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "none",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data1,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 


dat_13 <- subset(dat, trt %in% c("No dose change", "Dose increase"))
dat_13$trt <- droplevels(dat_13$trt)

prep_13 <- bal.tab(
  fmla,
  data = dat_13,
  weights = dat_13$wt_trunc,
  estimand = "ATE",
  m.threshold = 0.1,
  s.d.denom = "pooled",
  un = TRUE,
  stats = "mean.diffs"
)




plot_rprev_3 <- love.plot(
  prep_13,
  stat = "mean.diffs",
  star = "raw",
  stars = "none",
  thresholds = c(m, 0.1),
  abs = FALSE,
  var.names = new.names,
  sample.names = c("Before propensity score weighting", "After propensity score weighting"),
  position = "bottom",
  line = FALSE,
  colors = c("black", "red"),
  shapes = c(16, 16),
  title = ""  # Add this line
) +
  geom_rect(
    data = bg_data,
    aes(ymin = var_index - 0.5, ymax = var_index + 0.5,
        xmin = -Inf, xmax = Inf, fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  scale_x_continuous(
    limits = c(x_low, x_high),   # set the range of the x-axis
    breaks = x_breaks,  # customize tick marks
    expand = c(0, 0)         # avoid extra padding
  ) +
  scale_fill_identity() +
  labs(
    x = "SMD",  # <-- customize x-axis title here
    color = "", 
    shape = ""
  ) + 
  theme(
    title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.text.x =  element_text(color = "black", size = 14),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black")
  ) 

library(patchwork)


row1 <- plot_incid_1 + plot_spacer() + plot_incid_2 + plot_spacer() + plot_incid_3 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

row3 <- plot_rprev_1 + plot_spacer() + plot_rprev_2 + plot_spacer() + plot_rprev_3 +
  plot_layout(widths = c(1, 0.07, 1, 0.07, 1))

plot_5 = plot_spacer() / 
  row1/
  plot_spacer() /
  row3 + 
  plot_layout(
    heights = c(0.06, 1, 0.04, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))

library(cowplot)


plot_5s <- ggdraw(plot_5) + 
  draw_label("A. Buprenorphin/naloxone interruption of 6-14 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphin/naloxone interruption of 6-14 days among prevalent new users",  
             x = 0.02, y = 0.485, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("* <=4mg buprenorphine/naloxone\nNote: Confounders with SMD >= 0.1 were included in the \noutcome models as covariates.", 
             x = 0.73, y = 0.01, hjust = 0, vjust = 0, 
             size = 13, color = "black") 

name = "smd_plot_tr5_Jan27_iptw_MK"
figure=paste0(ps_folder, name,".png")
ppi = 300
png(figure, width=6200, height=6000, res=300)
plot_5s
dev.off()


