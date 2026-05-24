plot_dose_flip <- function(df, title_lab, 
                           y_lab="Methadone dose (mg)",
                           x_lab = NULL,
                           y_min=0, 
                           y_max=222, 
                           y_points=c(0, 40, 80, 120, 160, 200),
                           legend_option=TRUE){
  ggplot(df, aes(x = trt, y = dose, fill = dose_type)) +
    geom_boxplot(
      outlier.shape = 16,
      outlier.size = 1.5,
      alpha = 0.7,
      position = position_dodge(width = 0.8)
    ) +
    # mean diamonds — now grouped by dose_type
    stat_summary(
      aes(group = dose_type),
      fun = mean,
      geom = "point",
      shape = 23,             # diamond
      size = 3,
      color = "black",
      fill = "white",
      position = position_dodge(width = 0.8)
    ) +
    facet_wrap(~ ret_length_cat, ncol = 3) +
    scale_fill_brewer(palette = "Set2", name = "Dose Type") +
    scale_y_continuous(breaks = y_points, limits = c(y_min, y_max)) +
    labs(x = x_lab, y = y_lab, title = title_lab) +
    theme_minimal(base_size = 13) +
    coord_flip() +
    theme(
      title = element_text(color = "black", size = 14),
      #panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      axis.title.y = element_text(size = 14),
      legend.position = if(legend_option) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 13)
    )
}

comparison_colors_v3 <- c(
  "No dose change"              = "darkblue", 
  "Dose reduction"              = "#000000",        
  "Dose increase"               = "#D62728",  
  "Reductions of <50%"            = "#2CA02C",  
  "Reduction of >=50%/initial dose*"  = "#FF7F0E"  
)


# Boxplot for ddose  ###########################
met <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_m.sas7bdat") %>%
  mutate(trt_s456 = case_when(
    trt_s456 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s456
  ),
  ret_length_cat = factor(case_when(
    ret_length <= 30 ~ "0-30",
    ret_length > 30 & ret_length <= 90 ~ "31-90",
    ret_length > 90 ~ ">90"
  ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
  high_ddose = factor(case_when(
    ddose_prev >= 80 ~ ">=80", 
    ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  # filter(trial_cat == "1-3" & ret_length_cat == "0-90") %>%
  # rename(trt = !!trt_name) %>%
  # filter(trt %in% trt_strategy) %>%
  # mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()


bnx <- read_sas("E:/SAS Data/Missed dose protocol_v2/TT_analytic_rprev_sens_nw_b.sas7bdat") %>%
  mutate(trt_s5 = case_when(
    trt_s5 %in% c("1-<25%", "25-<50%") ~ "1-<50%",
    TRUE ~ trt_s5), 
    ret_length_cat = factor(case_when(
      ret_length <= 30 ~ "0-30",
      ret_length > 30 & ret_length <= 90 ~ "31-90",
      ret_length > 90 ~ ">90"
    ), levels = c("0-30", "31-90", ">90"), ordered = TRUE), 
    high_ddose = factor(case_when(
      ddose_prev >= 80 ~ ">=80", 
      ddose_prev < 80 ~ "<80"), levels = c("<80", ">=80"), ordered = TRUE)) %>%
  # filter(trial_cat == "1-3" & ret_length_cat == "0-90") %>%
  # rename(trt = !!trt_name) %>%
  # filter(trt %in% trt_strategy) %>%
  # mutate(trt = factor(trt, levels = trt_strategy)) %>% 
  group_by(studyid_rev, episode) %>%
  mutate(
    prior_gap = cumsum(gap_days),  # Cumulative sum including current row
    prior_gap_num = row_number()   # Row number = order of gap_days within group
  ) %>%
  ungroup()



# Trial 1 ####
library(tidyverse)

trt_strategy <- c("no reversion", "reversion", "increase")
trt_name <- "trt_p1"
#dat1 %>% count(trt)

dat1 <- met %>%
  filter(trial_gap %in% 1:2) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  # compute 1st and 99th percentiles for both doses (ignoring NAs)
  mutate(
    ddose_prev_cap_low = quantile(ddose_prev, 0.01, na.rm = TRUE),
    ddose_prev_cap_high = quantile(ddose_prev, 0.99, na.rm = TRUE),
    ddose_cap_low = quantile(ddose, 0.01, na.rm = TRUE),
    ddose_cap_high = quantile(ddose, 0.99, na.rm = TRUE)
  ) %>%
  # cap both variables
  mutate(
    ddose_prev = pmin(pmax(ddose_prev, ddose_prev_cap_low), ddose_prev_cap_high),
    ddose = pmin(pmax(ddose, ddose_cap_low), ddose_cap_high)
  ) %>%
  select(ret_length_cat, trt, ddose_prev, ddose)

# Reshape data to long format
dat1_long <- dat1 %>%
  pivot_longer(
    cols = c(ddose_prev, ddose),
    names_to = "dose_type",
    values_to = "dose"
  ) %>%
  mutate(
    dose_type = recode(dose_type,
                       ddose_prev = "Before interruption",
                       ddose = "After interruption"),
    dose_type = factor(dose_type, levels = c("Before interruption", "After interruption")),
    ret_length_cat = recode(ret_length_cat,
                            "0-30" = "Within 30 days of OAT initiation",
                 "31-90" = "Within 31 to 90 days of OAT initiation",
                 ">90" = "After 90 days of OAT initiation"),
    ret_length_cat = factor(ret_length_cat,
                            levels = c("Within 30 days of OAT initiation", 
                                       "Within 31 to 90 days of OAT initiation", 
                                       "After 90 days of OAT initiation")),
    trt = recode(trt,
                 "no reversion" = "No dose change",
                 "reversion" = "Dose reduction",
                 "increase" = "Dose increase"),
        trt = factor(trt,
                     levels = c("No dose change",
                                "Dose reduction",
                                "Dose increase"))
        
      ) 

#dat1_long %>% count(dose_type)

dat1_long <- dat1_long %>% mutate(trt = factor(trt, levels = rev(levels(dat1_long$trt))),
                                  dose_type = factor(dose_type,
                                                     levels = rev(levels(dat1_long$dose_type))))





plot_tr1 = plot_dose_flip(df=dat1_long, title_lab="Trial 1: Methadone interruption of 1-2 days", 
                     y_lab="Methadone dose (mg)",
                     x_lab = NULL,
                     y_min=0, 
                     y_max=222, 
                     y_points=c(0, 40, 80, 120, 160, 200),
                     legend_option=FALSE)


# Trial 2 ##########

trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat2 <- met %>%
  filter(trial_gap %in% 3:4) %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  # compute 1st and 99th percentiles for both doses (ignoring NAs)
  mutate(
    ddose_prev_cap_low = quantile(ddose_prev, 0.01, na.rm = TRUE),
    ddose_prev_cap_high = quantile(ddose_prev, 0.99, na.rm = TRUE),
    ddose_cap_low = quantile(ddose, 0.01, na.rm = TRUE),
    ddose_cap_high = quantile(ddose, 0.99, na.rm = TRUE)
  ) %>%
  # cap both variables
  mutate(
    ddose_prev = pmin(pmax(ddose_prev, ddose_prev_cap_low), ddose_prev_cap_high),
    ddose = pmin(pmax(ddose, ddose_cap_low), ddose_cap_high)
  ) %>%
  select(ret_length_cat, trt, ddose_prev, ddose)

# Reshape data to long format
dat2_long <- dat2 %>%
  pivot_longer(
    cols = c(ddose_prev, ddose),
    names_to = "dose_type",
    values_to = "dose"
  ) %>%
  mutate(
    dose_type = recode(dose_type,
                       ddose_prev = "Before interruption",
                       ddose = "After interruption"),
    dose_type = factor(dose_type, levels = c("Before interruption", "After interruption")),
    ret_length_cat = recode(ret_length_cat,
                            "0-30" = "Within 30 days of OAT initiation",
                            "31-90" = "Within 31 to 90 days of OAT initiation",
                            ">90" = "After 90 days of OAT initiation"),
    ret_length_cat = factor(ret_length_cat,
                            levels = c("Within 30 days of OAT initiation", 
                                       "Within 31 to 90 days of OAT initiation", 
                                       "After 90 days of OAT initiation")),
    trt = recode(trt,
                 "no reversion" = "No dose change",
                 "1-<50%" = "Reductions of <50%",
                 ">=50%/<=40mg" = "Reduction of >=50%/initial dose*",
                 "increase" = "Dose increase"),
    trt = factor(trt,
                 levels = c("No dose change",
                            "Reductions of <50%",
                            "Reduction of >=50%/initial dose*",
                            "Dose increase"))
    
  ) 





dat2_long <- dat2_long %>% mutate(trt = factor(trt, levels = rev(levels(dat2_long$trt))), 
                                  dose_type = factor(dose_type,
                                                     levels = rev(levels(dat2_long$dose_type))))

plot_tr2 = plot_dose_flip(df=dat2_long, title_lab="Tria 2: Methadone interruption of 3-4 days", 
                     y_lab="Methadone dose (mg)",
                     x_lab = NULL,
                     y_min=0, 
                     y_max=222, 
                     y_points=c(0, 40, 80, 120, 160, 200),
                     legend_option=FALSE)


# Trial 3 ##########
trt_strategy = c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
trt_name = "trt_s456"

dat3 <- met %>%
  filter(trial_cat == "5-14") %>%
  rename(trt = !!trt_name) %>%
  filter(trt %in% trt_strategy) %>%
  mutate(trt = factor(trt, levels = trt_strategy)) %>%
  # compute 1st and 99th percentiles for both doses (ignoring NAs)
  mutate(
    ddose_prev_cap_low = quantile(ddose_prev, 0.01, na.rm = TRUE),
    ddose_prev_cap_high = quantile(ddose_prev, 0.99, na.rm = TRUE),
    ddose_cap_low = quantile(ddose, 0.01, na.rm = TRUE),
    ddose_cap_high = quantile(ddose, 0.99, na.rm = TRUE)
  ) %>%
  # cap both variables
  mutate(
    ddose_prev = pmin(pmax(ddose_prev, ddose_prev_cap_low), ddose_prev_cap_high),
    ddose = pmin(pmax(ddose, ddose_cap_low), ddose_cap_high)
  ) %>%
  select(ret_length_cat, trt, ddose_prev, ddose)

# Reshape data to long format
dat3_long <- dat3 %>%
  pivot_longer(
    cols = c(ddose_prev, ddose),
    names_to = "dose_type",
    values_to = "dose"
  ) %>%
  mutate(
    dose_type = recode(dose_type,
                       ddose_prev = "Before interruption",
                       ddose = "After interruption"),
    dose_type = factor(dose_type, levels = c("Before interruption", "After interruption")),
    ret_length_cat = recode(ret_length_cat,
                            "0-30" = "Within 30 days of OAT initiation",
                            "31-90" = "Within 31 to 90 days of OAT initiation",
                            ">90" = "After 90 days of OAT initiation"),
    ret_length_cat = factor(ret_length_cat,
                            levels = c("Within 30 days of OAT initiation", 
                                       "Within 31 to 90 days of OAT initiation", 
                                       "After 90 days of OAT initiation")),
    trt = recode(trt,
                 "no reversion" = "No dose change",
                 "1-<50%" = "Reductions of <50%",
                 ">=50%/<=40mg" = "Reduction of >=50%/initial dose*",
                 "increase" = "Dose increase"),
    trt = factor(trt,
                 levels = c("No dose change",
                            "Reductions of <50%",
                            "Reduction of >=50%/initial dose*",
                            "Dose increase"))
    
  ) 




dat3_long <- dat3_long %>% mutate(trt = factor(trt, levels = rev(levels(dat3_long$trt))),
                                  dose_type = factor(dose_type,
                                                     levels = rev(levels(dat3_long$dose_type))))


plot_tr3 = plot_dose_flip(df=dat3_long, title_lab="Trial 3: Methadone interruption of 5-14 days", 
                     y_lab="Methadone dose (mg)",
                     x_lab = NULL,
                     y_min=0, 
                     y_max=222, 
                     y_points=c(0, 40, 80, 120, 160, 200),
                     legend_option=FALSE)


## creating plot ###########
library(patchwork)
library(cowplot)

plot_met_dose <- 
  plot_tr1/ 
  plot_spacer()/
  plot_tr2/
  plot_spacer()/
  plot_tr3 + 
  plot_layout(
    heights = c(1, 0.01, 1, 0.01, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_met_dose1 <- ggdraw(plot_met_dose) + 
  draw_label("* <=40mg methadone", 
             x = 0.63, y = 0.01, hjust = 0, vjust = 0, 
             size = 13, color = "black")


sa_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/prior_dose/"
name = "trial_1_3_flipped_MK"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=4500, height=6500, res=300)
plot_met_dose1
dev.off()

