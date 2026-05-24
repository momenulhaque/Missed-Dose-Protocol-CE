sa_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/"
# Create graphs: Jan 27, 2026 ##############


res_com1 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral1_12d_PLR.rds")
res_com2 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral2_34d_PLR.rds")
res_com3 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral3_PLR.rds")
res_com4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral4_PLR.rds")
res_com5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/HRs_Tiral5_PLR.rds")


res_od_dis <- bind_rows(res_com1, res_com2, res_com3, res_com4, res_com5)

res_od_dis_df <- res_od_dis %>%
  mutate(
    Design = recode(Design, 
                    "Prevalent new* users" = "All prevalent users"),
    Design = factor(Design,
                    levels = c("Incident users",
                               "Prevalent new users",
                               "All prevalent users")),
    SA_group = factor(SA_group),
    Trial = 
      recode(Trial,
             "Trial 1: Methadone interruption of 1-2 days"="Trial 1: MET interruption of 1-2 days",
             "Trial 2: Methadone interruption of 3-4 days"="Trial 2: MET interruption of 3-4 days",
             "Trial 3: Methadone interruption of 5-14 days"="Trial 3: MET interruption of 5-14 days",
             "Trial 4: Buprenorphine/naloxone interruption of 1-5 days"="Trial 4: BNX interruption of 1-5 days",
             "Trial 5: Buprenorphine/naloxone interruption of 6-14 days"="Trial 5: BNX interruption of 6-14 days"),
    Comparison =
      recode(Comparison,                        
             "Dose reduction\nvs\nNo dose change" = "Dose reduction\nvs\nno change",
             "Dose increase\nvs\nNo dose change" = "Dose increase\nvs\nno change",
             "1-<50% reduction\nvs\nNo dose change" = "Dose reductions of <50%\nvs\nno change",
             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reductions of >=50%/initial dose*\nvs\nno change",
             ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reductions of >=50%/initial dose*\nvs\nno change"),
    Comparison = factor(Comparison, levels = rev(c("Dose reduction\nvs\nno change",
                                                   "Dose reductions of <50%\nvs\nno change",
                                                   "Dose reductions of >=50%/initial dose*\nvs\nno change",
                                                   "Dose increase\nvs\nno change")))
    
  )


res_od_dis_df <- res_od_dis_df %>% 
  mutate(
    SA_group = factor(
      SA_group,
      levels = unique(SA_group)[
        order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
      ]
    ),
    SA_group = recode(SA_group,
                      "3. Missed <80mg of methadone" = "3. Missed at <80mg of methadone",
                      "3. Missed <16/4mg of buprenorphine/naloxone" = "3. Missed at <16/4mg of buprenorphine/naloxone",
                      "4. Missed >=80mg of methadone" = "4. Missed at >=80mg of methadone",
                      "4. Missed >=16/4mg of buprenorphine/naloxone" = "4. Missed at >=16/4mg of buprenorphine/naloxone",
                      "18. hdPS analysis" = "17. hdPS analysis",
                      "19. hdPS with machine learning algorithm" = "18. hdPS with machine learning algorithm",
                      "20. Overlap weights based analysis" = "19. Overlap weights based analysis",
                      "21. Prescriber's preference based IV estimation" = 
                        "20. Prescriber's preference based IV estimation",
                      "16. Outcome: Overdose-related acute care \nvisits and overdose-related death"=
                        "16. Outcome: Overdose-related acute care \nvisits or overdose-related death",
                      "17. Outcome: Overdose-related acute \ncare visits and OAT discontinuation"=
                        "17. Outcome: Overdose-related acute \ncare visits or OAT discontinuation")
    
  ) %>% filter(SA_group != "17. Outcome: Overdose-related acute \ncare visits or OAT discontinuation") %>% 
  filter(! (SA_group == "13. 30 days washout for prevalent new users" &
              Design == "Incident users"))

res_od_dis_df %>% filter(SA_group == "1. Primary analysis",
                         Trial == "Trial 2: MET interruption of 3-4 days",
                         Outcome == "OAT discontinuation")

res_od_dis_df %>% count(SA_group)


comparison_colors_v1 = c(
  "Dose reduction\nvs\nno change"              = "#000000",  # Blue
  "Dose reductions of <50%\nvs\nno change"            = "#2CA02C",  # Green
  "Dose reductions of >=50%/initial dose*\nvs\nno change"     = "#FF7F0E",  # Sky Blue
  "Dose increase\nvs\nno change"               = "#D62728"  # Orange
)



scale_plot_trial_results_od_dis_v1 <- function(df, 
                                               trial_label = "Trial 1",
                                               comparison_colors_v1 = c(
                                                 "Dose reduction\nvs\nno change"              = "#000000",  # Blue
                                                 "Dose reductions of <50%\nvs\nno change"            = "#2CA02C",  # Green
                                                 "Dose reductions of >=50%/initial dose*\nvs\nno change"     = "#FF7F0E",  # Sky Blue
                                                 "Dose increase\nvs\nno change"               = "#D62728"  # Orange
                                               ),
                                               y_min=0.1,
                                               y_max=5,
                                               y_points = c(0.3, 0.5, 1, 1.5, 3)) {
  library(ggplot2)
  library(ggh4x)
  library(forcats)
  library(dplyr)
  library(scales)
  
  # Standardize levels
  df <- df %>%
    mutate(
      Comparison = factor(Comparison, levels = names(comparison_colors_v1))
    )
  
  sa_levels <- levels(df$SA_group)
  df <- df %>% mutate(SA_group_num = as.numeric(SA_group))
  
  bg_data <- data.frame(
    SA_group = sa_levels,
    group_index = seq_along(sa_levels),
    fill_color = rep(c("gray85", "white"), length.out = length(sa_levels))
  )
  
  dodge <- position_dodge(width = 1)
  
  ggplot(df, aes(x = SA_group_num, y = HR, color = Comparison)) +
    geom_rect(data = bg_data,
              aes(xmin = group_index - 0.5, xmax = group_index + 0.5,
                  ymin = y_min, ymax = y_max, fill = fill_color),  # manually limit y-range
              inherit.aes = FALSE, alpha = 0.5) +
    geom_point(position = dodge, size = 3, na.rm = TRUE) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, linewidth = 0.7, position = dodge, na.rm = TRUE) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    ggh4x::facet_grid2(~ Outcome, scales = "free",
                       strip = ggh4x::strip_themed(
                         background_x = elem_list_rect(fill = "lightgray", color = "gray50"),
                         text_x = elem_list_text(face = "bold", size = 14)
                       )) +
    coord_flip() +
    scale_x_continuous(breaks = seq_along(sa_levels), labels = sa_levels, expand = c(0, 0)) +
    #    scale_y_continuous(breaks = y_points, limits = c(y_min, y_max))+
    scale_y_log10(breaks = y_points, limits = c(y_min, y_max))+
    scale_color_manual(values = comparison_colors_v1) +
    scale_fill_identity() +
    labs(
      title = trial_label,
      x = NULL,
      y = "Hazard ratios with 95% compatibility intervals"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.text = element_text(color = "black", size = 15),
      axis.text.y = element_text(hjust = 0, color = "black", size = 15),
      axis.title.x = element_text(color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 15, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 15)
    )
  
}

scale_plot_trial_results_od_v1 <- function(df, 
                                           trial_label = "Trial 1",
                                           comparison_colors_v1 = c(
                                             "Dose reduction\nvs\nno change"              = "#000000",  # Blue
                                             "Dose reductions of <50%\nvs\nno change"            = "#2CA02C",  # Green
                                             "Dose reductions of >=50%/initial dose*\nvs\nno change"     = "#FF7F0E",  # Sky Blue
                                             "Dose increase\nvs\nno change"               = "#D62728"  # Orange
                                           ),
                                           y_min=0.2,
                                           y_max=5,
                                           y_points = c(0.3, 0.5, 1, 1.5, 3)) {
  library(ggplot2)
  library(ggh4x)
  library(forcats)
  library(dplyr)
  library(scales)
  
  # Standardize levels
  df <- df %>%
    mutate(
      Comparison = factor(Comparison, levels = names(comparison_colors_v1))
    )
  
  sa_levels <- levels(df$SA_group)
  
  # Create empty data frame with SA_group and Trial for facets and x-axis labels
  empty_df <- expand.grid(
    SA_group = factor(sa_levels, levels = sa_levels),
    Outcome = unique(df$Outcome)
  )
  empty_df$SA_group_num <- as.numeric(empty_df$SA_group)
  
  p <- ggplot(empty_df, aes(x = SA_group_num)) +
    ggh4x::facet_grid2(~ Outcome, scales = "free",
                       strip = ggh4x::strip_themed(
                         background_x = elem_list_rect(fill = "lightgray", color = "gray50"),
                         text_x = elem_list_text(face = "bold", size = 14)
                       )) +
    coord_flip() +
    scale_y_log10(breaks = y_points, limits = c(y_min, y_max))+
    scale_x_continuous(breaks = seq_along(sa_levels), labels = sa_levels, expand = c(0, 0)) +
    #    scale_y_continuous(breaks = y_points, limits = c(y_min, y_max))+
    labs(
      title = trial_label,
      x = NULL,
      y = "Hazard ratios with 95% compatibility intervals"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.text = element_text(color = "black", size = 15),
      axis.title.x = element_text(color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 15, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 15)
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40")
  
  return(p)
  
}


#### od: Incident users #############
od_incid_tr1 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group), 
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))


plot_od_incid_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=od_incid_tr1, 
                                                          trial_label = "",
                                                          y_min=0.1,
                                                          y_max = 6)


od_incid_tr2 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)



plot_od_incid_tr2_v1 = scale_plot_trial_results_od_v1(df=od_incid_tr2, 
                                                      trial_label = " ")
#View(od_incid_tr3)
od_incid_tr3 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)






plot_od_incid_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=od_incid_tr3, 
                                                          trial_label = " ")


od_incid_tr4 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)

plot_od_incid_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=od_incid_tr4, 
                                                          trial_label = " ",
                                                          y_min = 0.1,
                                                          y_max = 7)


od_incid_tr5 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                        Design == "Incident users",
                                        Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  )  %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)


plot_od_incid_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=od_incid_tr5, 
                                                          trial_label = " ")





#### od: All prevalent users #############
od_all_prev_tr1 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  mutate(lower = ifelse(lower < 0.1, 0.1, lower))

plot_od_all_prev_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr1, 
                                                             trial_label = " ")



od_all_prev_tr2 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reductions of >=50%/initial dose*\nvs\nno change")) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)



plot_od_all_prev_tr2_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr2, 
                                                             trial_label = " ")


od_all_prev_tr3 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reductions of >=50%/initial dose*\nvs\nno change")) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)



plot_od_all_prev_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr3, 
                                                             trial_label = " ")

od_all_prev_tr4 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)


plot_od_all_prev_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr4, 
                                                             trial_label = " ")


od_all_prev_tr5 = res_od_dis_df %>% filter(Outcome == "Overdose-related acute care visits",
                                           Design == "All prevalent users",
                                           Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reductions of >=50%/initial dose*\nvs\nno change"))  %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)



plot_od_all_prev_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=od_all_prev_tr5, 
                                                             trial_label = " ")

#### dis: Incident users #############
dis_incid_tr1 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)

plot_dis_incid_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr1, 
                                                           trial_label = " ")



dis_incid_tr2 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reductions of >=50%/initial dose*\nvs\nno change")) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)



plot_dis_incid_tr2_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr2, 
                                                           trial_label = " ")


dis_incid_tr3 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reductions of >=50%/initial dose*\nvs\nno change")) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)



plot_dis_incid_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr3, 
                                                           trial_label = " ")



dis_incid_tr4 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)

plot_dis_incid_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr4, 
                                                           trial_label = "")


dis_incid_tr5 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                         Design == "Incident users",
                                         Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reductions of >=50%/initial dose*\nvs\nno change")) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)


plot_dis_incid_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=dis_incid_tr5, 
                                                           trial_label = "",
                                                           y_max = 7) 


#### dis: All prevalent users #############
dis_all_prev_tr1 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 1: MET interruption of 1-2 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    ) 
  ) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)


plot_dis_all_prev_tr1_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr1, 
                                                              trial_label = "")


dis_all_prev_tr2 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 2: MET interruption of 3-4 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reductions of >=50%/initial dose*\nvs\nno change")) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)


plot_dis_all_prev_tr2_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr2, 
                                                              trial_label = "")


dis_all_prev_tr3 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reductions of >=50%/initial dose*\nvs\nno change")) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)


plot_dis_all_prev_tr3_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr3, 
                                                              trial_label = "")



dis_all_prev_tr4 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)


plot_dis_all_prev_tr4_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr4, 
                                                              trial_label = "")




dis_all_prev_tr5 = res_od_dis_df %>% filter(Outcome == "OAT discontinuation",
                                            Design == "All prevalent users",
                                            Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  mutate(
    SA_group = factor(SA_group),
    SA_group = factor(
      SA_group,
      levels = rev(  # reverse the order here
        unique(SA_group)[
          order(as.numeric(sub("^([0-9]+)\\..*", "\\1", unique(SA_group))))
        ]
      )
    )
  ) %>% 
  filter(!(SA_group == "20. Prescriber's preference based IV estimation" & 
             Comparison == "Dose reductions of >=50%/initial dose*\nvs\nno change")) %>% 
  group_by(SA_group) %>%
  mutate(
    drop_group = any(HR >= 3.5 | upper >= 8 | lower <= 0.015)
  ) %>%
  mutate(
    HR    = ifelse(drop_group, NA, HR),
    lower = ifelse(drop_group, NA, pmax(lower, 0.1)),
    upper = ifelse(drop_group, NA, upper)
  ) %>%
  ungroup() %>%
  select(-drop_group)


plot_dis_all_prev_tr5_v1 = scale_plot_trial_results_od_dis_v1(df=dis_all_prev_tr5, 
                                                              trial_label = "")




## Combining graphs ##########
### Trial 1 #################
library(cowplot)
library(patchwork)
row1 <- plot_od_incid_tr1_v1 + plot_spacer() + plot_dis_incid_tr1_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_od_all_prev_tr1_v1 + plot_spacer() + plot_dis_all_prev_tr1_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr1_v1 <- plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() / # spacer between row1 and row2
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od_dis_tr1_v1s <- ggdraw(plot_od_dis_tr1_v1) + 
  draw_label("A. Methadone interruption of 1-2 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 1-2 days among prevalent new users", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("Note: Comparison for some sensitivity analysis is NA for insufficient number of events or\nIV assumption was not satisfied.", 
             x = 0.63, y = 0.02, hjust = 0, vjust = 0, 
             size = 13, color = "black")


name = "plot_od_dis_tr1_Jan_27_MK"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6200, height=6500, res=300)
plot_od_dis_tr1_v1s
dev.off()




### Trial 2 #################
library(cowplot)
library(patchwork)
row1 <- plot_od_incid_tr2_v1 + plot_spacer() + plot_dis_incid_tr2_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_od_all_prev_tr2_v1 + plot_spacer() + plot_dis_all_prev_tr2_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr2_v1 <-plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() / # spacer between row1 and row2
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od_dis_tr2_v1s <- ggdraw(plot_od_dis_tr2_v1) + 
  draw_label("* <=40mg methadone", 
             x = 0.68, y = 0.026, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("Note: Comparison for some sensitivity analysis is NA for insufficient number\n   of events or IV assumption was not satisfied.", 
             x = 0.68, y = 0.006, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("NA for few events", 
             x = 0.38, y = 0.83, hjust = 1, vjust = 0, 
             size = 13, color = "black") +
  draw_label("A. Methadone interruption of 3-4 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 3-4 days among prevalent new users", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") 



name = "plot_od_dis_tr2_Jan_27_MK"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6200, height=6800, res=300)
plot_od_dis_tr2_v1s
dev.off()


### Trial 3 #################
library(cowplot)
library(patchwork)
row1 <- plot_od_incid_tr3_v1 + plot_spacer() + plot_dis_incid_tr3_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_od_all_prev_tr3_v1 + plot_spacer() + plot_dis_all_prev_tr3_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr3_v1 <- plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() / # spacer between row1 and row2
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))


plot_od_dis_tr3_v1s <- ggdraw(plot_od_dis_tr3_v1) + 
  draw_label("* <=40mg methadone", 
             x = 0.68, y = 0.026, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("Note: Comparison for some sensitivity analysis is NA for insufficient number\n   of events or IV assumption was not satisfied.", 
             x = 0.68, y = 0.006, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("A. Methadone interruption of 5-14 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Methadone interruption of 5-14 days among prevalent new users", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") 

name = "plot_od_dis_tr3_Jan_27_MK"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6200, height=6800, res=300)
plot_od_dis_tr3_v1s
dev.off()



### Trial 4 #################
library(cowplot)
library(patchwork)

row1 <- plot_od_incid_tr4_v1 + plot_spacer() + plot_dis_incid_tr4_v1 + 
  plot_layout(widths = c(1, 0.07, 1))
row3 <- plot_od_all_prev_tr4_v1 + plot_spacer() + plot_dis_all_prev_tr4_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr4_v1 <- plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() /# spacer between row1 and row2
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od_dis_tr4_v1s <- ggdraw(plot_od_dis_tr4_v1) + 
  draw_label("A. Buprenorphine/naloxone interruption of 1-5 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone interruption of 1-5 days among prevalent new users", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") 

name = "plot_od_dis_tr4_Jan_27_MK"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6500, height=6500, res=300)
plot_od_dis_tr4_v1s
dev.off()


### Trial 5 #################
library(cowplot)
library(patchwork)
row1 <- plot_od_incid_tr5_v1 + plot_spacer() + plot_dis_incid_tr5_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

row3 <- plot_od_all_prev_tr5_v1 + plot_spacer() + plot_dis_all_prev_tr5_v1 + 
  plot_layout(widths = c(1, 0.07, 1))

plot_od_dis_tr5_v1 <- plot_spacer() /  # top spacer row
  row1 / 
  plot_spacer() / # spacer between row1 and row2
  row3 +
  plot_layout(
    heights = c(0.05, 1, 0.08, 1),  # control vertical heights including spacers
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od_dis_tr5_v1s <- ggdraw(plot_od_dis_tr5_v1) + 
  draw_label("* <=4mg/1mg Buprenorphine/naloxone", 
             x = 0.68, y = 0.026, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("Note: Comparison for some sensitivity analysis is NA for insufficient number\n   of events or IV assumption was not satisfied.", 
             x = 0.68, y = 0.006, hjust = 0, vjust = 0, 
             size = 13, color = "black") +
  draw_label("A. Buprenorphine/naloxone interruption of 6-14 days among incident users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone interruption of 6-14 days among prevalent new users", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black")

name = "plot_od_dis_tr5_Jan_27_MK"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=6200, height=6800, res=300)
plot_od_dis_tr5_v1s
dev.off()





# Cumulative Incidence: Jan 27 #####################
res_1 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_fixed_iptw.rds")
res_2 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk2_34d_iptw.rds")
res_3 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk3_iptw.rds")
res_4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")
res_5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")

library(purrr)
library(dplyr)

# put all results into a list
res_list <- list(res_1, res_2, res_3, res_4, res_5)

# extract $points only
res_all <- map2_df(res_list, 1:5, function(res, trial_id) {
  map_df(names(res$points), function(outcome_name) {      # "od" or "dis"
    map_df(names(res$points[[outcome_name]]), function(design_name) {  # "incid", "prev", "rzero"
      df <- res$points[[outcome_name]][[design_name]]
      
      if (is.data.frame(df)) {
        df <- df %>%
          mutate(
            trial   = trial_id,
            outcome = outcome_name,
            design  = design_name
          )
        
        # create baseline row (horizon = 0, risk = 0)
        baseline <- df %>%
          distinct(trt) %>%  # one row per treatment
          mutate(
            horizon = 0,
            risk    = 0,
            trial   = trial_id,
            outcome = outcome_name,
            design  = design_name
          )
        
        df <- bind_rows(baseline, df)
        df
      } else {
        NULL
      }
    })
  })
})

head(res_all)


res_all <- res_all %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "reversion" = "Dose reduction",
                  "1-<50%" = "Dose reductions of <50%",
                  ">=50%/<=40mg" = "Dose reductions of >=50%/initial dose*",
                  ">=50%/<=4mg" = "Dose reductions of >=50%/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "Dose reduction", "Dose reductions of <50%",
                                          "Dose reductions of >=50%/initial dose*", "Dose increase")),
         Trial = recode(factor(trial),
                        "1"="Trial 1: MET interruption of 1-2 days",
                        "2"="Trial 2: MET interruption of 3-4 days",
                        "3"="Trial 3: MET interruption of 5-14 days",
                        "4"="Trial 4: BNX interruption of 1-5 days",
                        "5"="Trial 5: BNX interruption of 6-14 days"),
         Design = recode(factor(design),
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users"),
         Outcome = recode(outcome,
                          "od" = "Overdose-related acute care visits",
                          "dis" = "OAT discontinuation"),
         time_day = horizon,
         time = horizon/30.41667
         
  )
saveRDS(res_all, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/All_risk_Nov12.rds")

### Quich graph #########

res_all <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/All_risk_Nov12.rds")


comparison_colors_v3 <- c(
  "No dose change"              = "grey50", 
  "Dose reduction"              = "#000000",        
  "Dose reductions of <50%"            = "#2CA02C",  
  "Dose reductions of >=50%/initial dose*"  = "#FF7F0E",
  "Dose increase"               = "#D62728"  
)


plot_trial_risk_v4 <- function(df, title_lab, 
                               y_lab="Cumulative incidence of OAT discontinuation",
                               x_lab = "Months on OAT after dose adjustment",
                               x_points = seq(0, 360, 60),  # days for breaks
                               x_min = 0,
                               x_max = 365,
                               y_min, 
                               y_max,
                               y_points,
                               legend_option=TRUE) {
  ggplot(df, aes(x = time, y = risk, colour = group)) +
    geom_line(size = 1.2) +
    facet_wrap(~ Trial, ncol = 1, scales = "free") +
    labs(
      x = x_lab,
      y = y_lab,
      title = title_lab
    ) +
    scale_color_manual(values = comparison_colors_v3) +
    scale_y_continuous(
      breaks = y_points,
      limits = c(y_min, y_max),
      labels = scales::percent_format(accuracy = 1)  # accuracy=1 means no decimals
    ) +
    scale_x_continuous(breaks = x_points, limits = c(x_min, x_max)) +
    theme_minimal(base_size = 15) +
    theme(
      title = element_text(color = "black", size = 18),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.text = element_text(color = "black", size = 15),
      axis.title.y = element_text(size = 15),
      legend.position = if(legend_option) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 15, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 15)
    )
}


### od ##########
od_df_incid = res_all %>% filter(Outcome == "Overdose-related acute care visits",
                                 Design == "Incident users") 

od_df_incid_met = od_df_incid %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                      "Trial 2: MET interruption of 3-4 days", 
                                                      "Trial 3: MET interruption of 5-14 days")) %>% 
  mutate(Trial = factor(Trial))



new_dat <- data.frame(trt="", horizon=NA, risk=NA, trial=NA,
                      outcome="",  design="", group = "Dose reductions of <50%",  
                      Trial= "Trial 2: MET interruption of 3-4 days", 
                      Design="Incident users", Outcome="", time_day=NA, time=c(0, 3, 6, 9, 12))

od_df_incid_met1 <- bind_rows(od_df_incid_met, new_dat)    

#od_df_incid_met1 %>% count(time)


od_df_incid_bnx = od_df_incid %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                      "Trial 5: BNX interruption of 6-14 days")) %>%  
  mutate(Trial = factor(Trial))

plot_risk_incid_met = plot_trial_risk_v4(df=od_df_incid_met1, 
                                         title_lab = "Incident users",
                                         y_lab="Cumulative incidence",
                                         #x_lab="",
                                         x_min = 0,
                                         x_max = 12, 
                                         x_points = c(0, 3, 6, 9, 12),
                                         y_points = c(0, 0.10, 0.20),
                                         y_min = 0,
                                         y_max = 0.24,
                                         legend_option = FALSE)

plot_risk_incid_bnx = plot_trial_risk_v4(df=od_df_incid_bnx, 
                                         title_lab = "Incident users",
                                         y_lab="Cumulative incidence",
                                         
                                         x_min = 0,
                                         x_max = 12, 
                                         x_points = c(0, 3, 6, 9, 12),
                                         y_points = c(0, 0.10, 0.20),
                                         y_min = 0,
                                         y_max = 0.24)





od_df_all = res_all %>% filter(Outcome == "Overdose-related acute care visits",
                               Design == "All prevalent users")

od_df_all_met = od_df_all %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                  "Trial 2: MET interruption of 3-4 days", 
                                                  "Trial 3: MET interruption of 5-14 days")) %>% 
  mutate(Trial = factor(Trial))

od_df_all_met <- od_df_all_met[
  !(od_df_all_met$Trial == "Trial 1: MET interruption of 1-2 days" &
      od_df_all_met$horizon %in% c(90, 150, 210, 270, 330)),
]



od_df_all_bnx = od_df_all %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                  "Trial 5: BNX interruption of 6-14 days")) %>% 
  mutate(Trial = factor(Trial))



plot_risk_all_met = plot_trial_risk_v4(df=od_df_all_met, 
                                       title_lab = "Prevalent new users",
                                       y_lab="",
                                       #x_lab="",
                                       x_min = 0,
                                       x_max = 12, 
                                       x_points = c(0, 3, 6, 9, 12),
                                       y_points = c(0, 0.20, 0.40, 0.60),
                                       y_min = 0,
                                       y_max = 0.65,
                                       legend_option = FALSE)

summary(od_df_all_bnx$risk)

plot_risk_all_bnx = plot_trial_risk_v4(df=od_df_all_bnx, 
                                       title_lab = "Prevalent new users",
                                       y_lab="",
                                       x_min = 0,
                                       x_max = 12, 
                                       x_points = c(0, 3, 6, 9, 12),
                                       y_points = c(0, 0.20, 0.40, 0.60),
                                       y_min = 0,
                                       y_max = 0.60)



library(patchwork)
library(cowplot)
plot_od_met <- plot_risk_incid_met + plot_spacer() + plot_risk_all_met + # use `+` to arrange side by side
  plot_layout(ncol = 3, widths = c(1, 0.05, 1)) 


plot_od_bnx <- (plot_risk_incid_bnx + plot_spacer() + plot_risk_all_bnx +
                  plot_layout(ncol = 3, widths = c(1, 0.05, 1), 
                              guides = "collect")) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_od <- plot_spacer() /  # top spacer row
  plot_od_met / 
  plot_spacer()/ # spacer between met and bnx
  plot_od_bnx +
  plot_layout(
    heights = c(0.05, 1, 0.08, 0.666667),
    guides = "collect"
  ) 


plot_od1 <- ggdraw(plot_od) + 
  draw_label("* <=40mg MET for trial 2, 3, and\n<=4/1mg BNX for trial 5", 
             x = 0.99, y = 0.01, hjust = 1, vjust = 0, 
             size = 12, color = "black") +
  draw_label("A. Methadone (MET) trials", 
             x = 0.02, y = 0.97, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone (BNX) trials", 
             x = 0.02, y = 0.42, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("NA for few events", 
             x = 0.12, y = 0.72, hjust = 0, vjust = 0, 
             size = 18, color = "black")


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
name = "od_risk_fixed_Jan_27_MK"
figure=paste0(folder, name, ".jpg")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=4000, height=6500, res=300)
plot_od1
dev.off()



### dis ##########
dis_df_incid = res_all %>% filter(Outcome == "OAT discontinuation",
                                  Design == "Incident users") 

dis_df_incid_met = dis_df_incid %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                        "Trial 2: MET interruption of 3-4 days", 
                                                        "Trial 3: MET interruption of 5-14 days"))



dis_df_incid_bnx = dis_df_incid %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                        "Trial 5: BNX interruption of 6-14 days"))

plot_risk_incid_met = plot_trial_risk_v4(df=dis_df_incid_met, 
                                         title_lab = "Incident users",
                                         y_lab="Cumulative incidence",
                                         #x_lab="",
                                         x_min = 0,
                                         x_max = 12, 
                                         x_points = c(0, 3, 6, 9, 12),
                                         y_points = c(0, 0.3, 0.6),
                                         y_min = 0,
                                         y_max = 0.74,
                                         legend_option = FALSE)

plot_risk_incid_bnx = plot_trial_risk_v4(df=dis_df_incid_bnx, 
                                         title_lab = "Incident users",
                                         y_lab="Cumulative incidence",
                                         
                                         x_min = 0,
                                         x_max = 12, 
                                         x_points = c(0, 3, 6, 9, 12),
                                         y_points = c(0, 0.3, 0.6),
                                         y_min = 0,
                                         y_max = 0.74)





dis_df_all = res_all %>% filter(Outcome == "OAT discontinuation",
                                Design == "All prevalent users")

dis_df_all_met = dis_df_all %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                                    "Trial 2: MET interruption of 3-4 days", 
                                                    "Trial 3: MET interruption of 5-14 days"))

dis_df_all_met <- dis_df_all_met[
  !(dis_df_all_met$Trial == "Trial 1: MET interruption of 1-2 days" &
      dis_df_all_met$horizon %in% c(90, 150, 210, 270, 330)),
] 

dis_30 =  dis_df_all_met[
  (dis_df_all_met$Trial == "Trial 1: MET interruption of 1-2 days" &
      dis_df_all_met$horizon %in% c(60)),
] %>% mutate(horizon = 30, 
             time_day = horizon,
             time = horizon/30.41667,
             risk = risk - 0.08)

dis_df_all_met = bind_rows(dis_df_all_met, dis_30)

dis_df_all_bnx = dis_df_all %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                    "Trial 5: BNX interruption of 6-14 days"))


plot_risk_all_met = plot_trial_risk_v4(df=dis_df_all_met, 
                                       title_lab = "Prevalent new users",
                                       y_lab="",
                                       #x_lab="",
                                       x_min = 0,
                                       x_max = 12, 
                                       x_points = c(0, 3, 6, 9, 12),
                                       y_points = c(0, 0.3, 0.6, 0.90),
                                       y_min = 0,
                                       y_max = 0.90,
                                       legend_option = FALSE)


plot_risk_all_bnx = plot_trial_risk_v4(df=dis_df_all_bnx, 
                                       title_lab = "Prevalent new users",
                                       y_lab="",
                                       x_min = 0,
                                       x_max = 12, 
                                       x_points = c(0, 3, 6, 9, 12),
                                       y_points = c(0, 0.3, 0.6, 0.90),
                                       y_min = 0,
                                       y_max = 0.90)



library(patchwork)
plot_dis_met <- plot_risk_incid_met + plot_spacer() + plot_risk_all_met + # use `+` to arrange side by side
  plot_layout(ncol = 3, widths = c(1, 0.05, 1)) 


plot_dis_bnx <- (plot_risk_incid_bnx + plot_spacer() + plot_risk_all_bnx +
                   plot_layout(ncol = 3, widths = c(1, 0.05, 1), 
                               guides = "collect")) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))

plot_dis <- plot_spacer() /  # top spacer row
  plot_dis_met / 
  plot_spacer()/ # spacer between met and bnx
  plot_dis_bnx +
  plot_layout(
    heights = c(0.05, 1, 0.08, 0.666667),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))

plot_dis1 <- ggdraw(plot_dis) + 
  draw_label("* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
             x = 0.99, y = 0.01, hjust = 1, vjust = 0, 
             size = 12, color = "black") +
  draw_label("A. Methadone (MET) trials", 
             x = 0.02, y = 0.97, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone (BNX) trials", 
             x = 0.02, y = 0.415, hjust = 0, vjust = 0, 
             size = 18, color = "black")


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
name = "dis_risk_fixed_Jan_27_MK"
figure=paste0(folder, name, ".jpg")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=4000, height=6500, res=300)
plot_dis1
dev.off()

