library(tidyverse)
library(survival)
library(haven)
library(Publish)
library(ggplot2)
library(survminer)
require(tableone)
require(survey)
library(lubridate)
library(naniar) # For missing value
library(cobalt)
library(cowplot)


# Results ##################
res11 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res1_iptw.rds")
res12 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res2_iptw.rds")
res13 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res3_iptw.rds")
res14 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res4_iptw.rds")
res15 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res5_iptw.rds")

res21 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res1_ow.rds")
res22 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res2_ow.rds")
res23 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res3_ow.rds")
res24 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res4_ow.rds")
res25 <- readRDS(file="H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/res5_ow.rds")


extract_all_results <- function(res_obj, trial_label = "Trial 1", model_type = "IPTW") {
  designs <- c("incid" = "Incident users", 
               "prev" = "Prevalent new users", 
               "rzero" = "Prevalent new* users")
  outcomes <- c("od" = "Overdose-related acute care visits", 
                "dis" = "OAT discontinuation")
  
  bind_rows(
    lapply(names(outcomes), function(outcome) {
      lapply(names(designs), function(design) {
        result <- res_obj[[outcome]][[design]]
        if (is.null(result)) return(NULL)
        
        data.frame(
          Comparison = result$names,
          HR = result$OR,
          lower = result$lower,
          upper = result$upper,
          Design = designs[[design]],
          Outcome = outcomes[[outcome]],
          Model = model_type,
          Trial = trial_label,
          stringsAsFactors = FALSE
        )
      })
    }) |> unlist(recursive = FALSE)
  )
}



comparison_colors_v1 <- c(
  "Dose reduction\nvs\nNo dose change"              = "#000000",  # Blue
  "1-<50% reduction\nvs\nNo dose change"            = "#2CA02C",  # Green
  ">=50% reduction/initial dose*\nvs\nNo dose change"        = "#FF7F0E",  # Sky Blue
  "Dose increase\nvs\nNo dose change"               = "#D62728"  # Orange
  
)

plot_trial_results_wot_v1 <- function(df, title_lab, y_lab="Hazard ratios with 95% compatibility (‘confidence’) intervals") {
  ggplot(df, aes(x = Comparison, y = HR, color = Comparison)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.07, size = 0.7, position = position_dodge(width = 0.5)) +
    facet_wrap(~ Trial, ncol = 1, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    labs(
      title = title_lab,
      x = NULL,
      y = y_lab,
      color = " "
    ) +
    scale_color_manual(values = comparison_colors_v1) +
    # Optionally: scale_y_continuous(breaks = seq(y_min, y_max, by = y_by), limits = c(y_min, y_max)) +
    theme_minimal() +
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.x = element_line(linewidth = 1, color = "black"),
      axis.ticks.length  = unit(0.3, "cm"),
      axis.text.x = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 14),
      axis.title.y = element_text(size = 14),       # <-- added this line for y-axis label size
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 13)
    ) 
}




## Trial 1 #################################
res_df1_ipw <- extract_all_results(res11, 
                               trial_label = "Trial 1: MET interruption of 1-3 days",
                               model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    ))
  )

res_df1_ow <- extract_all_results(res21, 
                               trial_label = "Trial 1: MET interruption of 1-3 days",
                               model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    ))
  )

## Trial 2 #################################
res_df2_ipw <- extract_all_results(res12, 
                               trial_label = "Trial 2: MET interruption of 4 days",
                               model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change"))
  )

res_df2_ow <- extract_all_results(res22, 
                               trial_label = "Trial 2: MET interruption of 4 days",
                               model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change"))
  )

## Trial 3 #################################
res_df3_ipw <- extract_all_results(res13, 
                                   trial_label = "Trial 3: MET interruption of 5-14 days",
                                   model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change"))
    
  )

res_df3_ow <- extract_all_results(res23, 
                                   trial_label = "Trial 3: MET interruption of 5-14 days",
                                   model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change"))
    
  )


## Trial 4 #################################
res_df4_ipw <- extract_all_results(res14, 
                                   trial_label = "Trial 4: BNX interruption of 1-5 days",
                                   model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    ))
  )

res_df4_ow <- extract_all_results(res24, 
                                   trial_label = "Trial 4: BNX interruption of 1-5 days",
                                   model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    ))
  )

## Trial 5 #################################
res_df5_ipw <- extract_all_results(res15, 
                               trial_label = "Trial 5: BNX interruption of 6-14 days",
                               model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change"))
    
  )

res_df5_ow <- extract_all_results(res25, 
                                   trial_label = "Trial 5: BNX interruption of 6-14 days",
                                   model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change"))
    
  )



res_df <- bind_rows(res_df1_ipw, res_df2_ipw, res_df3_ipw, res_df4_ipw, res_df5_ipw,
                    res_df1_ow, res_df2_ow, res_df3_ow, res_df4_ow, res_df5_ow)
res_df <- res_df %>%
  mutate(
    Design = recode(Design, 
                    "Prevalent new* users" = "All prevalent users"),
    Design = factor(Design,
                    levels = c("Incident users",
                               "Prevalent new users",
                               "All prevalent users")
    ))


saveRDS(res_df, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_ow.rds")



# Quick graph ###################
comparison_colors_v4 <- c(
  "Dose reduction\nvs\nNo dose change"              = "#000000",  # Blue
  "1-<50% reduction\nvs\nNo dose change"            = "#2CA02C",  # Green
  ">=50% reduction/initial dose*\nvs\nNo dose change" = "#FF7F0E", 
  "Dose increase\nvs\nNo dose change"               = "#D62728"
)


plot_trial_results_wot_v4 <- function(df, title_lab, 
                                      y_lab="Hazard ratios with 95% CI",
                                      y_min, 
                                      y_max, 
                                      y_points,
                                      legend_option=TRUE) {
  
  df <- df %>% mutate(lower = ifelse(lower < y_min, y_min, lower),
                      upper = ifelse(upper > y_max, y_max, upper))
  
  ggplot(df, aes(x = Comparison, y = HR, color = Comparison)) +
    geom_point(position = position_dodge(width = 0.8), size = 3.2, stroke = 1) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0, size = 0.8, position = position_dodge(width = 0.5)) +
    facet_wrap(~ Trial, ncol = 1, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    labs(
      title = title_lab,
      x = NULL,
      y = y_lab,
      color = " "
    ) +
    scale_color_manual(values = comparison_colors_v4) +
    scale_y_log10(breaks = y_points, limits = c(y_min, y_max))+
    theme_minimal() +
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.text.x = element_blank(),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.text = element_text(color = "black", size = 14),
      axis.title.y = element_text(size = 14),       # <-- added this line for y-axis label size
      legend.position = if(legend_option) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 13)
    ) 
}

## iptw: od ###############

res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_ow.rds")

res_df <- res_df %>% mutate(Trial = 
                              recode(Trial,
                                     "Trial 1: Methadone interruption of 1-3 days"="Trial 1: MET interruption of 1-3 days",
                                     "Trial 2: Methadone interruption of 4 days"="Trial 2: MET interruption of 4 days",
                                     "Trial 3: Methadone interruption of 5-14 days"="Trial 3: MET interruption of 5-14 days",
                                     "Trial 4: Buprenorphine/naloxone interruption of 1-5 days"="Trial 4: BNX interruption of 1-5 days",
                                     "Trial 5: Buprenorphine/naloxone interruption of 6-14 days"="Trial 5: BNX interruption of 6-14 days"),
                            Comparison=
                              recode(Comparison,
                                     ">=50% reduction/<=40mg\nvs\nNo dose change"=">=50% reduction/initial dose*\nvs\nNo dose change",
                                     ">=50% reduction/<=4mg\nvs\nNo dose change"=">=50% reduction/initial dose*\nvs\nNo dose change")
)


od_df_incid = res_df %>% filter(Outcome == "Overdose-related acute care visits",
                                Model == "IPTW", 
                                Design == "Incident users") 

od_df_incid[c(od_df_incid$Trial == "Trial 2: MET interruption of 4 days" & od_df_incid$Design == "Incident users"), c("HR", "lower", "upper")] <- NA



od_df_incid_met = od_df_incid %>% filter(Trial %in% c("Trial 1: MET interruption of 1-3 days", 
                                                      "Trial 2: MET interruption of 4 days", 
                                                      "Trial 3: MET interruption of 5-14 days"))

od_df_incid_bnx = od_df_incid %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                      "Trial 5: BNX interruption of 6-14 days"))

od_plot_incid_met = plot_trial_results_wot_v4(df=od_df_incid_met, 
                                              title_lab = "Incident users",
                                              y_min = 0.45,
                                              y_max = 3.2, 
                                              y_points = c(0.5, 1, 1.5, 3),
                                              legend_option=FALSE)

od_plot_incid_bnx = plot_trial_results_wot_v4(df=od_df_incid_bnx, 
                                              title_lab = "Incident users",
                                              y_min = 0.45,
                                              y_max = 3.8, 
                                              y_points = c(0.5, 1, 1.5, 3))

od_df_prev = res_df %>% filter(Outcome == "Overdose-related acute care visits",
                               Model == "IPTW", 
                               Design == "Prevalent new users")




od_df_prev_met = od_df_prev %>% filter(Trial %in% c("Trial 1: MET interruption of 1-3 days", 
                                                    "Trial 2: MET interruption of 4 days", 
                                                    "Trial 3: MET interruption of 5-14 days"))

od_df_prev_bnx = od_df_prev %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                    "Trial 5: BNX interruption of 6-14 days"))

od_plot_prev_met = plot_trial_results_wot_v4(df=od_df_prev_met, 
                                             title_lab = "Prevalent new users",
                                             y_lab="",
                                             y_min = 0.45,
                                             y_max = 3.2, 
                                             y_points = c(0.5, 1, 1.5, 3),
                                             legend_option=FALSE)

od_plot_prev_bnx = plot_trial_results_wot_v4(df=od_df_prev_bnx, 
                                             title_lab = "Prevalent new users",
                                             y_lab="",
                                             y_min = 0.45,
                                             y_max = 3.2, 
                                             y_points = c(0.5, 1, 1.5, 3))



od_df_all = res_df %>% filter(Outcome == "Overdose-related acute care visits",
                              Model == "IPTW", 
                              Design == "All prevalent users")



od_df_all_met = od_df_all %>% filter(Trial %in% c("Trial 1: MET interruption of 1-3 days", 
                                                  "Trial 2: MET interruption of 4 days", 
                                                  "Trial 3: MET interruption of 5-14 days"))

od_df_all_bnx = od_df_all %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                  "Trial 5: BNX interruption of 6-14 days"))

od_plot_all_met = plot_trial_results_wot_v4(df=od_df_all_met, 
                                            title_lab = "All prevalent users",
                                            y_lab="",
                                            y_min = 0.45,
                                            y_max = 3.2, 
                                            y_points = c(0.5, 1, 1.5, 3),
                                            legend_option=FALSE)

od_plot_all_bnx = plot_trial_results_wot_v4(df=od_df_all_bnx, 
                                            title_lab = "All prevalent users",
                                            y_lab="",
                                            y_min = 0.45,
                                            y_max = 3.2, y_points = c(0.5, 1, 1.5, 3))




library(patchwork)
plot_od_met <- od_plot_incid_met + plot_spacer() + od_plot_prev_met + plot_spacer() + od_plot_all_met + # use `+` to arrange side by side
  plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1, 0.1)) 

plot_od_bnx <- od_plot_incid_bnx + plot_spacer() + od_plot_prev_bnx + plot_spacer() + od_plot_all_bnx + # use `+` to arrange side by side
  plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1, 0.1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_od <- plot_spacer() /  # top spacer row
  plot_od_met / 
  plot_spacer()/ # spacer between met and bnx
  plot_od_bnx +
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
  guides(color = guide_legend(nrow = 1))



plot_od1 <- ggdraw(plot_od) + 
  draw_label("* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
             x = 0.99, y = 0.01, hjust = 1, vjust = 0, 
             size = 12, color = "black") +
  draw_label("A. Methadone (MET) trials", 
             x = 0.02, y = 0.97, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone (BNX) trials", 
             x = 0.02, y = 0.436, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("NA for few events", 
             x = 0.13, y = 0.695, hjust = 0, vjust = 0, 
             size = 12, color = "black")



ggsave("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/od_plot_Nov_12_PLR.pdf", 
       plot = plot_od1, width = 13, height = 12, units = "in", dpi = 600)



name = "od_plot_Nov_12_PLR"
figure=paste0("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=4000, height=4200, res=300)
plot_od1
dev.off()





## iptw: dis ###############

dis_df_incid = res_df %>% filter(Outcome == "OAT discontinuation",
                                 Model == "IPTW", 
                                 Design == "Incident users")


dis_df_incid_met = dis_df_incid %>% filter(Trial %in% c("Trial 1: MET interruption of 1-3 days", 
                                                        "Trial 2: MET interruption of 4 days", 
                                                        "Trial 3: MET interruption of 5-14 days"))

dis_df_incid_bnx = dis_df_incid %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                        "Trial 5: BNX interruption of 6-14 days"))

dis_plot_incid_met = plot_trial_results_wot_v4(df=dis_df_incid_met, 
                                               title_lab = "Incident users", 
                                               y_min = 0.45,
                                               y_max = 2.5, 
                                               y_points = c(0.5, 1, 1.5),
                                               legend_option=FALSE)

dis_plot_incid_bnx = plot_trial_results_wot_v4(df=dis_df_incid_bnx, 
                                               title_lab = "Incident users",
                                               y_min = 0.45,
                                               y_max = 2.5, 
                                               y_points = c(0.5, 1, 1.5))







dis_df_prev = res_df %>% filter(Outcome == "OAT discontinuation",
                                Model == "IPTW", 
                                Design == "Prevalent new users")


dis_df_prev_met = dis_df_prev %>% filter(Trial %in% c("Trial 1: MET interruption of 1-3 days", 
                                                      "Trial 2: MET interruption of 4 days", 
                                                      "Trial 3: MET interruption of 5-14 days"))

dis_df_prev_bnx = dis_df_prev %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                      "Trial 5: BNX interruption of 6-14 days"))

dis_plot_prev_met = plot_trial_results_wot_v4(df=dis_df_prev_met, 
                                              title_lab = "Prevalent new users", 
                                              y_min = 0.45,
                                              y_max = 2.5, 
                                              y_points = c(0.5, 1, 1.5),
                                              legend_option=FALSE)

dis_plot_prev_bnx = plot_trial_results_wot_v4(df=dis_df_incid_bnx, 
                                              title_lab = "Prevalent new users", 
                                              y_min = 0.45,
                                              y_max = 2.5, 
                                              y_points = c(0.5, 1, 1.5))


dis_df_all = res_df %>% filter(Outcome == "OAT discontinuation",
                               Model == "IPTW", 
                               Design == "All prevalent users")


dis_df_all_met = dis_df_all %>% filter(Trial %in% c("Trial 1: MET interruption of 1-3 days", 
                                                    "Trial 2: MET interruption of 4 days", 
                                                    "Trial 3: MET interruption of 5-14 days"))

dis_df_all_bnx = dis_df_all %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                                    "Trial 5: BNX interruption of 6-14 days"))

dis_plot_all_met = plot_trial_results_wot_v4(df=dis_df_all_met, 
                                             title_lab = "All prevalent users", 
                                             y_min = 0.45,
                                             y_max = 2.5, 
                                             y_points = c(0.5, 1, 1.5),
                                             legend_option=FALSE)

dis_plot_all_bnx = plot_trial_results_wot_v4(df=dis_df_all_bnx, 
                                             title_lab = "All prevalent users", 
                                             y_min = 0.45,
                                             y_max = 2.5, 
                                             y_points = c(0.5, 1, 1.5))



library(patchwork)
plot_dis_met <- dis_plot_incid_met + plot_spacer() + dis_plot_prev_met + plot_spacer() + dis_plot_all_met + # use `+` to arrange side by side
  plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1, 0.1)) 

plot_dis_bnx <- dis_plot_incid_bnx + plot_spacer() + dis_plot_prev_bnx + plot_spacer() + dis_plot_all_bnx + # use `+` to arrange side by side
  plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1, 0.1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

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
  guides(color = guide_legend(nrow = 1))

plot_dis1 <- ggdraw(plot_dis) + 
  draw_label("* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
             x = 0.99, y = 0.01, hjust = 1, vjust = 0, 
             size = 12, color = "black") +
  draw_label("A. Methadone (MET) trials", 
             x = 0.02, y = 0.97, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Buprenorphine/naloxone (BNX) trials", 
             x = 0.02, y = 0.436, hjust = 0, vjust = 0, 
             size = 18, color = "black")



ggsave("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/dis_plot_Nov_12_PLR.pdf", 
       plot = plot_dis1, width = 13, height = 12, units = "in", dpi = 600)



name = "dis_plot_Nov_12_PLR"
figure=paste0("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=4000, height=4200, res=300)
plot_dis1
dev.off()







## Result section (HR) #############

res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/all_HRs.rds")
od_df1 = res_df %>% filter(Outcome == "Overdose-related acute care visits",
                           Model == "IPTW")
dis_df1 = res_df %>% filter(Outcome == "OAT discontinuation",
                            Model == "IPTW")

### od: dose reduction ################

(harm = od_df1 %>% filter(Comparison != "Dose increase\nvs\nNo dose change") %>% 
   filter(lower > 1) %>% select(Comparison, HR, lower, upper, Design, Trial))

(benefit = od_df1 %>% filter(Comparison != "Dose increase\nvs\nNo dose change") %>% 
    filter(upper < 1) %>% select(Comparison, HR, lower, upper, Design, Trial))

(noeffect = od_df1 %>% filter(Comparison != "Dose increase\nvs\nNo dose change") %>% 
    filter(lower< 1, upper > 1) %>% select(Comparison, HR, lower, upper, Design, Trial))


### od: dose increase ##############

(harm = od_df1 %>% filter(Comparison == "Dose increase\nvs\nNo dose change") %>% 
   filter(lower > 1) %>% select(Comparison, HR, lower, upper, Design, Trial))

write.csv(harm, "harm.csv")

(benefit = od_df1 %>% filter(Comparison == "Dose increase\nvs\nNo dose change") %>% 
    filter(upper < 1) %>% select(Comparison, HR, lower, upper, Design, Trial))
write.csv(harm, "harm.csv")
write.csv(benefit, "benefit.csv")

(noeffect = od_df1 %>% filter(Comparison == "Dose increase\nvs\nNo dose change") %>% 
    filter(lower< 1, upper > 1) %>% select(Comparison, HR, lower, upper, Design, Trial))

write.csv(harm, "harm.csv")
write.csv(benefit, "benefit.csv")
write.csv(noeffect, "noeffect.csv")


od_df1 = res_df %>% filter(Outcome == "Overdose-related acute care visits",
                           Model == "IPTW")
dis_df1 = res_df %>% filter(Outcome == "OAT discontinuation",
                            Model == "IPTW")

### dis: dose reduction ################

(harm = dis_df1 %>% filter(Comparison != "Dose increase\nvs\nNo dose change") %>% 
   filter(lower > 1) %>% select(Comparison, HR, lower, upper, Design, Trial))

(benefit = dis_df1 %>% filter(Comparison != "Dose increase\nvs\nNo dose change") %>% 
    filter(upper < 1) %>% select(Comparison, HR, lower, upper, Design, Trial))

(noeffect = dis_df1 %>% filter(Comparison != "Dose increase\nvs\nNo dose change") %>% 
    filter(lower< 1, upper > 1) %>% select(Comparison, HR, lower, upper, Design, Trial))


write.csv(harm, "harm.csv")
write.csv(benefit, "benefit.csv")
write.csv(noeffect, "noeffect.csv")


### dis: dose increase ##############

(harm = dis_df1 %>% filter(Comparison == "Dose increase\nvs\nNo dose change") %>% 
   filter(lower > 1) %>% select(Comparison, HR, lower, upper, Design, Trial))

(benefit = dis_df1 %>% filter(Comparison == "Dose increase\nvs\nNo dose change") %>% 
    filter(upper < 1) %>% select(Comparison, HR, lower, upper, Design, Trial))

(noeffect = dis_df1 %>% filter(Comparison == "Dose increase\nvs\nNo dose change") %>% 
    filter(lower< 1, upper > 1) %>% select(Comparison, HR, lower, upper, Design, Trial))

write.csv(harm, "harm.csv")
write.csv(benefit, "benefit.csv")
write.csv(noeffect, "noeffect.csv")




## Result section (CID) ##############
risk_all <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/All_risk.rds")
names(risk_all)
View(risk_all)
risk_all %>% count(group)

risk_all %>% filter(Outcome == "Overdose-related acute care visits",
                    time_day %in% c(90, 365),
                    Design == "Prevalent new users",
                    Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  select(time_day, group, risk, lower, upper)

risk_all %>% filter(Outcome == "Overdose-related acute care visits",
                    time_day %in% c(365),
                    #Design == "Prevalent new users",
                    Trial == "Trial 4: BNX interruption of 1-5 days") %>% 
  select(time_day, group, risk, lower, upper, Design)

risk_all %>% filter(Outcome == "Overdose-related acute care visits",
                    time_day %in% c(365),
                    #Design == "Prevalent new users",
                    Trial == "Trial 5: BNX interruption of 6-14 days") %>% 
  select(time_day, group, risk, lower, upper, Design)



risk_all %>% filter(Outcome == "OAT discontinuation",
                    time_day %in% c(365),
                    #Design == "Prevalent new users",
                    Trial == "Trial 1: MET interruption of 1-3 days") %>% 
  select(time_day, group, risk, lower, upper, Design)



risk_all %>% filter(Outcome == "OAT discontinuation",
                    time_day %in% c(365),
                    #Design == "Prevalent new users",
                    Trial == "Trial 3: MET interruption of 5-14 days") %>% 
  select(time_day, group, risk, lower, upper, Design)





# Forestplotter ##################

extract_all_results <- function(res_obj, trial_label = "Trial 1", model_type = "IPTW") {
  designs <- c("incid" = "Incident users", 
               "prev" = "Prevalent new users", 
               "rzero" = "Prevalent new* users")
  outcomes <- c("od" = "Overdose-related acute care visits", 
                "dis" = "OAT discontinuation")
  
  bind_rows(
    lapply(names(outcomes), function(outcome) {
      lapply(names(designs), function(design) {
        result <- res_obj[[outcome]][[design]]
        if (is.null(result)) return(NULL)
        
        data.frame(
          Comparison = result$names,
          HR = result$OR,
          lower = result$lower,
          upper = result$upper,
          Design = designs[[design]],
          Outcome = outcomes[[outcome]],
          Model = model_type,
          Trial = trial_label,
          stringsAsFactors = FALSE
        )
      })
    }) |> unlist(recursive = FALSE)
  )
}
extract_all_results <- function(res_obj, trial_label = "Trial 1", model_type = "IPTW") {
  designs <- c("incid" = "Incident users", 
               "prev" = "Prevalent new users", 
               "rzero" = "Prevalent new* users")
  outcomes <- c("od" = "Overdose-related acute care visits", 
                "dis" = "OAT discontinuation")
  
  safe_extract <- function(res) {
    # Return NULL if res is NULL or empty list
    if (is.null(res) || length(res) == 0) return(NULL)
    
    # Ensure required fields exist
    required <- c("names", "OR", "lower", "upper")
    if (!all(required %in% names(res))) return(NULL)
    
    # Build output row
    data.frame(
      Comparison = res$names,
      HR = res$OR,
      lower = res$lower,
      upper = res$upper,
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(
    lapply(names(outcomes), function(outcome) {
      lapply(names(designs), function(design) {
        result <- safe_extract(res_obj[[outcome]][[design]])
        if (is.null(result)) return(NULL)
        
        cbind(
          result,
          Design = designs[[design]],
          Outcome = outcomes[[outcome]],
          Model = model_type,
          Trial = trial_label,
          stringsAsFactors = FALSE
        )
      })
    }) |> unlist(recursive = FALSE)
  )
}


## Read results Nov 12 #############
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"
res1_lt30d_iptw <- readRDS(file=paste0(folder, "res1_12d_iptw.rds"))
res2_lt30d_iptw <- readRDS(file=paste0(folder, "res2_34d_iptw.rds"))
res3_lt30d_iptw <- readRDS(file=paste0(folder, "res3_iptw.rds"))
res4_lt30d_iptw <- readRDS(file=paste0(folder, "res4_iptw.rds"))
res5_lt30d_iptw <- readRDS(file=paste0(folder, "res5_iptw.rds"))

res1_lt30d_ow <- readRDS(file=paste0(folder, "res1_12d_ow.rds"))
res2_lt30d_ow <- readRDS(file=paste0(folder, "res2_34d_ow.rds"))
res3_lt30d_ow <- readRDS(file=paste0(folder, "res3_ow.rds"))
res4_lt30d_ow <- readRDS(file=paste0(folder, "res4_ow.rds"))
res5_lt30d_ow <- readRDS(file=paste0(folder, "res5_ow.rds"))

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_30_lt90d/"
res1_30_lt90d_iptw <- readRDS(file=paste0(folder, "res1_12d_iptw.rds"))
res2_30_lt90d_iptw <- readRDS(file=paste0(folder, "res2_34d_iptw.rds"))
res3_30_lt90d_iptw <- readRDS(file=paste0(folder, "res3_iptw.rds"))
res4_30_lt90d_iptw <- readRDS(file=paste0(folder, "res4_iptw.rds"))
res5_30_lt90d_iptw <- readRDS(file=paste0(folder, "res5_iptw.rds"))

res1_30_lt90d_ow <- readRDS(file=paste0(folder, "res1_12d_ow.rds"))
res2_30_lt90d_ow <- readRDS(file=paste0(folder, "res2_34d_ow.rds"))
res3_30_lt90d_ow <- readRDS(file=paste0(folder, "res3_ow.rds"))
res4_30_lt90d_ow <- readRDS(file=paste0(folder, "res4_ow.rds"))
res5_30_lt90d_ow <- readRDS(file=paste0(folder, "res5_ow.rds"))


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_gt90d/"
res1_gt90d_iptw <- readRDS(file=paste0(folder, "res1_12d_iptw.rds"))
res2_gt90d_iptw <- readRDS(file=paste0(folder, "res2_34d_iptw.rds"))
res3_gt90d_iptw <- readRDS(file=paste0(folder, "res3_iptw.rds"))
res4_gt90d_iptw <- readRDS(file=paste0(folder, "res4_iptw.rds"))
res5_gt90d_iptw <- readRDS(file=paste0(folder, "res5_iptw.rds"))

res1_gt90d_ow <- readRDS(file=paste0(folder, "res1_12d_ow.rds"))
res2_gt90d_ow <- readRDS(file=paste0(folder, "res2_34d_ow.rds"))
res3_gt90d_ow <- readRDS(file=paste0(folder, "res3_ow.rds"))
res4_gt90d_ow <- readRDS(file=paste0(folder, "res4_ow.rds"))
res5_gt90d_ow <- readRDS(file=paste0(folder, "res5_ow.rds"))


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
res1_full_iptw <- readRDS(file=paste0(folder, "res1_12d_iptw.rds"))
res2_full_iptw <- readRDS(file=paste0(folder, "res2_34d_iptw.rds"))
res3_full_iptw <- readRDS(file=paste0(folder, "res3_iptw.rds"))
res4_full_iptw <- readRDS(file=paste0(folder, "res4_iptw.rds"))
res5_full_iptw <- readRDS(file=paste0(folder, "res5_iptw.rds"))

res1_full_ow <- readRDS(file=paste0(folder, "res1_12d_ow.rds"))
res2_full_ow <- readRDS(file=paste0(folder, "res2_34d_ow.rds"))
res3_full_ow <- readRDS(file=paste0(folder, "res3_ow.rds"))
res4_full_ow <- readRDS(file=paste0(folder, "res4_ow.rds"))
res5_full_ow <- readRDS(file=paste0(folder, "res5_ow.rds"))





comparison_colors_v1 <- c(
  "Dose reduction\nvs\nNo dose change"              = "#000000",  # Blue
  "1-<50% reduction\nvs\nNo dose change"            = "#2CA02C",  # Green
  ">=50% reduction/initial dose*\nvs\nNo dose change"        = "#FF7F0E",  # Sky Blue
  "Dose increase\nvs\nNo dose change"               = "#D62728"  # Orange
  
)

plot_trial_results_wot_v1 <- function(df, title_lab, y_lab="Hazard ratios with 95% compatibility (‘confidence’) intervals") {
  ggplot(df, aes(x = Comparison, y = HR, color = Comparison)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.07, size = 0.7, position = position_dodge(width = 0.5)) +
    facet_wrap(~ Trial, ncol = 1, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    labs(
      title = title_lab,
      x = NULL,
      y = y_lab,
      color = " "
    ) +
    scale_color_manual(values = comparison_colors_v1) +
    # Optionally: scale_y_continuous(breaks = seq(y_min, y_max, by = y_by), limits = c(y_min, y_max)) +
    theme_minimal() +
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.x = element_line(linewidth = 1, color = "black"),
      axis.ticks.length  = unit(0.3, "cm"),
      axis.text.x = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 14),
      axis.title.y = element_text(size = 14),       # <-- added this line for y-axis label size
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 13)
    ) 
}


### Trial 1 #################################
res_df1_lt30d_iptw <- extract_all_results(res1_lt30d_iptw, 
                                          trial_label = "Trial 1: MET interruption of 1-2 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Within 30 days of OAT initiation"
  )



res_df1_lt30d_ow <- extract_all_results(res1_lt30d_ow, 
                                        trial_label = "Trial 1: MET interruption of 1-2 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Within 30 days of OAT initiation"
  )



res_df1_30_lt90d_iptw <- extract_all_results(res1_30_lt90d_iptw, 
                                             trial_label = "Trial 1: MET interruption of 1-2 days",
                                             model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Within 31-90 days of OAT initiation"
  )

res_df1_30_lt90d_ow <- extract_all_results(res1_30_lt90d_ow, 
                                           trial_label = "Trial 1: MET interruption of 1-2 days",
                                           model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Within 31-90 days of OAT initiation"
  )

res_df1_gt90d_iptw <- extract_all_results(res1_gt90d_iptw, 
                                          trial_label = "Trial 1: MET interruption of 1-2 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "After 90 days of OAT initiation"
  )

res_df1_gt90d_ow <- extract_all_results(res1_gt90d_ow, 
                                        trial_label = "Trial 1: MET interruption of 1-2 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "After 90 days of OAT initiation"
  )



res_df1_full_iptw <- extract_all_results(res1_full_iptw, 
                                         trial_label = "Trial 1: MET interruption of 1-2 days",
                                         model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Primary analysis"
  )

res_df1_full_ow <- extract_all_results(res1_full_ow, 
                                       trial_label = "Trial 1: MET interruption of 1-2 days",
                                       model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Primary analysis"
  )


### Trial 2 #################################
res_df2_lt30d_iptw <- extract_all_results(res2_lt30d_iptw, 
                                          trial_label = "Trial 2: MET interruption of 3-4 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 30 days of OAT initiation"
  )

res_df2_lt30d_ow <- extract_all_results(res2_lt30d_ow, 
                                        trial_label = "Trial 2: MET interruption of 3-4 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 30 days of OAT initiation"
  )

res_df2_30_lt90d_iptw <- extract_all_results(res2_30_lt90d_iptw, 
                                             trial_label = "Trial 2: MET interruption of 3-4 days",
                                             model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 31-90 days of OAT initiation"
  )

res_df2_30_lt90d_ow <- extract_all_results(res2_30_lt90d_ow, 
                                           trial_label = "Trial 2: MET interruption of 3-4 days",
                                           model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 31-90 days of OAT initiation"
  )



res_df2_gt90d_iptw <- extract_all_results(res2_gt90d_iptw, 
                                          trial_label = "Trial 2: MET interruption of 3-4 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "After 90 days of OAT initiation"
  )

res_df2_gt90d_ow <- extract_all_results(res2_gt90d_ow, 
                                        trial_label = "Trial 2: MET interruption of 3-4 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "After 90 days of OAT initiation"
  )



res_df2_full_iptw <- extract_all_results(res2_full_iptw, 
                                         trial_label = "Trial 2: MET interruption of 3-4 days",
                                         model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Primary analysis"
  )

res_df2_full_ow <- extract_all_results(res2_full_ow, 
                                       trial_label = "Trial 2: MET interruption of 3-4 days",
                                       model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Primary analysis"
  )



### Trial 3 #################################
res_df3_lt30d_iptw <- extract_all_results(res3_lt30d_iptw, 
                                          trial_label = "Trial 3: MET interruption of 5-14 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 30 days of OAT initiation"
  )

res_df3_lt30d_ow <- extract_all_results(res3_lt30d_ow, 
                                        trial_label = "Trial 3: MET interruption of 5-14 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 30 days of OAT initiation"
  )

res_df3_30_lt90d_iptw <- extract_all_results(res3_30_lt90d_iptw, 
                                             trial_label = "Trial 3: MET interruption of 5-14 days",
                                             model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 31-90 days of OAT initiation"
  )

res_df3_30_lt90d_ow <- extract_all_results(res3_30_lt90d_ow, 
                                           trial_label = "Trial 3: MET interruption of 5-14 days",
                                           model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 31-90 days of OAT initiation"
  )



res_df3_gt90d_iptw <- extract_all_results(res3_gt90d_iptw, 
                                          trial_label = "Trial 3: MET interruption of 5-14 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "After 90 days of OAT initiation"
  )

res_df3_gt90d_ow <- extract_all_results(res3_gt90d_ow, 
                                        trial_label = "Trial 3: MET interruption of 5-14 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "After 90 days of OAT initiation"
  )



res_df3_full_iptw <- extract_all_results(res3_full_iptw, 
                                         trial_label = "Trial 3: MET interruption of 5-14 days",
                                         model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Primary analysis"
  )

res_df3_full_ow <- extract_all_results(res3_full_ow, 
                                       trial_label = "Trial 3: MET interruption of 5-14 days",
                                       model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=40mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=40mg" = ">=50% reduction/<=40mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=40mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Primary analysis"
  )



### Trial 4 #################################
res_df4_lt30d_iptw <- extract_all_results(res4_lt30d_iptw, 
                                          trial_label = "Trial 4: BNX interruption of 1-5 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Within 30 days of OAT initiation"
  )

res_df4_lt30d_ow <- extract_all_results(res4_lt30d_ow, 
                                        trial_label = "Trial 4: BNX interruption of 1-5 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Within 30 days of OAT initiation"
  )

res_df4_30_lt90d_iptw <- extract_all_results(res4_30_lt90d_iptw, 
                                             trial_label = "Trial 4: BNX interruption of 1-5 days",
                                             model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Within 31-90 days of OAT initiation"
  )

res_df4_30_lt90d_ow <- extract_all_results(res4_30_lt90d_ow, 
                                           trial_label = "Trial 4: BNX interruption of 1-5 days",
                                           model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Within 31-90 days of OAT initiation"
  )



res_df4_gt90d_iptw <- extract_all_results(res4_gt90d_iptw, 
                                          trial_label = "Trial 4: BNX interruption of 1-5 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "After 90 days of OAT initiation"
  )

res_df4_gt90d_ow <- extract_all_results(res4_gt90d_ow, 
                                        trial_label = "Trial 4: BNX interruption of 1-5 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "After 90 days of OAT initiation"
  )



res_df4_full_iptw <- extract_all_results(res4_full_iptw, 
                                         trial_label = "Trial 4: BNX interruption of 1-5 days",
                                         model_type = "IPTW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Primary analysis"
  )

res_df4_full_ow <- extract_all_results(res4_full_ow, 
                                       trial_label = "Trial 4: BNX interruption of 1-5 days",
                                       model_type = "OW") %>%
  filter(Comparison %in% c("trtreversion", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trtreversion" = "Dose reduction\nvs\nNo dose change",
                        "trtincrease"  = "Dose increase\nvs\nNo dose change"
    ),
    Comparison = factor(Comparison, levels = c(
      "Dose reduction\nvs\nNo dose change",
      "Dose increase\nvs\nNo dose change"
    )),
    res_type = "Primary analysis"
  )


### Trial 5 #################################
res_df5_lt30d_iptw <- extract_all_results(res5_lt30d_iptw, 
                                          trial_label = "Trial 5: BNX interruption of 6-14 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 30 days of OAT initiation"
  )

res_df5_lt30d_ow <- extract_all_results(res5_lt30d_ow, 
                                        trial_label = "Trial 5: BNX interruption of 6-14 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 30 days of OAT initiation"
  )

res_df5_30_lt90d_iptw <- extract_all_results(res5_30_lt90d_iptw, 
                                             trial_label = "Trial 5: BNX interruption of 6-14 days",
                                             model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 31-90 days of OAT initiation"
  )

res_df5_30_lt90d_ow <- extract_all_results(res5_30_lt90d_ow, 
                                           trial_label = "Trial 5: BNX interruption of 6-14 days",
                                           model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Within 31-90 days of OAT initiation"
  )



res_df5_gt90d_iptw <- extract_all_results(res5_gt90d_iptw, 
                                          trial_label = "Trial 5: BNX interruption of 6-14 days",
                                          model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "After 90 days of OAT initiation"
  )

res_df5_gt90d_ow <- extract_all_results(res5_gt90d_ow, 
                                        trial_label = "Trial 5: BNX interruption of 6-14 days",
                                        model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "After 90 days of OAT initiation"
  )



res_df5_full_iptw <- extract_all_results(res5_full_iptw, 
                                         trial_label = "Trial 5: BNX interruption of 6-14 days",
                                         model_type = "IPTW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Primary analysis"
  )

res_df5_full_ow <- extract_all_results(res5_full_ow, 
                                       trial_label = "Trial 5: BNX interruption of 6-14 days",
                                       model_type = "OW") %>%
  filter(Comparison %in% c("trt1-<50%", "trt>=50%/<=4mg", "trtincrease")) %>%
  mutate(
    Comparison = recode(Comparison,
                        "trt1-<50%" = "1-<50% reduction\nvs\nNo dose change",
                        "trt>=50%/<=4mg" = ">=50% reduction/<=4mg\nvs\nNo dose change",
                        "trtincrease" = "Dose increase\nvs\nNo dose change"),
    Comparison = factor(Comparison,
                        levels = c("1-<50% reduction\nvs\nNo dose change",
                                   ">=50% reduction/<=4mg\nvs\nNo dose change",
                                   "Dose increase\nvs\nNo dose change")),
    res_type = "Primary analysis"
  )


res_df <- bind_rows(res_df1_full_ow, res_df2_full_ow, res_df3_full_ow, res_df4_full_ow, res_df5_full_ow,
                    res_df1_full_iptw, res_df2_full_iptw, res_df3_full_iptw, res_df4_full_iptw, res_df5_full_iptw,
                    res_df1_lt30d_ow, res_df2_lt30d_ow, res_df3_lt30d_ow, res_df4_lt30d_ow, res_df5_lt30d_ow,
                    res_df1_lt30d_iptw, res_df2_lt30d_iptw, res_df3_lt30d_iptw, res_df4_lt30d_iptw, res_df5_lt30d_iptw,
                    res_df1_30_lt90d_ow, res_df2_30_lt90d_ow, res_df3_30_lt90d_ow, res_df4_30_lt90d_ow, res_df5_30_lt90d_ow,
                    res_df1_30_lt90d_iptw, res_df2_30_lt90d_iptw, res_df3_30_lt90d_iptw, res_df4_30_lt90d_iptw, res_df5_30_lt90d_iptw,
                    res_df1_gt90d_ow, res_df2_gt90d_ow, res_df3_gt90d_ow, res_df4_gt90d_ow, res_df5_gt90d_ow,
                    res_df1_gt90d_iptw, res_df2_gt90d_iptw, res_df3_gt90d_iptw, res_df4_gt90d_iptw, res_df5_gt90d_iptw)

res_df <- res_df %>%
  mutate(
    Design = recode(Design, 
                    "Prevalent new* users" = "All prevalent users"),
    Design = factor(Design,
                    levels = c("Incident users",
                               "Prevalent new users",
                               "All prevalent users")
    )) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design)



saveRDS(res_df, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change")
  ) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design)

### od: primary ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"

#### Trial 1: od: ow #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_ow_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                         Model == "OW",
                                         Outcome == "Overdose-related acute care visits",
                                         Trial == "Trial 1: MET interruption of 1-2 days",
                                         Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_ow_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                          Model == "OW",
                                          Outcome == "Overdose-related acute care visits",
                                          Trial == "Trial 1: MET interruption of 1-2 days",
                                          Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_ow_od <- bind_rows(header1, tr1_full_ow_od_red, tr1_full_ow_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: ow #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_ow_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "OW",
                                             Outcome == "Overdose-related acute care visits",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_ow_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "OW",
                                             Outcome == "Overdose-related acute care visits",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_ow_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                          Model == "OW",
                                          Outcome == "Overdose-related acute care visits",
                                          Trial == "Trial 2: MET interruption of 3-4 days",
                                          Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr2_full_ow_od <- bind_rows(header2, tr2_full_ow_od_redlt50, tr2_full_ow_od_redgt50, tr2_full_ow_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: ow #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_ow_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "OW",
                                             Outcome == "Overdose-related acute care visits",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_ow_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "OW",
                                             Outcome == "Overdose-related acute care visits",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_ow_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                          Model == "OW",
                                          Outcome == "Overdose-related acute care visits",
                                          Trial == "Trial 3: MET interruption of 5-14 days",
                                          Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_ow_od <- bind_rows(header3, tr3_full_ow_od_redlt50, tr3_full_ow_od_redgt50, tr3_full_ow_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 4: od: ow #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_ow_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                         Model == "OW",
                                         Outcome == "Overdose-related acute care visits",
                                         Trial == "Trial 4: BNX interruption of 1-5 days",
                                         Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_ow_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                          Model == "OW",
                                          Outcome == "Overdose-related acute care visits",
                                          Trial == "Trial 4: BNX interruption of 1-5 days",
                                          Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_ow_od <- bind_rows(header4, tr4_full_ow_od_red, tr4_full_ow_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: ow #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_ow_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "OW",
                                             Outcome == "Overdose-related acute care visits",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_ow_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "OW",
                                             Outcome == "Overdose-related acute care visits",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_ow_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                          Model == "OW",
                                          Outcome == "Overdose-related acute care visits",
                                          Trial == "Trial 5: BNX interruption of 6-14 days",
                                          Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_ow_od <- bind_rows(header5, tr5_full_ow_od_redlt50, tr5_full_ow_od_redgt50, tr5_full_ow_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
full_ow_od <- bind_rows(tr1_full_ow_od,
                        tr2_full_ow_od,
                        tr3_full_ow_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_ow_od$`Overdose-related acute care visits` <- "" 
full_ow_od$`HR(95% CI)` <- full_ow_od$HR_display
full_ow_od$`Dose adjustment after interruption` <- full_ow_od$Comparison
full_ow_od$`User design` <- full_ow_od$Design


library(ggplot2)
library(forestploter)
library(tidyverse)
p_od_met <- forest(
  full_ow_od[, c("Dose adjustment after interruption", 
                 "User design", "HR(95% CI)", "Overdose-related acute care visits")],
  est   = full_ow_od$HR,
  lower = full_ow_od$lower,
  upper = full_ow_od$upper,
  ci_column = 4,
  ref_line = 1,
  arrow_lab = c("Beneficial", "Harmful"),
  x_trans = "log",
  xlim = c(0.45, 3),
  ticks_at = c(0.5, 1, 1.5),
  #xlab = "Adjusted hazard ratio",
  title = "A. Methadone (MET) trials",
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 16,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


#### bnx forestplot ##########################
full_ow_od <- bind_rows(tr4_full_ow_od,
                        tr5_full_ow_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_ow_od$`Overdose-related acute care visits` <- "" 
full_ow_od$`HR(95% CI)` <- full_ow_od$HR_display
full_ow_od$`Dose adjustment after interruption` <- full_ow_od$Comparison
full_ow_od$`User design` <- full_ow_od$Design

p_od_bnx <- forest(
  full_ow_od[, c("Dose adjustment after interruption", 
                 "User design", "HR(95% CI)", "Overdose-related acute care visits")],
  est   = full_ow_od$HR,
  lower = full_ow_od$lower,
  upper = full_ow_od$upper,
  ci_column = 4,
  ref_line = 1,
  x_trans = "log",
  xlim = c(0.45, 3),
  ticks_at = c(0.5, 1, 1.5),
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5\n** NA for few events", 
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 16,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))





library(patchwork)
plot_od <- plot_spacer() /  # top spacer row
  p_od_met / 
  plot_spacer() /
  p_od_bnx +
  plot_layout(
    heights = c(0.09, 1, 0.02, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))




name = "full_od_ow_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3000, height=3900, res=300)
plot_od
dev.off()


### dis: primary ######

#### Trial 1: od: ow #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_ow_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                          Model == "OW",
                                          Outcome == "OAT discontinuation",
                                          Trial == "Trial 1: MET interruption of 1-2 days",
                                          Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr1_full_ow_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "OW",
                                           Outcome == "OAT discontinuation",
                                           Trial == "Trial 1: MET interruption of 1-2 days",
                                           Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_ow_dis <- bind_rows(header1, tr1_full_ow_dis_red, tr1_full_ow_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: ow #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_ow_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                              Model == "OW",
                                              Outcome == "OAT discontinuation",
                                              Trial == "Trial 2: MET interruption of 3-4 days",
                                              Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_ow_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                              Model == "OW",
                                              Outcome == "OAT discontinuation",
                                              Trial == "Trial 2: MET interruption of 3-4 days",
                                              Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_ow_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "OW",
                                           Outcome == "OAT discontinuation",
                                           Trial == "Trial 2: MET interruption of 3-4 days",
                                           Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_ow_dis <- bind_rows(header2, tr2_full_ow_dis_redlt50, tr2_full_ow_dis_redgt50, tr2_full_ow_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

#### Trial 3: od: ow #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_ow_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                              Model == "OW",
                                              Outcome == "OAT discontinuation",
                                              Trial == "Trial 3: MET interruption of 5-14 days",
                                              Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_ow_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                              Model == "OW",
                                              Outcome == "OAT discontinuation",
                                              Trial == "Trial 3: MET interruption of 5-14 days",
                                              Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_ow_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "OW",
                                           Outcome == "OAT discontinuation",
                                           Trial == "Trial 3: MET interruption of 5-14 days",
                                           Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_ow_dis <- bind_rows(header3, tr3_full_ow_dis_redlt50, tr3_full_ow_dis_redgt50, tr3_full_ow_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 4: od: ow #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_ow_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                          Model == "OW",
                                          Outcome == "OAT discontinuation",
                                          Trial == "Trial 4: BNX interruption of 1-5 days",
                                          Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_ow_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "OW",
                                           Outcome == "OAT discontinuation",
                                           Trial == "Trial 4: BNX interruption of 1-5 days",
                                           Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_ow_dis <- bind_rows(header4, tr4_full_ow_dis_red, tr4_full_ow_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: ow #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_ow_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                              Model == "OW",
                                              Outcome == "OAT discontinuation",
                                              Trial == "Trial 5: BNX interruption of 6-14 days",
                                              Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_ow_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                              Model == "OW",
                                              Outcome == "OAT discontinuation",
                                              Trial == "Trial 5: BNX interruption of 6-14 days",
                                              Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_ow_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "OW",
                                           Outcome == "OAT discontinuation",
                                           Trial == "Trial 5: BNX interruption of 6-14 days",
                                           Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_ow_dis <- bind_rows(header5, tr5_full_ow_dis_redlt50, tr5_full_ow_dis_redgt50, tr5_full_ow_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
full_ow_dis <- bind_rows(tr1_full_ow_dis,
                         tr2_full_ow_dis,
                         tr3_full_ow_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_ow_dis$`OAT discontinuation` <- "" 
full_ow_dis$`HR(95% CI)` <- full_ow_dis$HR_display
full_ow_dis$`Dose adjustment after interruption` <- full_ow_dis$Comparison
full_ow_dis$`User design` <- full_ow_dis$Design


library(ggplot2)
library(forestploter)
library(tidyverse)
p_dis_met <- forest(
  full_ow_dis[, c("Dose adjustment after interruption", 
                  "User design", "HR(95% CI)", "OAT discontinuation")],
  est   = full_ow_dis$HR,
  lower = full_ow_dis$lower,
  upper = full_ow_dis$upper,
  ci_column = 4,
  ref_line = 1,
  arrow_lab = c("Beneficial", "Harmful"),
  x_trans = "log",
  xlim = c(0.45, 3),
  ticks_at = c(0.5, 1, 1.5),
  #xlab = "Adjusted hazard ratio",
  title = "A. Methadone (MET) trials",
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 16,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


#### bnx forestplot ##########################
full_ow_dis <- bind_rows(tr4_full_ow_dis,
                         tr5_full_ow_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_ow_dis$`OAT discontinuation` <- "" 
full_ow_dis$`HR(95% CI)` <- full_ow_dis$HR_display
full_ow_dis$`Dose adjustment after interruption` <- full_ow_dis$Comparison
full_ow_dis$`User design` <- full_ow_dis$Design

p_dis_bnx <- forest(
  full_ow_dis[, c("Dose adjustment after interruption", 
                  "User design", "HR(95% CI)", "OAT discontinuation")],
  est   = full_ow_dis$HR,
  lower = full_ow_dis$lower,
  upper = full_ow_dis$upper,
  ci_column = 4,
  ref_line = 1,
  x_trans = "log",
  xlim = c(0.45, 3),
  ticks_at = c(0.5, 1, 1.5),
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 16,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))






library(patchwork)
plot_dis <- plot_spacer() /  # top spacer row
  p_dis_met / 
  plot_spacer() /
  p_dis_bnx +
  plot_layout(
    heights = c(0.09, 1, 0.02, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))




name = "full_dis_ow_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3000, height=3900, res=300)
plot_dis
dev.off()





# Re-graph from June 24 ###########
## Read results files #############
folder = "H:/SAS-Work/Missed doses Paper v2/Rate_Folder_PLR/"
res1 <- readRDS(file = paste0(folder, "HRs_Tiral1_PLR.rds"))
res2 <- readRDS(file = paste0(folder, "HRs_Tiral2_PLR.rds"))
res3 <- readRDS(file = paste0(folder, "HRs_Tiral3_PLR.rds"))
res4 <- readRDS(file = paste0(folder, "HRs_Tiral4_PLR.rds"))
res5 <- readRDS(file = paste0(folder, "HRs_Tiral5_PLR.rds"))

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"
res1_lt30d_ow <- readRDS(file=paste0(folder, "res1_ow.rds"))
res2_lt30d_ow <- readRDS(file=paste0(folder, "res2_ow.rds"))
res3_lt30d_ow <- readRDS(file=paste0(folder, "res3_ow.rds"))
res4_lt30d_ow <- readRDS(file=paste0(folder, "res4_ow.rds"))
res5_lt30d_ow <- readRDS(file=paste0(folder, "res5_ow.rds"))

folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_30_lt90d/"
res1_30_lt90d_ow <- readRDS(file=paste0(folder, "res1_ow.rds"))
res2_30_lt90d_ow <- readRDS(file=paste0(folder, "res2_ow.rds"))
res3_30_lt90d_ow <- readRDS(file=paste0(folder, "res3_ow.rds"))
res4_30_lt90d_ow <- readRDS(file=paste0(folder, "res4_ow.rds"))
res5_30_lt90d_ow <- readRDS(file=paste0(folder, "res5_ow.rds"))


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_gt90d/"
res1_gt90d_ow <- readRDS(file=paste0(folder, "res1_ow.rds"))
res2_gt90d_ow <- readRDS(file=paste0(folder, "res2_ow.rds"))
res3_gt90d_ow <- readRDS(file=paste0(folder, "res3_ow.rds"))
res4_gt90d_ow <- readRDS(file=paste0(folder, "res4_ow.rds"))
res5_gt90d_ow <- readRDS(file=paste0(folder, "res5_ow.rds"))


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
res1_full_ow <- readRDS(file=paste0(folder, "res1_ow.rds"))
res2_full_ow <- readRDS(file=paste0(folder, "res2_ow.rds"))
res3_full_ow <- readRDS(file=paste0(folder, "res3_ow.rds"))
res4_full_ow <- readRDS(file=paste0(folder, "res4_ow.rds"))
res5_full_ow <- readRDS(file=paste0(folder, "res5_ow.rds"))


comparison_colors_v1 <- c(
  "Dose reduction\nvs\nNo dose change"              = "#000000",  # Blue
  "1-<50% reduction\nvs\nNo dose change"            = "#2CA02C",  # Green
  ">=50% reduction/initial dose*\nvs\nNo dose change"        = "#FF7F0E",  # Sky Blue
  "Dose increase\nvs\nNo dose change"               = "#D62728"  # Orange
  
)

plot_trial_results_wot_v1 <- function(df, title_lab, y_lab="Hazard ratios with 95% compatibility (‘confidence’) intervals") {
  ggplot(df, aes(x = Comparison, y = HR, color = Comparison)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.07, size = 0.7, position = position_dodge(width = 0.5)) +
    facet_wrap(~ Trial, ncol = 1, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    labs(
      title = title_lab,
      x = NULL,
      y = y_lab,
      color = " "
    ) +
    scale_color_manual(values = comparison_colors_v1) +
    # Optionally: scale_y_continuous(breaks = seq(y_min, y_max, by = y_by), limits = c(y_min, y_max)) +
    theme_minimal() +
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.x = element_line(linewidth = 1, color = "black"),
      axis.ticks.length  = unit(0.3, "cm"),
      axis.text.x = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 14),
      axis.title.y = element_text(size = 14),       # <-- added this line for y-axis label size
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 13)
    ) 
}


### Trial 1 #################################
res_df1_lt30d_iptw <- res1 %>% filter(SA_group == "5. Days since OAT initiation: <=30 days") %>% 
  mutate(res_type = "Within 30 days of OAT initiation",
         Comparison = recode(Comparison,
                             "Dose reduction\nvs\nNo dose change"="Dose reduction vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 1: MET interruption of 1-3 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)

res1 %>% count(SA_group)

res_df1_30_lt90d_iptw <- res1 %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days") %>% 
  mutate(res_type = "Within 31-90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "Dose reduction\nvs\nNo dose change"="Dose reduction vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 1: MET interruption of 1-3 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)

res_df1_gt90d_iptw <- res1 %>% filter(SA_group == "7. Days since OAT initiation: >90 days") %>% 
  mutate(res_type = "After 90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "Dose reduction\nvs\nNo dose change"="Dose reduction vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 1: MET interruption of 1-3 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


res_df1_full_iptw <- res1 %>% filter(SA_group == "1. Primary analysis") %>% 
  mutate(res_type = "Primary analysis",
         Comparison = recode(Comparison,
                             "Dose reduction\nvs\nNo dose change"="Dose reduction vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 1: MET interruption of 1-3 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)



### Trial 2 #################################
res_df2_lt30d_iptw <- res2 %>% filter(SA_group == "5. Days since OAT initiation: <=30 days") %>% 
  mutate(res_type = "Within 30 days of OAT initiation",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 2: MET interruption of 4 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


res_df2_30_lt90d_iptw <- res2 %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days") %>% 
  mutate(res_type = "Within 31-90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 2: MET interruption of 4 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)

res_df2_gt90d_iptw <- res2 %>% filter(SA_group == "7. Days since OAT initiation: >90 days") %>% 
  mutate(res_type = "After 90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 2: MET interruption of 4 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


res_df2_full_iptw <- res2 %>% filter(SA_group == "1. Primary analysis") %>% 
  mutate(res_type = "Primary analysis",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 2: MET interruption of 4 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


### Trial 3 #################################
res_df3_lt30d_iptw <- res3 %>% filter(SA_group == "5. Days since OAT initiation: <=30 days") %>% 
  mutate(res_type = "Within 30 days of OAT initiation",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 3: MET interruption of 5-14 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


res_df3_30_lt90d_iptw <- res3 %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days") %>% 
  mutate(res_type = "Within 31-90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 3: MET interruption of 5-14 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)

res_df3_gt90d_iptw <- res3 %>% filter(SA_group == "7. Days since OAT initiation: >90 days") %>% 
  mutate(res_type = "After 90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 3: MET interruption of 5-14 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


res_df3_full_iptw <- res3 %>% filter(SA_group == "1. Primary analysis") %>% 
  mutate(res_type = "Primary analysis",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 3: MET interruption of 5-14 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)



### Trial 4 #################################
res_df4_lt30d_iptw <- res4 %>% filter(SA_group == "5. Days since OAT initiation: <=30 days") %>% 
  mutate(res_type = "Within 30 days of OAT initiation",
         Comparison = recode(Comparison,
                             "Dose reduction\nvs\nNo dose change"="Dose reduction vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 4: BNX interruption of 1-5 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)

res4 %>% count(SA_group)

res_df4_30_lt90d_iptw <- res4 %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days") %>% 
  mutate(res_type = "Within 31-90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "Dose reduction\nvs\nNo dose change"="Dose reduction vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 4: BNX interruption of 1-5 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)

res_df4_gt90d_iptw <- res4 %>% filter(SA_group == "7. Days since OAT initiation: >90 days") %>% 
  mutate(res_type = "After 90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "Dose reduction\nvs\nNo dose change"="Dose reduction vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 4: BNX interruption of 1-5 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


res_df4_full_iptw <- res4 %>% filter(SA_group == "1. Primary analysis") %>% 
  mutate(res_type = "Primary analysis",
         Comparison = recode(Comparison,
                             "Dose reduction\nvs\nNo dose change"="Dose reduction vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 4: BNX interruption of 1-5 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


### Trial 5 #################################
res_df5_lt30d_iptw <- res5 %>% filter(SA_group == "5. Days since OAT initiation: <=30 days") %>% 
  mutate(res_type = "Within 30 days of OAT initiation",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 5: BNX interruption of 6-14 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


res_df5_30_lt90d_iptw <- res5 %>% filter(SA_group == "6. Days since OAT initiation: 31-90 days") %>% 
  mutate(res_type = "Within 31-90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 5: BNX interruption of 6-14 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)

res_df5_gt90d_iptw <- res5 %>% filter(SA_group == "7. Days since OAT initiation: >90 days") %>% 
  mutate(res_type = "After 90 days of OAT initiation",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 5: BNX interruption of 6-14 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


res_df5_full_iptw <- res5 %>% filter(SA_group == "1. Primary analysis") %>% 
  mutate(res_type = "Primary analysis",
         Comparison = recode(Comparison,
                             "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                             ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                             "Dose increase\nvs\nNo dose change" = "Dose increase vs no change"),
         Comparison = factor(Comparison, levels = c("Dose reduction to <50% vs no change",
                                                    "Dose reduction to >=50%/initial dose* vs no change",
                                                    "Dose increase vs no change")),
         Design = recode(Design,
                         "Prevalent new* users" = "All prevalent users"),
         Trial = "Trial 5: BNX interruption of 6-14 days") %>% 
  select(Outcome, res_type, Trial, Comparison, Design, HR, lower, upper)


res_df <- bind_rows(res_df1_full_iptw, res_df2_full_iptw, res_df3_full_iptw, res_df4_full_iptw, res_df5_full_iptw,
                    res_df1_lt30d_iptw, res_df2_lt30d_iptw, res_df3_lt30d_iptw, res_df4_lt30d_iptw, res_df5_lt30d_iptw,
                    res_df1_30_lt90d_iptw, res_df2_30_lt90d_iptw, res_df3_30_lt90d_iptw, res_df4_30_lt90d_iptw, res_df5_30_lt90d_iptw,
                    res_df1_gt90d_iptw, res_df2_gt90d_iptw, res_df3_gt90d_iptw, res_df4_gt90d_iptw, res_df5_gt90d_iptw)
names(res_df)

res_df <- res_df %>%
  mutate(
    Design = factor(Design,
                    levels = c("Incident users",
                               "Prevalent new users",
                               "All prevalent users")
    )) %>% arrange(res_type, Outcome, Trial, Comparison, Design)


saveRDS(res_df, file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw.rds")





res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = factor(Comparison, levels = c("Dose reduction vs no change", 
                                               "Dose reduction to <50% vs no change", 
                                               "Dose reduction to >=50%/initial dose* vs no change",
                                               "Dose increase vs no change"))
  ) %>% arrange(res_type, Outcome, Trial, Comparison, Design)



### od: primary, iptw ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"

#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 1: MET interruption of 1-2 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 2: MET interruption of 3-4 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 4: BNX interruption of 1-5 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
full_iptw_od <- bind_rows(tr1_full_iptw_od,
                          tr2_full_iptw_od,
                          tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_od$`Overdose-related acute care visits` <- "" 
full_iptw_od$`HR(95% CI)` <- full_iptw_od$HR_display
full_iptw_od$`Dose adjustment after interruption` <- full_iptw_od$Comparison
full_iptw_od$`User design` <- full_iptw_od$Design


library(ggplot2)
library(forestploter)
library(tidyverse)

p_od_met <- forest(
  full_iptw_od[, c("Dose adjustment after interruption", 
                   "User design", "HR(95% CI)", "Overdose-related acute care visits")],
  est   = full_iptw_od$HR,
  lower = full_iptw_od$lower,
  upper = full_iptw_od$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  arrow_lab = c("Beneficial", "Harmful"),
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5, 3),
  #xlab = "Adjusted hazard ratio",
  title = "A. Methadone (MET) trials",
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


#### bnx forestplot ##########################
full_iptw_od <- bind_rows(tr4_full_iptw_od,
                          tr5_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_od$`Overdose-related acute care visits` <- "" 
full_iptw_od$`HR(95% CI)` <- full_iptw_od$HR_display
full_iptw_od$`Dose adjustment after interruption` <- full_iptw_od$Comparison
full_iptw_od$`User design` <- full_iptw_od$Design

p_od_bnx <- forest(
  full_iptw_od[, c("Dose adjustment after interruption", 
                   "User design", "HR(95% CI)", "Overdose-related acute care visits")],
  est   = full_iptw_od$HR,
  lower = full_iptw_od$lower,
  upper = full_iptw_od$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5, 3), 
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5\n** NA for few events", 
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))





library(patchwork)
plot_od <- plot_spacer() /  # top spacer row
  p_od_met / 
  plot_spacer() /
  p_od_bnx +
  plot_layout(
    heights = c(0.09, 1, 0.02, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))




name = "full_od_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3000, height=3900, res=300)
plot_od
dev.off()


### dis: primary, iptw ######

#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr1_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 1: MET interruption of 1-2 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 4: BNX interruption of 1-5 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                           tr2_full_iptw_dis,
                           tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_dis$`OAT discontinuation` <- "" 
full_iptw_dis$`HR(95% CI)` <- full_iptw_dis$HR_display
full_iptw_dis$`Dose adjustment after interruption` <- full_iptw_dis$Comparison
full_iptw_dis$`User design` <- full_iptw_dis$Design


library(ggplot2)
library(forestploter)
library(tidyverse)
p_dis_met <- forest(
  full_iptw_dis[, c("Dose adjustment after interruption", 
                    "User design", "HR(95% CI)", "OAT discontinuation")],
  est   = full_iptw_dis$HR,
  lower = full_iptw_dis$lower,
  upper = full_iptw_dis$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  arrow_lab = c("Beneficial", "Harmful"),
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5, 3),
  #xlab = "Adjusted hazard ratio",
  title = "A. Methadone (MET) trials",
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


#### bnx forestplot ##########################
full_iptw_dis <- bind_rows(tr4_full_iptw_dis,
                           tr5_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_dis$`OAT discontinuation` <- "" 
full_iptw_dis$`HR(95% CI)` <- full_iptw_dis$HR_display
full_iptw_dis$`Dose adjustment after interruption` <- full_iptw_dis$Comparison
full_iptw_dis$`User design` <- full_iptw_dis$Design

p_dis_bnx <- forest(
  full_iptw_dis[, c("Dose adjustment after interruption", 
                    "User design", "HR(95% CI)", "OAT discontinuation")],
  est   = full_iptw_dis$HR,
  lower = full_iptw_dis$lower,
  upper = full_iptw_dis$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5, 3),
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))






library(patchwork)
plot_dis <- plot_spacer() /  # top spacer row
  p_dis_met / 
  plot_spacer() /
  p_dis_bnx +
  plot_layout(
    heights = c(0.09, 1, 0.02, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))




name = "full_dis_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3000, height=3900, res=300)
plot_dis
dev.off()


### od & dis: primary, iptw ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"

#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 1: MET interruption of 1-2 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 2: MET interruption of 3-4 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 4: BNX interruption of 1-5 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                              tr2_full_iptw_od,
                              tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


met_full_iptw_od$`OD-related acute care visits` <- "" 
met_full_iptw_od$`HR(95% CI)` <- met_full_iptw_od$HR_display
met_full_iptw_od$`Dose adjustment after interruption` <- met_full_iptw_od$Comparison
met_full_iptw_od$`User design` <- met_full_iptw_od$Design

# 
# library(ggplot2)
# library(forestploter)
# library(tidyverse)
# 
# p_od_met <- forest(
#   full_iptw_od[, c("Dose adjustment after interruption", 
#                    "User design", "HR(95% CI)", "Overdose-related acute care visits")],
#   est   = full_iptw_od$HR,
#   lower = full_iptw_od$lower,
#   upper = full_iptw_od$upper,
#   ci_column = 4,
#   ref_line = 1,
#   sizes = .6,
#   arrow_lab = c("Beneficial", "Harmful"),
#   x_trans = "log",
#   xlim = c(0.45, 3),
#   ticks_at = c(0.5, 1, 1.5),
#   #xlab = "Adjusted hazard ratio",
#   title = "A. Methadone (MET) trials",
#   nudge_y = 0.2,
#   theme = 
#     forest_theme(
#       base_size = 11,
#       ci_pch = 15,
#       ci_lwd = 1.7,
#       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
#       footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
#                          fontface="plain", col = "black"),
#       arrow_length = 0.06,
#       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
#       xlab_adjust = c("refline", "center")
#     ))
# 

#### bnx forestplot ##########################
bnx_full_iptw_od <- bind_rows(tr4_full_iptw_od,
                              tr5_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


bnx_full_iptw_od$`OD-related acute care visits` <- "" 
bnx_full_iptw_od$`HR(95% CI)` <- bnx_full_iptw_od$HR_display
bnx_full_iptw_od$`Dose adjustment after interruption` <- bnx_full_iptw_od$Comparison
bnx_full_iptw_od$`User design` <- bnx_full_iptw_od$Design

# p_od_bnx <- forest(
#   full_iptw_od[, c("Dose adjustment after interruption", 
#                    "User design", "HR(95% CI)", "Overdose-related acute care visits")],
#   est   = full_iptw_od$HR,
#   lower = full_iptw_od$lower,
#   upper = full_iptw_od$upper,
#   ci_column = 4,
#   ref_line = 1,
#   sizes = .6,
#   x_trans = "log",
#   xlim = c(0.45, 3),
#   ticks_at = c(0.5, 1, 1.5),
#   xlab = "Adjusted hazard ratio",
#   arrow_lab = c("Beneficial", "Harmful"),
#   title = "B. Buprenorphine/naloxone (BNX) trials",
#   footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5\n** NA for few events", 
#   nudge_y = 0.2,
#   theme = 
#     forest_theme(
#       base_size = 11,
#       ci_pch = 15,
#       ci_lwd = 1.7,
#       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
#       footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
#                          fontface="plain", col = "black"),
#       arrow_length = 0.06,
#       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
#       xlab_adjust = c("refline", "center")
#     ))
# 
# 
# 
# 
# 
# library(patchwork)
# plot_od <- plot_spacer() /  # top spacer row
#   p_od_met / 
#   plot_spacer() /
#   p_od_bnx +
#   plot_layout(
#     heights = c(0.09, 1, 0.02, 1),
#     guides = "collect"
#   ) &
#   theme(
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     legend.justification = "left",
#     legend.box.just = "left"
#   ) &
#   guides(color = guide_legend(nrow = 1))
# 
# 
# 
# 
# name = "full_od_iptw_Nov12"
# figure=paste0(folder, name, ".png")
# ppi = 300
# png(figure, width=3000, height=3900, res=300)
# plot_od
# dev.off()


#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr1_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 1: MET interruption of 1-2 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 4: BNX interruption of 1-5 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
met_full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                               tr2_full_iptw_dis,
                               tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


met_full_iptw_dis$`OAT discontinuation` <- "" 
met_full_iptw_dis$`HR(95% CI)` <- met_full_iptw_dis$HR_display
met_full_iptw_dis$`Dose adjustment after interruption` <- met_full_iptw_dis$Comparison
met_full_iptw_dis$`User design` <- met_full_iptw_dis$Design

# 
# library(ggplot2)
# library(forestploter)
# library(tidyverse)
# p_dis_met <- forest(
#   full_iptw_dis[, c("HR(95% CI)", "OAT discontinuation")],
#   est   = full_iptw_dis$HR,
#   lower = full_iptw_dis$lower,
#   upper = full_iptw_dis$upper,
#   ci_column = 2,
#   ref_line = 1,
#   sizes = .6,
#   arrow_lab = c("Beneficial", "Harmful"),
#   x_trans = "log",
#   xlim = c(0.45, 3),
#   ticks_at = c(0.5, 1, 1.5),
#   #xlab = "Adjusted hazard ratio",
#   title = NULL,
#   nudge_y = 0.2,
#   theme = 
#     forest_theme(
#       base_size = 11,
#       ci_pch = 15,
#       ci_lwd = 1.7,
#       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
#       footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
#                          fontface="plain", col = "black"),
#       arrow_length = 0.06,
#       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
#       xlab_adjust = c("refline", "center")
#     ))


#### bnx forestplot ##########################
bnx_full_iptw_dis <- bind_rows(tr4_full_iptw_dis,
                               tr5_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


bnx_full_iptw_dis$`OAT discontinuation` <- "" 
bnx_full_iptw_dis$`HR(95% CI)` <- bnx_full_iptw_dis$HR_display
bnx_full_iptw_dis$`Dose adjustment after interruption` <- bnx_full_iptw_dis$Comparison
bnx_full_iptw_dis$`User design` <- bnx_full_iptw_dis$Design



# p_dis_bnx <- forest(
#   full_iptw_dis[, c("HR(95% CI)", "OAT discontinuation")],
#   est   = full_iptw_dis$HR,
#   lower = full_iptw_dis$lower,
#   upper = full_iptw_dis$upper,
#   ci_column = 2,
#   ref_line = 1,
#   sizes = .6,
#   x_trans = "log",
#   xlim = c(0.45, 3),
#   ticks_at = c(0.5, 1, 1.5),
#   xlab = "Adjusted hazard ratio",
#   arrow_lab = c("Beneficial", "Harmful"),
#   title = NULL,
#   footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
#   nudge_y = 0.2,
#   theme = 
#     forest_theme(
#       base_size = 11,
#       ci_pch = 15,
#       ci_lwd = 1.7,
#       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
#       footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
#                          fontface="plain", col = "black"),
#       arrow_length = 0.06,
#       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
#       xlab_adjust = c("refline", "center")
#     ))

#### creating the graph ######

met_full_iptw_dis_v1 <- met_full_iptw_dis %>% select(`HR(95% CI)`, `OAT discontinuation`)
met_full_iptw <- bind_cols(met_full_iptw_od, met_full_iptw_dis_v1)

names(met_full_iptw)[c(8, 11)] <- c("HR(95% CI)", "    HR(95% CI)")
met_full_iptw$`    HR(95% CI)` <- paste0("    ", met_full_iptw$`    HR(95% CI)`)

p_od_dis_met <- forest(
  met_full_iptw[, c(9,10,8,7,11,12)],
  est   = list(met_full_iptw_od$HR, met_full_iptw_dis$HR),
  lower = list(met_full_iptw_od$lower, met_full_iptw_dis$lower),
  upper = list(met_full_iptw_od$upper, met_full_iptw_dis$upper),
  ci_column = c(4, 6),
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5, 3),
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "A. Methadone (MET) trials",
  footnote = NULL,
  nudge_y = 0.2,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))



bnx_full_iptw_dis_v1 <- bnx_full_iptw_dis %>% select(`HR(95% CI)`, `OAT discontinuation`)
bnx_full_iptw <- bind_cols(bnx_full_iptw_od, bnx_full_iptw_dis_v1)
names(bnx_full_iptw)
names(bnx_full_iptw)[c(8, 11)] <- c("HR(95% CI)", "    HR(95% CI)")
bnx_full_iptw$`    HR(95% CI)` <- paste0("    ", bnx_full_iptw$`    HR(95% CI)`)

p_od_dis_bnx <- forest(
  bnx_full_iptw[, c(9,10,8,7,11,12)],
  est   = list(bnx_full_iptw_od$HR, bnx_full_iptw_dis$HR),
  lower = list(bnx_full_iptw_od$lower, bnx_full_iptw_dis$lower),
  upper = list(bnx_full_iptw_od$upper, bnx_full_iptw_dis$upper),
  ci_column = c(4, 6),
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5, 3),
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5\n** NA for few events; OD. overdose",
  nudge_y = 0.2,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      footnote_parse = FALSE,
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))




library(patchwork)
plot_od_dis <- plot_spacer() /  # top spacer row
  p_od_dis_met / 
  plot_spacer() /
  p_od_dis_bnx/plot_spacer() +
  plot_layout(
    heights = c(0.12, 1, 0.08, 1, 0.01),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))



name = "full_od_dis_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=4700, height=4600, res=300)
plot_od_dis
dev.off()



### od: Less than 30 days, iptw ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"

res_df2 %>% count(res_type)

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

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


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

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
full_iptw_od <- bind_rows(tr1_full_iptw_od,
                          tr2_full_iptw_od,
                          tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_od$`Overdose-related acute care visits` <- "" 
full_iptw_od$`HR(95% CI)` <- full_iptw_od$HR_display
full_iptw_od$`Dose adjustment after interruption` <- full_iptw_od$Comparison
full_iptw_od$`User design` <- full_iptw_od$Design


library(ggplot2)
library(forestploter)
library(tidyverse)

p_od_met <- forest(
  full_iptw_od[, c("Dose adjustment after interruption", 
                   "User design", "HR(95% CI)", "Overdose-related acute care visits")],
  est   = full_iptw_od$HR,
  lower = full_iptw_od$lower,
  upper = full_iptw_od$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  arrow_lab = c("Beneficial", "Harmful"),
  x_trans = "log",
  xlim = c(0.45, 3),
  ticks_at = c(0.5, 1, 1.5),
  #xlab = "Adjusted hazard ratio",
  title = "A. Methadone (MET) trials",
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


#### bnx forestplot ##########################
full_iptw_od <- bind_rows(tr4_full_iptw_od,
                          tr5_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_od$`Overdose-related acute care visits` <- "" 
full_iptw_od$`HR(95% CI)` <- full_iptw_od$HR_display
full_iptw_od$`Dose adjustment after interruption` <- full_iptw_od$Comparison
full_iptw_od$`User design` <- full_iptw_od$Design

p_od_bnx <- forest(
  full_iptw_od[, c("Dose adjustment after interruption", 
                   "User design", "HR(95% CI)", "Overdose-related acute care visits")],
  est   = full_iptw_od$HR,
  lower = full_iptw_od$lower,
  upper = full_iptw_od$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.45, 3),
  ticks_at = c(0.5, 1, 1.5),
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5\n** NA for few events", 
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))





library(patchwork)
lt30d_plot_od <- plot_spacer() /  # top spacer row
  p_od_met / 
  plot_spacer() /
  p_od_bnx +
  plot_layout(
    heights = c(0.09, 1, 0.02, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))




name = "lt30d_od_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3000, height=3900, res=300)
lt30d_plot_od
dev.off()


### dis: Less than 30 days, iptw ######

#### Trial 1: od: iptw #######
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

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 4: od: iptw #######
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

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                           tr2_full_iptw_dis,
                           tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_dis$`OAT discontinuation` <- "" 
full_iptw_dis$`HR(95% CI)` <- full_iptw_dis$HR_display
full_iptw_dis$`Dose adjustment after interruption` <- full_iptw_dis$Comparison
full_iptw_dis$`User design` <- full_iptw_dis$Design


library(ggplot2)
library(forestploter)
library(tidyverse)
p_dis_met <- forest(
  full_iptw_dis[, c("Dose adjustment after interruption", 
                    "User design", "HR(95% CI)", "OAT discontinuation")],
  est   = full_iptw_dis$HR,
  lower = full_iptw_dis$lower,
  upper = full_iptw_dis$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  arrow_lab = c("Beneficial", "Harmful"),
  x_trans = "log",
  xlim = c(0.45, 3),
  ticks_at = c(0.5, 1, 1.5),
  #xlab = "Adjusted hazard ratio",
  title = "A. Methadone (MET) trials",
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


#### bnx forestplot ##########################
full_iptw_dis <- bind_rows(tr4_full_iptw_dis,
                           tr5_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_dis$`OAT discontinuation` <- "" 
full_iptw_dis$`HR(95% CI)` <- full_iptw_dis$HR_display
full_iptw_dis$`Dose adjustment after interruption` <- full_iptw_dis$Comparison
full_iptw_dis$`User design` <- full_iptw_dis$Design

p_dis_bnx <- forest(
  full_iptw_dis[, c("Dose adjustment after interruption", 
                    "User design", "HR(95% CI)", "OAT discontinuation")],
  est   = full_iptw_dis$HR,
  lower = full_iptw_dis$lower,
  upper = full_iptw_dis$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.45, 3),
  ticks_at = c(0.5, 1, 1.5),
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))






library(patchwork)
lt30d_plot_dis <- plot_spacer() /  # top spacer row
  p_dis_met / 
  plot_spacer() /
  p_dis_bnx +
  plot_layout(
    heights = c(0.09, 1, 0.02, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))




name = "lt30d_dis_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3000, height=3900, res=300)
lt30d_plot_dis
dev.off()



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

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()


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

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),    
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()




#### met forestplot ##########################
full_iptw_od <- bind_rows(tr1_full_iptw_od,
                          tr2_full_iptw_od,
                          tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_od$`Overdose-related acute care visits` <- "" 
full_iptw_od$`HR(95% CI)` <- full_iptw_od$HR_display
full_iptw_od$`Dose adjustment after interruption` <- full_iptw_od$Comparison
full_iptw_od$`User design` <- full_iptw_od$Design


library(ggplot2)
library(forestploter)
library(tidyverse)

p_od_met <- forest(
  full_iptw_od[, c("Dose adjustment after interruption", 
                   "User design", "HR(95% CI)", "Overdose-related acute care visits")],
  est   = full_iptw_od$HR,
  lower = full_iptw_od$lower,
  upper = full_iptw_od$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  arrow_lab = c("Beneficial", "Harmful"),
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5),
  #xlab = "Adjusted hazard ratio",
  title = "A. Methadone (MET) trials",
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


#### bnx forestplot ##########################
full_iptw_od <- bind_rows(tr4_full_iptw_od,
                          tr5_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_od$`Overdose-related acute care visits` <- "" 
full_iptw_od$`HR(95% CI)` <- full_iptw_od$HR_display
full_iptw_od$`Dose adjustment after interruption` <- full_iptw_od$Comparison
full_iptw_od$`User design` <- full_iptw_od$Design

p_od_bnx <- forest(
  full_iptw_od[, c("Dose adjustment after interruption", 
                   "User design", "HR(95% CI)", "Overdose-related acute care visits")],
  est   = full_iptw_od$HR,
  lower = full_iptw_od$lower,
  upper = full_iptw_od$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5),
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5\n** NA for few events", 
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))





library(patchwork)
f31to90d_plot_od <- plot_spacer() /  # top spacer row
  p_od_met / 
  plot_spacer() /
  p_od_bnx +
  plot_layout(
    heights = c(0.09, 1, 0.02, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))




name = "f31to90d_od_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3000, height=3900, res=300)
f31to90d_plot_od
dev.off()


### dis: 31 to 90 days, iptw ######

#### Trial 1: od: iptw #######
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

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()


#### Trial 4: od: iptw #######
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

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()




#### met forestplot ##########################
full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                           tr2_full_iptw_dis,
                           tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_dis$`OAT discontinuation` <- "" 
full_iptw_dis$`HR(95% CI)` <- full_iptw_dis$HR_display
full_iptw_dis$`Dose adjustment after interruption` <- full_iptw_dis$Comparison
full_iptw_dis$`User design` <- full_iptw_dis$Design


library(ggplot2)
library(forestploter)
library(tidyverse)
p_dis_met <- forest(
  full_iptw_dis[, c("Dose adjustment after interruption", 
                    "User design", "HR(95% CI)", "OAT discontinuation")],
  est   = full_iptw_dis$HR,
  lower = full_iptw_dis$lower,
  upper = full_iptw_dis$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  arrow_lab = c("Beneficial", "Harmful"),
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5),
  #xlab = "Adjusted hazard ratio",
  title = "A. Methadone (MET) trials",
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


#### bnx forestplot ##########################
full_iptw_dis <- bind_rows(tr4_full_iptw_dis,
                           tr5_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


full_iptw_dis$`OAT discontinuation` <- "" 
full_iptw_dis$`HR(95% CI)` <- full_iptw_dis$HR_display
full_iptw_dis$`Dose adjustment after interruption` <- full_iptw_dis$Comparison
full_iptw_dis$`User design` <- full_iptw_dis$Design

p_dis_bnx <- forest(
  full_iptw_dis[, c("Dose adjustment after interruption", 
                    "User design", "HR(95% CI)", "OAT discontinuation")],
  est   = full_iptw_dis$HR,
  lower = full_iptw_dis$lower,
  upper = full_iptw_dis$upper,
  ci_column = 4,
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5),
  xlab = "Adjusted hazard ratio",
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
  nudge_y = 0.2,
  theme = 
    forest_theme(
      base_size = 11,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
                         fontface="plain", col = "black"),
      arrow_length = 0.06,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))






library(patchwork)
f31to90d_plot_dis <- plot_spacer() /  # top spacer row
  p_dis_met / 
  plot_spacer() /
  p_dis_bnx +
  plot_layout(
    heights = c(0.09, 1, 0.02, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))




name = "f31to90d_dis_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3000, height=3900, res=300)
f31to90d_plot_dis
dev.off()

# Sens results: forestplot #####################
res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change")
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

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


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

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
lt30d_met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                                    tr2_full_iptw_od,
                                    tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


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
  )


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

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: dis: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

#### Trial 3: dis: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


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

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: dis: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
lt30d_met_full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                                     tr2_full_iptw_dis,
                                     tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


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
  )


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

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()


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

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),    
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()




#### met forestplot ##########################
f31to90_met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                                      tr2_full_iptw_od,
                                      tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


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
  )


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

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: dis: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

#### Trial 3: dis: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()


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

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: dis: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()




#### met forestplot ##########################
f31to90_met_full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                                       tr2_full_iptw_dis,
                                       tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


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
  )


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

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()


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

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),    
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()




#### met forestplot ##########################
af90d_met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                                    tr2_full_iptw_od,
                                    tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


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
  )


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

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: dis: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

#### Trial 3: dis: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()


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

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: dis: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)) %>%
  ungroup()




#### met forestplot ##########################
af90d_met_full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                                     tr2_full_iptw_dis,
                                     tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


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
  )


af90d_bnx_full_iptw_dis$`OAT discontinuation` <- "" 
af90d_bnx_full_iptw_dis$`HR(95% CI)` <- af90d_bnx_full_iptw_dis$HR_display
af90d_bnx_full_iptw_dis$`Dose adjustment after interruption` <- af90d_bnx_full_iptw_dis$Comparison
af90d_bnx_full_iptw_dis$`User design` <- af90d_bnx_full_iptw_dis$Design



#### creating forestplot ############


tm <- forest_theme(ci_lty = "solid",
                   ci_pch = c(15, 15, 15),
                   base_size = 14,
                   #ci_pch = 15,
                   ci_lwd = 1.8,
                   ci_col = c("black", "springgreen3", "hotpink"),
                   legend_name = NULL,
                   legend_value = c(" Within 30 day of OAT initiation ", 
                                    " 31 to 90 day of OAT initiation ", 
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


txt <- "<sup>*</sup> <=40mg methadone<br><sup>a</sup> NA for few events<br><sup>b</sup> adjusted hazard ratio with 95% compatibility intervals"



met_full_iptw_od <- lt30d_met_full_iptw_od %>% select(`Dose adjustment after interruption`,
                                                      `User design`)

met_full_iptw_od <- met_full_iptw_od %>%  mutate(aHR=ifelse(`User design` %in% c(""), "", 
                                                            paste0(lt30d_met_full_iptw_od$`HR(95% CI)`, "\n",
                                                                   f31to90_met_full_iptw_od$`HR(95% CI)`, "\n",
                                                                   af90d_met_full_iptw_od$`HR(95% CI)`)),
                                                 `User design` = paste0(" ", `User design`)
                                                 
)
met_full_iptw_od$`aHR(95% CI)^b` <-  met_full_iptw_od$aHR
met_full_iptw_od$`OD-related acute care visits` <- ""


p_od_met <- forest(
  met_full_iptw_od[, c(1, 2, 4, 5)],
  est   = list(lt30d_met_full_iptw_od$HR, 
               f31to90_met_full_iptw_od$HR, 
               af90d_met_full_iptw_od$HR),
  lower = list(lt30d_met_full_iptw_od$lower, 
               f31to90_met_full_iptw_od$lower, 
               af90d_met_full_iptw_od$lower),
  upper = list(lt30d_met_full_iptw_od$upper, 
               f31to90_met_full_iptw_od$upper, 
               af90d_met_full_iptw_od$upper),
  ci_column = c(4),
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 6),
  ticks_at = c(0.5, 1, 1.5, 3),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "A. Methadone (MET) trials",
  footnote = NULL,
  nudge_y = 0.3,
  theme = tm)
p_od_met <- add_grob(p_od_met, 
                     row = 28, 
                     col = 1:2,
                     order = "background",
                     gb_fn = gridtext::richtext_grob,
                     text = txt,
                     gp = gpar(fontsize = 11),
                     hjust = 0, vjust = 1, halign = 0, valign = 1,
                     x = unit(0, "npc"), y = unit(1, "npc"))



folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"
name = "str_full_od_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3800, height=6700, res=300)
p_od_met
dev.off()


#length(af90d_met_full_iptw_dis$HR)

met_full_iptw_dis <- lt30d_met_full_iptw_dis %>% select(`Dose adjustment after interruption`,
                                                        `User design`)

met_full_iptw_dis <- met_full_iptw_dis %>%  mutate(aHR=ifelse(`User design` %in% c("", "Incident users**"), "", 
                                                              paste0(lt30d_met_full_iptw_dis$`HR(95% CI)`, "\n",
                                                                     f31to90_met_full_iptw_dis$`HR(95% CI)`, "\n",
                                                                     af90d_met_full_iptw_dis$`HR(95% CI)`)),
                                                   `User design` = paste0(" ", `User design`)
                                                   
)
met_full_iptw_dis$`aHR(95% CI)^b` <-  met_full_iptw_od$aHR

met_full_iptw_dis$`OAT discontinuation` <- ""


p_dis_met <- forest(
  met_full_iptw_dis[, c(1, 2, 4, 5)],
  est   = list(lt30d_met_full_iptw_dis$HR, 
               f31to90_met_full_iptw_dis$HR, 
               af90d_met_full_iptw_dis$HR),
  lower = list(lt30d_met_full_iptw_dis$lower, 
               f31to90_met_full_iptw_dis$lower, 
               af90d_met_full_iptw_dis$lower),
  upper = list(lt30d_met_full_iptw_dis$upper, 
               f31to90_met_full_iptw_dis$upper, 
               af90d_met_full_iptw_dis$upper),
  ci_column = c(4),
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 6),
  ticks_at = c(0.5, 1, 1.5, 3),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "A. Methadone (MET) trials",
  footnote = NULL,
  nudge_y = 0.3,
  theme = tm)

p_dis_met <- add_grob(p_dis_met, 
                      row = 28, 
                      col = 1:2,
                      order = "background",
                      gb_fn = gridtext::richtext_grob,
                      text = txt,
                      gp = gpar(fontsize = 11),
                      hjust = 0, vjust = 1, halign = 0, valign = 1,
                      x = unit(0, "npc"), y = unit(1, "npc"))




folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/ret_lt30d/"
name = "str_full_dis_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3800, height=6700, res=300)
p_dis_met
dev.off()





# Figure 3 #################

res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change")
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
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 30 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Within 31-90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "After 90 days of OAT initiation",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
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


tm_met <- forest_theme(ci_lty = "solid",
                       ci_pch = c(15, 15, 15),
                       base_size = 14,
                       #ci_pch = 15,
                       ci_lwd = 1.8,
                       ci_col = c("black", "springgreen3", "hotpink"),
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
                   ci_pch = c(15, 15, 15),
                   base_size = 14,
                   #ci_pch = 15,
                   ci_lwd = 1.8,
                   ci_col = c("black", "springgreen3", "hotpink"),
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



met_full_iptw %>% count(`Dose adjustment after interruption`)

met_full_iptw <- lt30d_met_full_iptw_od %>% select(`Dose adjustment after interruption`,
                                                   `User design`) %>% 
  mutate(`Dose adjustment after interruption`=recode(`Dose adjustment after interruption`,
                                                     "   Dose reduction to >=50%/initial dose* vs no change"=
                                                       "   Dose reduction to >=50%/initial dose* vs\n   no change"))

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
  title = "A. Methadone (MET) trials among all prevalent users",
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
                                                     "   Dose reduction to >=50%/initial dose* vs no change"=
                                                       "   Dose reduction to >=50%/initial dose* vs\n   no change"))

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
  title = "A. Buprenorphine/naloxone (BNX) trials among all prevalent users",
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
name = "Str_all_od_dis_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=3700, height=4700, res=300)
plot_od_dis
dev.off()



## Figure 2: iptw ################

res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change")
  ) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design)


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"

#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 1: MET interruption of 1-2 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users^a"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users^a"))

tr2_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 2: MET interruption of 3-4 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users^a"))



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users^a", NA, HR),
    lower = ifelse(Design == "Incident users^a", NA, lower),
    upper = ifelse(Design == "Incident users^a", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 4: BNX interruption of 1-5 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, tr5_full_iptw_od_redgt50, tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                              tr2_full_iptw_od,
                              tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


met_full_iptw_od$`Overdose-related acute care visits` <- "" 
met_full_iptw_od$`HR(95% CI)**` <- met_full_iptw_od$HR_display
met_full_iptw_od$`Dose adjustment after interruption` <- met_full_iptw_od$Comparison
met_full_iptw_od$`User design` <- met_full_iptw_od$Design

# 
# library(ggplot2)
# library(forestploter)
# library(tidyverse)
# 
# p_od_met <- forest(
#   full_iptw_od[, c("Dose adjustment after interruption", 
#                    "User design", "HR(95% CI)", "Overdose-related acute care visits")],
#   est   = full_iptw_od$HR,
#   lower = full_iptw_od$lower,
#   upper = full_iptw_od$upper,
#   ci_column = 4,
#   ref_line = 1,
#   sizes = .6,
#   arrow_lab = c("Beneficial", "Harmful"),
#   x_trans = "log",
#   xlim = c(0.45, 3),
#   ticks_at = c(0.5, 1, 1.5),
#   #xlab = "Adjusted hazard ratio",
#   title = "A. Methadone (MET) trials",
#   nudge_y = 0.2,
#   theme = 
#     forest_theme(
#       base_size = 11,
#       ci_pch = 15,
#       ci_lwd = 1.7,
#       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
#       footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
#                          fontface="plain", col = "black"),
#       arrow_length = 0.06,
#       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
#       xlab_adjust = c("refline", "center")
#     ))
# 

#### bnx forestplot ##########################
bnx_full_iptw_od <- bind_rows(tr4_full_iptw_od,
                              tr5_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


bnx_full_iptw_od$`Overdose-related acute care visits` <- "" 
bnx_full_iptw_od$`HR(95% CI)` <- bnx_full_iptw_od$HR_display
bnx_full_iptw_od$`Dose adjustment after interruption` <- bnx_full_iptw_od$Comparison
bnx_full_iptw_od$`User design` <- bnx_full_iptw_od$Design

# p_od_bnx <- forest(
#   full_iptw_od[, c("Dose adjustment after interruption", 
#                    "User design", "HR(95% CI)", "Overdose-related acute care visits")],
#   est   = full_iptw_od$HR,
#   lower = full_iptw_od$lower,
#   upper = full_iptw_od$upper,
#   ci_column = 4,
#   ref_line = 1,
#   sizes = .6,
#   x_trans = "log",
#   xlim = c(0.45, 3),
#   ticks_at = c(0.5, 1, 1.5),
#   xlab = "Adjusted hazard ratio",
#   arrow_lab = c("Beneficial", "Harmful"),
#   title = "B. Buprenorphine/naloxone (BNX) trials",
#   footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5\n** NA for few events", 
#   nudge_y = 0.2,
#   theme = 
#     forest_theme(
#       base_size = 11,
#       ci_pch = 15,
#       ci_lwd = 1.7,
#       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
#       footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
#                          fontface="plain", col = "black"),
#       arrow_length = 0.06,
#       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
#       xlab_adjust = c("refline", "center")
#     ))
# 
# 
# 
# 
# 
# library(patchwork)
# plot_od <- plot_spacer() /  # top spacer row
#   p_od_met / 
#   plot_spacer() /
#   p_od_bnx +
#   plot_layout(
#     heights = c(0.09, 1, 0.02, 1),
#     guides = "collect"
#   ) &
#   theme(
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     legend.justification = "left",
#     legend.box.just = "left"
#   ) &
#   guides(color = guide_legend(nrow = 1))
# 
# 
# 
# 
# name = "full_od_iptw_Nov12"
# figure=paste0(folder, name, ".png")
# ppi = 300
# png(figure, width=3000, height=3900, res=300)
# plot_od
# dev.off()


#### Trial 1: dis: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr1_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 1: MET interruption of 1-2 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: dis: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

#### Trial 3: dis: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 4: dis: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr4_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 4: BNX interruption of 1-5 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: dis: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, tr5_full_iptw_dis_redgt50, tr5_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### met forestplot ##########################
met_full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                               tr2_full_iptw_dis,
                               tr3_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


met_full_iptw_dis$`OAT discontinuation` <- "" 
met_full_iptw_dis$`HR(95% CI)` <- met_full_iptw_dis$HR_display
met_full_iptw_dis$`Dose adjustment after interruption` <- met_full_iptw_dis$Comparison
met_full_iptw_dis$`User design` <- met_full_iptw_dis$Design

# 
# library(ggplot2)
# library(forestploter)
# library(tidyverse)
# p_dis_met <- forest(
#   full_iptw_dis[, c("HR(95% CI)", "OAT discontinuation")],
#   est   = full_iptw_dis$HR,
#   lower = full_iptw_dis$lower,
#   upper = full_iptw_dis$upper,
#   ci_column = 2,
#   ref_line = 1,
#   sizes = .6,
#   arrow_lab = c("Beneficial", "Harmful"),
#   x_trans = "log",
#   xlim = c(0.45, 3),
#   ticks_at = c(0.5, 1, 1.5),
#   #xlab = "Adjusted hazard ratio",
#   title = NULL,
#   nudge_y = 0.2,
#   theme = 
#     forest_theme(
#       base_size = 11,
#       ci_pch = 15,
#       ci_lwd = 1.7,
#       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
#       footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
#                          fontface="plain", col = "black"),
#       arrow_length = 0.06,
#       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
#       xlab_adjust = c("refline", "center")
#     ))


#### bnx forestplot ##########################
bnx_full_iptw_dis <- bind_rows(tr4_full_iptw_dis,
                               tr5_full_iptw_dis) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


bnx_full_iptw_dis$`OAT discontinuation` <- "" 
bnx_full_iptw_dis$`HR(95% CI)` <- bnx_full_iptw_dis$HR_display
bnx_full_iptw_dis$`Dose adjustment after interruption` <- bnx_full_iptw_dis$Comparison
bnx_full_iptw_dis$`User design` <- bnx_full_iptw_dis$Design



# p_dis_bnx <- forest(
#   full_iptw_dis[, c("HR(95% CI)", "OAT discontinuation")],
#   est   = full_iptw_dis$HR,
#   lower = full_iptw_dis$lower,
#   upper = full_iptw_dis$upper,
#   ci_column = 2,
#   ref_line = 1,
#   sizes = .6,
#   x_trans = "log",
#   xlim = c(0.45, 3),
#   ticks_at = c(0.5, 1, 1.5),
#   xlab = "Adjusted hazard ratio",
#   arrow_lab = c("Beneficial", "Harmful"),
#   title = NULL,
#   footnote = "* <=40mg MET for trial 2, 3, and <=4/1mg BNX for trial 5", 
#   nudge_y = 0.2,
#   theme = 
#     forest_theme(
#       base_size = 11,
#       ci_pch = 15,
#       ci_lwd = 1.7,
#       refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
#       footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8, 
#                          fontface="plain", col = "black"),
#       arrow_length = 0.06,
#       arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
#       xlab_adjust = c("refline", "center")
#     ))

#### creating the graph ######

met_full_iptw_dis_v1 <- met_full_iptw_dis %>% select(`HR(95% CI)`, `OAT discontinuation`)
met_full_iptw <- bind_cols(met_full_iptw_od, met_full_iptw_dis_v1)

names(met_full_iptw)[c(8, 11)] <- c("aHR(95% CI)**", "    aHR(95% CI)")
met_full_iptw$`    aHR(95% CI)` <- paste0("    ", met_full_iptw$`    aHR(95% CI)`)

met_full_iptw <- met_full_iptw %>% mutate(`Dose adjustment after interruption`=recode(`Dose adjustment after interruption`,
                                                   "   Dose reduction to >=50%/initial dose* vs no change"=
                                                     "   Dose reduction to >=50%/initial dose* vs\n   no change"))

p_od_dis_met <- forest(
  met_full_iptw[, c(9,10,8,7,11,12)],
  est   = list(met_full_iptw_od$HR, met_full_iptw_dis$HR),
  lower = list(met_full_iptw_od$lower, met_full_iptw_dis$lower),
  upper = list(met_full_iptw_od$upper, met_full_iptw_dis$upper),
  ci_column = c(4, 6),
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5, 3),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "A. Methadone (MET) trials",
  footnote = NULL,
  nudge_y = 0.02,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))



bnx_full_iptw_dis_v1 <- bnx_full_iptw_dis %>% select(`HR(95% CI)`, `OAT discontinuation`)
bnx_full_iptw <- bind_cols(bnx_full_iptw_od, bnx_full_iptw_dis_v1)
names(bnx_full_iptw)
names(bnx_full_iptw)[c(8, 11)] <- c("aHR(95% CI)", "    aHR(95% CI)")
bnx_full_iptw$`    aHR(95% CI)` <- paste0("    ", bnx_full_iptw$`    aHR(95% CI)`)

bnx_full_iptw <- bnx_full_iptw %>% mutate(`Dose adjustment after interruption`=recode(`Dose adjustment after interruption`,
                                                                                      "   Dose reduction to >=50%/initial dose* vs no change"=
                                                                                        "   Dose reduction to >=50%/initial dose* vs\n   no change"))



p_od_dis_bnx <- forest(
  bnx_full_iptw[, c(9,10,8,7,11,12)],
  est   = list(bnx_full_iptw_od$HR, bnx_full_iptw_dis$HR),
  lower = list(bnx_full_iptw_od$lower, bnx_full_iptw_dis$lower),
  upper = list(bnx_full_iptw_od$upper, bnx_full_iptw_dis$upper),
  ci_column = c(4, 6),
  ref_line = 1,
  sizes = .6,
  x_trans = "log",
  xlim = c(0.35, 4),
  ticks_at = c(0.5, 1, 1.5, 3),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = NULL,
  nudge_y = 0.02,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      footnote_parse = FALSE,
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


txt <- "<sup>*</sup> <=40mg MET for Trial 2, 3 and <=4/1mg BNX for Trial 5<br><sup>**</sup> adjusted hazard ratio with 95% compatibility intervals<br><sup>a</sup> NA for few events"

p_od_dis_bnx <- add_grob(p_od_dis_bnx, 
                         row = 18, 
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
  p_od_dis_bnx/plot_spacer() +
  plot_layout(
    heights = c(0.12, 1, 0.08, 1, 0.01),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))



name = "Figure2_od_dis_iptw_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=4700, height=4700, res=300)
plot_od_dis
dev.off()



# Figure 3: Ret-spline  #####################

res1_od <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr1_12d_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 1: MET interruption of 1-2 days", 
         Outcome = "Overdose-related acute care visits")

res2_od <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr2_34d_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 2: MET interruption of 3-4 days", 
         Outcome = "Overdose-related acute care visits")

res3_od <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr3_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 3: MET interruption of 5-14 days", 
         Outcome = "Overdose-related acute care visits")

res4_od <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr4_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 4: BNX interruption of 1-5 days", 
         Outcome = "Overdose-related acute care visits")

res5_od <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr5_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 5: BNX interruption of 6-14 days", 
         Outcome = "Overdose-related acute care visits")



res1_dis <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr1_12d_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 1: MET interruption of 1-2 days", 
         Outcome = "OAT discontinuation")

res2_dis <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_od_tr2_34d_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 2: MET interruption of 3-4 days", 
         Outcome = "OAT discontinuation")

res3_dis <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr3_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 3: MET interruption of 5-14 days", 
         Outcome = "OAT discontinuation")

res4_dis <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr4_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 4: BNX interruption of 1-5 days", 
         Outcome = "OAT discontinuation")

res5_dis <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Sens_spline/ret_dis_tr5_rprev_ow.rds") %>% 
  select(ret_length, contrast, HR, HR_low, HR_high) %>% 
  mutate(Trial = "Trial 5: BNX interruption of 6-14 days", 
         Outcome = "OAT discontinuation")


res5_dis %>% count(contrast)


comparison_colors_v1 = c(
  "Dose reduction\nvs\nno change"              = "#000000",  # Blue
  "Dose reduction to <50%\nvs\nno change"            = "#2CA02C",  # Green
  "Dose reduction to >=50%/initial dose*\nvs\nno change"     = "#FF7F0E",  # Sky Blue
  "Dose increase\nvs\nno change"               = "#D62728"  # Orange
)


plot_trial_results_ret_v4 <- function(df, title_lab, 
                                      y_lab="Adjusted hazard ratios",
                                      x_lab = "Duration of OAT episode prior to the index interruption (in days)",
                                      y_min, y_max, y_points,
                                      legend_option=TRUE) {
  ggplot(df, aes(x = ret_length, y = HR, color = contrast)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 30, linetype = "dashed", color = "red4") +
    facet_wrap(~ Trial, ncol = 1, scales = "free") +  # Use facet_wrap to avoid axis clash
    labs(
      x = x_lab,
      y = y_lab,
      title = title_lab
    ) +
    scale_color_manual(values = comparison_colors_v1) +
    scale_x_continuous(breaks = c(0, 30, 90, 180, 270, 365)) + 
    scale_y_log10(breaks = y_points, limits = c(y_min, y_max))+
    theme_minimal() +
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.x = element_line(linewidth = 1.1, color = "black"),
      axis.ticks.y = element_line(linewidth = 1.1, color = "black"),
      axis.ticks.length  = unit(0.15, "cm"),
      axis.text = element_text(color = "black", size = 14),
      axis.title.x = element_text(color = "black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      legend.position = if(legend_option) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 14)
    )
  
}


### OD #########

res_od <- bind_rows(res1_od, res2_od, res3_od, res4_od, res5_od) %>% 
  filter(ret_length <= 365)

res_od_met = res_od %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                            "Trial 2: MET interruption of 3-4 days", 
                                            "Trial 3: MET interruption of 5-14 days"))

res_od_bnx = res_od %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                            "Trial 5: BNX interruption of 6-14 days"))


od_plot_met = plot_trial_results_ret_v4(df=res_od_met, 
                                        title_lab = "Overdose-related acute care visits",
                                        #x_lab = "",
                                        y_min = 0.5,
                                        y_max = 2, 
                                        y_points = c(0.5, 1, 1.5, 2),
                                        legend_option = FALSE) 


od_plot_bnx = plot_trial_results_ret_v4(df=res_od_bnx, 
                                        title_lab = "Overdose-related acute care visits",
                                        y_min = 0.5,
                                        y_max = 2, 
                                        y_points = c(0.5, 1, 1.5, 2),
                                        legend_option = FALSE) 



### DIS #########
res_dis <- bind_rows(res1_dis, res2_dis, res3_dis, res4_dis, res5_dis) %>% 
  filter(ret_length <= 365)


res_dis_met = res_dis %>% filter(Trial %in% c("Trial 1: MET interruption of 1-2 days", 
                                              "Trial 2: MET interruption of 3-4 days", 
                                              "Trial 3: MET interruption of 5-14 days"))

res_dis_bnx = res_dis %>% filter(Trial %in% c("Trial 4: BNX interruption of 1-5 days", 
                                              "Trial 5: BNX interruption of 6-14 days"))


dis_plot_met = plot_trial_results_ret_v4(df=res_dis_met, 
                                         title_lab = "OAT discontinuation",
                                         #x_lab = "",
                                         y_lab = "",
                                         y_min = 0.5,
                                         y_max = 2, 
                                         y_points = c(0.5, 1, 1.5, 2),
                                         legend_option = FALSE) 



dis_plot_bnx = plot_trial_results_ret_v4(df=res_dis_bnx, 
                                         title_lab = "OAT discontinuation",
                                         
                                         y_lab = "",
                                         y_min = 0.5,
                                         y_max = 2, 
                                         y_points = c(0.5, 1, 1.5, 2)) 




library(patchwork)
plot_met <- od_plot_met + plot_spacer() + dis_plot_met + # use `+` to arrange side by side
  plot_layout(ncol = 3, widths = c(1, 0.05, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_met1 <- plot_spacer() /  # top spacer row
  plot_met +
  plot_layout(
    heights = c(0.03, 1), guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))



plot_met2 <- ggdraw(plot_met1) + 
  draw_label("* <=40mg MET for Trial 2 and 3", 
             x = 0.99, y = 0.01, hjust = 1, vjust = 0, 
             size = 12, color = "black") +
  draw_label("A. Methadone (MET) trials among all prevalent users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
name = "Figure3_Panel_A_Nov_12"
figure=paste0(folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=4000, height=4800, res=300)
plot_met2
dev.off()





plot_bnx <- od_plot_bnx + plot_spacer() + dis_plot_bnx + # use `+` to arrange side by side
  plot_layout(ncol = 3, widths = c(1, 0.05, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))


plot_bnx1 <- plot_spacer() /  # top spacer row
  plot_bnx +
  plot_layout(
    heights = c(0.03, 1), guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))




plot_bnx2 <- ggdraw(plot_bnx1) + 
  draw_label("* <=4/1mg BNX for Trial 5", 
             x = 0.99, y = 0.01, hjust = 1, vjust = 0, 
             size = 12, color = "black") +
  draw_label("B. Buprenorphine/naloxone (BNX) trials among all prevalent users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
name = "Figure3_Panel_B_Nov_12"
figure=paste0(folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=4000, height=3400, res=300)
plot_bnx2
dev.off()







plot_bnx <- od_plot_bnx + plot_spacer() + dis_plot_bnx + # use `+` to arrange side by side
  plot_layout(ncol = 3, widths = c(1, 0.05, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))

plot_bnx1 <- plot_spacer() /  # top spacer row
  plot_bnx +
  plot_layout(
    heights = c(0.02, 1), guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 1))



plot_bnx2 <- ggdraw(plot_bnx1) + 
  draw_label("* <=4mg buprenorphine/naloxone", 
             x = 0.99, y = 0.01, hjust = 1, vjust = 0, 
             size = 12, color = "black") +
  draw_label("B. Buprenorphine/naloxone (BNX) trials", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black")




folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
name = "Figure3_Panel_B_Nov12_v1"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=4800, height=4200, res=300)
plot_bnx2
dev.off()




# Figure 2: Panel A ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change")
  ) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design)


res_df2 %>% count(res_type)

#### HRs- met forestplot ##########################

#### Trial 1: od: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 1: MET interruption of 1-2 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_od <- bind_rows(header1, tr1_full_iptw_od_red, tr1_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: od: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 2: MET interruption of 3-4 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))

tr2_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 2: MET interruption of 3-4 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison),
         Design = recode(Design,
                         "Incident users"="Incident users**"))



tr2_full_iptw_od <- bind_rows(header2, tr2_full_iptw_od_redlt50, tr2_full_iptw_od_redgt50, tr2_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " "),
    HR    = ifelse(Design == "Incident users**", NA, HR),
    lower = ifelse(Design == "Incident users**", NA, lower),
    upper = ifelse(Design == "Incident users**", NA, upper)
  ) %>%
  ungroup()

#### Trial 3: od: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 3: MET interruption of 5-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 3: MET interruption of 5-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_od <- bind_rows(header3, tr3_full_iptw_od_redlt50, tr3_full_iptw_od_redgt50, tr3_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()



met_full_iptw_od <- bind_rows(tr1_full_iptw_od,
                              tr2_full_iptw_od,
                              tr3_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


met_full_iptw_od$`overdose-related acute care visits` <- "" 
met_full_iptw_od$`Adjusted hazard\nratio (95% CI)` <- met_full_iptw_od$HR_display
met_full_iptw_od$`Dose adjustment after interruption` <- met_full_iptw_od$Comparison
met_full_iptw_od$`User design` <- met_full_iptw_od$Design


#### Risk- met forestplot ###############
extract_365 <- function(res_obj, type = c("od", "dis")) {
  type <- match.arg(type)
  
  pops <- names(res_obj$points_ci[[type]])
  
  out <- lapply(pops, function(p) {
    pop_obj <- res_obj$points_ci[[type]][[p]]
    
    # skip if NULL, NA, or not a list
    if (is.null(pop_obj) || all(is.na(pop_obj)) || !is.list(pop_obj))
      return(NULL)
    
    # skip if contrasts is missing or NULL
    if (is.null(pop_obj$contrasts) || all(is.na(pop_obj$contrasts)))
      return(NULL)
    
    df <- pop_obj$contrasts
    
    # skip if empty or horizon column missing
    if (nrow(df) == 0 || !("horizon" %in% names(df)))
      return(NULL)
    
    df %>%
      dplyr::filter(horizon == 365) %>%
      dplyr::mutate(pop = p, type = type)
  })
  
  # bind only non-null results
  dplyr::bind_rows(out)
}


res_1 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_iptw.rds")
res_2 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk2_34d_iptw.rds")
res_3 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk3_iptw.rds")
res_4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")
res_5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")

#### Trial 1: 1-2 days ##############
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)

tr1_risk_iptw_od_red <- extract_365(res_1, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_risk_iptw_od_incr <- extract_365(res_1, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_risk_iptw_od <- bind_rows(header1, tr1_risk_iptw_od_red, tr1_risk_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### Trial 2: 3-4 days #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)
lt50_risk <- extract_365(res_2, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper)

tr2_risk_iptw_od_redlt50 <- bind_rows(lt50_risk, data.frame(
  Comparison = c("Dose reduction to <50% vs no change",
                 "Dose reduction to >=50%/initial dose* vs no change",
                 "Dose increase vs no change"),
  Design = "Incident users**",
  RD = NA, RD_lower = NA, RD_upper = NA
  
)) %>% mutate(Design = factor(Design, levels = c("Incident users**", "Prevalent new users",
                                                 "All prevalent users"))) %>% 
  arrange(Design) %>% 
  filter(Comparison == "Dose reduction to <50% vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr2_risk_iptw_od_redgt50 <- bind_rows(lt50_risk, data.frame(
  Comparison = c("Dose reduction to <50% vs no change",
                 "Dose reduction to >=50%/initial dose* vs no change",
                 "Dose increase vs no change"),
  Design = "Incident users**",
  RD = NA, RD_lower = NA, RD_upper = NA
  
)) %>% mutate(Design = factor(Design, levels = c("Incident users**", "Prevalent new users",
                                                 "All prevalent users"))) %>% 
  arrange(Design) %>% 
  filter(Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr2_risk_iptw_od_incr <- bind_rows(lt50_risk, data.frame(
  Comparison = c("Dose reduction to <50% vs no change",
                 "Dose reduction to >=50%/initial dose* vs no change",
                 "Dose increase vs no change"),
  Design = "Incident users**",
  RD = NA, RD_lower = NA, RD_upper = NA
  
)) %>% mutate(Design = factor(Design, levels = c("Incident users**", "Prevalent new users",
                                                 "All prevalent users"))) %>% 
  arrange(Design) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr2_risk_iptw_od <- bind_rows(header2, 
                              tr2_risk_iptw_od_redlt50, 
                              tr2_risk_iptw_od_redgt50, 
                              tr2_risk_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()

#### Trial 3  #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)
tr3_risk_iptw_od_redlt50 <- extract_365(res_3, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction to <50% vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr3_risk_iptw_od_redgt50 <- extract_365(res_3, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr3_risk_iptw_od_incr <- extract_365(res_3, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr3_risk_iptw_od <- bind_rows(header3, 
                              tr3_risk_iptw_od_redlt50, 
                              tr3_risk_iptw_od_redgt50, 
                              tr3_risk_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()


names(met_risk_iptw_od)
met_risk_iptw_od <- bind_rows(tr1_risk_iptw_od,
                              tr2_risk_iptw_od,
                              tr3_risk_iptw_od) %>%
  mutate(RD_display = ifelse(
    is.na(RD),
    " ",
    sprintf("%.2f(%.2f,%.2f)", RD*100, RD_lower*100, RD_upper*100))
  )


met_risk_iptw_od$`overdose-related acute care visits` <- "" 
met_risk_iptw_od$`Cumulative incidence difference (95% CI)%` <- met_risk_iptw_od$RD_display
met_risk_iptw_od$`Dose adjustment after interruption` <- met_risk_iptw_od$Comparison
met_risk_iptw_od$`User design` <- met_risk_iptw_od$Design

#### creating the graph ######

met_risk_iptw_od_v1 <- met_risk_iptw_od %>% select(`Cumulative incidence difference (95% CI)%`, 
                                                   `overdose-related acute care visits`) %>% 
  rename(`Difference in cumulative\nincidence (95% CI)%`=`Cumulative incidence difference (95% CI)%`)
met_full_risk_iptw <- bind_cols(met_full_iptw_od, met_risk_iptw_od_v1)

names(met_full_risk_iptw)[c(7, 12)] <- c("                             ", 
                                         "                         ")

met_full_risk_iptw <- met_full_risk_iptw %>% 
  mutate(`Dose adjustment after interruption`= 
           recode(`Dose adjustment after interruption`,
                  "   Dose reduction to >=50%/initial dose* vs no change"=
                    "   Dose reduction to >=50%/initial dose* vs no change "))



p_od_dis_met <- forest(
  met_full_risk_iptw[, c(9,10,8,7,11,12)],
  est   = list(met_full_risk_iptw$HR, met_risk_iptw_od$RD*100),
  lower = list(met_full_iptw_od$lower, met_risk_iptw_od$RD_lower*100),
  upper = list(met_full_iptw_od$upper, met_risk_iptw_od$RD_upper*100),
  ci_column = c(4, 6),
  ref_line = c(1,0),
  sizes = .6,
  x_trans = c("log", "none"),
  xlim = list(c(0.35, 4), c(-5, 5)),
  ticks_at = list(c(0.5, 1, 1.5, 3), c(-2.5, 0, 2.5)),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "A. Methadone (MET) trials",
  footnote = NULL,
  nudge_y = 0.02,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


txt <- "<sup>*</sup> <=40mg MET for Trial 2 and 3 <br><sup>**</sup> NA for few events"

p_od_dis_met1 <- add_grob(p_od_dis_met, 
                          row = 28, 
                          col = 1:2,
                          order = "background",
                          gb_fn = gridtext::richtext_grob,
                          text = txt,
                          gp = gpar(fontsize = 11),
                          hjust = 0, vjust = 1, halign = 0, valign = 1,
                          x = unit(0, "npc"), y = unit(1, "npc"))


name = "Figure2_od_risk_Panel_A_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=4300, height=2700, res=300)
p_od_dis_met1
dev.off()




### Figure 2: Panel B ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change")
  ) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design)


res_df2 %>% count(res_type)

#### HRs- bnx forestplot ##########################

#### Trial 4: od: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_od_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                           Model == "IPTW",
                                           Outcome == "Overdose-related acute care visits",
                                           Trial == "Trial 4: BNX interruption of 1-5 days",
                                           Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr4_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_od <- bind_rows(header4, tr4_full_iptw_od_red, tr4_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: od: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_od_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                               Model == "IPTW",
                                               Outcome == "Overdose-related acute care visits",
                                               Trial == "Trial 5: BNX interruption of 6-14 days",
                                               Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_od_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "Overdose-related acute care visits",
                                            Trial == "Trial 5: BNX interruption of 6-14 days",
                                            Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_od <- bind_rows(header5, tr5_full_iptw_od_redlt50, 
                              tr5_full_iptw_od_redgt50, 
                              tr5_full_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()



bnx_full_iptw_od <- bind_rows(tr4_full_iptw_od,
                              tr5_full_iptw_od) %>%
  mutate(HR_display = ifelse(
    is.na(HR),
    " ",
    sprintf("%.2f(%.2f,%.2f)", HR, lower, upper))
  )


bnx_full_iptw_od$`overdose-related acute care visits` <- "" 
bnx_full_iptw_od$`Adjusted hazard\nratio (95% CI)` <- bnx_full_iptw_od$HR_display
bnx_full_iptw_od$`Dose adjustment after interruption` <- bnx_full_iptw_od$Comparison
bnx_full_iptw_od$`User design` <- bnx_full_iptw_od$Design


#### Risk- bnx forestplot ###############
extract_365 <- function(res_obj, type = c("od", "dis")) {
  type <- match.arg(type)
  
  pops <- names(res_obj$points_ci[[type]])
  
  out <- lapply(pops, function(p) {
    pop_obj <- res_obj$points_ci[[type]][[p]]
    
    # skip if NULL, NA, or not a list
    if (is.null(pop_obj) || all(is.na(pop_obj)) || !is.list(pop_obj))
      return(NULL)
    
    # skip if contrasts is missing or NULL
    if (is.null(pop_obj$contrasts) || all(is.na(pop_obj$contrasts)))
      return(NULL)
    
    df <- pop_obj$contrasts
    
    # skip if empty or horizon column missing
    if (nrow(df) == 0 || !("horizon" %in% names(df)))
      return(NULL)
    
    df %>%
      dplyr::filter(horizon == 365) %>%
      dplyr::mutate(pop = p, type = type)
  })
  
  # bind only non-null results
  dplyr::bind_rows(out)
}

res_4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")
res_5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")

#### Trial 1: 1-2 days ##############
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)

tr4_risk_iptw_od_red <- extract_365(res_4, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr4_risk_iptw_od_incr <- extract_365(res_4, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr4_risk_iptw_od <- bind_rows(header4, tr4_risk_iptw_od_red, tr4_risk_iptw_od_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### Trial 5: 6-14 days #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)
tr5_risk_od_redlt50 <- extract_365(res_5, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=4mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction to <50% vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr5_risk_od_redgt50 <- extract_365(res_5, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=4mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>%  
  filter(Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_risk_od_incr <- extract_365(res_5, "od") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=4mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "od" = "Overdose-related acute care visits"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>%  
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr5_risk_iptw_od <- bind_rows(header5, 
                              tr5_risk_od_redlt50, 
                              tr5_risk_od_redgt50, 
                              tr5_risk_od_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()

names(bnx_risk_iptw_od)
bnx_risk_iptw_od <- bind_rows(tr4_risk_iptw_od,
                              tr5_risk_iptw_od) %>%
  mutate(RD_display = ifelse(
    is.na(RD),
    " ",
    sprintf("%.2f(%.2f,%.2f)", RD*100, RD_lower*100, RD_upper*100))
  )


bnx_risk_iptw_od$`overdose-related acute care visits` <- "" 
bnx_risk_iptw_od$`Cumulative incidence difference (95% CI)%` <- bnx_risk_iptw_od$RD_display
bnx_risk_iptw_od$`Dose adjustment after interruption` <- bnx_risk_iptw_od$Comparison
bnx_risk_iptw_od$`User design` <- bnx_risk_iptw_od$Design

#### creating the graph ######

bnx_risk_iptw_od_v1 <- bnx_risk_iptw_od %>% select(`Cumulative incidence difference (95% CI)%`, 
                                                   `overdose-related acute care visits`) %>% 
  rename(`Difference in cumulative\nincidence (95% CI)%`=`Cumulative incidence difference (95% CI)%`)
bnx_full_risk_iptw <- bind_cols(bnx_full_iptw_od, bnx_risk_iptw_od_v1)

names(bnx_full_risk_iptw)[c(7, 12)] <- c("                             ", 
                                         "                         ")

bnx_full_risk_iptw <- bnx_full_risk_iptw %>% 
  mutate(`Dose adjustment after interruption`= 
           recode(`Dose adjustment after interruption`,
                  "   Dose reduction to >=50%/initial dose* vs no change"=
                    "   Dose reduction to >=50%/initial dose* vs no change "))



p_od_bnx <- forest(
  bnx_full_risk_iptw[, c(9,10,8,7,11,12)],
  est   = list(bnx_full_risk_iptw$HR, bnx_risk_iptw_od$RD*100),
  lower = list(bnx_full_iptw_od$lower, bnx_risk_iptw_od$RD_lower*100),
  upper = list(bnx_full_iptw_od$upper, bnx_risk_iptw_od$RD_upper*100),
  ci_column = c(4, 6),
  ref_line = c(1,0),
  sizes = .6,
  x_trans = c("log", "none"),
  xlim = list(c(0.35, 4), c(-5, 5)),
  ticks_at = list(c(0.5, 1, 1.5, 3), c(-2.5, 0, 2.5)),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = NULL,
  nudge_y = 0.02,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


txt <- "<sup>*</sup> <=40mg BNX for Trial 5"

p_od_bnx1 <- add_grob(p_od_bnx, 
                      row = 18, 
                      col = 1:2,
                      order = "background",
                      gb_fn = gridtext::richtext_grob,
                      text = txt,
                      gp = gpar(fontsize = 11),
                      hjust = 0, vjust = 1, halign = 0, valign = 1,
                      x = unit(0, "npc"), y = unit(1, "npc"))


name = "Figure2_od_risk_Panel_B_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=4300, height=1900, res=300)
p_od_bnx1
dev.off()









# Figure 4: Panel A ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change")
  ) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design)


#### HRs- met forestplot ##########################

#### Trial 1: dis: iptw #######
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr1_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 1: MET interruption of 1-2 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 1: MET interruption of 1-2 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr1_full_iptw_dis <- bind_rows(header1, tr1_full_iptw_dis_red, tr1_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 2: dis: iptw #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr2_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 2: MET interruption of 3-4 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr2_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 2: MET interruption of 3-4 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr2_full_iptw_dis <- bind_rows(header2, tr2_full_iptw_dis_redlt50, tr2_full_iptw_dis_redgt50, tr2_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()

#### Trial 3: dis: iptw #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr3_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 3: MET interruption of 5-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr3_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 3: MET interruption of 5-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr3_full_iptw_dis <- bind_rows(header3, tr3_full_iptw_dis_redlt50, tr3_full_iptw_dis_redgt50, tr3_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()

format_num <- function(x) {
  ifelse(
    abs(x) >= 10,
    sprintf("%.1f", x),
    sprintf("%.2f", x)
  )
}

met_full_iptw_dis <- bind_rows(tr1_full_iptw_dis,
                               tr2_full_iptw_dis,
                               tr3_full_iptw_dis) %>%
  mutate(
    HR_display = dplyr::case_when(
      is.na(HR) ~ " ",
      TRUE ~ sprintf(
        "%s (%s, %s)",
        format_num(HR),
        format_num(lower),
        format_num(upper)
      )
    )
  )


met_full_iptw_dis$`OAT discontinuation` <- "" 
met_full_iptw_dis$`Adjusted hazard\nratio (95% CI)` <- met_full_iptw_dis$HR_display
met_full_iptw_dis$`Dose adjustment after interruption` <- met_full_iptw_dis$Comparison
met_full_iptw_dis$`User design` <- met_full_iptw_dis$Design


#### Risk- met forestplot ###############
extract_365 <- function(res_obj, type = c("od", "dis")) {
  type <- match.arg(type)
  
  pops <- names(res_obj$points_ci[[type]])
  
  out <- lapply(pops, function(p) {
    pop_obj <- res_obj$points_ci[[type]][[p]]
    
    # skip if NULL, NA, or not a list
    if (is.null(pop_obj) || all(is.na(pop_obj)) || !is.list(pop_obj))
      return(NULL)
    
    # skip if contrasts is missing or NULL
    if (is.null(pop_obj$contrasts) || all(is.na(pop_obj$contrasts)))
      return(NULL)
    
    df <- pop_obj$contrasts
    
    # skip if empty or horizon column missing
    if (nrow(df) == 0 || !("horizon" %in% names(df)))
      return(NULL)
    
    df %>%
      dplyr::filter(horizon == 365) %>%
      dplyr::mutate(pop = p, type = type)
  })
  
  # bind only non-null results
  dplyr::bind_rows(out)
}


res_1 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk1_12d_iptw.rds")
res_2 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk2_34d_iptw.rds")
res_3 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk3_iptw.rds")

#### Trial 1: 1-2 days ##############
header1 <- data.frame(Comparison = "Trial 1: After MET interruption of 1-2 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)

tr1_risk_iptw_dis_red <- extract_365(res_1, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_risk_iptw_dis_incr <- extract_365(res_1, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr1_risk_iptw_dis <- bind_rows(header1, tr1_risk_iptw_dis_red, tr1_risk_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### Trial 2: 3-4 days #######
header2 <- data.frame(Comparison = "Trial 2: After MET interruption of 3-4 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)

tr2_risk_iptw_dis_redlt50 <- extract_365(res_2, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction to <50% vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr2_risk_iptw_dis_redgt50 <- extract_365(res_2, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr2_risk_iptw_dis_incr <- extract_365(res_2, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr2_risk_iptw_dis <- bind_rows(header2, 
                               tr2_risk_iptw_dis_redlt50, 
                               tr2_risk_iptw_dis_redgt50, 
                               tr2_risk_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()

#### Trial 3  #######
header3 <- data.frame(Comparison = "Trial 3: After MET interruption of 5-14 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)

tr3_risk_iptw_dis_redlt50 <- extract_365(res_3, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction to <50% vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr3_risk_iptw_dis_redgt50 <- extract_365(res_3, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr3_risk_iptw_dis_incr <- extract_365(res_3, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=40mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr3_risk_iptw_dis <- bind_rows(header3, 
                               tr3_risk_iptw_dis_redlt50, 
                               tr3_risk_iptw_dis_redgt50, 
                               tr3_risk_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()

format_num <- function(x) {
  ifelse(
    abs(x) >= 10,
    sprintf("%.1f", x),
    sprintf("%.2f", x)
  )
}
met_risk_iptw_dis <- bind_rows(
  tr1_risk_iptw_dis,
  tr2_risk_iptw_dis,
  tr3_risk_iptw_dis
) %>%
  mutate(
    RD_display = dplyr::case_when(
      is.na(RD) ~ " ",
      TRUE ~ sprintf(
        "%s (%s, %s)",
        format_num(RD * 100),
        format_num(RD_lower * 100),
        format_num(RD_upper * 100)
      )
    )
  )




met_risk_iptw_dis$`OAT discontinuation` <- "" 
met_risk_iptw_dis$`Cumulative incidence difference (95% CI)%` <- met_risk_iptw_dis$RD_display
met_risk_iptw_dis$`Dose adjustment after interruption` <- met_risk_iptw_dis$Comparison
met_risk_iptw_dis$`User design` <- met_risk_iptw_dis$Design

#### creating the graph ######

met_risk_iptw_dis_v1 <- met_risk_iptw_dis %>% select(`Cumulative incidence difference (95% CI)%`, 
                                                     `OAT discontinuation`) %>% 
  rename(`Difference in cumulative\nincidence (95% CI)%`=`Cumulative incidence difference (95% CI)%`)
met_full_risk_iptw <- bind_cols(met_full_iptw_dis, met_risk_iptw_dis_v1)

names(met_full_risk_iptw)[c(7, 12)] <- c("                             ", 
                                         "                                ")

met_full_risk_iptw <- met_full_risk_iptw %>% 
  mutate(`Dose adjustment after interruption`= 
           recode(`Dose adjustment after interruption`,
                  "   Dose reduction to >=50%/initial dose* vs no change"=
                    "   Dose reduction to >=50%/initial dose* vs no change "))

library(forestploter)

p_dis_met <- forest(
  met_full_risk_iptw[, c(9,10,8,7,11,12)],
  est   = list(met_full_risk_iptw$HR, met_risk_iptw_dis$RD*100),
  lower = list(met_full_iptw_dis$lower, met_risk_iptw_dis$RD_lower*100),
  upper = list(met_full_iptw_dis$upper, met_risk_iptw_dis$RD_upper*100),
  ci_column = c(4, 6),
  ref_line = c(1,0),
  sizes = .6,
  x_trans = c("log", "none"),
  xlim = list(c(0.35, 4), c(-17, 17)),
  ticks_at = list(c(0.5, 1, 1.5, 3), c(-10, -5, 0, 5, 10)),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "A. Methadone (MET) trials",
  footnote = NULL,
  nudge_y = 0.02,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


txt <- "<sup>*</sup> <=40mg MET for Trial 2 and 3"

p_dis_met1 <- add_grob(p_dis_met, 
                       row = 28, 
                       col = 1:2,
                       order = "background",
                       gb_fn = gridtext::richtext_grob,
                       text = txt,
                       gp = gpar(fontsize = 11),
                       hjust = 0, vjust = 1, halign = 0, valign = 1,
                       x = unit(0, "npc"), y = unit(1, "npc"))


name = "Figure4_dis_risk_Panel_A_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=4350, height=2700, res=300)
p_dis_met1
dev.off()




### Figure 2: Panel B ######
folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
res_df <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/all_HRs_iptw_ow.rds")

res_df2 <- res_df %>%
  mutate(
    Comparison = recode(Comparison, 
                        "Dose reduction\nvs\nNo dose change" = "Dose reduction vs no change",
                        "Dose increase\nvs\nNo dose change" = "Dose increase vs no change",
                        "1-<50% reduction\nvs\nNo dose change" = "Dose reduction to <50% vs no change",
                        ">=50% reduction/<=40mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change",
                        ">=50% reduction/<=4mg\nvs\nNo dose change" = "Dose reduction to >=50%/initial dose* vs no change")
  ) %>% arrange(res_type, Model, Outcome, Trial, Comparison, Design)


res_df2 %>% count(res_type)

#### HRs- bnx forestplot ##########################

#### Trial 4: dis: iptw #######
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr4_full_iptw_dis_red <- res_df2 %>% filter(res_type == "Primary analysis",
                                            Model == "IPTW",
                                            Outcome == "OAT discontinuation",
                                            Trial == "Trial 4: BNX interruption of 1-5 days",
                                            Comparison == "Dose reduction vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr4_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 4: BNX interruption of 1-5 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr4_full_iptw_dis <- bind_rows(header4, tr4_full_iptw_dis_red, tr4_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()


#### Trial 5: dis: iptw #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      HR = NA, lower = NA, upper = NA)
tr5_full_iptw_dis_redlt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to <50% vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_redgt50 <- res_df2 %>% filter(res_type == "Primary analysis",
                                                Model == "IPTW",
                                                Outcome == "OAT discontinuation",
                                                Trial == "Trial 5: BNX interruption of 6-14 days",
                                                Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_full_iptw_dis_incr <- res_df2 %>% filter(res_type == "Primary analysis",
                                             Model == "IPTW",
                                             Outcome == "OAT discontinuation",
                                             Trial == "Trial 5: BNX interruption of 6-14 days",
                                             Comparison == "Dose increase vs no change") %>% 
  select(Comparison, Design, HR, lower, upper) %>% 
  mutate(Comparison = paste0("   ", Comparison))



tr5_full_iptw_dis <- bind_rows(header5, tr5_full_iptw_dis_redlt50, 
                               tr5_full_iptw_dis_redgt50, 
                               tr5_full_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()



format_num <- function(x) {
  ifelse(
    abs(x) >= 10,
    sprintf("%.1f", x),
    sprintf("%.2f", x)
  )
}


bnx_full_iptw_dis <- bind_rows(tr4_full_iptw_dis,
                               tr5_full_iptw_dis) %>%
  mutate(
    HR_display = dplyr::case_when(
      is.na(HR) ~ " ",
      TRUE ~ sprintf(
        "%s (%s, %s)",
        format_num(HR),
        format_num(lower),
        format_num(upper)
      )
    )
  )


bnx_full_iptw_dis$`OAT discontinuation` <- "" 
bnx_full_iptw_dis$`Adjusted hazard\nratio (95% CI)` <- bnx_full_iptw_dis$HR_display
bnx_full_iptw_dis$`Dose adjustment after interruption` <- bnx_full_iptw_dis$Comparison
bnx_full_iptw_dis$`User design` <- bnx_full_iptw_dis$Design


#### Risk- bnx forestplot ###############
extract_365 <- function(res_obj, type = c("od", "dis")) {
  type <- match.arg(type)
  
  pops <- names(res_obj$points_ci[[type]])
  
  out <- lapply(pops, function(p) {
    pop_obj <- res_obj$points_ci[[type]][[p]]
    
    # skip if NULL, NA, or not a list
    if (is.null(pop_obj) || all(is.na(pop_obj)) || !is.list(pop_obj))
      return(NULL)
    
    # skip if contrasts is missing or NULL
    if (is.null(pop_obj$contrasts) || all(is.na(pop_obj$contrasts)))
      return(NULL)
    
    df <- pop_obj$contrasts
    
    # skip if empty or horizon column missing
    if (nrow(df) == 0 || !("horizon" %in% names(df)))
      return(NULL)
    
    df %>%
      dplyr::filter(horizon == 365) %>%
      dplyr::mutate(pop = p, type = type)
  })
  
  # bind only non-null results
  dplyr::bind_rows(out)
}

res_4 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk4_iptw.rds")
res_5 <- readRDS(file = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/risk5_iptw.rds")

#### Trial 1: 1-2 days ##############
header4 <- data.frame(Comparison = "Trial 4: After BNX interruption of 1-5 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)

tr4_risk_iptw_dis_red <- extract_365(res_4, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr4_risk_iptw_dis_incr <- extract_365(res_4, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "reversion" = "Dose reduction vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr4_risk_iptw_dis <- bind_rows(header4, tr4_risk_iptw_dis_red, tr4_risk_iptw_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(Comparison = ifelse(row_number() == 1, Comparison, " ")) %>%
  ungroup()




#### Trial 5: 6-14 days #######
header5 <- data.frame(Comparison = "Trial 5: After BNX interruption of 6-14 days",
                      Design = "",
                      RD = NA, RD_lower = NA, RD_upper = NA)
tr5_risk_dis_redlt50 <- extract_365(res_5, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=4mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>% 
  filter(Comparison == "Dose reduction to <50% vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr5_risk_dis_redgt50 <- extract_365(res_5, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=4mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>%  
  filter(Comparison == "Dose reduction to >=50%/initial dose* vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))

tr5_risk_dis_incr <- extract_365(res_5, "dis") %>% 
  mutate(Comparison = trt, 
         Comparison = recode(trt,
                             "1-<50%" = "Dose reduction to <50% vs no change",
                             ">=50%/<=4mg" = "Dose reduction to >=50%/initial dose* vs no change",
                             "increase" = "Dose increase vs no change"),
         Outcome = recode(type,
                          "dis" = "OAT discontinuation"),
         Design = recode(pop,
                         "incid" = "Incident users",
                         "prev" = "Prevalent new users",
                         "rzero" = "All prevalent users")) %>% 
  select(Comparison, Design, RD, RD_lower, RD_upper) %>%  
  filter(Comparison == "Dose increase vs no change") %>% 
  mutate(Comparison = paste0("   ", Comparison))


tr5_risk_iptw_dis <- bind_rows(header5, 
                               tr5_risk_dis_redlt50, 
                               tr5_risk_dis_redgt50, 
                               tr5_risk_dis_incr) %>%
  group_by(Comparison) %>%
  mutate(
    Comparison = ifelse(row_number() == 1, Comparison, " ")
  ) %>%
  ungroup()


format_num <- function(x) {
  ifelse(
    abs(x) >= 10,
    sprintf("%.1f", x),
    sprintf("%.2f", x)
  )
}

bnx_risk_iptw_dis <- bind_rows(tr4_risk_iptw_dis,
                               tr5_risk_iptw_dis) %>%
  mutate(
    RD_display = dplyr::case_when(
      is.na(RD) ~ " ",
      TRUE ~ sprintf(
        "%s (%s, %s)",
        format_num(RD * 100),
        format_num(RD_lower * 100),
        format_num(RD_upper * 100)
      )
    )
  )



bnx_risk_iptw_dis$`OAT discontinuation` <- "" 
bnx_risk_iptw_dis$`Cumulative incidence difference (95% CI)%` <- bnx_risk_iptw_dis$RD_display
bnx_risk_iptw_dis$`Dose adjustment after interruption` <- bnx_risk_iptw_dis$Comparison
bnx_risk_iptw_dis$`User design` <- bnx_risk_iptw_dis$Design

#### creating the graph ######

bnx_risk_iptw_dis_v1 <- bnx_risk_iptw_dis %>% select(`Cumulative incidence difference (95% CI)%`, 
                                                     `OAT discontinuation`) %>% 
  rename(`Difference in cumulative\nincidence (95% CI)%`=`Cumulative incidence difference (95% CI)%`)
bnx_full_risk_iptw <- bind_cols(bnx_full_iptw_dis, bnx_risk_iptw_dis_v1)

names(bnx_full_risk_iptw)[c(7, 12)] <- c("                             ", 
                                         "                                ")


bnx_full_risk_iptw <- bnx_full_risk_iptw %>% 
  mutate(`Dose adjustment after interruption`= 
           recode(`Dose adjustment after interruption`,
                  "   Dose reduction to >=50%/initial dose* vs no change"=
                    "   Dose reduction to >=50%/initial dose* vs no change "))



p_dis_bnx <- forest(
  bnx_full_risk_iptw[, c(9,10,8,7,11,12)],
  est   = list(bnx_full_risk_iptw$HR, bnx_risk_iptw_dis$RD*100),
  lower = list(bnx_full_iptw_dis$lower, bnx_risk_iptw_dis$RD_lower*100),
  upper = list(bnx_full_iptw_dis$upper, bnx_risk_iptw_dis$RD_upper*100),
  ci_column = c(4, 6),
  ref_line = c(1,0),
  sizes = .6,
  x_trans = c("log", "none"),
  xlim = list(c(0.35, 4), c(-17, 17)),
  ticks_at = list(c(0.5, 1, 1.5, 3), c(-10, -5, 0, 5, 10)),
  xlab = NULL,
  arrow_lab = c("Beneficial", "Harmful"),
  title = "B. Buprenorphine/naloxone (BNX) trials",
  footnote = NULL,
  nudge_y = 0.02,
  theme =
    forest_theme(
      base_size = 14,
      ci_pch = 15,
      ci_lwd = 1.7,
      refline_gp = gpar(lwd = 1, lty = "dashed", col = "red"),
      footnote_gp = gpar(fontsize = 11, fontfamily = "", cex = 0.8,
                         fontface="plain", col = "black"),
      arrow_length = 0.05,
      arrow_gp = gpar(fontsize = 11, fontfamily = "", lwd = 0.7),
      xlab_adjust = c("refline", "center")
    ))


txt <- "<sup>*</sup> <=40mg BNX for Trial 5"

p_dis_bnx1 <- add_grob(p_dis_bnx, 
                       row = 18, 
                       col = 1:2,
                       order = "background",
                       gb_fn = gridtext::richtext_grob,
                       text = txt,
                       gp = gpar(fontsize = 11),
                       hjust = 0, vjust = 1, halign = 0, valign = 1,
                       x = unit(0, "npc"), y = unit(1, "npc"))


name = "Figure4_dis_risk_Panel_B_Nov12"
figure=paste0(folder, name, ".png")
ppi = 300
png(figure, width=4300, height=1900, res=300)
p_dis_bnx1
dev.off()



