
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


comparison_colors_v1 = c(
  "Dose reduction\nvs\nno change"              = "#000000",  # Blue
  "Dose reduction of <50%\nvs\nno change"            = "#2CA02C",  # Green
  "Dose reduction of >=50%/initial dose*\nvs\nno change"     = "#FF7F0E",  # Sky Blue
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


### od #########

res_od <- bind_rows(res1_od, res2_od, res3_od, res4_od, res5_od) %>% 
  filter(ret_length <= 365) %>% 
  mutate(contrast = recode(contrast,
                           "Dose reduction to <50%\nvs\nno change"=
                             "Dose reduction of <50%\nvs\nno change",
                           "Dose reduction to >=50%/initial dose*\nvs\nno change"=
                             "Dose reduction of >=50%/initial dose*\nvs\nno change"))

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



### dis #########
res_dis <- bind_rows(res1_dis, res2_dis, res3_dis, res4_dis, res5_dis) %>% 
  filter(ret_length <= 365) %>% 
  mutate(contrast = recode(contrast,
                           "Dose reduction to <50%\nvs\nno change"=
                             "Dose reduction of <50%\nvs\nno change",
                           "Dose reduction to >=50%/initial dose*\nvs\nno change"=
                             "Dose reduction of >=50%/initial dose*\nvs\nno change"))



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
  draw_label("A. Methadone (MET) trials among prevalent new users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
name = "Figure3_Panel_A_Jan_27_MK"
figure=paste0(folder, name, ".jpg")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
jpeg(figure, width=4000, height=4800, res=300)
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
  draw_label("B. Buprenorphine/naloxone (BNX) trials among prevalent new users", 
             x = 0.02, y = 0.98, hjust = 0, vjust = 0, 
             size = 18, color = "black") 


folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/"
name = "Figure3_Panel_B_Jan_27_MK"
figure=paste0(folder, name, ".jpg")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
jpeg(figure, width=4000, height=3400, res=300)
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
name = "Figure3_Panel_B_Nov12_v1_MK"
figure=paste0(folder, name, ".jpg")
ppi = 300
jpeg(figure, width=4800, height=4200, res=300)
plot_bnx2
dev.off()

