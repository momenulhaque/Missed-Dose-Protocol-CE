
## histogram ###############
### Trial 1: 1-2 days ############
pres_1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr1_12d_ow.rds")
#### pres #############

pres_incid_1 = plot_iv_logistic(fit = pres_1$grp_red$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")


pres_rprev_1 = plot_iv_logistic(fit = pres_1$grp_red$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Prevalent new users")



pres_incid_2 = plot_iv_logistic(fit = pres_1$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


pres_rprev_2 = plot_iv_logistic(fit = pres_1$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_1 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr1_12d_ow.rds")

faci_incid_1 = plot_iv_logistic(fit = faci_1$grp_red$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)
faci_rprev_1 = plot_iv_logistic(fit = faci_1$grp_red$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_1$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


faci_rprev_2 = plot_iv_logistic(fit = faci_1$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))

row2 <- pres_incid_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))


row3 <- faci_incid_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))


row4 <- faci_incid_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))


plot_1 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer() + 
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.16, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_1s <- ggdraw(plot_1) + 
  draw_label("A. Prescriber's facicribing preference for dose reduction versus no dose change after methadone interruption of 1-2 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose increase versus no dose change after methadone interruption of 1-2 days",  
             x = 0.02, y = 0.73, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Facility-level prescribing preference for dose reduction versus no dose change after methadone interruption of 1-2 days", 
             x = 0.02, y = 0.48, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("D. Facility-level prescribing preference for dose increase versus no dose change after methadone interruption of 1-2 days",  
             x = 0.02, y = 0.24, hjust = 0, vjust = 0, 
             size = 18, color = "black")  




name = "histogram_1_12d_MK_Jan_27"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=5000, height=7000, res=300)
plot_1s
dev.off()


### Trial 2: 3-4 days ############
pres_2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr2_34d_ow.rds")
#### pres #############
pres_incid_0 = plot_iv_logistic(fit = pres_2$grp_lt50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_rprev_0 = plot_iv_logistic(fit = pres_2$grp_lt50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2), 
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Prevalent new users")

pres_incid_1 = plot_iv_logistic(fit = pres_2$grp_ge50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_rprev_1 = plot_iv_logistic(fit = pres_2$grp_ge50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



pres_incid_2 = plot_iv_logistic(fit = pres_2$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_2$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_2 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr2_34d_ow.rds")


faci_incid_0 = plot_iv_logistic(fit = faci_2$grp_lt50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_rprev_0 = plot_iv_logistic(fit = faci_2$grp_lt50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


faci_incid_1 = plot_iv_logistic(fit = faci_2$grp_ge50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_2$grp_ge50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_2$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_2$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_0 + plot_spacer() + pres_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1))

row2 <- pres_incid_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- pres_incid_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row4 <- faci_incid_0 + plot_spacer() + faci_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1))

row5 <- faci_incid_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))

row6 <- faci_incid_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_2 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer()/
  row5/
  plot_spacer()/
  row6/
  plot_spacer() +
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_2s <- ggdraw(plot_2) + 
  draw_label("A. Prescriber's facicribing preference for dose reductions of <50% versus no dose change after methadone interruption of 3-4 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose reductions of >=50%/initial dose* versus no dose change after methadone interruption of 3-4 days",  
             x = 0.02, y = 0.815, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Prescriber's facicribing preference for dose increase versus no dose change after methadone interruption of 3-4 days",  
             x = 0.02, y = 0.656, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("D. Facility-level prescribing preference for dose reductions of <50% versus no dose change after methadone interruption of 3-4 days", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("E. Facility-level prescribing preference for dose reductions of >=50%/initial dose* versus no dose change after methadone interruption of 3-4 days",  
             x = 0.02, y = 0.325, hjust = 0, vjust = 0, 
             size = 18, color = "black") + 
  draw_label("F. Facility-level prescribing preference for dose increase versus no dose change after methadone interruption of 3-4 days",  
             x = 0.02, y = 0.168, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("* <=40mg methadone", 
             x = 0.88, y = 0.003, hjust = 0, vjust = 0, 
             size = 15, color = "black")



name = "histogram_2_34d_MK_Jan27"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=5000, height=8500, res=300)
plot_2s
dev.off()




### Trial 3 ############
pres_3 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr3.rds")
#### pres #############
pres_incid_0 = plot_iv_logistic(fit = pres_3$grp_lt50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_rprev_0 = plot_iv_logistic(fit = pres_3$grp_lt50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2), 
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Prevalent new users")

pres_incid_1 = plot_iv_logistic(fit = pres_3$grp_ge50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_rprev_1 = plot_iv_logistic(fit = pres_3$grp_ge50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



pres_incid_2 = plot_iv_logistic(fit = pres_3$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_3$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_3 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr3.rds")


faci_incid_0 = plot_iv_logistic(fit = faci_3$grp_lt50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_rprev_0 = plot_iv_logistic(fit = faci_3$grp_lt50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


faci_incid_1 = plot_iv_logistic(fit = faci_3$grp_ge50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_3$grp_ge50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_3$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_3$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_0 + plot_spacer() + pres_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1))

row2 <- pres_incid_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- pres_incid_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row4 <- faci_incid_0 + plot_spacer() + faci_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1))

row5 <- faci_incid_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))

row6 <- faci_incid_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_3 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer()/
  row5/
  plot_spacer()/
  row6/
  plot_spacer() +
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.04),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_3s <- ggdraw(plot_3) + 
  draw_label("A. Prescriber's facicribing preference for dose reductions of <50% versus no dose change after methadone interruption of 5-14 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose reductions of >=50%/initial dose* versus no dose change after methadone interruption of 5-14 days",  
             x = 0.02, y = 0.815, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Prescriber's facicribing preference for dose increase versus no dose change after methadone interruption of 5-14 days",  
             x = 0.02, y = 0.652, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("D. Facility-level prescribing preference for dose reductions of <50% versus no dose change after methadone interruption of 5-14 days", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("E. Facility-level prescribing preference for dose reductions of >=50%/initial dose* versus no dose change after methadone interruption of 5-14 days",  
             x = 0.02, y = 0.325, hjust = 0, vjust = 0, 
             size = 18, color = "black") + 
  draw_label("F. Facility-level prescribing preference for dose increase versus no dose change after methadone interruption of 5-14 days",  
             x = 0.02, y = 0.165, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("* <=40mg methadone", 
             x = 0.88, y = 0.003, hjust = 0, vjust = 0, 
             size = 15, color = "black")



name = "histogram_3_MK_Jan27"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=5000, height=8500, res=300)
plot_3s
dev.off()



### Trial 4 ############
pres_4 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr4.rds")
#### pres #############

pres_incid_1 = plot_iv_logistic(fit = pres_4$grp_red$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_prev_1 = plot_iv_logistic(fit = pres_4$grp_red$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1, 
                               plot_title = "Prevalent new users")

pres_rprev_1 = plot_iv_logistic(fit = pres_4$grp_red$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "All prevalent users")



pres_incid_2 = plot_iv_logistic(fit = pres_4$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_2 = plot_iv_logistic(fit = pres_4$grp_incr$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_4$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_4 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr4.rds")

faci_incid_1 = plot_iv_logistic(fit = faci_4$grp_red$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_1 = plot_iv_logistic(fit = faci_4$grp_red$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_4$grp_red$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose reduction",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_4$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_2 = plot_iv_logistic(fit = faci_4$grp_incr$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_4$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose increase",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))

row2 <- pres_incid_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))


row3 <- faci_incid_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))


row4 <- faci_incid_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))


plot_4 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer() + 
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.02),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_4s <- ggdraw(plot_4) + 
  draw_label("A. Prescriber's facicribing preference for dose reduction versus no dose change after buprenorphine/naloxone interruption of 1-5 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose increase versus no dose change after buprenorphine/naloxone interruption of 1-5 days",  
             x = 0.02, y = 0.73, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Facility-level prescribing preference for dose reduction versus no dose change after buprenorphine/naloxone interruption of 1-5 days", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("D. Facility-level prescribing preference for dose increase versus no dose change after buprenorphine/naloxone interruption of 1-5 days",  
             x = 0.02, y = 0.24, hjust = 0, vjust = 0, 
             size = 18, color = "black")  


name = "histogram_4_MK_Jan27"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=5000, height=7000, res=300)
plot_4s
dev.off()



### Trial 5 ############
pres_5 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/pred_p_tr5.rds")
#### pres #############
pres_incid_0 = plot_iv_logistic(fit = pres_5$grp_lt50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "Incident users")

pres_prev_0 = plot_iv_logistic(fit = pres_5$grp_lt50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1, 
                               plot_title = "Prevalent new users")

pres_rprev_0 = plot_iv_logistic(fit = pres_5$grp_lt50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1, 
                                plot_title = "All prevalent users")

pres_incid_1 = plot_iv_logistic(fit = pres_5$grp_ge50$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_1 = plot_iv_logistic(fit = pres_5$grp_ge50$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_1 = plot_iv_logistic(fit = pres_5$grp_ge50$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



pres_incid_2 = plot_iv_logistic(fit = pres_5$grp_incr$incid,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

pres_prev_2 = plot_iv_logistic(fit = pres_5$grp_incr$prev,  
                               xlab_text = "Provider's prescribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

pres_rprev_2 = plot_iv_logistic(fit = pres_5$grp_incr$rzero,  
                                xlab_text = "Provider's prescribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


#### faci #############
faci_5 <- readRDS("H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/SA Results/sens_IV/faci_pred_p_tr5.rds")


faci_incid_0 = plot_iv_logistic(fit = faci_5$grp_lt50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_0 = plot_iv_logistic(fit = faci_5$grp_lt50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_0 = plot_iv_logistic(fit = faci_5$grp_lt50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)


faci_incid_1 = plot_iv_logistic(fit = faci_5$grp_ge50$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_1 = plot_iv_logistic(fit = faci_5$grp_ge50$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_1 = plot_iv_logistic(fit = faci_5$grp_ge50$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = "Probability of assignment to dose adjustment strategy",
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)



faci_incid_2 = plot_iv_logistic(fit = faci_5$grp_incr$incid,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = "Proportion of sample",
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

faci_prev_2 = plot_iv_logistic(fit = faci_5$grp_incr$prev,  
                               xlab_text = "Facility-level facicribing preference (%)",
                               ylab_text = NULL,
                               sec_y_text = NULL,
                               y_points1 = seq(0, 1, 0.2),
                               y_points2 = seq(0, 1, 0.2),
                               y_min = 0, y_max = 1)

faci_rprev_2 = plot_iv_logistic(fit = faci_5$grp_incr$rzero,  
                                xlab_text = "Facility-level facicribing preference (%)",
                                ylab_text = NULL,
                                sec_y_text = NULL,
                                y_points1 = seq(0, 1, 0.2),
                                y_points2 = seq(0, 1, 0.2),
                                y_min = 0, y_max = 1)

#### combine #############
library(patchwork)
library(cowplot)

row1 <- pres_incid_0 + plot_spacer() + pres_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1))

row2 <- pres_incid_1 + plot_spacer() + pres_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))

row3 <- pres_incid_2 + plot_spacer() + pres_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

row4 <- faci_incid_0 + plot_spacer() + faci_rprev_0 +
  plot_layout(widths = c(1, 0.07, 1))

row5 <- faci_incid_1 + plot_spacer() + faci_rprev_1 +
  plot_layout(widths = c(1, 0.07, 1))

row6 <- faci_incid_2 + plot_spacer() + faci_rprev_2 +
  plot_layout(widths = c(1, 0.07, 1))

plot_5 = plot_spacer() / 
  row1/
  plot_spacer() /
  row2/
  plot_spacer() /
  row3/
  plot_spacer()/
  row4/
  plot_spacer()/
  row5/
  plot_spacer()/
  row6/
  plot_spacer() +
  plot_layout(
    heights = c(0.06, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.11, 1, 0.04),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.box.just = "left"
  ) &
  guides(color = guide_legend(nrow = 2))


plot_5s <- ggdraw(plot_5) + 
  draw_label("A. Prescriber's facicribing preference for dose reductions of <50% versus no dose change after buprenorphine/naloxone interruption of 6-14 days", 
             x = 0.02, y = 0.986, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("B. Prescriber's facicribing preference for dose reductions of >=50%/initial dose* versus no dose change after buprenorphine/naloxone interruption of 6-14 days",  
             x = 0.02, y = 0.815, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("C. Prescriber's facicribing preference for dose increase versus no dose change after buprenorphine/naloxone interruption of 6-14 days",  
             x = 0.02, y = 0.656, hjust = 0, vjust = 0, 
             size = 18, color = "black")  +
  draw_label("D. Facility-level prescribing preference for dose reductions of <50% versus no dose change after buprenorphine/naloxone interruption of 6-14 days", 
             x = 0.02, y = 0.49, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("E. Facility-level prescribing preference for dose reductions of >=50%/initial dose* versus no dose change after buprenorphine/naloxone interruption of 6-14 days",  
             x = 0.02, y = 0.325, hjust = 0, vjust = 0, 
             size = 18, color = "black") + 
  draw_label("F. Facility-level prescribing preference for dose increase versus no dose change after buprenorphine/naloxone interruption of 6-14 days",  
             x = 0.02, y = 0.168, hjust = 0, vjust = 0, 
             size = 18, color = "black") +
  draw_label("* <=4/1mg buprenorphine/naloxone", 
             x = 0.80, y = 0.003, hjust = 0, vjust = 0, 
             size = 15, color = "black")



name = "histogram_5_MK_Jan27"
figure=paste0("H:\\SAS-Work\\Missed doses Paper v2\\Result_Robust_PLR\\SA Results\\sens_IV\\plots\\",name,".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=5000, height=8500, res=300)
plot_5s
dev.off()

