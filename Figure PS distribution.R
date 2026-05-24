

comparison_colors_v3 <- c(
  "No dose change"              = "darkblue", 
  "Dose reduction"              = "#000000",        
  "Dose increase"               = "#D62728",  
  "1-<50% reduction"            = "#2CA02C",  
  ">=50% reduction/initial dose*"  = "#FF7F0E"  
)


ps_score_plot_c <- function(dat,
                            title_name = "Incident users",
                            x_lab = "Propensity Score",
                            y_lab = "Density",
                            legend_option = FALSE,
                            x_points = c(0, .25, .5, .75, 1)){
  ggplot(dat, aes(x = ps, colour = group, fill = group)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.2, bins = 30) +
    geom_density(alpha = 0.05, size = 0.5) +
    scale_colour_manual(values = comparison_colors_v3) +
    scale_fill_manual(values = comparison_colors_v3) +
    facet_wrap(~Trial+Design, nrow = 5, scales = "free") +
    # facet_grid(
    #   rows = vars(Trial),
    #   cols = vars(Design),
    #   scales = "free"
    # ) +
    labs(
      title = NULL,
      x = x_lab,
      y = y_lab
    ) +
    theme_minimal(base_size = 13) +
    scale_x_continuous(breaks = x_points) +
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
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

ps_score_plot_c3 <- function(dat,
                             title_name = "Incident users",
                             x_lab = "Propensity Score",
                             y_lab = "Density",
                             legend_position = c(0.02, 0.02),   # bottom-left inside plot
                             legend_justification = c(0, 0),
                             x_points = c(0, .25, .5, .75, 1)
) {
  
  ggplot(dat, aes(x = ps, colour = group, fill = group)) +
    geom_histogram(aes(y = ..density..),
                   position = "identity",
                   alpha = 0.2,
                   bins = 30) +
    geom_density(alpha = 0.05, size = 0.5) +
    scale_colour_manual(values = comparison_colors_v3) +
    scale_fill_manual(values = comparison_colors_v3) +
    facet_wrap(~Trial_design, nrow = 5, scales = "free") +
    labs(
      title = NULL,
      x = x_lab,
      y = y_lab
    ) +
    theme_minimal(base_size = 13) +
    scale_x_continuous(breaks = x_points) +
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
      axis.ticks.x = element_line(linewidth = 1.2, color = "black"),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      axis.title.y = element_text(size = 14),
      
      # legend: bottom-left (inside plot area)
      legend.position = legend_position,
      legend.justification = legend_justification,
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text = element_text(face = "bold", size = 13)
    )
}




ps_score_plot <- function(dat,
                          title_name = "Incident users",
                          x_lab = "Propensity Score",
                          y_lab = "Density",
                          legend_option = FALSE,
                          x_points = c(0, .25, .5, .75, 1)){
  ggplot(dat, aes(x = ps, colour = group, fill = group)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.2, bins = 30) +
    geom_density(alpha = 0.05, size = 0.5) +
    scale_colour_manual(values = comparison_colors_v3) +
    scale_fill_manual(values = comparison_colors_v3) +
    facet_wrap(~Design, nrow = 1, scales = "free") +
    labs(
      title = title_name,
      x = x_lab,
      y = y_lab
    ) +
    theme_minimal(base_size = 13) +
    scale_x_continuous(breaks = x_points) +
    theme(
      title = element_text(color = "black", size = 14),
      panel.grid = element_blank(),
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


# PS histogram #############
ps_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust/"

## Trial 1 #################
tr1 <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_ow.rds"))$data_weights
tr1_prev <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_prev_ow.rds"))$data_weights
tr1_rprev <- readRDS(file = paste0(ps_folder, "ps_weight_tr1_12d_rprev_ow.rds"))$data_weights

tr1 <- tr1 %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "reversion" = "Dose reduction",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "Dose reduction", "Dose increase")),
         Trial = "Trial 1: MET interruption of 1-2 days",
         Design = "Incident users"
  ) %>% select(group, ps, Trial, Design)

tr1_prev <- tr1_prev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "reversion" = "Dose reduction",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "Dose reduction", "Dose increase")),
         Trial = "Trial 1: MET interruption of 1-2 days",
         Design = "Prevalent new users"
  ) %>% select(group, ps, Trial, Design)

tr1_rprev <- tr1_rprev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "reversion" = "Dose reduction",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "Dose reduction", "Dose increase")),
         Trial = "Trial 1: MET interruption of 1-2 days",
         Design = "All prevalent users"
  ) %>% select(group, ps, Trial, Design)

dat_tr1 = bind_rows(tr1, tr1_prev, tr1_rprev) %>% mutate(
  Design = factor(Design, levels = c("Incident users", "Prevalent new users", "All prevalent users"))
)


## Trial 2 #######
tr2 <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_ow.rds"))$data_weights
tr2_prev <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_prev_ow.rds"))$data_weights
tr2_rprev <- readRDS(file = paste0(ps_folder, "ps_weight_tr2_34d_rprev_ow.rds"))$data_weights



tr2 <- tr2 %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "1-<50%" = "1-<50% reduction",
                  ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "1-<50% reduction",
                                          ">=50% reduction/initial dose*", "Dose increase")),
         Trial = "Trial 2: MET interruption of 3-4 days",
         Design = "Incident users"
  )  %>% select(group, ps, Trial, Design)

tr2_prev <- tr2_prev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "1-<50%" = "1-<50% reduction",
                  ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "1-<50% reduction",
                                          ">=50% reduction/initial dose*", "Dose increase")),
         Trial = "Trial 2: MET interruption of 3-4 days",
         Design = "Prevalent new users"
  )  %>% select(group, ps, Trial, Design)

tr2_rprev <- tr2_rprev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "1-<50%" = "1-<50% reduction",
                  ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "1-<50% reduction",
                                          ">=50% reduction/initial dose*", "Dose increase")),
         Trial = "Trial 2: MET interruption of 3-4 days",
         Design = "All prevalent users"
  )  %>% select(group, ps, Trial, Design)

dat_tr2 = bind_rows(tr2, tr2_prev, tr2_rprev) %>% mutate(
  Design = factor(Design, levels = c("Incident users", "Prevalent new users", "All prevalent users"))
)

## Trial 3 ########
tr3 <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_ow.rds"))$data_weights
tr3_prev <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_prev_ow.rds"))$data_weights
tr3_rprev <- readRDS(file = paste0(ps_folder, "ps_weight_tr3_rprev_ow.rds"))$data_weights

tr3 <- tr3 %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "1-<50%" = "1-<50% reduction",
                  ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "1-<50% reduction",
                                          ">=50% reduction/initial dose*", "Dose increase")),
         Trial = "Trial 3: MET interruption of 5-14 days",
         Design = "Incident users"
  )  %>% select(group, ps, Trial, Design)

tr3_prev <- tr3_prev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "1-<50%" = "1-<50% reduction",
                  ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "1-<50% reduction",
                                          ">=50% reduction/initial dose*", "Dose increase")),
         Trial = "Trial 3: MET interruption of 5-14 days",
         Design = "Prevalent new users"
  ) %>% select(group, ps, Trial, Design)

tr3_rprev <- tr3_rprev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "1-<50%" = "1-<50% reduction",
                  ">=50%/<=40mg" = ">=50% reduction/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "1-<50% reduction",
                                          ">=50% reduction/initial dose*", "Dose increase")),
         Trial = "Trial 3: MET interruption of 5-14 days",
         Design = "All prevalent users"
  ) %>% select(group, ps, Trial, Design)

dat_tr3 = bind_rows(tr3, tr3_prev, tr3_rprev) %>% mutate(
  Design = factor(Design, levels = c("Incident users", "Prevalent new users", "All prevalent users"))
)


## Trial 4 ####################
tr4 <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_ow.rds"))$data_weights
tr4_prev <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_prev_ow.rds"))$data_weights
tr4_rprev <- readRDS(file = paste0(ps_folder, "ps_weight_tr4_rprev_ow.rds"))$data_weights

tr4 <- tr4 %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "reversion" = "Dose reduction",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "Dose reduction", "Dose increase")),
         Trial = "Trial 4: BNX interruption of 1-5 days",
         Design = "Incident users"
  ) %>% select(group, ps, Trial, Design)

tr4_prev <- tr4_prev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "reversion" = "Dose reduction",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "Dose reduction", "Dose increase")),
         Trial = "Trial 4: BNX interruption of 1-5 days",
         Design = "Prevalent new users"
  ) %>% select(group, ps, Trial, Design)

tr4_rprev <- tr4_rprev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "reversion" = "Dose reduction",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "Dose reduction", "Dose increase")),
         Trial = "Trial 4: BNX interruption of 1-5 days",
         Design = "All prevalent users"
  ) %>% select(group, ps, Trial, Design)

dat_tr4 = bind_rows(tr4, tr4_prev, tr4_rprev) %>% mutate(
  Design = factor(Design, levels = c("Incident users", "Prevalent new users", "All prevalent users"))
)

## Trial 5 ####################
tr5 <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_ow.rds"))$data_weights
tr5_prev <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_prev_ow.rds"))$data_weights
tr5_rprev <- readRDS(file = paste0(ps_folder, "ps_weight_tr5_rprev_ow.rds"))$data_weights


tr5 <- tr5 %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "1-<50%" = "1-<50% reduction",
                  ">=50%/<=4mg" = ">=50% reduction/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "1-<50% reduction",
                                          ">=50% reduction/initial dose*", "Dose increase")),
         Trial = "Trial 5: MET interruption of 6-14 days",
         Design = "Incident users"
  ) %>% select(group, ps, Trial, Design)
tr5_prev <- tr5_prev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "1-<50%" = "1-<50% reduction",
                  ">=50%/<=4mg" = ">=50% reduction/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "1-<50% reduction",
                                          ">=50% reduction/initial dose*", "Dose increase")),
         Trial = "Trial 5: BNX interruption of 6-14 days",
         Design = "Prevalent new users"
  ) %>% select(group, ps, Trial, Design)

tr5_rprev <- tr5_rprev %>% 
  mutate(group =
           recode(trt,
                  "no reversion" = "No dose change",
                  "1-<50%" = "1-<50% reduction",
                  ">=50%/<=4mg" = ">=50% reduction/initial dose*",
                  "increase" = "Dose increase"),
         group = factor(group, levels = c("No dose change", "1-<50% reduction",
                                          ">=50% reduction/initial dose*", "Dose increase")),
         Trial = "Trial 5: BNX interruption of 6-14 days",
         Design = "All prevalent users"
  ) %>% select(group, ps, Trial, Design)


dat_tr5 = bind_rows(tr5, tr5_prev, tr5_rprev) %>% mutate(
  Design = factor(Design, levels = c("Incident users", "Prevalent new users", "All prevalent users"))
)



## Combined #########

levels(dat_tr5$Design)
dat_all <- bind_rows(dat_tr1, dat_tr2, dat_tr3, dat_tr4, dat_tr5)

dat_all <- dat_all %>% mutate(
  Trial_design = paste0(Trial, "\n", Design)
) %>% 
  filter(!(Trial_design %in% c("Trial 1: MET interruption of 1-2 days\nPrevalent new users",
                               "Trial 2: MET interruption of 3-4 days\nPrevalent new users",
                               "Trial 3: MET interruption of 5-14 days\nPrevalent new users",
                               "Trial 4: BNX interruption of 1-5 days\nPrevalent new users",
                               "Trial 5: BNX interruption of 6-14 days\nPrevalent new users"))) %>% 
  mutate(Trial_design = recode(Trial_design, 
                               "Trial 1: MET interruption of 1-2 days\nAll prevalent users"=
                                 "Trial 1: MET interruption of 1-2 days\nPrevalent new users", 
                               "Trial 2: MET interruption of 3-4 days\nAll prevalent users"=
                                 "Trial 2: MET interruption of 3-4 days\nPrevalent new users",
                               "Trial 3: MET interruption of 5-14 days\nAll prevalent users"=
                                 "Trial 3: MET interruption of 5-14 days\nPrevalent new users",
                               "Trial 4: BNX interruption of 1-5 days\nAll prevalent users"=
                                 "Trial 4: BNX interruption of 1-5 days\nPrevalent new users",
                               "Trial 5: BNX interruption of 6-14 days\nAll prevalent users"=
                                 "Trial 5: BNX interruption of 6-14 days\nPrevalent new users"),
         Trial_design = factor(Trial_design, levels = c("Trial 1: MET interruption of 1-2 days\nIncident users",
                                                        "Trial 1: MET interruption of 1-2 days\nPrevalent new users",
                                                        
                                                        "Trial 2: MET interruption of 3-4 days\nIncident users",
                                                        "Trial 2: MET interruption of 3-4 days\nPrevalent new users",
                                                        
                                                        "Trial 3: MET interruption of 5-14 days\nIncident users",
                                                        "Trial 3: MET interruption of 5-14 days\nPrevalent new users",
                                                        
                                                        "Trial 4: BNX interruption of 1-5 days\nIncident users",
                                                        "Trial 4: BNX interruption of 1-5 days\nPrevalent new users",
                                                       
                                                        "Trial 5: MET interruption of 6-14 days\nIncident users",
                                                        "Trial 5: BNX interruption of 6-14 days\nPrevalent new users"
                                                        )))


library(patchwork)
levels(dat_all$Trial_design)

plot2 = ps_score_plot_c3(dat=dat_all,
                         #x_lab = NULL,
                         #y_lab = NULL,
                         legend_position = "bottom") + 
  guides(color = guide_legend(nrow = 2))


library(cowplot)
plot2_dis1 <- ggdraw(plot2) +
  draw_label(
    "* <=40mg methadone for Trial 2, 3, and\n <=4/1mg buprenorphine/naloxone for Trial 5",
    x = 0.99, y = 0.01, hjust = 1, vjust = 0,
    size = 12, color = "black"
  )



sa_folder = "H:/SAS-Work/Missed doses Paper v2/Result_Robust_PLR/prior_dose/"
name = "ps_dist_plot_all_Jan27"
figure=paste0(sa_folder, name, ".png")
ppi = 300
# png(figure,width=9.61,height=3.95, units = "in", res=ppi)
png(figure, width=4400, height=7500, res=300)
plot2_dis1
dev.off()

