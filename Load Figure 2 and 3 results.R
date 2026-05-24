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


## Read res Nov 12 #############
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
