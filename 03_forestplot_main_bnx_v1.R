library(tidyverse)
library(forestploter)
library(grid)
library(gridtext)


# 1. PATHS AND DATA LOADING ##########

RES_DIR <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/BNX/"

FIG_DIR <- paste0(
  RES_DIR,
  "Figures/"
)

dir.create(
  FIG_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# Risk-result objects
risk_tr4 <- readRDS(
  paste0(
    RES_DIR,
    "Risk/risk_tr4_rmst.rds"
  )
)

risk_tr5 <- readRDS(
  paste0(
    RES_DIR,
    "Risk/risk_tr5_rmst.rds"
  )
)


# Propensity-score weighted analytic datasets
dat_tr4_all <- readRDS(
  paste0(
    RES_DIR,
    "PScore/ps_weight_tr4_all.rds"
  )
)

dat_tr5_all <- readRDS(
  paste0(
    RES_DIR,
    "PScore/ps_weight_tr5_all.rds"
  )
)


# 2. SHARED TREATMENT COLORS ###########

trt_colors <- c(
  "No dose change"         = "#000000",
  "<50% reduction"         = "#2CA02C",
  "≥50% reduction/≤4mg"    = "#FF7F0E",
  "Dose increase"          = "#D62728"
)


# 3. SUMMARY STATISTICS WITH RAW EVENT COUNTS ##########

compute_summary_stats_timewin <- function(
    Dat,
    time_var = "time_od_acute",
    event_var = "sens_od_acute",
    weight_var = "wt_trunc",
    trt_var = "trt",
    max_time = 90
) {
  
  required_vars <- c(
    time_var,
    event_var,
    weight_var,
    trt_var
  )
  
  missing_vars <- setdiff(
    required_vars,
    names(Dat)
  )
  
  if (length(missing_vars) > 0) {
    stop(
      "Missing variable(s) in analytic dataset: ",
      paste(
        missing_vars,
        collapse = ", "
      )
    )
  }
  
  Dat %>%
    mutate(
      # Raw observed event occurring by the selected horizon.
      #
      # An observation is counted as an event only if:
      #   1. the event indicator equals 1; and
      #   2. the recorded event time is at or before max_time.
      event_in_window = case_when(
        is.na(.data[[event_var]]) ~ 0,
        
        .data[[event_var]] == 1 &
          !is.na(.data[[time_var]]) &
          .data[[time_var]] <= max_time ~ 1,
        
        TRUE ~ 0
      )
    ) %>%
    group_by(
      .data[[trt_var]]
    ) %>%
    summarise(
      # Raw number of episodes at reinitiation
      Episode = n(),
      
      # Raw, unweighted number of observed events by the horizon
      Event = sum(
        event_in_window,
        na.rm = TRUE
      ),
      
      # Retained for internal checking, but not displayed
      Weighted_Episode = sum(
        .data[[weight_var]],
        na.rm = TRUE
      ),
      
      Weighted_event = sum(
        .data[[weight_var]] * event_in_window,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    rename(
      trt_group = 1
    ) %>%
    select(
      trt_group,
      Episode,
      Event,
      Weighted_Episode,
      Weighted_event
    )
}


# 4. EXTRACT RESULTS #######

extract_rd_data <- function(
    risk_obj,
    dat_weights,
    max_time,
    outcome = "od",
    ref_name = "no reversion",
    int_label
) {
  
  time_v <- if (
    outcome == "dis"
  ) {
    "time_discont"
  } else {
    "time_od_acute"
  }
  
  event_v <- if (
    outcome == "dis"
  ) {
    "sens_discont"
  } else {
    "sens_od_acute"
  }
  
  
  
  counts <- compute_summary_stats_timewin(
    Dat = dat_weights,
    time_var = time_v,
    event_var = event_v,
    weight_var = "wt_trunc",
    trt_var = "trt",
    max_time = max_time
  ) %>%
    select(
      trt_group,
      Episode,
      Event
    )
  
  
  
  risks <- risk_obj$risks %>%
    filter(
      horizon == max_time
    ) %>%
    mutate(
      risk_pct =
        risk * 100,
      
      risk_lower_pct =
        risk_lower * 100,
      
      risk_upper_pct =
        risk_upper * 100,
      
      abs_risk_disp = ifelse(
        abs(risk_pct) >= 10,
        
        sprintf(
          "%.1f (%.1f, %.1f)",
          risk_pct,
          risk_lower_pct,
          risk_upper_pct
        ),
        
        sprintf(
          "%.2f (%.2f, %.2f)",
          risk_pct,
          risk_lower_pct,
          risk_upper_pct
        )
      )
    ) %>%
    select(
      trt_group = trt,
      abs_risk_disp
    )
  
  
  
  contrasts <- risk_obj$contrasts %>%
    filter(
      horizon == max_time
    ) %>%
    mutate(
      rd_pct =
        RD * 100,
      
      rd_lower_pct =
        RD_lower * 100,
      
      rd_upper_pct =
        RD_upper * 100,
      
      rd_disp = ifelse(
        abs(rd_pct) >= 10,
        
        sprintf(
          "%.1f (%.1f, %.1f)",
          rd_pct,
          rd_lower_pct,
          rd_upper_pct
        ),
        
        sprintf(
          "%.2f (%.2f, %.2f)",
          rd_pct,
          rd_lower_pct,
          rd_upper_pct
        )
      )
    ) %>%
    select(
      trt_group = trt,
      rd_pct,
      rd_lower_pct,
      rd_upper_pct,
      rd_disp
    )
  
  
  
  counts %>%
    left_join(
      risks,
      by = "trt_group"
    ) %>%
    left_join(
      contrasts,
      by = "trt_group"
    ) %>%
    mutate(
      rd_disp = if_else(
        trt_group == ref_name,
        "Reference",
        rd_disp
      ),
      
      rd_pct = if_else(
        trt_group == ref_name,
        0,
        rd_pct
      ),
      
      rd_lower_pct = if_else(
        trt_group == ref_name,
        0,
        rd_lower_pct
      ),
      
      rd_upper_pct = if_else(
        trt_group == ref_name,
        0,
        rd_upper_pct
      ),
      
      Interruption =
        int_label
    )
}


# 5. FIGURE 2: OVERDOSE RISK DIFFERENCE AT 30 DAYS #######
risk_at_od <- 30


# Trial 4: interruption of 3-5 days
tr4_od_rd <- extract_rd_data(
  risk_obj = risk_tr4$all$od$ci,
  dat_weights = dat_tr4_all$data_weights,
  max_time = risk_at_od,
  outcome = "od",
  int_label = "3-5 days"
)


# Trial 5: interruption of 6-14 days
tr5_od_rd <- extract_rd_data(
  risk_obj = risk_tr5$all$od$ci,
  dat_weights = dat_tr5_all$data_weights,
  max_time = risk_at_od,
  outcome = "od",
  int_label = "6-14 days"
)


# Prepare overdose forest-plot table ########

plot_data_4a <- bind_rows(
  tr4_od_rd,
  tr5_od_rd
) %>%
  mutate(
    trt_group = recode(
      trt_group,
      "no reversion" = "No dose change",
      "1-<50%" = "<50% reduction",
      ">=50%/<=4mg" = "≥50% reduction/≤4mg",
      "increase" = "Dose increase"
    ),
    
    rd_pct = ifelse(
      trt_group == "No dose change",
      NA_real_,
      rd_pct
    ),
    
    rd_lower_pct = ifelse(
      trt_group == "No dose change",
      NA_real_,
      rd_lower_pct
    ),
    
    rd_upper_pct = ifelse(
      trt_group == "No dose change",
      NA_real_,
      rd_upper_pct
    ),
    
    Interruption = ifelse(
      duplicated(Interruption),
      "",
      Interruption
    ),
    
    # Display 1-9 events as "<10"; retain 0 and exact counts >=10
    Event_display = case_when(
      is.na(Event)           ~ "",
      Event > 0 & Event < 10 ~ "<10",
      TRUE                   ~ scales::comma(Event, accuracy = 1)
    ),
    
    `% Difference in 1-month Risk (95% CI)` =
      paste(rep(" ", 16), collapse = " ")
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = trt_group,
    `Episodes at risk\nat reinitiation, n` = Episode,
    `Events by\n1 month, n` = Event_display,
    `% risk at 1 month\n(95% CI)` = abs_risk_disp,
    ` ` = rd_disp
  )

# Create overdose forest plot ########

p_4a <- forest(
  plot_data_4a[
    ,
    c(
      "Interruption",
      "Dose adjustment\nafter interruption",
      "Episodes at risk\nat reinitiation, n",
      "Events by\n1 month, n",
      "% risk at 1 month\n(95% CI)",
      " ",
      "% Difference in 1-month Risk (95% CI)"
    )
  ],
  
  est =
    plot_data_4a$rd_pct,
  
  lower =
    plot_data_4a$rd_lower_pct,
  
  upper =
    plot_data_4a$rd_upper_pct,
  
  # The forest canvas is now the seventh column.
  ci_column = 7,
  
  ref_line = 0,
  
  x_trans = "none",
  
  xlim = c(
    -2.5,
    2.5
  ),
  
  ticks_at = c(
    -2,
    -1,
    0,
    1,
    2
  ),
  
  xlab =
    "Risk difference (%)",
  
  title =
    "",
  
  theme = forest_theme(
    base_size = 14,
    
    ci_pch = 15,
    
    ci_lwd = 1.7,
    
    ci_Theight = 0.2,
    
    refline_gp = gpar(
      lwd = 1,
      lty = "dashed",
      col = "red"
    ),
    
    xlab_adjust = c(
      "refline",
      "center"
    )
  )
)


# Align overdose headers and table contents ######

p_4a <- edit_plot(
  p_4a,
  col = 1:2,
  part = "header",
  hjust = 0,
  x = unit(
    0,
    "npc"
  )
)

p_4a <- edit_plot(
  p_4a,
  col = 1:2,
  part = "body",
  hjust = 0,
  x = unit(
    0,
    "npc"
  )
)

p_4a <- edit_plot(
  p_4a,
  col = 3:7,
  part = "header",
  hjust = 0.5,
  x = unit(
    0.5,
    "npc"
  )
)

p_4a <- edit_plot(
  p_4a,
  col = 3:6,
  part = "body",
  hjust = 0.5,
  x = unit(
    0.5,
    "npc"
  )
)

p_4a <- edit_plot(
  p_4a,
  col = 7,
  part = "header",
  hjust = 0.5,
  x = unit(
    0.15,
    "npc"
  )
)


# Apply dose-adjustment colors to overdose risk-difference estimates #######

ci_col_4a <- trt_colors[
  plot_data_4a$`Dose adjustment\nafter interruption`
]

for (
  i in seq_len(
    nrow(
      plot_data_4a
    )
  )
) {
  
  if (
    !is.na(
      plot_data_4a$rd_pct[i]
    )
  ) {
    
    current_color <- ci_col_4a[i]
    
    p_4a <- edit_plot(
      p_4a,
      
      row = i,
      
      col = 7,
      
      which = "ci",
      
      gp = gpar(
        col = current_color,
        fill = current_color,
        cex = 1.2
      )
    )
  }
}



p_4a_final <- add_grob(
  p_4a,
  
  row =
    nrow(plot_data_4a) + 2,
  
  col =
    1:2,
  
  order =
    "background",
  
  gb_fn =
    gridtext::richtext_grob,
  
  text =
    paste0(
      "Risk: cumulative incidence<br>",
      "CI: compatibility ('confidence') intervals"
    ),
  
  gp = gpar(
    fontsize = 11
  ),
  
  hjust = 0,
  
  vjust = 1,
  
  halign = 0,
  
  valign = 1,
  
  x = unit(
    0,
    "npc"
  ),
  
  y = unit(
    1,
    "npc"
  )
)


# Save overdose forest plot ##########

jpeg(
  filename = paste0(
    FIG_DIR,
    "Figure2_OD_RD_v1.jpg"
  ),
  
  width = 5500,
  
  height = 1700,
  
  res = 400,
  
  quality = 100,
  
  bg = "white"
)

print(
  p_4a_final
)

dev.off()


# 6. FIGURE 3: DISCONTINUATION RISK DIFFERENCE AT 365 DAYS #########

risk_at_dis <- 365


# Trial 4: interruption of 3-5 days
tr4_dis_rd <- extract_rd_data(
  risk_obj = risk_tr4$all$dis$ci,
  dat_weights = dat_tr4_all$data_weights,
  max_time = risk_at_dis,
  outcome = "dis",
  int_label = "3-5 days"
)


# Trial 5: interruption of 6-14 days
tr5_dis_rd <- extract_rd_data(
  risk_obj = risk_tr5$all$dis$ci,
  dat_weights = dat_tr5_all$data_weights,
  max_time = risk_at_dis,
  outcome = "dis",
  int_label = "6-14 days"
)


# Prepare discontinuation forest-plot table #######

plot_data_4b <- bind_rows(
  tr4_dis_rd,
  tr5_dis_rd
) %>%
  mutate(
    trt_group = recode(
      trt_group,
      "no reversion" = "No dose change",
      "1-<50%" = "<50% reduction",
      ">=50%/<=4mg" = "≥50% reduction/≤4mg",
      "increase" = "Dose increase"
    ),
    
    rd_pct = ifelse(
      trt_group == "No dose change",
      NA_real_,
      rd_pct
    ),
    
    rd_lower_pct = ifelse(
      trt_group == "No dose change",
      NA_real_,
      rd_lower_pct
    ),
    
    rd_upper_pct = ifelse(
      trt_group == "No dose change",
      NA_real_,
      rd_upper_pct
    ),
    
    Interruption = ifelse(
      duplicated(Interruption),
      "",
      Interruption
    ),
    
    # Display 1-9 events as "<10"; retain 0 and exact counts >=10
    Event_display = case_when(
      is.na(Event)           ~ "",
      Event > 0 & Event < 10 ~ "<10",
      TRUE                   ~ scales::comma(Event, accuracy = 1)
    ),
    
    `% Difference in 12-month Risk (95% CI)` =
      paste(rep(" ", 16), collapse = " ")
  ) %>%
  rename(
    `Dose adjustment\nafter interruption` = trt_group,
    `Episodes at risk\nat reinitiation, n` = Episode,
    `Events by\n12 months, n` = Event_display,
    `% risk at 12 months\n(95% CI)` = abs_risk_disp,
    ` ` = rd_disp
  )


# Create discontinuation forest plot ########

p_4b <- forest(
  plot_data_4b[
    ,
    c(
      "Interruption",
      "Dose adjustment\nafter interruption",
      "Episodes at risk\nat reinitiation, n",
      "Events by\n12 months, n",
      "% risk at 12 months\n(95% CI)",
      " ",
      "% Difference in 12-month Risk (95% CI)"
    )
  ],
  
  est =
    plot_data_4b$rd_pct,
  
  lower =
    plot_data_4b$rd_lower_pct,
  
  upper =
    plot_data_4b$rd_upper_pct,
  
  # Forest canvas is now the seventh column.
  ci_column = 7,
  
  ref_line = 0,
  
  x_trans = "none",
  
  xlim = c(
    -12.5,
    12.5
  ),
  
  ticks_at = c(
    -10,
    -5,
    0,
    5,
    10
  ),
  
  xlab =
    "Risk difference (%)",
  
  title =
    "",
  
  theme = forest_theme(
    base_size = 14,
    
    ci_pch = 15,
    
    ci_lwd = 1.7,
    
    ci_Theight = 0.2,
    
    refline_gp = gpar(
      lwd = 1,
      lty = "dashed",
      col = "red"
    ),
    
    xlab_adjust = c(
      "refline",
      "center"
    )
  )
)


# Align discontinuation headers and table contents ######

p_4b <- edit_plot(
  p_4b,
  col = 1:2,
  part = "header",
  hjust = 0,
  x = unit(
    0,
    "npc"
  )
)

p_4b <- edit_plot(
  p_4b,
  col = 1:2,
  part = "body",
  hjust = 0,
  x = unit(
    0,
    "npc"
  )
)

p_4b <- edit_plot(
  p_4b,
  col = 3:7,
  part = "header",
  hjust = 0.5,
  x = unit(
    0.5,
    "npc"
  )
)

p_4b <- edit_plot(
  p_4b,
  col = 3:6,
  part = "body",
  hjust = 0.5,
  x = unit(
    0.5,
    "npc"
  )
)

p_4b <- edit_plot(
  p_4b,
  col = 7,
  part = "header",
  hjust = 0.5,
  x = unit(
    0.15,
    "npc"
  )
)


# Apply dose-adjustment colors ########

ci_col_4b <- trt_colors[
  plot_data_4b$`Dose adjustment\nafter interruption`
]

for (
  i in seq_len(
    nrow(
      plot_data_4b
    )
  )
) {
  
  if (
    !is.na(
      plot_data_4b$rd_pct[i]
    )
  ) {
    
    current_color <- ci_col_4b[i]
    
    p_4b <- edit_plot(
      p_4b,
      
      row = i,
      
      col = 7,
      
      which = "ci",
      
      gp = gpar(
        col = current_color,
        fill = current_color,
        cex = 1.2
      )
    )
  }
}


# Add discontinuation footnote #######

p_4b_final <- add_grob(
  p_4b,
  
  row =
    nrow(plot_data_4b) + 2,
  
  col =
    1:2,
  
  order =
    "background",
  
  gb_fn =
    gridtext::richtext_grob,
  
  text =
    paste0(
      "Risk: cumulative incidence<br>",
      "CI: compatibility ('confidence') intervals"
    ),
  
  gp = gpar(
    fontsize = 11
  ),
  
  hjust = 0,
  
  vjust = 1,
  
  halign = 0,
  
  valign = 1,
  
  x = unit(
    0,
    "npc"
  ),
  
  y = unit(
    1,
    "npc"
  )
)


# Save discontinuation forest plot #####

jpeg(
  filename = paste0(
    FIG_DIR,
    "Figure3_Discont_RD_v1.jpg"
  ),
  
  width = 5700,
  
  height = 1700,
  
  res = 400,
  
  quality = 100,
  
  bg = "white"
)

print(
  p_4b_final
)

dev.off()

