# Overall: MET ##########

library(tidyverse)
library(haven)
library(survival)
library(patchwork)

# ── 1. Paths & Setup ────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/LE30D/"

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
TRT_LABELS <- c("No dose change", "<50% reduction", "\u226550% reduction/\u226440mg", "Dose increase")

PLOT_COLORS <- c("No dose change" = "#333333", "<50% reduction" = "#4DAF4A", 
                 "\u226550% reduction/\u226440mg" = "#FF7F00", "Dose increase" = "#E41A1C")

# ── 2. Load and Clean Raw Data (Any Episode) ────────────────────────────────
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))


dat_clean <- raw_all |>
  mutate(
    trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)
  ) |>
  rename(trt = trt_s456) |>
  filter(trt %in% TRT_LEVELS) |>
  mutate(
    trt = factor(trt, levels = TRT_LEVELS),
    # Create the retention stratification variable
    ret_strata = if_else(ret_length <= 30, "Interruption within first \u226430 days", "Interruption after >30 days"),
    ret_strata = factor(ret_strata, levels = c("Interruption within first \u226430 days", "Interruption after >30 days"))
  )

# Split into Trials
dat_tr2 <- dat_clean |> filter(trial_gap %in% 3:4)
dat_tr3 <- dat_clean |> filter(trial_cat == "5-14")

rm(raw_all, dat_clean); gc()

# ── 3. Helper Function: Unadjusted RMST via Kaplan-Meier ────────────────────
get_unadjusted_rmst <- function(df, horizon = 365) {
  # Fit unadjusted Kaplan-Meier
  fit <- survfit(Surv(time_discont, sens_discont) ~ trt, data = df)
  
  # Extract Restricted Mean (rmean) at specified horizon
  sum_fit <- summary(fit, rmean = horizon)$table
  
  # Parse rownames and calculate 95% CIs (Removed the asterisks)
  tibble(
    trt_raw = rownames(sum_fit),
    RMST    = sum_fit[, "rmean"],      
    se      = sum_fit[, "se(rmean)"]    
  ) |>
    mutate(
      trt        = str_remove(trt_raw, "trt="),
      trt        = factor(trt, levels = TRT_LEVELS, labels = TRT_LABELS),
      RMST_lower = RMST - (1.96 * se),
      RMST_upper = RMST + (1.96 * se)
    ) |>
    select(trt, RMST, RMST_lower, RMST_upper)
}


# ── 4. Calculate Stratified RMST ────────────────────────────────────────────

# Trial 4 (3-4 days)
rmst_tr2_le30 <- get_unadjusted_rmst(dat_tr2 |> filter(ret_strata == "Interruption within first \u226430 days")) |> 
  mutate(panel = "Interruption of 3-4 days", ret_strata = "Interruption within first \u226430 days")

rmst_tr2_gt30 <- get_unadjusted_rmst(dat_tr2 |> filter(ret_strata == "Interruption after >30 days")) |> 
  mutate(panel = "Interruption of 3-4 days", ret_strata = "Interruption after >30 days")

# Trial 5 (5-14 days)
rmst_tr3_le30 <- get_unadjusted_rmst(dat_tr3 |> filter(ret_strata == "Interruption within first \u226430 days")) |> 
  mutate(panel = "Interruption of 5-14 days", ret_strata = "Interruption within first \u226430 days")

rmst_tr3_gt30 <- get_unadjusted_rmst(dat_tr3 |> filter(ret_strata == "Interruption after >30 days")) |> 
  mutate(panel = "Interruption of 5-14 days", ret_strata = "Interruption after >30 days")

# Combine all unadjusted results
rmst_abs <- bind_rows(rmst_tr2_le30, rmst_tr2_gt30, rmst_tr3_le30, rmst_tr3_gt30) |>
  mutate(ret_strata = factor(ret_strata, levels = c("Interruption within first \u226430 days", "Interruption after >30 days")))


# ── 5. Visualization ────────────────────────────────────────────────────────
plot_bar <- ggplot(rmst_abs, aes(x = trt, y = RMST, fill = trt)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = RMST_lower, ymax = RMST_upper), width = 0.2, linewidth = 0.7) +
  geom_text(aes(label = paste0(round(RMST), " days"), y = RMST / 2), 
            color = "white", angle = 90, size = 4, fontface = "bold") +
  facet_grid(ret_strata~panel) +
  scale_fill_manual(values = PLOT_COLORS) +
  scale_y_continuous(limits = c(0, 365), expand = c(0, 0)) +
  labs(
    title = NULL, 
    x = NULL, 
    y = "Unadjusted restricted mean survival days"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey85", color = "grey85"),
    strip.text = element_text(face = "bold", size = 14, color = "black"),
    axis.ticks        = element_line(color = "black", linewidth = 0.8), axis.ticks.length = unit(0.2, "cm"),
    axis.text.y         = element_text(color = "black", size = 13), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13, color = "black"),
    axis.title = element_text(color = "black", size = 14), 
    panel.grid = element_blank()
  )

ggsave(paste0(RES_DIR, "Figures/eFigure6_Unadjusted_BarPlot_MET.jpg"), plot_bar, width = 12, height = 10, dpi = 300)


# LE30D vs GT30D: MET ##########

library(tidyverse)
library(haven)
library(survival)
library(patchwork)

# ── 1. Paths & Setup ────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/LE30D/"
FIG_DIR   <- paste0(RES_DIR, "Figures/")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
TRT_LABELS <- c("No dose change", "<50% reduction",
                "\u226550% reduction/\u226440mg", "Dose increase")

PLOT_COLORS <- c("No dose change"                 = "#333333",
                 "<50% reduction"                 = "#4DAF4A",
                 "\u226550% reduction/\u226440mg" = "#FF7F00",
                 "Dose increase"                  = "#E41A1C")

# ── 2. Load and Clean Raw Data (Any Episode) ────────────────────────────────
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_clean <- raw_all |>
  mutate(
    trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)
  ) |>
  rename(trt = trt_s456) |>
  filter(trt %in% TRT_LEVELS) |>
  mutate(
    trt = factor(trt, levels = TRT_LEVELS),
    ret_strata = if_else(ret_length <= 30,
                         "Interruption within first \u226430 days",
                         "Interruption after >30 days"),
    ret_strata = factor(ret_strata,
                        levels = c("Interruption within first \u226430 days",
                                   "Interruption after >30 days"))
  )

dat_tr2 <- dat_clean |> filter(trial_gap %in% 3:4)
dat_tr3 <- dat_clean |> filter(trial_cat == "5-14")

rm(raw_all, dat_clean); gc()

# ── 3. Helper Function: Unadjusted RMST via Kaplan-Meier ────────────────────
get_unadjusted_rmst <- function(df, horizon = 365) {
  fit <- survfit(Surv(time_discont, sens_discont) ~ trt, data = df)
  sum_fit <- summary(fit, rmean = horizon)$table
  
  tibble(
    trt_raw = rownames(sum_fit),
    RMST    = sum_fit[, "rmean"],
    se      = sum_fit[, "se(rmean)"]
  ) |>
    mutate(
      trt        = str_remove(trt_raw, "trt="),
      trt        = factor(trt, levels = TRT_LEVELS, labels = TRT_LABELS),
      RMST_lower = RMST - 1.96 * se,
      RMST_upper = RMST + 1.96 * se
    ) |>
    select(trt, RMST, RMST_lower, RMST_upper)
}

# ── 4. Calculate Stratified RMST per Trial ──────────────────────────────────
# TR2: 3-4 days
rmst_tr2 <- bind_rows(
  get_unadjusted_rmst(dat_tr2 |> filter(ret_strata == "Interruption within first \u226430 days")) |>
    mutate(ret_strata = "Interruption within first \u226430 days"),
  get_unadjusted_rmst(dat_tr2 |> filter(ret_strata == "Interruption after >30 days")) |>
    mutate(ret_strata = "Interruption after >30 days")
) |>
  mutate(ret_strata = factor(ret_strata,
                             levels = c("Interruption within first \u226430 days",
                                        "Interruption after >30 days")))

# TR3: 5-14 days
rmst_tr3 <- bind_rows(
  get_unadjusted_rmst(dat_tr3 |> filter(ret_strata == "Interruption within first \u226430 days")) |>
    mutate(ret_strata = "Interruption within first \u226430 days"),
  get_unadjusted_rmst(dat_tr3 |> filter(ret_strata == "Interruption after >30 days")) |>
    mutate(ret_strata = "Interruption after >30 days")
) |>
  mutate(ret_strata = factor(ret_strata,
                             levels = c("Interruption within first \u226430 days",
                                        "Interruption after >30 days")))

# ── 5. Plot Builder ─────────────────────────────────────────────────────────
build_rmst_barplot <- function(df) {
  ggplot(df, aes(x = trt, y = RMST, fill = trt)) +
    geom_col(color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = RMST_lower, ymax = RMST_upper),
                  width = 0.2, linewidth = 0.7) +
    geom_text(aes(label = paste0(round(RMST), " days"), y = RMST / 2),
              color = "white", angle = 90, size = 4, fontface = "bold") +
    facet_wrap(~ ret_strata, ncol = 2) +
    scale_fill_manual(values = PLOT_COLORS) +
    scale_y_continuous(limits = c(0, 365), expand = c(0, 0)) +
    labs(title = NULL, x = NULL,
         y = "Unadjusted restricted mean survival days") +
    theme_bw() +
    theme(
      legend.position   = "none",
      strip.background  = element_rect(fill = "grey85", color = "grey85"),
      strip.text        = element_text(face = "bold", size = 14, color = "black"),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text.y       = element_text(color = "black", size = 13),
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 13, color = "black"),
      axis.title        = element_text(color = "black", size = 14),
      panel.grid        = element_blank()
    )
}

# ── 6. Build and Save Two Separate Figures ──────────────────────────────────
plot_tr2 <- build_rmst_barplot(rmst_tr2)
plot_tr3 <- build_rmst_barplot(rmst_tr3)

ggsave(paste0(FIG_DIR, "eFigure6_Unadjusted_BarPlot_MET_TR2_34d.jpg"),
       plot_tr2, width = 12, height = 6, dpi = 300)

ggsave(paste0(FIG_DIR, "eFigure6_Unadjusted_BarPlot_MET_TR3_514d.jpg"),
       plot_tr3, width = 12, height = 6, dpi = 300)

cat("\n\u2713 Both figures saved to:\n  ", FIG_DIR, "\n", sep = "")





# SUMMARY: TEXT ##############

#  SUMMARY STATISTICS FOR MANUSCRIPT PARAGRAPH
#  Compares retention (RMST) and 12-month discontinuation across strata
library(tidyverse)
library(haven)
library(survival)

# ── Paths & Setup ────────────────────────────────────────────────────────────
DATA_DIR <- "E:/SAS Data/Missed dose protocol_v2/"
OUT_DIR  <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/LE30D/Tables/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
TRT_LABELS <- c("No dose change", "<50% reduction",
                "\u226550% reduction/\u226440mg", "Dose increase")

# ── Load & Clean ─────────────────────────────────────────────────────────────
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_m.sas7bdat"))

dat_clean <- raw_all |>
  mutate(trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)) |>
  rename(trt = trt_s456) |>
  filter(trt %in% TRT_LEVELS) |>
  mutate(
    trt = factor(trt, levels = TRT_LEVELS),
    ret_strata = if_else(ret_length <= 30, "LE30D", "GT30D"),
    ret_strata = factor(ret_strata, levels = c("LE30D", "GT30D"))
  )

dat_tr2 <- dat_clean |> filter(trial_gap %in% 3:4)
dat_tr3 <- dat_clean |> filter(trial_cat == "5-14")

rm(raw_all, dat_clean); gc()

# ── Helper 1: Per-arm RMST + 12-month cumulative incidence ─────────────────
summarise_arm <- function(df, horizon = 365) {
  fit <- survfit(Surv(time_discont, sens_discont) ~ trt, data = df)
  
  # RMST
  rmst_tbl <- summary(fit, rmean = horizon)$table |>
    as.data.frame() |>
    rownames_to_column("trt_raw") |>
    transmute(
      trt        = factor(str_remove(trt_raw, "trt="),
                          levels = TRT_LEVELS, labels = TRT_LABELS),
      N          = as.integer(records),
      Events     = as.integer(events),
      RMST       = round(rmean, 0),
      RMST_lo    = round(rmean - 1.96 * `se(rmean)`, 0),
      RMST_hi    = round(rmean + 1.96 * `se(rmean)`, 0)
    )
  
  # 12-month cumulative incidence per arm
  ci12 <- summary(fit, times = horizon, extend = TRUE)
  cum12 <- tibble(
    trt_raw   = ci12$strata,
    surv12    = ci12$surv,
    cum12_pct = round((1 - ci12$surv) * 100, 1)
  ) |>
    mutate(trt = factor(str_remove(as.character(trt_raw), "trt_label=|trt="),
                        levels = TRT_LEVELS, labels = TRT_LABELS)) |>
    select(trt, cum12_pct)
  
  rmst_tbl |> left_join(cum12, by = "trt")
}

# ── Helper 2: Overall (pooled) discontinuation at 12 months ────────────────
summarise_overall <- function(df, horizon = 365, label = "") {
  fit <- survfit(Surv(time_discont, sens_discont) ~ 1, data = df)
  s12 <- summary(fit, times = horizon, extend = TRUE)
  rmean <- summary(fit, rmean = horizon)$table
  
  tibble(
    Cohort         = label,
    N_episodes     = as.integer(rmean["records"]),
    N_events       = as.integer(rmean["events"]),
    Cum_disc_12mo  = round((1 - s12$surv) * 100, 1),
    RMST_days      = round(rmean["rmean"], 0)
  )
}

# ── 1. OVERALL discontinuation by trial (pooling both strata) ───────────────
overall_stats <- bind_rows(
  summarise_overall(dat_tr2, label = "TR2 (3-4 days, all)"),
  summarise_overall(dat_tr2 |> filter(ret_strata == "LE30D"),
                    label = "TR2 (3-4 days, LE30D)"),
  summarise_overall(dat_tr2 |> filter(ret_strata == "GT30D"),
                    label = "TR2 (3-4 days, GT30D)"),
  summarise_overall(dat_tr3, label = "TR3 (5-14 days, all)"),
  summarise_overall(dat_tr3 |> filter(ret_strata == "LE30D"),
                    label = "TR3 (5-14 days, LE30D)"),
  summarise_overall(dat_tr3 |> filter(ret_strata == "GT30D"),
                    label = "TR3 (5-14 days, GT30D)")
)

# ── 2. PER-ARM stats by trial × stratum ─────────────────────────────────────
per_arm_stats <- bind_rows(
  summarise_arm(dat_tr2)                                   |> mutate(cohort = "TR2 (3-4 days, all)"),
  summarise_arm(dat_tr2 |> filter(ret_strata == "LE30D"))  |> mutate(cohort = "TR2 (3-4 days, LE30D)"),
  summarise_arm(dat_tr2 |> filter(ret_strata == "GT30D"))  |> mutate(cohort = "TR2 (3-4 days, GT30D)"),
  summarise_arm(dat_tr3)                                   |> mutate(cohort = "TR3 (5-14 days, all)"),
  summarise_arm(dat_tr3 |> filter(ret_strata == "LE30D"))  |> mutate(cohort = "TR3 (5-14 days, LE30D)"),
  summarise_arm(dat_tr3 |> filter(ret_strata == "GT30D"))  |> mutate(cohort = "TR3 (5-14 days, GT30D)")
) |>
  select(cohort, trt, N, Events, cum12_pct, RMST, RMST_lo, RMST_hi) |>
  mutate(RMST_CI = sprintf("%d (%d-%d)", RMST, RMST_lo, RMST_hi))

# ── 3. LE30D vs GT30D RMST DIFFERENCE per arm & trial ──────────────────────
rmst_diff <- per_arm_stats |>
  filter(cohort %in% c("TR2 (3-4 days, LE30D)", "TR2 (3-4 days, GT30D)",
                       "TR3 (5-14 days, LE30D)", "TR3 (5-14 days, GT30D)")) |>
  mutate(
    trial   = if_else(str_detect(cohort, "TR2"), "3-4 days", "5-14 days"),
    stratum = if_else(str_detect(cohort, "LE30D"), "LE30D", "GT30D")
  ) |>
  select(trial, stratum, trt, RMST, cum12_pct) |>
  pivot_wider(names_from = stratum, values_from = c(RMST, cum12_pct)) |>
  mutate(
    RMST_diff_GT_minus_LE = RMST_GT30D - RMST_LE30D,
    Disc_diff_GT_minus_LE = cum12_pct_GT30D - cum12_pct_LE30D
  )

# ── 4. PRINT & SAVE ─────────────────────────────────────────────────────────
cat("\n============================================================\n")
cat(" OVERALL 12-MONTH DISCONTINUATION (pooled across arms) \n")
cat("============================================================\n")
print(overall_stats, n = Inf)

cat("\n============================================================\n")
cat(" PER-ARM SUMMARY (RMST + 12-mo discontinuation %) \n")
cat("============================================================\n")
print(per_arm_stats |> select(cohort, trt, N, Events, cum12_pct, RMST_CI))

cat("\n============================================================\n")
cat(" LE30D vs GT30D COMPARISON (per arm, per trial) \n")
cat("============================================================\n")
print(rmst_diff, n = Inf)

write_csv(overall_stats,  paste0(OUT_DIR, "T1_overall_discontinuation.csv"))
write_csv(per_arm_stats,  paste0(OUT_DIR, "T2_per_arm_RMST.csv"))
write_csv(rmst_diff,      paste0(OUT_DIR, "T3_LE30_vs_GT30_comparison.csv"))

cat(sprintf("\n\u2713 CSVs saved to: %s\n", OUT_DIR))



apply(rmst_diff[rmst_diff$trial == "3-4 days", 3:8], 2, mean)  
apply(rmst_diff[rmst_diff$trial == "5-14 days", 3:8], 2, mean)  

head(per_arm_stats)

apply(per_arm_stats[per_arm_stats$cohort == "TR2 (3-4 days, LE30D)", 4:7], 2, mean)
apply(per_arm_stats[per_arm_stats$cohort == "TR2 (3-4 days, GT30D)", 4:7], 2, mean)

apply(per_arm_stats[per_arm_stats$cohort == "TR3 (5-14 days, LE30D)", 4:7], 2, mean)
apply(per_arm_stats[per_arm_stats$cohort == "TR3 (5-14 days, GT30D)", 4:7], 2, mean)


apply(rmst_diff[rmst_diff$trial == "5-14 days", 3:8], 2, mean)  
