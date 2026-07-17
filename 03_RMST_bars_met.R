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