library(tidyverse)
library(haven)
library(survival)
library(patchwork)

# ── 1. Paths & Setup ────────────────────────────────────────────────────────
DATA_DIR  <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/BNX/"

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=4mg", "increase")
TRT_LABELS <- c("No dose change", "<50% reduction", "\u226550% reduction/\u22644mg", "Dose increase")

PLOT_COLORS <- c("No dose change" = "#333333", "<50% reduction" = "#4DAF4A", 
                 "\u226550% reduction/\u22644mg" = "#FF7F00", "Dose increase" = "#E41A1C")

# ── 2. Load and Clean Raw Data (Any Episode) ────────────────────────────────
raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_sens_nw_b.sas7bdat"))

dat_clean <- raw_all |>
  mutate(
    trt_s5 = if_else(trt_s5 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s5)
  ) |>
  rename(trt = trt_s5) |>
  filter(trt %in% TRT_LEVELS) |>
  mutate(
    trt = factor(trt, levels = TRT_LEVELS),
    # Create the retention stratification variable
    ret_strata = if_else(ret_length <= 30, "\u226430 days", ">30 days"),
    ret_strata = factor(ret_strata, levels = c("\u226430 days", ">30 days"))
  )

# Split into Trials
dat_tr4 <- dat_clean |> filter(trial_gap %in% 3:5)
dat_tr5 <- dat_clean |> filter(trial_cat == "6-14")

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
rmst_tr4_le30 <- get_unadjusted_rmst(dat_tr4 |> filter(ret_strata == "\u226430 days")) |> 
  mutate(panel = "Interruption of 3-5 days", ret_strata = "\u226430 days")

rmst_tr4_gt30 <- get_unadjusted_rmst(dat_tr4 |> filter(ret_strata == ">30 days")) |> 
  mutate(panel = "Interruption of 3-5 days", ret_strata = ">30 days")

# Trial 5 (5-14 days)
rmst_tr5_le30 <- get_unadjusted_rmst(dat_tr5 |> filter(ret_strata == "\u226430 days")) |> 
  mutate(panel = "Interruption of 6-14 days", ret_strata = "\u226430 days")

rmst_tr5_gt30 <- get_unadjusted_rmst(dat_tr5 |> filter(ret_strata == ">30 days")) |> 
  mutate(panel = "Interruption of 6-14 days", ret_strata = ">30 days")

# Combine all unadjusted results
rmst_abs <- bind_rows(rmst_tr4_le30, rmst_tr4_gt30, rmst_tr5_le30, rmst_tr5_gt30) |>
  mutate(ret_strata = factor(ret_strata, levels = c("\u226430 days", ">30 days")))


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

ggsave(paste0(RES_DIR, "Figures/eFigure6_Unadjusted_BarPlot_BNX.jpg"), plot_bar, width = 12, height = 10, dpi = 300)







get_unadjusted_rmst(dat_tr5)
get_unadjusted_rmst(dat_tr4)
