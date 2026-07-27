# =====================================================================
#  SA POSITIVITY & EVENT-COUNT AUDIT --- TRIAL 3 (5-14 DAYS)
#  Corrected file paths to match actual Trial 3 manifest.
# =====================================================================
library(tidyverse)
library(fs)

# ---- Paths -----------------------------------------------------------
RES_DIR   <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/MET/LE30D/"
SA_DIR    <- paste0(RES_DIR, "Sensitivities/")
PS_DIR    <- paste0(RES_DIR, "PScore/")
AUDIT_OUT <- paste0(SA_DIR, "SA_positivity_audit_TR3.csv")

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")

# ---- Outcome variable map --------------------------------------------
outcome_map <- tribble(
  ~sa_tag,              ~outcome_key, ~time_var,              ~event_var,
  # Defaults
  "default",            "od",         "time_od_acute",        "sens_od_acute",
  "default",            "dis",        "time_discont",         "sens_discont",
  "default",            "all",        "time_all_acute",       "sens_all_acute",
  # SA12 - outcome swap (no truncation)
  "SA12_NoTruncation",  "od",         "t_od_acute",           "s_od_acute",
  "SA12_NoTruncation",  "dis",        "t_discont",            "s_discont",
  # SA16 - composite OD
  "SA16_Composite_OD",  "od",         "time_od_comp_outcome", "sens_od_comp_outcome",
  "SA16_Composite_OD",  "dis",        "time_discont",         "sens_discont",
  # 30-day discont
  "discont30",          "dis",        "time_discont30",       "sens_discont30",
  # Per-protocol
  "per_protocol",       "od",         "time_od_acute_pp",     "sens_od_acute_pp",
  "per_protocol",       "dis",        "time_discont_pp",      "sens_discont_pp"
)

# ---- Enumerate every Trial 3 PS file (based on your manifest) --------
# Subgroup SAs under Sensitivities/ with `ps_wt_tr3_514d_*.rds` pattern
sa_ps_wt <- tibble(
  sa_tag = c("SA03_Dose_LT80", "SA04_Dose_GE80", "SA09_HighFentArea",
             "SA11_CumOAT_GE1yr", "SA13_30dWashout",
             "SA10_PostFentEra", "SA16_Composite_OD",
             # Include the ones not in your current manifest but exist on disk:
             "SA02_Induction", "SA05_Ret_LE30", "SA06_Ret_GT30",
             "SA12_NoTruncation"),
  trial  = "tr3_514d",
  file   = paste0(SA_DIR, "ps_wt_tr3_514d_",
                  c("SA03_Dose_LT80", "SA04_Dose_GE80", "SA09_HighFentArea",
                    "SA11_CumOAT_GE1yr", "SA13_30dWashout",
                    "SA10_PostFentEra", "SA16_Composite_OD",
                    "SA02_Induction", "SA05_Ret_LE30", "SA06_Ret_GT30",
                    "SA12_NoTruncation"),
                  ".rds")
)

# Files under Sensitivities/ with `ps_weight_tr3_514d_*.rds` pattern
sa_ps_weight <- tribble(
  ~sa_tag,       ~trial,       ~file,
  "5cat",        "tr3_514d",   paste0(SA_DIR, "ps_weight_tr3_514d_5cat.rds"),
  "discont30",   "tr3_514d",   paste0(SA_DIR, "ps_weight_tr3_514d_discont30.rds")
)

# Main / alternative PS models under PScore/
main_ps <- tribble(
  ~sa_tag,          ~trial,       ~file,
  # Primary analysis (any-episode) and all-acute for tr3
  "primary_all",    "tr3_514d",   paste0(PS_DIR, "ps_weight_tr3_all.rds"),
  "all_acute_v1",   "tr3_514d",   paste0(PS_DIR, "ps_weight_tr3_514d_all.rds"),  # in manifest
  # Incident users (first)
  "primary_first",  "tr3_514d",   paste0(PS_DIR, "ps_weight_tr3_first.rds"),
  # hdPS variants
  "hdPS",           "tr3_514d",   paste0(PS_DIR, "ps_weight_514d_hdPS.rds"),
  "hdPS_ML",        "tr3_514d",   paste0(PS_DIR, "ps_weight_514d_hdPS_ML.rds"),
  # Overlap weights and provider-covariate PS
  "OW",             "tr3_514d",   paste0(PS_DIR, "ps_weight_tr3_all.rds"),
  "cl",             "tr3_514d",   paste0(PS_DIR, "ps_weight_tr3_all_cl.rds")
)

all_files <- bind_rows(sa_ps_wt, sa_ps_weight, main_ps) |> distinct()

# ---- Core audit function ---------------------------------------------
audit_one <- function(sa_tag, trial, file) {
  
  .sa <- sa_tag
  .tr <- trial
  .fp <- file
  
  if (!file.exists(.fp)) {
    return(tibble(sa_tag = .sa, trial = .tr, outcome = NA_character_,
                  trt = NA_character_, n_episodes = NA_integer_,
                  n_events = NA_integer_, status = "MISSING"))
  }
  
  ps_obj <- tryCatch(readRDS(.fp), error = function(e) NULL)
  if (is.null(ps_obj)) {
    return(tibble(sa_tag = .sa, trial = .tr, outcome = NA_character_,
                  trt = NA_character_, n_episodes = NA_integer_,
                  n_events = NA_integer_, status = "UNREADABLE"))
  }
  
  dw <- if ("data_weights" %in% names(ps_obj)) ps_obj$data_weights else ps_obj
  if (is.null(dw) || !"trt" %in% names(dw)) {
    return(tibble(sa_tag = .sa, trial = .tr, outcome = NA_character_,
                  trt = NA_character_, n_episodes = NA_integer_,
                  n_events = NA_integer_, status = "NO_TRT"))
  }
  
  specific <- outcome_map |> dplyr::filter(sa_tag == .sa)
  defaults <- outcome_map |> dplyr::filter(sa_tag == "default")
  
  out_this <- dplyr::bind_rows(
    specific,
    defaults |> dplyr::filter(!outcome_key %in% specific$outcome_key)
  )
  
  purrr::pmap_dfr(
    list(out_this$outcome_key, out_this$time_var, out_this$event_var),
    function(okey, tvar, evar) {
      
      if (!all(c(tvar, evar, "trt") %in% names(dw))) {
        return(tibble(sa_tag = .sa, trial = .tr, outcome = okey,
                      trt = NA_character_, n_episodes = NA_integer_,
                      n_events = NA_integer_, status = "MISSING_VAR"))
      }
      
      dw |>
        dplyr::group_by(trt) |>
        dplyr::summarise(
          n_episodes = dplyr::n(),
          n_events   = sum(.data[[evar]] == 1, na.rm = TRUE),
          .groups    = "drop"
        ) |>
        dplyr::mutate(
          sa_tag  = .sa,
          trial   = .tr,
          outcome = okey,
          status  = "OK"
        ) |>
        dplyr::select(sa_tag, trial, outcome, trt,
                      n_episodes, n_events, status)
    }
  )
}

# ---- Run audit -------------------------------------------------------
cat("Auditing", nrow(all_files), "Trial 3 SA/PS files...\n\n")

# Show which files are being checked (helpful for troubleshooting)
cat("--- Files enumerated ---\n")
all_files |> mutate(exists = file.exists(file)) |>
  select(sa_tag, exists, file) |>
  print(n = Inf)

audit_tbl <- purrr::pmap_dfr(all_files, audit_one)

# ---- Add positivity flags --------------------------------------------
audit_flagged <- audit_tbl |>
  mutate(
    positivity_violation = (n_events == 0),
    sparse_events        = (n_events > 0 & n_events < 5),
    display_rule = case_when(
      positivity_violation ~ "SUPPRESS_ESTIMATE_show_N_only",
      sparse_events        ~ "SHOW_WITH_CAUTION_flag",
      TRUE                 ~ "SHOW_ESTIMATE"
    )
  )

write_csv(audit_flagged, AUDIT_OUT)

violation_summary <- audit_flagged |>
  group_by(sa_tag, trial, outcome) |>
  summarise(
    any_zero_arm   = any(positivity_violation, na.rm = TRUE),
    any_sparse_arm = any(sparse_events, na.rm = TRUE),
    zero_arms      = paste(trt[positivity_violation], collapse = "; "),
    sparse_arms    = paste(trt[sparse_events], collapse = "; "),
    min_n_events   = suppressWarnings(min(n_events, na.rm = TRUE)),
    total_N        = sum(n_episodes, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(any_zero_arm), sa_tag, trial, outcome)

cat("\n=== TRIAL 3: Positivity violations (any arm with 0 events) ===\n")
print(violation_summary |> filter(any_zero_arm), n = Inf)

cat("\n=== TRIAL 3: Sparse events (any arm with 1-4 events) ===\n")
print(violation_summary |> filter(any_sparse_arm & !any_zero_arm), n = Inf)

cat("\n=== TRIAL 3: Files that could not be loaded ===\n")
audit_flagged |>
  filter(status %in% c("MISSING", "UNREADABLE", "NO_TRT", "MISSING_VAR")) |>
  distinct(sa_tag, trial, status) |>
  print(n = Inf)

cat(sprintf("\nFull Trial 3 audit saved to: %s\n", AUDIT_OUT))