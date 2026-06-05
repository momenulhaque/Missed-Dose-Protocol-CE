# =============================================================================
# 07_sa_confounding_by_cluster.R
# -----------------------------------------------------------------------------
# Provider-level assessment of potential confounding by cluster, at TWO levels:
#   - Facility / provider-network  (com)            <- primary ("clinic" proxy)
#   - Prescriber                   (prac_studyid_rev) <- secondary
#
#   PART 1 - Quantify between-cluster variation in EXPOSURE and OUTCOMES
#            (ICC and median odds ratio, MOR) from empty random-intercept
#            logistic models, at both levels.
#
#   PART 2 - Marginal between-within sensitivity analysis (Sjolander 2021;
#            Li et al. 2013): add each cluster's mean exposure level to the
#            weighted pooled-logistic outcome model, cluster SEs at that level,
#            and compare 1-year cumulative incidence with primary.
#
# Inputs : Results/PScore/*.rds (from 01_main_analysis.R)
#          TT_bw_tr2 / TT_bw_tr3 (from 02_confounding_by_cluster_prep.sas)
# Output : Results/Tables/cluster_confounding.xlsx
#
# Note: neither com nor prac_studyid_rev is retained in
#       ps_model_wt()$data_weights, so both are pulled back in from the SAS
#       BW file inside prep_dat(), together with the cluster-mean terms.
# =============================================================================

library(tidyverse)
library(haven)
library(lme4)
library(openxlsx)

source("H:/GitHub/Missed-Dose-Protocol-CE/00_functions.R")  # expand_person_period(), predict_cuminc()

# -- Paths and constants ------------------------------------------------------
DATA_DIR   <- "E:/SAS Data/Missed dose protocol_v2/"
RES_DIR    <- "H:/GitHub/Missed-Dose-Protocol-CE/Results/"

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")
N_BOOT     <- 499      # lower (e.g. 199) for quick test runs
H          <- 365      # 1-year cumulative incidence

# Cluster levels to assess. "faci" first = primary (clinic/network proxy).
LEVELS <- list(
  faci = list(id = "com",              terms = c("pmean_red_faci", "pmean_big_faci", "pmean_incr_faci"), label = "BW (facility)"),
  prac = list(id = "prac_studyid_rev", terms = c("pmean_red_prac", "pmean_big_prac", "pmean_incr_prac"), label = "BW (prescriber)")
)

cfg <- list(
  tr2 = list(ps = "PScore/ps_weight_tr2_34d_all.rds", bw = "TT_bw_tr2.sas7bdat",
             covs = list(od = c("ddose_prev", "prior_episode"),
                         dis = c("ddose_prev", "prior_episode"))),
  tr3 = list(ps = "PScore/ps_weight_tr3_all.rds", bw = "TT_bw_tr3.sas7bdat",
             covs = list(od = c("ddose_prev", "year"),
                         dis = c("ddose_prev", "year")))
)

outcomes <- list(
  od  = list(time = "time_od_acute", event = "sens_od_acute", label = "Overdose"),
  dis = list(time = "time_discont",  event = "sens_discont",  label = "Discontinuation")
)


# -- Helpers ------------------------------------------------------------------

read_bw <- function(path) {
  read_sas(path) |>
    mutate(trt = if_else(trt_s456 %in% c("1-<25%", "25-<50%"), "1-<50%", trt_s456)) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(trt        = factor(trt, levels = TRT_LEVELS),
           any_change = as.integer(trt != "no reversion"))
}

icc_mor <- function(sigma2) {
  c(ICC = sigma2 / (sigma2 + pi^2 / 3),
    MOR = exp(sqrt(2 * sigma2) * qnorm(0.75)))
}

fit_empty <- function(dat, yvar, cluster) {
  f <- reformulate(sprintf("(1 | %s)", cluster), response = yvar)
  m <- glmer(f, data = dat, family = binomial,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  icc_mor(as.numeric(VarCorr(m)[[cluster]]))
}

std_risks <- function(x, h = H) {
  r  <- if (is.list(x) && !is.null(x$risks)) x$risks else x
  r  <- as.data.frame(r); nm <- names(r)
  hz <- intersect(c("horizon", "time", "day", "days", "t"), nm)
  if (length(hz) >= 1) r <- r[r[[hz[1]]] == h, , drop = FALSE]
  pick <- function(pat) nm[grep(pat, nm, ignore.case = TRUE)][1]
  tibble(
    trt  = r[[pick("^trt$|group|exposure")]],
    risk = r[[pick("^risk$|cuminc|incidence|^est")]],
    lo   = r[[pick("lower|lcl|low|ci_l")]],
    hi   = r[[pick("upper|ucl|high|ci_u")]]
  )
}

fmt_pct <- function(risk, lo, hi)
  sprintf("%.1f (%.1f, %.1f)", 100 * risk, 100 * lo, 100 * hi)


# =============================================================================
# PART 1 - ICC / MOR at both cluster levels
# =============================================================================
cat("\n=== PART 1: between-cluster variation (ICC, MOR) ===\n")

variation_rows <- list()
for (tr in names(cfg)) {
  d <- read_bw(file.path(DATA_DIR, cfg[[tr]]$bw))
  
  for (lv in names(LEVELS)) {
    id <- LEVELS[[lv]]$id
    dd <- d |> filter(!is.na(.data[[id]]))
    
    for (y in c("any_change", "sens_od_acute", "sens_discont")) {
      res <- tryCatch(fit_empty(dd, y, id),
                      error = function(e) c(ICC = NA_real_, MOR = NA_real_))
      variation_rows[[length(variation_rows) + 1]] <- tibble(
        Trial  = tr,
        Level  = if (lv == "faci") "Facility/network" else "Prescriber",
        Target = recode(y,
                        any_change    = "Exposure: any dose change",
                        sens_od_acute = "Outcome: overdose",
                        sens_discont  = "Outcome: discontinuation"),
        ICC = unname(res["ICC"]), MOR = unname(res["MOR"])
      )
    }
  }
  cat(sprintf("  %s done\n", tr))
}

variation_tab <- bind_rows(variation_rows) |>
  mutate(across(c(ICC, MOR), ~round(.x, 3)))
print(variation_tab)

> print(variation_tab)
# A tibble: 12 × 5
    Trial Level            Target                      ICC   MOR
     <chr> <chr>            <chr>                     <dbl> <dbl>
  1 tr2   Facility/network Exposure: any dose change 0      1   
  2 tr2   Facility/network Outcome: overdose         0.047  1.47
  3 tr2   Facility/network Outcome: discontinuation  0.008  1.16
  4 tr2   Prescriber       Exposure: any dose change 0.032  1.37
  5 tr2   Prescriber       Outcome: overdose         0.25   2.71
  6 tr2   Prescriber       Outcome: discontinuation  0.058  1.54
  7 tr3   Facility/network Exposure: any dose change 0.02   1.28
  8 tr3   Facility/network Outcome: overdose         0.037  1.40
  9 tr3   Facility/network Outcome: discontinuation  0.01   1.19
 10 tr3   Prescriber       Exposure: any dose change 0.059  1.54
 11 tr3   Prescriber       Outcome: overdose         0.236  2.61
 12 tr3   Prescriber       Outcome: discontinuation  0.062  1.56

 
# =============================================================================
# PART 2 - Between-within sensitivity analysis at both levels
# =============================================================================
cat("\n=== PART 2: between-within sensitivity analysis ===\n")

# Pull both cluster ids + both sets of mean-exposure terms onto the PS data.
prep_dat <- function(ps_rds, bw_path) {
  ps <- readRDS(file.path(RES_DIR, ps_rds))$data_weights
  bw <- read_sas(file.path(DATA_DIR, bw_path)) |>
    transmute(studyid_rev, time0, com, prac_studyid_rev,
              pmean_red_prac, pmean_big_prac, pmean_incr_prac,
              pmean_red_faci, pmean_big_faci, pmean_incr_faci)
  for (v in c("com", "prac_studyid_rev"))
    if (v %in% names(ps)) ps <- ps |> select(-all_of(v))
  out <- ps |> left_join(bw, by = c("studyid_rev", "time0"))
  n_na <- sum(is.na(out$com))
  if (n_na > 0)
    message(sprintf("  note: %d/%d episodes did not match a facility id (com)", n_na, nrow(out)))
  out
}

run_ci <- function(long, covs, cluster_var) {
  predict_cuminc(
    Dat_expanded = long, trt_var = "trt", covs = covs,
    horizons = H, cluster_var = cluster_var,
    samples = N_BOOT, conf_int = TRUE,
    ref_group = "no reversion", wt_type = "wt_trunc"
  ) |> std_risks()
}

bw_long <- list()
for (tr in names(cfg)) {
  cat(sprintf("\nTrial %s\n", tr))
  dat <- prep_dat(cfg[[tr]]$ps, cfg[[tr]]$bw)
  
  for (oc in names(outcomes)) {
    o    <- outcomes[[oc]]
    covs <- cfg[[tr]]$covs[[oc]]
    cat(sprintf("  outcome: %s\n", o$label))
    
    long <- expand_person_period(dat, time_var = o$time,
                                 event_var = o$event, interval_length = 14)
    
    # Primary (person-level clustering, as in 01_main_analysis.R)
    bw_long[[length(bw_long) + 1]] <-
      run_ci(long, covs, cluster_var = "studyid_rev") |>
      mutate(Trial = tr, Outcome = o$label, Model = "Primary")
    
    # BW-augmented at each cluster level
    for (lv in names(LEVELS)) {
      L <- LEVELS[[lv]]
      bw_long[[length(bw_long) + 1]] <-
        run_ci(long, c(covs, L$terms), cluster_var = L$id) |>
        mutate(Trial = tr, Outcome = o$label, Model = L$label)
    }
    rm(long); gc()
  }
}

bw_long <- bind_rows(bw_long) |>
  mutate(trt   = factor(trt, levels = TRT_LEVELS),
         Model = factor(Model, levels = c("Primary", "BW (facility)", "BW (prescriber)"))) |>
  arrange(Trial, Outcome, Model, trt)

make_wide <- function(df) {
  df |>
    mutate(`1-yr cumulative incidence, % (95% CI)` = fmt_pct(risk, lo, hi)) |>
    dplyr::select(Model, trt, `1-yr cumulative incidence, % (95% CI)`) |>
    pivot_wider(names_from = Model, values_from = `1-yr cumulative incidence, % (95% CI)`)
}

cat("\n=== 1-year cumulative incidence: primary vs between-within ===\n")
wide_tabs <- list()
for (tr in names(cfg)) for (oc in names(outcomes)) {
  key <- paste0(tr, "_", oc)
  wt  <- bw_long |> filter(Trial == tr, Outcome == outcomes[[oc]]$label) |> make_wide()
  wide_tabs[[key]] <- wt
  cat(sprintf("\nTrial %s - %s\n", tr, outcomes[[oc]]$label)); print(wt)
}


# =============================================================================
# SAVE
# =============================================================================
wb <- createWorkbook()
addWorksheet(wb, "ICC_MOR"); writeData(wb, "ICC_MOR", variation_tab)
addWorksheet(wb, "BW_long")
writeData(wb, "BW_long", bw_long |> mutate(across(c(risk, lo, hi), ~round(100 * .x, 2))))
for (key in names(wide_tabs)) { addWorksheet(wb, key); writeData(wb, key, wide_tabs[[key]]) }
dir.create(file.path(RES_DIR, "Tables"), showWarnings = FALSE)
saveWorkbook(wb, file.path(RES_DIR, "Tables", "cluster_confounding.xlsx"), overwrite = TRUE)
cat("\nSaved: Results/Tables/cluster_confounding.xlsx\n")
