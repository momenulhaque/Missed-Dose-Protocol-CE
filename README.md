# Missed Dose Protocols for Opioid Agonist Treatment: A Population-Based Target Trial Emulations

**Author:** Momenul Haque Mondol

This repository contains the SAS and R code used to emulate target trials evaluating alternative dose-adjustment strategies following interruptions (missed doses) of opioid agonist treatment (OAT), using population-based linked administrative health data from British Columbia (BC), Canada.

---

## About This Repository

This repository is built to support the analyses for **two companion papers** that share the same data infrastructure, analytic framework, and target trial emulation approach, but evaluate two different medications for opioid use disorder (OUD):

1. **Methadone (MET)**
2. **Buprenorphine/Naloxone (BNX)**

### File & folder conventions

To keep the two analyses organized within a single repository:

- Scripts and files ending in **`_met`** belong to the **methadone paper**.
- Scripts and files ending in **`_bnx`** belong to the **buprenorphine/naloxone paper**.
- **Results are removed**
- Shared helper functions and general pipeline scripts (where applicable) are used by both analyses.

Each paper's specific documentation is provided in its own section below.

---

# Paper 1 — Methadone (MET)

### Title
**Comparative effectiveness of missed dose management for methadone in British Columbia, Canada: population-based target trial emulation**

### Authors
Momenul Haque Mondol, Jeong Eun Min, Michelle Zanette, Md. Belal Hossain, Megan Kurz, Paxton Bach, Robert William Platt, Shaun Seaman, Sander Greenland, Lawrence C. McCandless, Mohammad Ehsanul Karim, Paul Gustafson, Maria Eugenia Socías, Jason M Sutherland, Bohdan Nosyk.

### Code (files ending in `_met`)
- `01. create_analytic_data_met` — build the methadone analytic dataset (SAS).
- `01. create_flow_chart_met_bnx` — cohort flow chart (SAS).
- `02_prepare_hdPS_analysis` — prepare data for hdPS (SAS).
- `00_functions` — helper functions (R).
- `02_main_analysis_rmst_met_LE30D, 02_main_analysis_rmst_met_GT30D` — primary analysis using propensity score weighting (R).
- `03_Timing of interruptions_met` — exploring the timing of interruptions (R).
- `02_sa_analysis_met_LE30D_1, 02_sa_analysis_met_LE30D_2, 02_sa_analysis_met_GT30D_1 and 02_sa_analysis_met_GT30D_2 ` — sensitivity analyses (R).
- `02_Table1` — create Tables of baseline characteristics (R).
- `03_forestplot_main_met` — create forestplot for primary results (R).
- `03_KM_curves_met, 03_risk_curves_met, and 03_RMST_bars_met` — Figures for KM estimates, risk, and restricted mean survival time respectively (R).
- `03_sa_Table_TRIAL 2_met_LE30D, 03_sa_Table_TRIAL 2_met_GT30D, 03_sa_Table_TRIAL 3_met_LE30D, and 03_sa_Table_TRIAL 3_met_GT30D` — eTables for sensitivity analysis (R).
- `04_assessing_covariates_dist_met, and 04_sa_assessing overlap_met` — Some diagnostic plots.
- `generate_BC_health_neighbour_csv` — Prepare geographical clusters part 1.

### Results
Stored in the **Results\>MET** folder (deleted now)

---

# Paper 2 — Buprenorphine/Naloxone (BNX)

### Title
**Comparative effectiveness of missed dose management for buprenorphine/naloxone in British Columbia, Canada: population-based target trial emulation**

### Authors
Momenul Haque Mondol, Kiera Tsakonas, Jeong Eun Min, Michelle Zanette, Md. Belal Hossain, Megan Kurz, Paxton Bach, Robert William Platt, Shaun Seaman, Sander Greenland, Lawrence C. McCandless, Mohammad Ehsanul Karim, Paul Gustafson, Maria Eugenia Socías, Jason M Sutherland, Bohdan Nosyk.

### Code (files ending in `_bnx`)
- `01. create_analytic_data_bnx` — build the BNX analytic dataset (SAS).
- `01. create_flow_chart_met_bnx` — cohort flow chart (SAS).
- `02_prepare_hdPS_analysis` — prepare data for hdPS (SAS).
- `00_functions` — helper functions (R).
- `02_main_analysis_rmst_bnx` — primary analysis using propensity score weighting (R).
- `03_Timing of interruptions_bnx` — exploring the timing of interruptions (R).
- `02_sa_analysis_bnx, 02_sa_analysis_ow_bnx, 02_sa_hdPS_bnx, 02_sa_pp_analysis_bnx` — sensitivity analyses (R).
- `03_forestplot_main_bnx` — create forestplot for primary results.
- `03_KM_curves_bnx, 03_risk_curves_bnx, and 03_RMST_bars_bnx` — Figures for KM estimates, risk, and restricted mean survival time respectively (R).
- `03_sa_forestplot_TRIAL 4_bnx, and 03_sa_forestplot_TRIAL 5_bnx` — Forestplot for sensitivity analyses (R).
- `04_assessing_covariates_dist_met, and 04_sa_assessing overlap_met` — Some diagnostic plots.
- 
### Results
Stored in the **Results>BNX** folder (deleted now)

---

## Data availability
The datasets of the study are not publicly available due to concerns on individual level administrative data and privacy restrictions.

---

## Funding

National Institutes of Health / National Institute on Drug Abuse, grant No. **R01-DA050629**.

---
