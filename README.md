# Comparative Effectiveness of Missed Dose Management for Methadone in British Columbia, Canada: A Population-Based Target Trial Emulation

**Author:** Momenul Haque Mondol

This repository contains the SAS and R code used to emulate target trials evaluating alternative methadone dose‑adjustment strategies following treatment interruptions (missed doses), using population‑based linked administrative health data from British Columbia (BC), Canada.

---

## Overview

Clinical guidelines for opioid use disorder (OUD) recommend adjusting the methadone dose after a patient misses doses, but the empirical evidence base is limited and relies heavily on expert consensus. Using linked administrative health records from BC (2010–2022), we emulated two target trials — one for methadone interruptions of **3–4 days** and one for **5–14 days** — to estimate how alternative reinitiation strategies affect two outcomes:

- **Overdose-related acute care visits** (hospital or emergency department encounters for poisoning; ICD-10 codes T40, T42.4, T43.6) — assessed at 1 month.
- **Treatment discontinuation** (no OAT dispensation for ≥15 consecutive days) — assessed at 12 months.

At reinitiation, the first dispensed dose was categorized relative to the previous dose as: **no change (reference)**, **<50% reduction**, **≥50% reduction** (or reduction to the guideline-recommended starting dose of ≤40 mg), or **dose increase**. Findings support **tailoring missed-dose protocols to treatment phase**: within the first 30 days of treatment, dose increases after interruptions were associated with lower risk of both outcomes, whereas later in treatment resuming the prior dose had similar or lower risk.

---

## Data

- **Source:** Linkage of provincial health administrative databases for all methadone recipients in BC between 01/01/2010 and 12/31/2022.
- **Access:** These are confidential, individual-level administrative data and **are not included in this repository**. Access is governed by the relevant BC data stewards and data-sharing agreements.
- **Design:** Retrospective cohort study emulating target trials, following the **TARGET** reporting checklist (Transparent Reporting of Observational Studies Emulating a Target trial).

---

## Repository Structure & Analysis Pipeline

The analysis proceeds from cohort/analytic-data construction in SAS to estimation, sensitivity analyses, and figure/table generation in R.

### 1. Analytic data creation (SAS)
- `3.0 analytic_data_met` — build the methadone analytic dataset.
- `5.0 Flow_chart` — cohort flow chart / attrition counts.
- `3.3 Missed_dose_hdPS` — high-dimensional propensity score (hdPS) data creation.

### 2. Core functions (R)
- `func_PLR` — helper functions for the pooled logistic regression (PLR) analysis.

### 3. Primary analysis (R)
- `Primary analysis_PSM` — propensity score / IPT weight estimation.
- `Primary analysis_with_PLR_HR` — pooled logistic regression models for hazard ratios (HRs).
- `Primary analysis_risk_with_PLR` — risk estimation with compatibility ("confidence") intervals.

### 4. Sensitivity analyses (R)
- `Sensitivity Analysis` and `Sensitivity Analysis_with_PLR` — risk-based sensitivity analyses.
- `Sensitivity Analysis_HR` — HR-based sensitivity analyses.

### 5. Diagnostics & subgroup analyses (R)
- `SMD plot` — standardized mean difference (balance) plots.
- `Analysis of Effect Modification` — effect-measure modification / subgroup analyses (e.g., by treatment-engagement duration and interruption length).

### 6. Figures & tables (R)
- `HR_graph_making codes_with_PLR` — main HR figures.
- `HR_graph_making codes SA_with_PLR` — sensitivity-analysis HR figures.
- `Figure 2`, `Figure 4`, `Figure 8` — manuscript figures.
- `Rate_Tables` — sample-size and event-rate tables.

> **Note:** File names in the repository may carry date stamps (e.g., `_Nov 12, 2025`) indicating the most recent version. Always use the latest dated version of each script.

---

## Reproducing the Analysis

1. Run the **SAS** scripts to build the analytic datasets and hdPS variables (Section 1 & 3.3 above).
2. Source `func_PLR` in **R**, then run the **primary analysis** scripts to estimate weights, HRs, and risks.
3. Run **sensitivity** and **effect-modification** scripts.
4. Generate **figures and tables** using the plotting and rate-table scripts.

**Software:** SAS (data construction) and R (estimation, inference, and visualization).

---

## Methods Summary

- **Confounding control:** A modified disjunctive cause criterion applied to an augmented covariate list from prior systematic reviews, including sociodemographic factors, disease severity/comorbidity, and OAT-related factors.
- **Estimand:** Risk difference (and hazard ratios) comparing each dose-adjustment strategy with resuming the prior dose (reference).
- **Inference:** 95% compatibility ("confidence") intervals.
- **Limitations:** Residual confounding may persist because clinical factors (e.g., intoxication or withdrawal severity) are unavailable in administrative data.

---

## Related Work

A parallel target trial emulation evaluates missed-dose management for **buprenorphine/naloxone** in BC using the same data infrastructure and analytic framework.

---

## Citation

If you use this code, please cite:

> Mondol MH, Min JE, Zanette M, Hossain MB, Kurz M, Bach P, Platt RW, Seaman S, Greenland S, McCandless LC, Karim ME, Gustafson P, Socías ME, Sutherland JM, Nosyk B. *Comparative effectiveness of missed dose management for methadone in British Columbia, Canada: population-based target trial emulation.*

---

## Funding

National Institutes of Health / National Institute on Drug Abuse, grant No. **R01-DA050629**.

---

## Contact

**Momenul Haque Mondol** — School of Population and Public Health, University of British Columbia, Vancouver, Canada.
