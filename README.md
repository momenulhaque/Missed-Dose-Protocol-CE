# Missed Dose Management for Opioid Agonist Treatment in British Columbia, Canada: Population-Based Target Trial Emulations

**Author:** Momenul Haque Mondol

This repository contains the SAS and R code used to emulate target trials evaluating alternative dose-adjustment strategies following opioid agonist treatment (OAT) interruptions (missed doses), using population-based linked administrative health data from British Columbia (BC), Canada.

---

## About This Repository

This repository is built to support the analyses for **two companion papers** that share the same data infrastructure, analytic framework, and target trial emulation approach, but evaluate two different medications for opioid use disorder (OUD):

1. **Methadone (MET)**
2. **Buprenorphine/Naloxone (BNX)**

### File & folder conventions

To keep the two analyses organized within a single repository:

- Scripts and files ending in **`_met`** belong to the **methadone paper**.
- Scripts and files ending in **`_bnx`** belong to the **buprenorphine/naloxone paper**.
- **Results are stored in separate folders** for each paper (e.g., a MET results folder and a BNX results folder), so outputs, figures, and tables do not overlap.
- Shared helper functions and general pipeline scripts (where applicable) are used by both analyses.

Each paper's specific documentation is provided in its own section below.

---

# Paper 1 — Methadone (MET)

### Title
**Comparative effectiveness of missed dose management for methadone in British Columbia, Canada: population-based target trial emulation**

### Authors
Momenul Haque Mondol, Jeong Eun Min, Michelle Zanette, Md. Belal Hossain, Megan Kurz, Paxton Bach, Robert William Platt, Shaun Seaman, Sander Greenland, Lawrence C. McCandless, Mohammad Ehsanul Karim, Paul Gustafson, Maria Eugenia Socías, Jason M Sutherland, Bohdan Nosyk.

### Overview
Using linked administrative health records from BC (2010–2022), we emulated two target trials — one for methadone interruptions of **3–4 days** and one for **5–14 days** — to estimate how alternative reinitiation strategies affect:

- **Overdose-related acute care visits** (hospital or ED encounters for poisoning; ICD-10 codes T40, T42.4, T43.6) — assessed at 1 month.
- **Treatment discontinuation** (no OAT dispensation for ≥15 consecutive days) — assessed at 12 months.

At reinitiation, the first dispensed dose was categorized relative to the previous dose as: **no change (reference)**, **<50% reduction**, **≥50% reduction** (or reduction to the guideline-recommended starting dose of ≤40 mg), or **dose increase**. Findings support **tailoring missed-dose protocols to treatment phase**: within the first 30 days of treatment, dose increases after interruptions were associated with lower risk of both outcomes, whereas later in treatment resuming the prior dose had similar or lower risk.

### Code (files ending in `_met`)
- `3.0 analytic_data_met` — build the methadone analytic dataset (SAS).
- `5.0 Flow_chart_..._met` — cohort flow chart / attrition counts (SAS).
- `3.3 Missed_dose_hdPS_met` — high-dimensional propensity score (hdPS) data creation (SAS).
- `func_PLR` — helper functions for pooled logistic regression (R; shared).
- `Primary analysis_..._PSM_met` — propensity score / IPT weight estimation (R).
- `Primary analysis_..._with_PLR_HR_met` — pooled logistic regression models for hazard ratios (R).
- `Primary analysis_..._risk_with_PLR_met` — risk estimation with compatibility intervals (R).
- `Sensitivity Analysis ..._met` — sensitivity analyses, risk- and HR-based (R).
- `SMD plot_..._met` — standardized mean difference (balance) plots (R).
- `Analysis of Effect Modification_met` — subgroup / effect-measure modification analyses (R).
- `HR_graph_making codes_..._met`, `Figure 2/4/8_met` — figures (R).
- `Rate_Tables_..._met` — sample-size and event-rate tables (R).

### Results
Stored in the **MET results folder**.

---

# Paper 2 — Buprenorphine/Naloxone (BNX)

### Title
**Comparative effectiveness of missed dose management for buprenorphine/naloxone in British Columbia, Canada: population-based target trial emulation**

### Authors
Momenul Haque Mondol, Kiera Tsakonas, Jeong Eun Min, Michelle Zanette, Md. Belal Hossain, Megan Kurz, Paxton Bach, Robert William Platt, Shaun Seaman, Sander Greenland, Lawrence C. McCandless, Mohammad Ehsanul Karim, Paul Gustafson, Maria Eugenia Socías, Jason M Sutherland, Bohdan Nosyk.

### Overview
This companion analysis applies the same target trial emulation framework to **buprenorphine/naloxone (BNX)** recipients in BC, evaluating how alternative dose-adjustment strategies following treatment interruptions affect overdose-related acute care visits and treatment discontinuation. BNX-specific clinical considerations (e.g., precipitated withdrawal risk) are addressed in the corresponding manuscript.

### Code (files ending in `_bnx`)
- `1.0 analytic_data_bnx` — build the BNX analytic dataset (SAS).
- `5.0 Flow_chart_..._bnx` — cohort flow chart / attrition counts (SAS).
- `..._Missed_dose_hdPS_bnx` — high-dimensional propensity score (hdPS) data creation (SAS).
- `func_PLR` — helper functions for pooled logistic regression (R; shared).
- `Primary analysis_..._PSM_bnx` — propensity score / IPT weight estimation (R).
- `Primary analysis_..._with_PLR_HR_bnx` — pooled logistic regression models for hazard ratios (R).
- `Primary analysis_..._risk_with_PLR_bnx` — risk estimation with compatibility intervals (R).
- `Sensitivity Analysis ..._bnx` — sensitivity analyses, risk- and HR-based (R).
- `SMD plot_..._bnx` — standardized mean difference (balance) plots (R).
- `Analysis of Effect Modification_bnx` — subgroup / effect-measure modification analyses (R).
- `HR_graph_making codes_..._bnx`, `Figure ..._bnx` — figures (R).
- `Rate_Tables_..._bnx` — sample-size and event-rate tables (R).

### Results
Stored in the **BNX results folder**.

---

## Data

- **Source:** Linkage of provincial health administrative databases for all methadone and BNX recipients in BC between 01/01/2010 and 12/31/2022.
- **Access:** These are confidential, individual-level administrative data and **are not included in this repository**. Access is governed by the relevant BC data stewards and data-sharing agreements.
- **Design:** Retrospective cohort studies emulating target trials, following the **TARGET** reporting checklist (Transparent Reporting of Observational Studies Emulating a Target trial).

---

## Reproducing the Analyses

For each paper (using the `_met` or `_bnx` scripts as appropriate):

1. Run the **SAS** scripts to build the analytic datasets and hdPS variables.
2. Source `func_PLR` in **R**, then run the **primary analysis** scripts to estimate weights, HRs, and risks.
3. Run **sensitivity** and **effect-modification** scripts.
4. Generate **figures and tables**; outputs are written to the paper's dedicated results folder.

**Software:** SAS (data construction) and R (estimation, inference, and visualization).

> **Note:** File names may carry date stamps (e.g., `_Nov 12, 2025`) indicating the most recent version. Always use the latest dated version of each script.

---

## Methods Summary

- **Confounding control:** A modified disjunctive cause criterion applied to an augmented covariate list from prior systematic reviews, including sociodemographic factors, disease severity/comorbidity, and OAT-related factors.
- **Estimand:** Risk difference (and hazard ratios) comparing each dose-adjustment strategy with resuming the prior dose (reference).
- **Inference:** 95% compatibility ("confidence") intervals.
- **Limitations:** Residual confounding may persist because clinical factors (e.g., intoxication or withdrawal severity) are unavailable in administrative data.

---

## Funding

National Institutes of Health / National Institute on Drug Abuse, grant No. **R01-DA050629**.

---

## Contact

**Momenul Haque Mondol** — School of Population and Public Health, University of British Columbia, Vancouver, Canada.
