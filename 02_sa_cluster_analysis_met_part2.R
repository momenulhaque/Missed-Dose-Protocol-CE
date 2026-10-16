# 1. PACKAGES AND PATHS ############################

library(tidyverse)
library(haven)
library(nnet)
library(splines)
library(sandwich)
library(mvtnorm)
library(tableone)
library(survey)
library(cobalt)
library(patchwork)
library(cowplot)

source("H:/GitHub/Missed-Dose-Protocols-PRI/00_functions.R")

DATA_DIR <- "E:/SAS Data/Missed dose protocol_v3/"
RES_DIR  <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/"

walk(c("PScore", "Risk", "Tables", "Figures"),
     ~ dir.create(file.path(RES_DIR, .x), recursive = TRUE, showWarnings = FALSE))





# 2. SHARED SETTING ############################

TRT_LEVELS <- c("no reversion", "1-<50%", ">=50%/<=40mg", "increase")

LHA_VAR  <- "clnt_lha"
HSDA_VAR <- "clnt_hsda_corrected"
HA_VAR   <- "clnt_ha"

MIN_PER_STRATEGY <- 1L
TRUNCATION_Q     <- 0.99
PS_EPSILON       <- 1e-8

LHA_NEIGHBOUR_FILE  <- file.path(DATA_DIR, "BC_LHA_nearby_within_HSDA.csv")
HSDA_NEIGHBOUR_FILE <- file.path(DATA_DIR, "BC_HSDA_nearby_within_HA.csv")
LHA_CROSSWALK_FILE  <- file.path(DATA_DIR, "BC_LHA_HSDA_HA_crosswalk.csv")
HSDA_CROSSWALK_FILE <- file.path(DATA_DIR, "BC_HSDA_HA_crosswalk.csv")

out_vars <- c(
  "sens_discont", "time_discont",
  "sens_od_acute", "time_od_acute"
)

clust_vars <- c(
  LHA_VAR, HSDA_VAR, HA_VAR,
  "studyid_rev", "episode", "time0"
)

baselinevars <- c(
  "sexM", "age_lt30", "age_30_39", "age_40_49", "age_ge50",
  "rural", "year", "po_preOUD", "phy_att_gr2", "social_assist",
  "ret_length", "home_grp_b", "incar_grp_bl", "drug_acute", "cci_gt0",
  "copd_asthma", "sed_psy_med", "po", "sud", "aud", "smd", "hcv", "tob",
  "udt", "ddose_prev", "clin_size_b", "virtual_care", "pss"
)

continuous_vars <- c("ddose_prev", "year", "ret_length")

STRATUM_LABELS <- c(
  dat2_early = "3-4 days, treatment duration <=30 days",
  dat2_late  = "3-4 days, treatment duration >30 days",
  dat3_early = "5-14 days, treatment duration <=30 days",
  dat3_late  = "5-14 days, treatment duration >30 days"
)




# 3. TREATMENT PREPROCESSING ############################

clean_trt_all <- function(dat) {
  
  dat |>
    mutate(
      trt_s456 = if_else(trt_s456 %in% c("1-<25%", "25-<50%"),
                         "1-<50%", trt_s456)
    ) |>
    rename(trt = trt_s456) |>
    filter(trt %in% TRT_LEVELS) |>
    mutate(
      trt = factor(trt, levels = TRT_LEVELS),
      prior_episode = as.numeric(episode - 1)
    ) |>
    group_by(studyid_rev, episode) |>
    mutate(
      prior_gap = cumsum(gap_days),
      prior_gap_num = row_number()
    ) |>
    ungroup()
}





# 4. LOAD AND VALIDATE CORRECTED GEOGRAPHY ############################

raw_all <- read_sas(file.path(DATA_DIR, "TT_analytic_prev_cl_sens_nw_m.sas7bdat"))

raw_all <- raw_all |> 
  mutate(
    clnt_hsda_original =
      clnt_hsda,
    
    hsda_recorded_numeric =
      suppressWarnings(
        as.integer(
          as.character(
            clnt_hsda_original
          )
        )
      ),
    
    hsda_from_lha =
      derive_hsda_from_lha(
        clnt_lha
      ),
    
    clnt_hsda_corrected = case_when(
      
      # Retain a valid recorded HSDA when it agrees with LHA.
      hsda_recorded_numeric %in%
        c(
          11L, 12L, 13L, 14L,
          21L, 22L, 23L,
          31L, 32L, 33L,
          41L, 42L, 43L,
          51L, 52L, 53L
        ) &
        hsda_recorded_numeric ==
        hsda_from_lha ~
        hsda_recorded_numeric,
      
      # Reconstruct HSDA from valid LHA when HSDA is 99 or missing.
      (
        hsda_recorded_numeric == 99L |
          is.na(hsda_recorded_numeric)
      ) &
        !is.na(hsda_from_lha) ~
        hsda_from_lha,
      
      # Preserve valid recorded HSDA if LHA is unavailable.
      hsda_recorded_numeric %in%
        c(
          11L, 12L, 13L, 14L,
          21L, 22L, 23L,
          31L, 32L, 33L,
          41L, 42L, 43L,
          51L, 52L, 53L
        ) &
        is.na(hsda_from_lha) ~
        hsda_recorded_numeric,
      
      # Do not silently resolve valid-code disagreements.
      TRUE ~
        NA_integer_
    ),
    
    hsda_correction_source = case_when(
      hsda_recorded_numeric %in%
        c(
          11L, 12L, 13L, 14L,
          21L, 22L, 23L,
          31L, 32L, 33L,
          41L, 42L, 43L,
          51L, 52L, 53L
        ) &
        hsda_recorded_numeric ==
        hsda_from_lha ~
        "Valid recorded HSDA",
      
      hsda_recorded_numeric == 99L &
        !is.na(hsda_from_lha) ~
        "HSDA 99 replaced using LHA hierarchy",
      
      is.na(hsda_recorded_numeric) &
        !is.na(hsda_from_lha) ~
        "Missing HSDA replaced using LHA hierarchy",
      
      is.na(hsda_from_lha) &
        hsda_recorded_numeric %in%
        c(
          11L, 12L, 13L, 14L,
          21L, 22L, 23L,
          31L, 32L, 33L,
          41L, 42L, 43L,
          51L, 52L, 53L
        ) ~
        "Valid recorded HSDA; LHA unavailable",
      
      !is.na(hsda_from_lha) &
        hsda_recorded_numeric !=
        hsda_from_lha &
        hsda_recorded_numeric != 99L ~
        "Unresolved disagreement",
      
      TRUE ~
        "Unresolved missing geography"
    )
  )

required_vars <- c(
  LHA_VAR, HSDA_VAR, HA_VAR,
  "trt_s456", "studyid_rev", "episode", "time0"
)

missing_vars <- setdiff(required_vars, names(raw_all))

if (length(missing_vars) > 0L) {
  stop(
    "Required variables are missing: ",
    paste(missing_vars, collapse = ", ")
  )
}


geography_validation <- raw_all |>
  summarise(
    Total_rows = n(),
    Missing_LHA_n = sum(is.na(.data[[LHA_VAR]])),
    Missing_HSDA_n = sum(is.na(.data[[HSDA_VAR]])),
    Missing_HSDA_pct = 100 * mean(is.na(.data[[HSDA_VAR]])),
    Missing_HA_n = sum(is.na(.data[[HA_VAR]])),
    Missing_HA_pct = 100 * mean(is.na(.data[[HA_VAR]]))
  )

print(geography_validation, width = Inf)

if (geography_validation$Missing_HSDA_n > 0L) {
  stop("The corrected HSDA variable still contains missing values.")
}

if (geography_validation$Missing_HA_n > 0L) {
  stop("The HA variable contains missing values.")
}





# 5. VALIDATE GEOGRAPHIC HIERARCHY ############################

geographic_hierarchy <- raw_all |>
  filter(
    !is.na(.data[[LHA_VAR]]),
    !is.na(.data[[HSDA_VAR]]),
    !is.na(.data[[HA_VAR]])
  ) |>
  transmute(
    lha_id = as.character(.data[[LHA_VAR]]),
    hsda_id = as.character(.data[[HSDA_VAR]]),
    ha_id = as.character(.data[[HA_VAR]])
  ) |>
  distinct()

invalid_lha_nesting <- geographic_hierarchy |>
  count(lha_id, name = "Number_of_HSDAs") |>
  filter(Number_of_HSDAs != 1L)

invalid_hsda_nesting <- geographic_hierarchy |>
  count(hsda_id, name = "Number_of_HAs") |>
  filter(Number_of_HAs != 1L)

if (nrow(invalid_lha_nesting) > 0L) {
  print(invalid_lha_nesting, n = Inf)
  stop("At least one LHA maps to more than one corrected HSDA.")
}

if (nrow(invalid_hsda_nesting) > 0L) {
  print(invalid_hsda_nesting, n = Inf)
  stop("At least one HSDA maps to more than one HA.")
}

write_csv(
  geographic_hierarchy,
  file.path(RES_DIR, "Tables", "validated_LHA_HSDA_HA_crosswalk.csv")
)



# 6. CREATE ANALYTIC STRATA ############################

analytic_all <- clean_trt_all(raw_all)

strata <- list(
  dat2_early = analytic_all |> filter(trial_gap %in% 3:4, ret_length <= 30),
  dat2_late  = analytic_all |> filter(trial_gap %in% 3:4, ret_length > 30),
  dat3_early = analytic_all |> filter(trial_cat == "5-14", ret_length <= 30),
  dat3_late  = analytic_all |> filter(trial_cat == "5-14", ret_length > 30)
)

rm(analytic_all)
gc()



# 7. PREPARE PS ANALYTIC SAMPLE ############################

get_valid_covariates <- function(dat) {
  
  candidates <- intersect(baselinevars, names(dat))
  
  candidates[
    vapply(
      dat[candidates],
      function(x) n_distinct(x[!is.na(x)]) > 1L,
      logical(1)
    )
  ]
}


prepare_ps_sample <- function(dat, stratum) {
  
  covs <- get_valid_covariates(dat)
  categorical_vars <- setdiff(covs, continuous_vars)
  
  d <- dat |>
    mutate(
      trt = factor(as.character(trt), levels = TRT_LEVELS),
      lha_id = as.character(.data[[LHA_VAR]]),
      hsda_id = as.character(.data[[HSDA_VAR]]),
      ha_id = as.character(.data[[HA_VAR]]),
      across(all_of(categorical_vars), factor)
    ) |>
    filter(
      !is.na(lha_id),
      !is.na(hsda_id),
      !is.na(ha_id),
      !is.na(trt),
      if_all(all_of(covs), ~ !is.na(.x))
    ) |>
    droplevels()
  
  sample_summary <- tibble(
    Stratum = stratum,
    Stratum_label = unname(STRATUM_LABELS[stratum]),
    Eligible_rows = nrow(dat),
    PS_complete_case_rows = nrow(d),
    PS_complete_case_pct = 100 * nrow(d) / nrow(dat),
    LHAs = n_distinct(d$lha_id),
    HSDAs = n_distinct(d$hsda_id),
    HAs = n_distinct(d$ha_id)
  )
  list(data = d, covariates = covs, summary = sample_summary)
}





# 8. READ GEOGRAPHIC CANDIDATE FILES ############################

required_geo_files <- c(
  LHA_NEIGHBOUR_FILE,
  HSDA_NEIGHBOUR_FILE,
  LHA_CROSSWALK_FILE,
  HSDA_CROSSWALK_FILE
)

missing_geo_files <- required_geo_files[!file.exists(required_geo_files)]

if (length(missing_geo_files) > 0L) {
  stop(
    "Missing geographic files: ",
    paste(missing_geo_files, collapse = ", ")
  )
}

lha_neighbours <- read_csv(
  LHA_NEIGHBOUR_FILE,
  col_types = cols(
    lha_id = col_character(),
    lha_name = col_character(),
    hsda_id = col_character(),
    nearby_rank = col_integer(),
    neighbor_lha_id = col_character(),
    neighbor_lha_name = col_character(),
    shares_boundary = col_logical(),
    boundary_distance_km = col_double()
  )
)

hsda_neighbours <- read_csv(
  HSDA_NEIGHBOUR_FILE,
  col_types = cols(
    hsda_id = col_character(),
    hsda_name = col_character(),
    ha_id = col_character(),
    nearby_rank = col_integer(),
    neighbor_hsda_id = col_character(),
    neighbor_hsda_name = col_character(),
    shares_boundary = col_logical(),
    boundary_distance_km = col_double()
  )
)

official_lha_crosswalk <- read_csv(
  LHA_CROSSWALK_FILE,
  col_types = cols(
    lha_id = col_character(),
    lha_name = col_character(),
    hsda_id = col_character(),
    ha_id = col_character()
  )
)

official_hsda_crosswalk <- read_csv(
  HSDA_CROSSWALK_FILE,
  col_types = cols(
    hsda_id = col_character(),
    hsda_name = col_character(),
    ha_id = col_character()
  )
)

stopifnot(
  all(lha_neighbours$lha_id != lha_neighbours$neighbor_lha_id),
  all(hsda_neighbours$hsda_id != hsda_neighbours$neighbor_hsda_id),
  !anyDuplicated(lha_neighbours[c("lha_id", "neighbor_lha_id")]),
  !anyDuplicated(hsda_neighbours[c("hsda_id", "neighbor_hsda_id")])
)




# 9. STANDARDIZE GEOGRAPHIC CANDIDATES ############################

lha_candidates <- lha_neighbours |>
  transmute(
    id = lha_id,
    neighbor = neighbor_lha_id,
    parent_id = hsda_id,
    nearby_rank,
    shares_boundary,
    distance_km = boundary_distance_km
  )

hsda_candidates <- hsda_neighbours |>
  transmute(
    id = hsda_id,
    neighbor = neighbor_hsda_id,
    parent_id = ha_id,
    nearby_rank,
    shares_boundary,
    distance_km = boundary_distance_km
  )



# 10. CALCULATE GROUP TREATMENT SUPPORT ############################

calculate_group_support <- function(members, strategy_cells, id_var,
                                    minimum = MIN_PER_STRATEGY) {
  
  counts <- strategy_cells |>
    filter(.data[[id_var]] %in% members) |>
    group_by(trt) |>
    summarise(n = sum(n), .groups = "drop") |>
    complete(
      trt = factor(TRT_LEVELS, levels = TRT_LEVELS),
      fill = list(n = 0L)
    )
  
  tibble(
    group_n = sum(counts$n),
    strategies_present = sum(counts$n > 0L),
    minimum_strategy_n = min(counts$n),
    adequate_support = all(counts$n >= minimum),
    no_reversion_n = counts$n[as.character(counts$trt) == "no reversion"][1],
    reduction_lt50_n = counts$n[as.character(counts$trt) == "1-<50%"][1],
    reduction_ge50_n = counts$n[as.character(counts$trt) == ">=50%/<=40mg"][1],
    increase_n = counts$n[as.character(counts$trt) == "increase"][1]
  )
}




# 11. CALCULATE BETWEEN-GROUP PROXIMITY ############################

calculate_group_proximity <- function(members_1, members_2, candidates) {
  
  pair_rows <- candidates |>
    filter(
      id %in% members_1 &
        neighbor %in% members_2
    )
  
  if (nrow(pair_rows) == 0L) {
    pair_rows <- candidates |>
      filter(
        id %in% members_2 &
          neighbor %in% members_1
      )
  }
  
  if (nrow(pair_rows) == 0L) {
    return(
      tibble(
        geographically_eligible = FALSE,
        shares_boundary = FALSE,
        minimum_distance_km = Inf,
        minimum_nearby_rank = Inf
      )
    )
  }
  
  tibble(
    geographically_eligible = TRUE,
    shares_boundary = any(pair_rows$shares_boundary %in% TRUE),
    minimum_distance_km = min(pair_rows$distance_km, na.rm = TRUE),
    minimum_nearby_rank = min(pair_rows$nearby_rank, na.rm = TRUE)
  )
}



# 12. SEQUENTIALLY MERGE GEOGRAPHIC UNITS ############################

make_empty_merge_history <- function() {
  
  tibble(
    parent_id = character(),
    geographic_level = character(),
    merge_step = integer(),
    group_1 = character(),
    group_2 = character(),
    new_group_id = character(),
    members_before_1 = character(),
    members_before_2 = character(),
    members_after = character(),
    members_n_after = integer(),
    group_n_after = integer(),
    strategies_present_after = integer(),
    minimum_strategy_n_after = integer(),
    adequate_support_after = logical(),
    added_episodes = integer(),
    shares_boundary = logical(),
    distance_km = double(),
    nearby_rank = integer()
  )
}


make_empty_merge_issue <- function() {
  
  tibble(
    parent_id = character(),
    geographic_level = character(),
    group_id = character(),
    members = character(),
    group_n = integer(),
    strategies_present = integer(),
    minimum_strategy_n = integer(),
    issue = character()
  )
}


merge_units_stepwise <- function(
    strategy_cells,
    candidates,
    id_var,
    parent_var,
    level,
    minimum = MIN_PER_STRATEGY) {
  
  required_cell_vars <- c(id_var, parent_var, "trt", "n")
  missing_cell_vars <- setdiff(required_cell_vars, names(strategy_cells))
  
  if (length(missing_cell_vars) > 0L) {
    stop(
      level,
      ": missing strategy-cell variables: ",
      paste(missing_cell_vars, collapse = ", ")
    )
  }
  
  required_candidate_vars <- c(
    "id", "neighbor", "parent_id", "nearby_rank",
    "shares_boundary", "distance_km"
  )
  
  missing_candidate_vars <- setdiff(
    required_candidate_vars,
    names(candidates)
  )
  
  if (length(missing_candidate_vars) > 0L) {
    stop(
      level,
      ": missing candidate variables: ",
      paste(missing_candidate_vars, collapse = ", ")
    )
  }
  
  parent_values <- strategy_cells |>
    distinct(.data[[parent_var]]) |>
    pull()
  
  if (length(parent_values) != 1L) {
    stop(
      level,
      ": strategy_cells must contain exactly one ",
      parent_var,
      "."
    )
  }
  
  parent_id <- as.character(parent_values[[1]])
  
  unit_ids <- strategy_cells |>
    distinct(.data[[id_var]]) |>
    pull() |>
    as.character() |>
    sort()
  
  if (length(unit_ids) == 0L) {
    stop(
      level,
      ": no geographic units were found for parent ",
      parent_id,
      "."
    )
  }
  
  parent_candidates <- candidates |>
    filter(
      parent_id == !!parent_id,
      id %in% unit_ids,
      neighbor %in% unit_ids,
      id != neighbor
    ) |>
    distinct(
      id,
      neighbor,
      .keep_all = TRUE
    )
  
  groups <- tibble(
    group_id = paste0(level, "_", unit_ids),
    members = map(unit_ids, ~ .x),
    merge_step = 0L
  )
  
  history <- make_empty_merge_history()
  issues <- make_empty_merge_issue()
  
  step <- 0L
  
  repeat {
    
    current_support <- map_dfr(
      seq_len(nrow(groups)),
      function(i) {
        
        calculate_group_support(
          members = groups$members[[i]],
          strategy_cells = strategy_cells,
          id_var = id_var,
          minimum = minimum
        ) |>
          mutate(
            group_id = groups$group_id[[i]],
            members = paste(
              sort(groups$members[[i]]),
              collapse = ";"
            ),
            n_members = length(
              groups$members[[i]]
            ),
            .before = 1
          )
      }
    )
    
    unsupported_groups <- current_support |>
      filter(!adequate_support) |>
      arrange(
        strategies_present,
        minimum_strategy_n,
        group_n,
        group_id
      )
    
    if (nrow(unsupported_groups) == 0L) {
      break
    }
    
    if (nrow(groups) == 1L) {
      
      issues <- bind_rows(
        issues,
        unsupported_groups |>
          transmute(
            parent_id = parent_id,
            geographic_level = level,
            group_id,
            members,
            group_n,
            strategies_present,
            minimum_strategy_n,
            issue =
              "All units in the parent area were combined but support remained inadequate"
          )
      )
      
      break
    }
    
    current_id <- unsupported_groups$group_id[[1]]
    current_i <- match(current_id, groups$group_id)
    current_members <- groups$members[[current_i]]
    
    other_indices <- setdiff(
      seq_len(nrow(groups)),
      current_i
    )
    
    candidate_scores <- map_dfr(
      other_indices,
      function(candidate_i) {
        
        candidate_id <- groups$group_id[[candidate_i]]
        candidate_members <- groups$members[[candidate_i]]
        
        proximity <- calculate_group_proximity(
          members_1 = current_members,
          members_2 = candidate_members,
          candidates = parent_candidates
        )
        
        if (!isTRUE(proximity$geographically_eligible[[1]])) {
          return(tibble())
        }
        
        merged_members <- sort(
          unique(
            c(
              current_members,
              candidate_members
            )
          )
        )
        
        merged_support <- calculate_group_support(
          members = merged_members,
          strategy_cells = strategy_cells,
          id_var = id_var,
          minimum = minimum
        )
        
        candidate_support <- calculate_group_support(
          members = candidate_members,
          strategy_cells = strategy_cells,
          id_var = id_var,
          minimum = minimum
        )
        
        tibble(
          candidate_group_id = candidate_id,
          candidate_index = candidate_i,
          candidate_members = list(candidate_members),
          merged_members = list(merged_members),
          adequate_support =
            merged_support$adequate_support,
          minimum_strategy_n =
            merged_support$minimum_strategy_n,
          strategies_present =
            merged_support$strategies_present,
          merged_group_n =
            merged_support$group_n,
          added_episodes =
            candidate_support$group_n,
          shares_boundary =
            proximity$shares_boundary,
          minimum_distance_km =
            proximity$minimum_distance_km,
          minimum_nearby_rank =
            proximity$minimum_nearby_rank
        )
      }
    )
    
    if (nrow(candidate_scores) == 0L) {
      
      current_row <- current_support |>
        filter(group_id == current_id)
      
      issues <- bind_rows(
        issues,
        current_row |>
          transmute(
            parent_id = parent_id,
            geographic_level = level,
            group_id,
            members,
            group_n,
            strategies_present,
            minimum_strategy_n,
            issue =
              "No eligible geographic neighbour was available within the parent area"
          )
      )
      
      break
    }
    
    selected <- candidate_scores |>
      arrange(
        desc(adequate_support),
        desc(minimum_strategy_n),
        desc(strategies_present),
        added_episodes,
        desc(shares_boundary),
        minimum_distance_km,
        minimum_nearby_rank,
        candidate_group_id
      ) |>
      slice(1)
    
    selected_id <- selected$candidate_group_id[[1]]
    selected_i <- selected$candidate_index[[1]]
    selected_members <- selected$candidate_members[[1]]
    merged_members <- selected$merged_members[[1]]
    
    step <- step + 1L
    
    new_group_id <- paste0(
      level,
      "_",
      parent_id,
      "_POOL_",
      str_pad(
        step,
        width = 2,
        pad = "0"
      )
    )
    
    history <- bind_rows(
      history,
      tibble(
        parent_id = parent_id,
        geographic_level = level,
        merge_step = step,
        group_1 = current_id,
        group_2 = selected_id,
        new_group_id = new_group_id,
        members_before_1 = paste(
          sort(current_members),
          collapse = ";"
        ),
        members_before_2 = paste(
          sort(selected_members),
          collapse = ";"
        ),
        members_after = paste(
          merged_members,
          collapse = ";"
        ),
        members_n_after =
          length(merged_members),
        group_n_after =
          selected$merged_group_n[[1]],
        strategies_present_after =
          selected$strategies_present[[1]],
        minimum_strategy_n_after =
          selected$minimum_strategy_n[[1]],
        adequate_support_after =
          selected$adequate_support[[1]],
        added_episodes =
          selected$added_episodes[[1]],
        shares_boundary =
          selected$shares_boundary[[1]],
        distance_km =
          selected$minimum_distance_km[[1]],
        nearby_rank =
          as.integer(
            selected$minimum_nearby_rank[[1]]
          )
      )
    )
    
    remove_indices <- sort(
      c(
        current_i,
        selected_i
      ),
      decreasing = TRUE
    )
    
    groups <- groups[
      -remove_indices,
      ,
      drop = FALSE
    ]
    
    groups <- bind_rows(
      groups,
      tibble(
        group_id = new_group_id,
        members = list(merged_members),
        merge_step = step
      )
    ) |>
      arrange(group_id)
  }
  
  final_support <- map_dfr(
    seq_len(nrow(groups)),
    function(i) {
      
      calculate_group_support(
        members = groups$members[[i]],
        strategy_cells = strategy_cells,
        id_var = id_var,
        minimum = minimum
      ) |>
        mutate(
          parent_id = parent_id,
          geographic_level = level,
          group_id = groups$group_id[[i]],
          members = paste(
            sort(groups$members[[i]]),
            collapse = ";"
          ),
          n_members =
            length(groups$members[[i]]),
          merge_step =
            groups$merge_step[[i]],
          .before = 1
        )
    }
  )
  
  crosswalk <- map_dfr(
    seq_len(nrow(groups)),
    function(i) {
      
      tibble(
        parent_id = parent_id,
        geographic_level = level,
        geographic_id =
          as.character(
            groups$members[[i]]
          ),
        group_id =
          groups$group_id[[i]],
        n_members =
          length(
            groups$members[[i]]
          ),
        merge_step =
          groups$merge_step[[i]]
      )
    }
  )
  
  unresolved <- final_support |>
    filter(!adequate_support)
  
  summary <- tibble(
    parent_id = parent_id,
    geographic_level = level,
    original_units = length(unit_ids),
    final_groups = nrow(groups),
    merges_completed = step,
    original_single_unit_groups =
      sum(final_support$n_members == 1L),
    pooled_groups =
      sum(final_support$n_members > 1L),
    unresolved_groups =
      nrow(unresolved),
    minimum_group_n =
      min(final_support$group_n),
    minimum_strategy_n =
      min(final_support$minimum_strategy_n),
    all_groups_supported =
      all(final_support$adequate_support)
  )
  
  list(
    groups = groups,
    crosswalk = crosswalk,
    support = final_support,
    history = history,
    issues = issues,
    unresolved = unresolved,
    summary = summary,
    all_supported =
      all(final_support$adequate_support)
  )
}



# 13. CONSTRUCT STRATEGY COUNT TABLES ############################

make_lha_strategy_cells <- function(dat) {
  
  dat |>
    count(ha_id, hsda_id, lha_id, trt, name = "n") |>
    complete(
      nesting(ha_id, hsda_id, lha_id),
      trt = factor(TRT_LEVELS, levels = TRT_LEVELS),
      fill = list(n = 0L)
    )
}

make_hsda_strategy_cells <- function(dat) {
  
  dat |>
    count(ha_id, hsda_id, trt, name = "n") |>
    complete(
      nesting(ha_id, hsda_id),
      trt = factor(TRT_LEVELS, levels = TRT_LEVELS),
      fill = list(n = 0L)
    )
}



# 14. CONSTRUCT HIERARCHICAL GEOGRAPHIC GROUPS ############################

construct_geographic_groups <- function(dat, stratum,
                                        minimum = MIN_PER_STRATEGY) {
  
  prepared <- prepare_ps_sample(dat, stratum)
  d <- prepared$data
  covs <- prepared$covariates
  
  lha_cells <- make_lha_strategy_cells(d)
  hsda_cells <- make_hsda_strategy_cells(d)
  
  lha_results <- lha_cells |>
    group_split(ha_id, hsda_id) |>
    map(
      ~ merge_units_stepwise(
        strategy_cells = .x,
        candidates = lha_candidates,
        id_var = "lha_id",
        parent_var = "hsda_id",
        level = "LHA",
        minimum = minimum
      )
    )
  
  lha_crosswalk <- map_dfr(lha_results, "crosswalk") |>
    transmute(
      hsda_id = parent_id,
      lha_id = geographic_id,
      lha_group_id = group_id
    ) |>
    left_join(
      d |>
        distinct(lha_id, hsda_id, ha_id),
      by = c("lha_id", "hsda_id")
    )
  
  lha_support <- map_dfr(lha_results, "support")
  lha_history <- map_dfr(lha_results, "history")
  
  unsupported_hsdas <- lha_support |>
    filter(!adequate_support) |>
    pull(parent_id) |>
    unique()
  
  hsda_results <- hsda_cells |>
    group_split(ha_id) |>
    map(
      ~ merge_units_stepwise(
        strategy_cells = .x,
        candidates = hsda_candidates,
        id_var = "hsda_id",
        parent_var = "ha_id",
        level = "HSDA",
        minimum = minimum
      )
    )
  
  hsda_crosswalk <- map_dfr(hsda_results, "crosswalk") |>
    transmute(
      ha_id = parent_id,
      hsda_id = geographic_id,
      hsda_group_id = group_id
    )
  
  hsda_support <- map_dfr(hsda_results, "support")
  hsda_history <- map_dfr(hsda_results, "history")
  
  hsda_group_size <- hsda_crosswalk |>
    count(hsda_group_id, name = "n_HSDAs")
  
  final_crosswalk <- lha_crosswalk |>
    left_join(hsda_crosswalk, by = c("ha_id", "hsda_id")) |>
    left_join(hsda_group_size, by = "hsda_group_id") |>
    mutate(
      final_geo_group = if_else(
        n_HSDAs > 1L,
        hsda_group_id,
        lha_group_id
      ),
      
      final_group_type = if_else(
        n_HSDAs > 1L,
        "Adjacent HSDAs pooled within HA",
        if_else(
          str_detect(lha_group_id, "_POOL_"),
          "Adjacent LHAs pooled within HSDA",
          "Original LHA retained"
        )
      )
    ) |>
    select(
      ha_id, hsda_id, lha_id,
      lha_group_id, hsda_group_id,
      final_geo_group, final_group_type
    )
  
  data_grouped <- d |>
    left_join(
      final_crosswalk,
      by = c("ha_id", "hsda_id", "lha_id")
    ) |>
    mutate(final_geo_group = factor(final_geo_group))
  
  if (any(is.na(data_grouped$final_geo_group))) {
    stop(stratum, ": missing geographic group after final crosswalk join.")
  }
  
  final_cells <- data_grouped |>
    count(final_geo_group, final_group_type, trt, name = "n") |>
    complete(
      nesting(final_geo_group, final_group_type),
      trt = factor(TRT_LEVELS, levels = TRT_LEVELS),
      fill = list(n = 0L)
    )
  
  final_support <- final_cells |>
    group_by(final_geo_group, final_group_type) |>
    summarise(
      final_group_n = sum(n),
      strategies_present = sum(n > 0L),
      minimum_strategy_n = min(n),
      adequate_support = all(n >= minimum),
      all_cells_ge2 = all(n >= 2L),
      all_cells_ge5 = all(n >= 5L),
      all_cells_ge10 = all(n >= 10L),
      no_reversion_n = n[as.character(trt) == "no reversion"][1],
      reduction_lt50_n = n[as.character(trt) == "1-<50%"][1],
      reduction_ge50_n = n[as.character(trt) == ">=50%/<=40mg"][1],
      increase_n = n[as.character(trt) == "increase"][1],
      .groups = "drop"
    ) |>
    mutate(
      Stratum = stratum,
      Minimum_required = minimum,
      .before = 1
    )
  
  unresolved <- final_support |>
    filter(!adequate_support)
  
  grouping_summary <- tibble(
    Stratum = stratum,
    Stratum_label = unname(STRATUM_LABELS[stratum]),
    Minimum_required = minimum,
    PS_rows = nrow(data_grouped),
    Original_LHAs = n_distinct(data_grouped$lha_id),
    HSDAs = n_distinct(data_grouped$hsda_id),
    HAs = n_distinct(data_grouped$ha_id),
    Final_geographic_groups = n_distinct(data_grouped$final_geo_group),
    
    Original_LHA_groups =
      sum(final_support$final_group_type == "Original LHA retained"),
    
    LHA_pooled_groups =
      sum(final_support$final_group_type ==
            "Adjacent LHAs pooled within HSDA"),
    
    HSDA_pooled_groups =
      sum(final_support$final_group_type ==
            "Adjacent HSDAs pooled within HA"),
    
    Groups_with_any_cell_lt2 =
      sum(!final_support$all_cells_ge2),
    
    Groups_with_any_cell_lt5 =
      sum(!final_support$all_cells_ge5),
    
    Groups_with_any_cell_lt10 =
      sum(!final_support$all_cells_ge10),
    
    Unresolved_groups = nrow(unresolved),
    Minimum_final_group_n = min(final_support$final_group_n),
    Minimum_final_strategy_cell_n = min(final_support$minimum_strategy_n),
    All_final_groups_supported = all(final_support$adequate_support)
  )
  
  list(
    data = data_grouped,
    covariates = covs,
    ps_sample_summary = prepared$summary,
    lha_crosswalk = lha_crosswalk,
    lha_support = lha_support,
    lha_merge_history = lha_history,
    hsda_crosswalk = hsda_crosswalk,
    hsda_support = hsda_support,
    hsda_merge_history = hsda_history,
    final_crosswalk = final_crosswalk,
    final_support_cells = final_cells,
    final_group_support = final_support,
    unresolved_groups = unresolved,
    grouping_summary = grouping_summary
  )
}



# 15. APPLY GROUPING TO ALL STRATA ############################

pooled_geography <- imap(
  strata,
  ~ construct_geographic_groups(
    dat = .x,
    stratum = .y,
    minimum = MIN_PER_STRATEGY
  )
)

View(strata_geo$dat2_early[, 150:197])


strata_geo$dat2_early |> count(lha_group_id)

strata_geo <- map(pooled_geography, "data")
grouping_summary <- map_dfr(pooled_geography, "grouping_summary")
ps_sample_summary <- map_dfr(pooled_geography, "ps_sample_summary")

print(grouping_summary, n = Inf, width = Inf)
print(ps_sample_summary, n = Inf, width = Inf)




# 16. VERIFY FINAL SUPPORT ############################

unresolved_geography <- imap_dfr(
  pooled_geography,
  ~ .x$unresolved_groups |>
    mutate(Stratum = .y, .before = 1)
)

if (nrow(unresolved_geography) > 0L) {
  
  write_csv(
    unresolved_geography,
    file.path(RES_DIR, "Tables", "unresolved_geographic_support.csv")
  )
  
  print(unresolved_geography, n = Inf, width = Inf)
  
  stop(
    "At least one geographic group remains unsupported after ",
    "sequential LHA and HSDA pooling."
  )
}

stopifnot(
  all(grouping_summary$All_final_groups_supported),
  all(map_lgl(strata_geo, ~ !any(is.na(.x$final_geo_group))))
)

cat(
  "\nSequential geographic grouping completed successfully.\n",
  "Every final geographic group contains all four treatment strategies.\n",
  sep = ""
)




# 17. SAVE GEOGRAPHIC GROUPING OUTPUTS ############################

final_crosswalk_all <- imap_dfr(
  pooled_geography,
  ~ .x$final_crosswalk |>
    mutate(Stratum = .y, .before = 1)
)

final_support_all <- imap_dfr(
  pooled_geography,
  ~ .x$final_group_support |>
    mutate(Stratum = .y, .before = 1)
)

final_cells_all <- imap_dfr(
  pooled_geography,
  ~ .x$final_support_cells |>
    mutate(Stratum = .y, .before = 1)
)

lha_merge_history_all <- imap_dfr(
  pooled_geography,
  ~ .x$lha_merge_history |>
    mutate(Stratum = .y, .before = 1)
)

hsda_merge_history_all <- imap_dfr(
  pooled_geography,
  ~ .x$hsda_merge_history |>
    mutate(Stratum = .y, .before = 1)
)

write_csv(
  ps_sample_summary,
  file.path(RES_DIR, "Tables", "geographic_PS_sample_summary.csv")
)

write_csv(
  grouping_summary,
  file.path(RES_DIR, "Tables", "geographic_grouping_summary.csv")
)

write_csv(
  final_crosswalk_all,
  file.path(RES_DIR, "Tables", "LHA_to_final_geographic_group_crosswalk.csv")
)

write_csv(
  final_support_all,
  file.path(RES_DIR, "Tables", "final_geographic_group_support.csv")
)

write_csv(
  final_cells_all,
  file.path(RES_DIR, "Tables", "final_geographic_group_strategy_counts.csv")
)

write_csv(
  lha_merge_history_all,
  file.path(RES_DIR, "Tables", "LHA_sequential_merge_history.csv")
)

write_csv(
  hsda_merge_history_all,
  file.path(RES_DIR, "Tables", "HSDA_sequential_merge_history.csv")
)

saveRDS(
  pooled_geography,
  file.path(RES_DIR, "PScore", "sequential_geographic_grouping_objects.rds")
)

saveRDS(
  strata_geo,
  file.path(RES_DIR, "PScore", "sequential_geographically_grouped_data.rds")
)







# 18. PS MODEL SETTING ############################

library(data.table)
library(Matrix)

PS_BASELINE_VARS <- unique(
  c(
    baselinevars,
    "prior_episode"
  )
)

PS_CONTINUOUS_VARS <- unique(
  c(
    continuous_vars,
    "prior_episode"
  )
)

PS_FIXED_EFFECT_VAR <- "final_geo_group"

PS_WEIGHT_TYPES <- c(
  "wt_stab",
  "wt_trunc",
  "ow_wt"
)

OUTCOME_COVARIATES <- c(
  "year",
  "ddose_prev",
  "prior_episode"
)

HORIZONS_CURVE <- c(
  seq(14, 350, by = 14),
  365
)

HORIZONS_CI <- sort(
  unique(
    c(
      HORIZONS_CURVE,
      30,
      60,
      90,
      180
    )
  )
)

BOOTSTRAP_SAMPLES <- 499L
REF_GROUP <- "no reversion"




OUTCOME_COVARIATES <- c(
  "year",
  "ddose_prev",
  "prior_episode"
)




# 19. SUPPORT AND UTILITY FUNCTIONS ############################

effective_sample_size <- function(w) {
  
  w <- w[is.finite(w)]
  
  if (length(w) == 0L || sum(w^2) == 0) {
    return(NA_real_)
  }
  
  sum(w)^2 / sum(w^2)
}


number_unique_nonmissing <- function(x) {
  
  n_distinct(
    x,
    na.rm = TRUE
  )
}


make_episode_id <- function(dat) {
  
  required_id_vars <- c(
    "studyid_rev",
    "episode",
    "time0"
  )
  
  missing_id_vars <- setdiff(
    required_id_vars,
    names(dat)
  )
  
  if (length(missing_id_vars) > 0L) {
    stop(
      "Missing episode-identifying variables: ",
      paste(missing_id_vars, collapse = ", ")
    )
  }
  
  interaction(
    dat$studyid_rev,
    dat$episode,
    dat$time0,
    drop = TRUE,
    lex.order = TRUE
  )
}





# 20. FIT FIXED-EFFECT PROPENSITY SCORE MODEL ############################

fit_fixed_geographic_ps <- function(
    grouped_result,
    stratum,
    exposure = "trt",
    baselinevars = PS_BASELINE_VARS,
    continuous_vars = PS_CONTINUOUS_VARS,
    fixed_effect_var = PS_FIXED_EFFECT_VAR,
    trunc_low = 0.01,
    trunc_high = 0.99,
    ps_epsilon = PS_EPSILON) {
  
  dat <- grouped_result$data
  
  required_variables <- unique(
    c(
      exposure,
      baselinevars,
      fixed_effect_var,
      out_vars,
      clust_vars
    )
  )
  
  missing_variables <- setdiff(
    required_variables,
    names(dat)
  )
  
  if (length(missing_variables) > 0L) {
    stop(
      stratum,
      ": missing variables: ",
      paste(missing_variables, collapse = ", ")
    )
  }
  
  invalid_continuous <- setdiff(
    continuous_vars,
    baselinevars
  )
  
  if (length(invalid_continuous) > 0L) {
    stop(
      stratum,
      ": continuous variables absent from baselinevars: ",
      paste(invalid_continuous, collapse = ", ")
    )
  }
  
  model_data <- dat |>
    select(
      all_of(required_variables)
    ) |>
    mutate(
      across(
        all_of(
          setdiff(
            baselinevars,
            continuous_vars
          )
        ),
        factor
      ),
      
      "{fixed_effect_var}" :=
        factor(.data[[fixed_effect_var]]),
      
      "{exposure}" :=
        factor(
          as.character(.data[[exposure]]),
          levels = TRT_LEVELS
        )
    ) |>
    filter(
      if_all(
        all_of(
          c(
            exposure,
            baselinevars,
            fixed_effect_var
          )
        ),
        ~ !is.na(.x)
      )
    ) |>
    droplevels()
  
  if (nrow(model_data) == 0L) {
    stop(
      stratum,
      ": no complete observations remained."
    )
  }
  
  if (!identical(levels(model_data[[exposure]]), TRT_LEVELS)) {
    stop(
      stratum,
      ": not all four treatment categories remained."
    )
  }
  
  valid_baselinevars <- baselinevars[
    vapply(
      model_data[baselinevars],
      number_unique_nonmissing,
      integer(1)
    ) > 1L
  ]
  
  removed_baselinevars <- setdiff(
    baselinevars,
    valid_baselinevars
  )
  
  if (length(removed_baselinevars) > 0L) {
    warning(
      stratum,
      ": nonvarying variables excluded: ",
      paste(removed_baselinevars, collapse = ", ")
    )
  }
  
  geographic_support <- model_data |>
    count(
      final_geo_group,
      .data[[exposure]],
      name = "n"
    ) |>
    complete(
      final_geo_group,
      !!sym(exposure) :=
        factor(
          TRT_LEVELS,
          levels = TRT_LEVELS
        ),
      fill = list(n = 0L)
    ) |>
    group_by(final_geo_group) |>
    mutate(
      strategies_present =
        sum(n > 0L),
      
      minimum_strategy_n =
        min(n),
      
      all_four_present =
        all(n >= 1L)
    ) |>
    ungroup()
  
  unsupported_groups <- geographic_support |>
    filter(!all_four_present)
  
  if (nrow(unsupported_groups) > 0L) {
    
    print(
      unsupported_groups,
      n = Inf,
      width = Inf
    )
    
    stop(
      stratum,
      ": complete-case PS data contain unsupported geographic groups."
    )
  }
  
  ps_covariates <- c(
    valid_baselinevars,
    fixed_effect_var
  )
  
  formula_ps <- reformulate(
    ps_covariates,
    response = exposure
  )
  
  max_weights <- max(
    100000L,
    10L *
      length(TRT_LEVELS) *
      (
        length(ps_covariates) +
          nlevels(model_data[[fixed_effect_var]])
      )
  )
  
  ps_fit <- multinom(
    formula_ps,
    data = model_data,
    trace = FALSE,
    Hess = TRUE,
    maxit = 2000,
    MaxNWts = max_weights
  )
  
  if (ps_fit$convergence != 0L) {
    warning(
      stratum,
      ": multinom returned convergence code ",
      ps_fit$convergence,
      "."
    )
  }
  
  ps_probs <- as.matrix(
    predict(
      ps_fit,
      newdata = model_data,
      type = "probs"
    )
  )
  
  if (is.null(dim(ps_probs))) {
    stop(
      stratum,
      ": multinomial probability matrix was not returned."
    )
  }
  
  if (is.null(colnames(ps_probs)) &&
      ncol(ps_probs) == length(TRT_LEVELS)) {
    colnames(ps_probs) <- TRT_LEVELS
  }
  
  missing_probability_columns <- setdiff(
    TRT_LEVELS,
    colnames(ps_probs)
  )
  
  if (length(missing_probability_columns) > 0L) {
    stop(
      stratum,
      ": missing probability columns: ",
      paste(missing_probability_columns, collapse = ", ")
    )
  }
  
  ps_probs <- ps_probs[
    ,
    TRT_LEVELS,
    drop = FALSE
  ]
  
  ps_probs <- pmax(
    ps_probs,
    ps_epsilon
  )
  
  ps_probs <- ps_probs /
    rowSums(ps_probs)
  
  received_column <- match(
    as.character(
      model_data[[exposure]]
    ),
    TRT_LEVELS
  )
  
  model_data$ps <- ps_probs[
    cbind(
      seq_len(nrow(model_data)),
      received_column
    )
  ]
  
  marginal_probabilities <- prop.table(
    table(
      model_data[[exposure]]
    )
  )
  
  model_data$stab_numerator <- as.numeric(
    marginal_probabilities[
      as.character(
        model_data[[exposure]]
      )
    ]
  )
  
  model_data$wt_stab <-
    model_data$stab_numerator /
    model_data$ps
  
  weight_limits <- quantile(
    model_data$wt_stab,
    probs = c(
      trunc_low,
      trunc_high
    ),
    na.rm = TRUE,
    names = FALSE
  )
  
  model_data$wt_trunc <- pmin(
    pmax(
      model_data$wt_stab,
      weight_limits[[1]]
    ),
    weight_limits[[2]]
  )
  
  harmonic_tilt <- 1 /
    rowSums(
      1 / ps_probs
    )
  
  model_data$ow_wt <-
    harmonic_tilt /
    model_data$ps
  
  invalid_weights <- vapply(
    model_data[PS_WEIGHT_TYPES],
    function(x) {
      any(!is.finite(x)) ||
        any(x < 0)
    },
    logical(1)
  )
  
  if (any(invalid_weights)) {
    stop(
      stratum,
      ": invalid weights detected in ",
      paste(
        names(invalid_weights)[invalid_weights],
        collapse = ", "
      )
    )
  }
  
  propensity_diagnostics <- model_data |>
    group_by(
      .data[[exposure]]
    ) |>
    summarise(
      Stratum = stratum,
      n = n(),
      minimum_ps = min(ps),
      q1_ps = quantile(ps, 0.25),
      median_ps = median(ps),
      q3_ps = quantile(ps, 0.75),
      maximum_ps = max(ps),
      ps_below_0_01_n = sum(ps < 0.01),
      ps_below_0_01_pct = 100 * mean(ps < 0.01),
      ps_below_0_05_n = sum(ps < 0.05),
      ps_below_0_05_pct = 100 * mean(ps < 0.05),
      .groups = "drop"
    )
  
  weight_diagnostics <- model_data |>
    select(
      all_of(
        c(
          exposure,
          PS_WEIGHT_TYPES
        )
      )
    ) |>
    pivot_longer(
      cols = all_of(PS_WEIGHT_TYPES),
      names_to = "weight_type",
      values_to = "weight"
    ) |>
    group_by(
      .data[[exposure]],
      weight_type
    ) |>
    summarise(
      Stratum = stratum,
      n = n(),
      minimum = min(weight),
      q1 = quantile(weight, 0.25),
      median = median(weight),
      mean = mean(weight),
      q3 = quantile(weight, 0.75),
      p95 = quantile(weight, 0.95),
      p99 = quantile(weight, 0.99),
      maximum = max(weight),
      sum_weights = sum(weight),
      effective_sample_size =
        effective_sample_size(weight),
      efficiency_pct =
        100 *
        effective_sample_size(weight) /
        n(),
      .groups = "drop"
    )
  
  overall_weight_diagnostics <- model_data |>
    select(
      all_of(PS_WEIGHT_TYPES)
    ) |>
    pivot_longer(
      cols = everything(),
      names_to = "weight_type",
      values_to = "weight"
    ) |>
    group_by(weight_type) |>
    summarise(
      Stratum = stratum,
      n = n(),
      minimum = min(weight),
      q1 = quantile(weight, 0.25),
      median = median(weight),
      mean = mean(weight),
      q3 = quantile(weight, 0.75),
      p95 = quantile(weight, 0.95),
      p99 = quantile(weight, 0.99),
      maximum = max(weight),
      sum_weights = sum(weight),
      effective_sample_size =
        effective_sample_size(weight),
      efficiency_pct =
        100 *
        effective_sample_size(weight) /
        n(),
      .groups = "drop"
    )
  
  list(
    stratum = stratum,
    model = ps_fit,
    formula = formula_ps,
    convergence = ps_fit$convergence,
    data_weights = model_data,
    ps_probabilities = ps_probs,
    propensity_diagnostics = propensity_diagnostics,
    weight_diagnostics = weight_diagnostics,
    overall_weight_diagnostics = overall_weight_diagnostics,
    geographic_support = geographic_support,
    valid_baselinevars = valid_baselinevars,
    removed_baselinevars = removed_baselinevars,
    stabilized_weight_truncation_limits = tibble(
      lower_percentile = trunc_low,
      upper_percentile = trunc_high,
      lower_weight_limit = weight_limits[[1]],
      upper_weight_limit = weight_limits[[2]]
    )
  )
}

# 21. FIT PS MODELS ############################

fixed_cluster_ps_fits <- imap(
  pooled_geography,
  ~ fit_fixed_geographic_ps(
    grouped_result = .x,
    stratum = .y
  )
)

walk2(
  fixed_cluster_ps_fits,
  names(fixed_cluster_ps_fits),
  ~ saveRDS(
    .x,
    file.path(
      RES_DIR,
      "PScore",
      paste0(
        "PS_fixed_geography_",
        .y,
        ".rds"
      )
    )
  )
)

saveRDS(
  fixed_cluster_ps_fits,
  file.path(
    RES_DIR,
    "PScore",
    "PS_fixed_geography_all_strata.rds"
  )
)






# 22. SAVE PS DIAGNOSTICS ############################

ps_diagnostics_all <- map_dfr(
  fixed_cluster_ps_fits,
  "propensity_diagnostics"
)

weight_diagnostics_all <- map_dfr(
  fixed_cluster_ps_fits,
  "weight_diagnostics"
)

overall_weight_diagnostics_all <- map_dfr(
  fixed_cluster_ps_fits,
  "overall_weight_diagnostics"
)

truncation_limits_all <- imap_dfr(
  fixed_cluster_ps_fits,
  ~ .x$stabilized_weight_truncation_limits |>
    mutate(
      Stratum = .y,
      .before = 1
    )
)

write_csv(
  ps_diagnostics_all,
  file.path(
    RES_DIR,
    "Tables",
    "PS_fixed_geography_propensity_diagnostics.csv"
  )
)

write_csv(
  weight_diagnostics_all,
  file.path(
    RES_DIR,
    "Tables",
    "PS_fixed_geography_weight_diagnostics.csv"
  )
)

write_csv(
  overall_weight_diagnostics_all,
  file.path(
    RES_DIR,
    "Tables",
    "PS_fixed_geography_overall_weight_diagnostics.csv"
  )
)

write_csv(
  truncation_limits_all,
  file.path(
    RES_DIR,
    "Tables",
    "PS_fixed_geography_truncation_limits.csv"
  )
)

print(
  ps_diagnostics_all,
  n = Inf,
  width = Inf
)

print(
  overall_weight_diagnostics_all,
  n = Inf,
  width = Inf
)






# 23. POSITIVE-DEFINITE COVARIANCE ############################

make_positive_definite <- function(V) {
  
  V <- as.matrix(V)
  V <- (V + t(V)) / 2
  
  eigenvalues <- eigen(
    V,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  if (
    all(is.finite(eigenvalues)) &&
    min(eigenvalues) > 0
  ) {
    return(V)
  }
  
  as.matrix(
    Matrix::nearPD(
      V,
      corr = FALSE,
      keepDiag = TRUE
    )$mat
  )
}



# 24. FIT POOLED LOGISTIC MODEL ############################

fit_pooled_logistic <- function(
    data,
    trt_var,
    covariates,
    weight_var,
    person_cluster_var,
    interval_df = 3L) {
  
  required_vars <- unique(
    c(
      "event",
      "interval",
      "stop",
      trt_var,
      covariates,
      weight_var,
      person_cluster_var,
      "episode_id"
    )
  )
  
  missing_vars <- setdiff(
    required_vars,
    names(data)
  )
  
  if (length(missing_vars) > 0L) {
    stop(
      "Missing pooled-logistic variables: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  valid_covariates <- intersect(
    covariates,
    names(data)
  )
  
  valid_covariates <- valid_covariates[
    vapply(
      data[valid_covariates],
      number_unique_nonmissing,
      integer(1)
    ) > 1L
  ]
  
  formula_terms <- c(
    trt_var,
    paste0(
      "ns(interval, df = ",
      interval_df,
      ")"
    ),
    valid_covariates
  )
  
  model_formula <- reformulate(
    formula_terms,
    response = "event"
  )
  
  model <- glm(
    model_formula,
    data = data,
    family = quasibinomial(),
    weights = data[[weight_var]],
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  
  covariance <- sandwich::vcovCL(
    model,
    cluster = data[[person_cluster_var]],
    type = "HC0",
    cadjust = TRUE,
    fix = FALSE
  )
  
  covariance <- make_positive_definite(
    covariance
  )
  
  list(
    model = model,
    covariance = covariance,
    coefficients = coef(model),
    valid_covariates = valid_covariates,
    formula = model_formula
  )
}




# 25. STANDARDIZE POOLED LOGISTIC RISKS ############################

standardize_pooled_risk <- function(
    data,
    model,
    coefficients,
    trt_var,
    treatment_levels,
    horizons) {
  
  dt <- as.data.table(copy(data))
  
  setorderv(
    dt,
    c("episode_id", "stop")
  )
  
  model_terms <- delete.response(
    terms(model)
  )
  
  model_matrices <- map(
    treatment_levels,
    function(level) {
      
      newdata <- copy(dt)
      
      newdata[[trt_var]] <- factor(
        level,
        levels = treatment_levels
      )
      
      model.matrix(
        model_terms,
        data = newdata,
        xlev = model$xlevels,
        contrasts.arg = model$contrasts
      )
    }
  )
  
  names(model_matrices) <- treatment_levels
  
  result <- map_dfr(
    treatment_levels,
    function(level) {
      
      linear_predictor <- as.numeric(
        model_matrices[[level]] %*%
          coefficients
      )
      
      dt[, hazard := plogis(linear_predictor)]
      
      dt[
        ,
        survival := cumprod(1 - hazard),
        by = episode_id
      ]
      
      map_dfr(
        horizons,
        function(horizon) {
          
          episode_risk <- dt[
            stop <= horizon,
            .(
              risk = 1 - survival[.N]
            ),
            by = episode_id
          ]
          
          tibble(
            trt = level,
            horizon = horizon,
            risk = mean(
              episode_risk$risk,
              na.rm = TRUE
            )
          )
        }
      )
    }
  )
  
  result
}



# 26. PARAMETRIC COEFFICIENT SIMULATION ############################

simulate_pooled_risks <- function(
    data,
    pooled_fit,
    trt_var = "trt",
    horizons,
    samples = BOOTSTRAP_SAMPLES,
    ref_group = REF_GROUP,
    seed = 20260922L,
    return_simulation_data = FALSE) {
  
  model <- pooled_fit$model
  coefficients <- pooled_fit$coefficients
  covariance <- pooled_fit$covariance
  
  treatment_levels <- levels(
    data[[trt_var]]
  )
  
  if (is.null(treatment_levels)) {
    treatment_levels <- unique(
      as.character(
        data[[trt_var]]
      )
    )
  }
  
  if (!ref_group %in% treatment_levels) {
    stop(
      "Reference treatment was not observed: ",
      ref_group
    )
  }
  
  missing_horizons <- setdiff(
    horizons,
    unique(data$stop)
  )
  
  if (length(missing_horizons) > 0L) {
    warning(
      "The following horizons do not correspond to exact stop times: ",
      paste(missing_horizons, collapse = ", "),
      ". Risks will use the latest interval ending on or before each horizon."
    )
  }
  
  point_risks <- standardize_pooled_risk(
    data = data,
    model = model,
    coefficients = coefficients,
    trt_var = trt_var,
    treatment_levels = treatment_levels,
    horizons = horizons
  )
  
  set.seed(seed)
  
  coefficient_draws <- mvtnorm::rmvnorm(
    n = samples,
    mean = coefficients,
    sigma = covariance
  )
  
  simulated_risks <- map_dfr(
    seq_len(samples),
    function(simulation) {
      
      standardize_pooled_risk(
        data = data,
        model = model,
        coefficients = coefficient_draws[simulation, ],
        trt_var = trt_var,
        treatment_levels = treatment_levels,
        horizons = horizons
      ) |>
        mutate(
          simulation = simulation,
          .before = 1
        )
    }
  )
  
  risk_intervals <- simulated_risks |>
    group_by(
      trt,
      horizon
    ) |>
    summarise(
      risk_lower = quantile(
        risk,
        0.025,
        na.rm = TRUE
      ),
      
      risk_upper = quantile(
        risk,
        0.975,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    left_join(
      point_risks,
      by = c(
        "trt",
        "horizon"
      )
    ) |>
    select(
      trt,
      horizon,
      risk,
      risk_lower,
      risk_upper
    )
  
  point_contrasts <- point_risks |>
    group_by(horizon) |>
    mutate(
      reference_risk = risk[
        trt == ref_group
      ][1]
    ) |>
    filter(
      trt != ref_group
    ) |>
    mutate(
      risk_difference =
        risk - reference_risk,
      
      risk_ratio = if_else(
        reference_risk > 0,
        risk / reference_risk,
        NA_real_
      )
    ) |>
    ungroup()
  
  simulated_contrasts <- simulated_risks |>
    group_by(
      simulation,
      horizon
    ) |>
    mutate(
      reference_risk = risk[
        trt == ref_group
      ][1]
    ) |>
    filter(
      trt != ref_group
    ) |>
    mutate(
      risk_difference =
        risk - reference_risk,
      
      risk_ratio = if_else(
        reference_risk > 0,
        risk / reference_risk,
        NA_real_
      )
    ) |>
    ungroup()
  
  contrast_intervals <- simulated_contrasts |>
    group_by(
      trt,
      horizon
    ) |>
    summarise(
      risk_difference_lower = quantile(
        risk_difference,
        0.025,
        na.rm = TRUE
      ),
      
      risk_difference_upper = quantile(
        risk_difference,
        0.975,
        na.rm = TRUE
      ),
      
      risk_ratio_lower = quantile(
        risk_ratio,
        0.025,
        na.rm = TRUE
      ),
      
      risk_ratio_upper = quantile(
        risk_ratio,
        0.975,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    left_join(
      point_contrasts |>
        select(
          trt,
          horizon,
          reference_risk,
          risk,
          risk_difference,
          risk_ratio
        ),
      by = c(
        "trt",
        "horizon"
      )
    ) |>
    select(
      trt,
      horizon,
      reference_risk,
      risk,
      risk_difference,
      risk_difference_lower,
      risk_difference_upper,
      risk_ratio,
      risk_ratio_lower,
      risk_ratio_upper
    )
  
  result <- list(
    risks = risk_intervals,
    contrasts = contrast_intervals,
    reference_group = ref_group,
    simulations = samples,
    seed = seed
  )
  
  if (return_simulation_data) {
    result$simulation_data <- simulated_risks
  }
  
  result
}





# 27. RUN ONE OUTCOME ANALYSIS ############################

run_one_pooled_outcome <- function(
    ps_fit,
    time_var,
    event_var,
    horizons,
    weight_var = "wt_trunc",
    samples = BOOTSTRAP_SAMPLES,
    seed = 20260922L) {
  
  subject_data <- ps_fit$data_weights |>
    mutate(
      episode_id =
        make_episode_id(
          ps_fit$data_weights
        )
    )
  
  long_data <- expand_person_period(
    subject_data,
    time_var = time_var,
    event_var = event_var,
    interval_length = 14
  ) |>
    mutate(
      trt = factor(
        as.character(trt),
        levels = TRT_LEVELS
      )
    )
  
  pooled_fit <- fit_pooled_logistic(
    data = long_data,
    trt_var = "trt",
    covariates = OUTCOME_COVARIATES,
    weight_var = weight_var,
    person_cluster_var =
      "studyid_rev",
    interval_df = 3L
  )
  
  estimates <- simulate_pooled_risks(
    data = long_data,
    pooled_fit = pooled_fit,
    trt_var = "trt",
    horizons = horizons,
    samples = samples,
    ref_group = REF_GROUP,
    seed = seed,
    return_simulation_data = FALSE
  )
  
  list(
    pooled_model =
      pooled_fit$model,
    
    pooled_covariance =
      pooled_fit$covariance,
    
    formula =
      pooled_fit$formula,
    
    risks =
      estimates$risks,
    
    contrasts =
      estimates$contrasts,
    
    reference_group =
      estimates$reference_group
  )
}




# 28. RUN ALL OUTCOMES WITHIN ONE STRATUM ############################

run_stratum_outcomes <- function(
    ps_fit,
    stratum,
    horizons = HORIZONS_CI,
    samples = BOOTSTRAP_SAMPLES) {
  
  outcome_definitions <- tribble(
    ~outcome, ~time_var, ~event_var,
    "overdose_acute", "time_od_acute", "sens_od_acute",
    "discontinuation", "time_discont", "sens_discont"
  )
  
  results <- pmap(
    outcome_definitions,
    function(outcome, time_var, event_var) {
      
      cat(
        "\n",
        stratum,
        ": ",
        outcome,
        "\n",
        sep = ""
      )
      
      run_one_pooled_outcome(
        ps_fit = ps_fit,
        time_var = time_var,
        event_var = event_var,
        horizons = horizons,
        weight_var = "wt_trunc",
        samples = samples,
        seed =
          20260922L +
          match(
            outcome,
            outcome_definitions$outcome
          )
      )
    }
  )
  
  names(results) <-
    outcome_definitions$outcome
  
  results
}







# 29. ESTIMATE STANDARDIZED RISKS ############################

risk_results <- imap(
  fixed_cluster_ps_fits,
  ~ run_stratum_outcomes(
    ps_fit = .x,
    stratum = .y,
    horizons = HORIZONS_CI,
    samples = BOOTSTRAP_SAMPLES
  )
)

walk2(
  risk_results,
  names(risk_results),
  ~ saveRDS(
    .x,
    file.path(
      RES_DIR,
      "Risk",
      paste0(
        "pooled_logistic_fixed_geography_",
        .y,
        ".rds"
      )
    )
  )
)

saveRDS(
  risk_results,
  file.path(
    RES_DIR,
    "Risk",
    "pooled_logistic_fixed_geography_all_strata.rds"
  )
)





dt = risk_results$dat2_early$overdose_acute$contrasts
dt[dt$horizon==30, ]


dt = risk_results$dat3_early$overdose_acute$contrasts
dt[dt$horizon==30, ]
dt = risk_results$dat3_early$overdose_acute$risks
dt[dt$horizon==30, ]






dt = risk_results$dat2_early$discontinuation$contrasts
dt[dt$horizon==365, ]


dt = risk_results$dat3_early$discontinuation$contrasts
dt[dt$horizon==365, ]
dt = risk_results$dat3_early$discontinuation$risks
dt[dt$horizon==365, ]




## 29.0 RISK TABLE ##################

CLUSTER_OUTCOMES <- tribble(
  ~Outcome, ~Outcome_label, ~Horizon_days,
  "overdose_acute", "Overdose-related acute care", 30L,
  "discontinuation", "Treatment discontinuation", 365L
)

TREATMENT_LABELS <- c(
  "no reversion" = "No change",
  "1-<50%" = "<50% reduction",
  ">=50%/<=40mg" = ">=50% reduction",
  "increase" = "Dose increase"
)

REF_GROUP <- "no reversion"

cluster_episode_counts <- imap_dfr(
  fixed_cluster_ps_fits,
  function(ps_fit, stratum) {
    
    ps_fit$data_weights |>
      count(
        trt,
        name = "Episodes_at_time_zero"
      ) |>
      transmute(
        Stratum = stratum,
        trt = as.character(trt),
        Episodes_at_time_zero
      )
  }
)

print(
  cluster_episode_counts,
  n = Inf,
  width = Inf
)

## 29.1 EXTRACT STANDARDIZED RISKS ############################

cluster_risks <- imap_dfr(
  risk_results,
  function(stratum_results, stratum) {
    
    pmap_dfr(
      CLUSTER_OUTCOMES,
      function(Outcome, Outcome_label, Horizon_days) {
        
        if (!Outcome %in% names(stratum_results)) {
          stop(
            stratum,
            ": outcome object not found: ",
            Outcome
          )
        }
        
        risk_rows <- stratum_results[[Outcome]]$risks |>
          filter(
            horizon == Horizon_days
          )
        
        if (nrow(risk_rows) == 0L) {
          stop(
            stratum,
            ": no ",
            Outcome,
            " risk estimates found at ",
            Horizon_days,
            " days."
          )
        }
        
        risk_rows |>
          transmute(
            Stratum = stratum,
            Outcome,
            Outcome_label,
            Horizon_days,
            trt = as.character(trt),
            Risk = risk,
            Risk_lower = risk_lower,
            Risk_upper = risk_upper
          )
      }
    )
  }
)

print(
  cluster_risks,
  n = Inf,
  width = Inf
)


## 29.2 EXTRACT RISK DIFFERENCES ############################

cluster_risk_differences <- imap_dfr(
  risk_results,
  function(stratum_results, stratum) {
    
    pmap_dfr(
      CLUSTER_OUTCOMES,
      function(Outcome, Outcome_label, Horizon_days) {
        
        contrast_rows <- stratum_results[[Outcome]]$contrasts |>
          filter(
            horizon == Horizon_days
          )
        
        if (nrow(contrast_rows) == 0L) {
          stop(
            stratum,
            ": no ",
            Outcome,
            " risk differences found at ",
            Horizon_days,
            " days."
          )
        }
        
        contrast_rows |>
          transmute(
            Stratum = stratum,
            Outcome,
            Horizon_days,
            trt = as.character(trt),
            Reference_risk = reference_risk,
            Risk_difference = risk_difference,
            Risk_difference_lower = risk_difference_lower,
            Risk_difference_upper = risk_difference_upper
          )
      }
    )
  }
)
## 29.3 ADD REFERENCE-TREATMENT ROWS ############################

cluster_reference_rows <- cluster_risks |>
  filter(
    trt == REF_GROUP
  ) |>
  transmute(
    Stratum,
    Outcome,
    Horizon_days,
    trt,
    Reference_risk = Risk,
    Risk_difference = 0,
    Risk_difference_lower = NA_real_,
    Risk_difference_upper = NA_real_
  )

cluster_risk_differences <- bind_rows(
  cluster_reference_rows,
  cluster_risk_differences
) |>
  distinct(
    Stratum,
    Outcome,
    Horizon_days,
    trt,
    .keep_all = TRUE
  )


## 29.4 COMBINE CLUSTER-ADJUSTED RESULTS ############################

cluster_risk_results <- cluster_risks |>
  left_join(
    cluster_risk_differences,
    by = c(
      "Stratum",
      "Outcome",
      "Horizon_days",
      "trt"
    )
  ) |>
  left_join(
    cluster_episode_counts,
    by = c(
      "Stratum",
      "trt"
    )
  ) |>
  mutate(
    Stratum_label = unname(
      STRATUM_LABELS[Stratum]
    ),
    
    Dose_adjustment = unname(
      TREATMENT_LABELS[trt]
    ),
    
    Treatment_order = match(
      trt,
      TRT_LEVELS
    ),
    
    Risk_pct = 100 * Risk,
    Risk_lower_pct = 100 * Risk_lower,
    Risk_upper_pct = 100 * Risk_upper,
    
    Risk_difference_pp =
      100 * Risk_difference,
    
    Risk_difference_lower_pp =
      100 * Risk_difference_lower,
    
    Risk_difference_upper_pp =
      100 * Risk_difference_upper
  ) |>
  arrange(
    factor(
      Stratum,
      levels = names(STRATUM_LABELS)
    ),
    factor(
      Outcome,
      levels = CLUSTER_OUTCOMES$Outcome
    ),
    Treatment_order
  )

print(
  cluster_risk_results,
  n = Inf,
  width = Inf
)


## 29.5 FORMAT CLUSTER-ADJUSTED RESULTS ############################

cluster_risk_table <- cluster_risk_results |>
  mutate(
    Reporting_time = case_when(
      Horizon_days == 30L ~ "1 month",
      Horizon_days == 365L ~ "12 months",
      TRUE ~ paste0(Horizon_days, " days")
    ),
    
    Risk_95CI = sprintf(
      "%.2f (%.2f, %.2f)",
      Risk_pct,
      Risk_lower_pct,
      Risk_upper_pct
    ),
    
    Risk_difference_95CI = case_when(
      trt == REF_GROUP ~
        "Reference",
      
      TRUE ~
        sprintf(
          "%.2f (%.2f, %.2f)",
          Risk_difference_pp,
          Risk_difference_lower_pp,
          Risk_difference_upper_pp
        )
    )
  ) |>
  transmute(
    Stratum,
    Stratum_label,
    Outcome,
    Outcome_label,
    Horizon_days,
    Reporting_time,
    Dose_adjustment,
    Episodes_at_time_zero,
    Risk = Risk,
    Risk_lower = Risk_lower,
    Risk_upper = Risk_upper,
    Risk_difference = Risk_difference,
    Risk_difference_lower = Risk_difference_lower,
    Risk_difference_upper = Risk_difference_upper,
    Risk_pct,
    Risk_lower_pct,
    Risk_upper_pct,
    Risk_difference_pp,
    Risk_difference_lower_pp,
    Risk_difference_upper_pp,
    Risk_95CI,
    Risk_difference_95CI
  )

print(
  cluster_risk_table,
  n = Inf,
  width = Inf
)

write_csv(
  cluster_risk_table,
  file.path(
    RES_DIR,
    "Tables",
    "cluster_adjusted_risk_results_v1.csv"
  )
)

## 29.6 CREATE PUBLICATION TABLE ############################

cluster_risk_table_publication <- cluster_risk_table |>
  transmute(
    `Analytic stratum` =
      Stratum_label,
    
    Outcome =
      Outcome_label,
    
    `Dose-adjustment strategy` =
      Dose_adjustment,
    
    `Episodes at time zero, n` =
      format(
        Episodes_at_time_zero,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      ),
    
    `Reporting time` =
      Reporting_time,
    
    `Standardized risk, % (95% CI)` =
      Risk_95CI,
    
    `Risk difference, percentage points (95% CI)` =
      Risk_difference_95CI
  )

print(
  cluster_risk_table_publication,
  n = Inf,
  width = Inf
)

write_csv(
  cluster_risk_table_publication,
  file.path(
    RES_DIR,
    "Tables",
    "cluster_adjusted_risk_results_publication_v1.csv"
  )
)


# 30. GROUPING SUMMARY SETTING ############################

REPORT_CELL_THRESHOLDS <- c(2L, 5L, 10L)

STRATEGY_LABELS <- c(
  "no reversion",
  "1-<50%",
  ">=50%/<=40mg",
  "increase"
)


# 31. SAFE SUMMARY FUNCTIONS ############################

safe_min <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  min(x)
}


safe_max <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  max(x)
}


safe_mean <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  mean(x)
}


safe_sd <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  
  sd(x)
}


safe_quantile <- function(x, probability) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  as.numeric(
    quantile(
      x,
      probs = probability,
      names = FALSE,
      type = 2
    )
  )
}


format_five_number <- function(
    minimum,
    q1,
    median,
    q3,
    maximum,
    digits = 0L) {
  
  format_string <- paste0(
    "%.",
    digits,
    "f [%.",
    digits,
    "f, %.",
    digits,
    "f, %.",
    digits,
    "f] %.",
    digits,
    "f"
  )
  
  sprintf(
    format_string,
    minimum,
    q1,
    median,
    q3,
    maximum
  )
}


format_n_pct <- function(
    n,
    denominator,
    digits = 1L) {
  
  if (is.na(denominator) || denominator == 0L) {
    return(NA_character_)
  }
  
  sprintf(
    paste0(
      "%s (%.",
      digits,
      "f%%)"
    ),
    format(
      n,
      big.mark = ",",
      scientific = FALSE
    ),
    100 * n / denominator
  )
}



# 32. ORIGINAL LHA CHARACTERISTICS ############################

summarize_original_lhas <- function(
    grouped_result,
    stratum) {
  
  d <- grouped_result$data
  
  original_lha <- d |>
    count(
      ha_id,
      hsda_id,
      lha_id,
      name = "Episodes"
    )
  
  treatment_support <- d |>
    count(
      ha_id,
      hsda_id,
      lha_id,
      trt,
      name = "n"
    ) |>
    complete(
      nesting(
        ha_id,
        hsda_id,
        lha_id
      ),
      
      trt = factor(
        TRT_LEVELS,
        levels = TRT_LEVELS
      ),
      
      fill = list(
        n = 0L
      )
    ) |>
    group_by(
      ha_id,
      hsda_id,
      lha_id
    ) |>
    summarise(
      Strategies_present =
        sum(n > 0L),
      
      Minimum_strategy_n =
        min(n),
      
      All_four_present =
        all(n >= 1L),
      
      All_cells_ge2 =
        all(n >= 2L),
      
      All_cells_ge5 =
        all(n >= 5L),
      
      All_cells_ge10 =
        all(n >= 10L),
      
      .groups = "drop"
    )
  
  original_lha <- original_lha |>
    left_join(
      treatment_support,
      by = c(
        "ha_id",
        "hsda_id",
        "lha_id"
      )
    )
  
  tibble(
    Stratum = stratum,
    Stratum_label =
      unname(
        STRATUM_LABELS[stratum]
      ),
    
    Episodes =
      nrow(d),
    
    Original_LHAs =
      nrow(original_lha),
    
    HSDAs =
      n_distinct(
        original_lha$hsda_id
      ),
    
    HAs =
      n_distinct(
        original_lha$ha_id
      ),
    
    Episodes_per_LHA_mean =
      mean(
        original_lha$Episodes
      ),
    
    Episodes_per_LHA_SD =
      sd(
        original_lha$Episodes
      ),
    
    Episodes_per_LHA_min =
      min(
        original_lha$Episodes
      ),
    
    Episodes_per_LHA_Q1 =
      safe_quantile(
        original_lha$Episodes,
        0.25
      ),
    
    Episodes_per_LHA_median =
      median(
        original_lha$Episodes
      ),
    
    Episodes_per_LHA_Q3 =
      safe_quantile(
        original_lha$Episodes,
        0.75
      ),
    
    Episodes_per_LHA_max =
      max(
        original_lha$Episodes
      ),
    
    Singleton_LHAs_n =
      sum(
        original_lha$Episodes == 1L
      ),
    
    Singleton_LHAs_pct =
      100 *
      mean(
        original_lha$Episodes == 1L
      ),
    
    Strategies_per_LHA_mean =
      mean(
        original_lha$Strategies_present
      ),
    
    Strategies_per_LHA_SD =
      sd(
        original_lha$Strategies_present
      ),
    
    Strategies_per_LHA_min =
      min(
        original_lha$Strategies_present
      ),
    
    Strategies_per_LHA_Q1 =
      safe_quantile(
        original_lha$Strategies_present,
        0.25
      ),
    
    Strategies_per_LHA_median =
      median(
        original_lha$Strategies_present
      ),
    
    Strategies_per_LHA_Q3 =
      safe_quantile(
        original_lha$Strategies_present,
        0.75
      ),
    
    Strategies_per_LHA_max =
      max(
        original_lha$Strategies_present
      ),
    
    LHAs_with_1_strategy_n =
      sum(
        original_lha$Strategies_present == 1L
      ),
    
    LHAs_with_2_strategies_n =
      sum(
        original_lha$Strategies_present == 2L
      ),
    
    LHAs_with_3_strategies_n =
      sum(
        original_lha$Strategies_present == 3L
      ),
    
    LHAs_with_4_strategies_n =
      sum(
        original_lha$Strategies_present == 4L
      ),
    
    LHAs_with_4_strategies_pct =
      100 *
      mean(
        original_lha$All_four_present
      ),
    
    LHAs_with_any_cell_lt2_n =
      sum(
        !original_lha$All_cells_ge2
      ),
    
    LHAs_with_any_cell_lt5_n =
      sum(
        !original_lha$All_cells_ge5
      ),
    
    LHAs_with_any_cell_lt10_n =
      sum(
        !original_lha$All_cells_ge10
      ),
    
    Episodes_in_all_four_strategy_LHAs_n =
      sum(
        original_lha$Episodes[
          original_lha$All_four_present
        ]
      ),
    
    Episodes_in_all_four_strategy_LHAs_pct =
      100 *
      sum(
        original_lha$Episodes[
          original_lha$All_four_present
        ]
      ) /
      sum(
        original_lha$Episodes
      )
  )
}


original_lha_summary <- imap_dfr(
  pooled_geography,
  summarize_original_lhas
)

print(
  original_lha_summary,
  n = Inf,
  width = Inf
)



# 33. FINAL GEOGRAPHIC GROUP CHARACTERISTICS ############################

summarize_final_groups <- function(
    grouped_result,
    stratum) {
  
  d <- grouped_result$data
  
  final_group_size <- d |>
    group_by(
      final_geo_group,
      final_group_type
    ) |>
    summarise(
      Episodes = n(),
      Original_LHAs =
        n_distinct(
          lha_id
        ),
      HSDAs =
        n_distinct(
          hsda_id
        ),
      HAs =
        n_distinct(
          ha_id
        ),
      .groups = "drop"
    )
  
  final_support <- grouped_result$
    final_group_support
  
  final_group_size <- final_group_size |>
    left_join(
      final_support |>
        select(
          final_geo_group,
          strategies_present,
          minimum_strategy_n,
          adequate_support,
          all_cells_ge2,
          all_cells_ge5,
          all_cells_ge10,
          no_reversion_n,
          reduction_lt50_n,
          reduction_ge50_n,
          increase_n
        ),
      by = "final_geo_group"
    )
  
  tibble(
    Stratum = stratum,
    Final_geographic_groups =
      nrow(final_group_size),
    
    Original_LHA_groups_n =
      sum(
        final_group_size$
          final_group_type ==
          "Original LHA retained"
      ),
    
    LHA_pooled_groups_n =
      sum(
        final_group_size$
          final_group_type ==
          "Adjacent LHAs pooled within HSDA"
      ),
    
    HSDA_pooled_groups_n =
      sum(
        final_group_size$
          final_group_type ==
          "Adjacent HSDAs pooled within HA"
      ),
    
    Original_LHAs_per_group_mean =
      mean(
        final_group_size$Original_LHAs
      ),
    
    Original_LHAs_per_group_SD =
      sd(
        final_group_size$Original_LHAs
      ),
    
    Original_LHAs_per_group_min =
      min(
        final_group_size$Original_LHAs
      ),
    
    Original_LHAs_per_group_Q1 =
      safe_quantile(
        final_group_size$Original_LHAs,
        0.25
      ),
    
    Original_LHAs_per_group_median =
      median(
        final_group_size$Original_LHAs
      ),
    
    Original_LHAs_per_group_Q3 =
      safe_quantile(
        final_group_size$Original_LHAs,
        0.75
      ),
    
    Original_LHAs_per_group_max =
      max(
        final_group_size$Original_LHAs
      ),
    
    Episodes_per_final_group_mean =
      mean(
        final_group_size$Episodes
      ),
    
    Episodes_per_final_group_SD =
      sd(
        final_group_size$Episodes
      ),
    
    Episodes_per_final_group_min =
      min(
        final_group_size$Episodes
      ),
    
    Episodes_per_final_group_Q1 =
      safe_quantile(
        final_group_size$Episodes,
        0.25
      ),
    
    Episodes_per_final_group_median =
      median(
        final_group_size$Episodes
      ),
    
    Episodes_per_final_group_Q3 =
      safe_quantile(
        final_group_size$Episodes,
        0.75
      ),
    
    Episodes_per_final_group_max =
      max(
        final_group_size$Episodes
      ),
    
    Minimum_strategy_cell_n =
      min(
        final_group_size$
          minimum_strategy_n
      ),
    
    Median_minimum_strategy_cell_n =
      median(
        final_group_size$
          minimum_strategy_n
      ),
    
    Maximum_minimum_strategy_cell_n =
      max(
        final_group_size$
          minimum_strategy_n
      ),
    
    Groups_with_any_cell_lt2_n =
      sum(
        !final_group_size$
          all_cells_ge2
      ),
    
    Groups_with_any_cell_lt2_pct =
      100 *
      mean(
        !final_group_size$
          all_cells_ge2
      ),
    
    Groups_with_any_cell_lt5_n =
      sum(
        !final_group_size$
          all_cells_ge5
      ),
    
    Groups_with_any_cell_lt5_pct =
      100 *
      mean(
        !final_group_size$
          all_cells_ge5
      ),
    
    Groups_with_any_cell_lt10_n =
      sum(
        !final_group_size$
          all_cells_ge10
      ),
    
    Groups_with_any_cell_lt10_pct =
      100 *
      mean(
        !final_group_size$
          all_cells_ge10
      ),
    
    All_groups_have_4_strategies =
      all(
        final_group_size$
          adequate_support
      )
  )
}


final_group_summary <- imap_dfr(
  pooled_geography,
  summarize_final_groups
)

print(
  final_group_summary,
  n = Inf,
  width = Inf
)




# 34. SEQUENTIAL MERGE PROCESS SUMMARY ############################

summarize_merge_process <- function(
    grouped_result,
    stratum) {
  
  lha_history <- grouped_result$
    lha_merge_history
  
  hsda_history <- grouped_result$
    hsda_merge_history
  
  tibble(
    Stratum = stratum,
    
    LHA_merge_steps =
      nrow(
        lha_history
      ),
    
    HSDAs_with_LHA_merges =
      n_distinct(
        lha_history$parent_id
      ),
    
    Original_LHAs_involved_in_merges =
      if (nrow(lha_history) == 0L) {
        0L
      } else {
        length(
          unique(
            unlist(
              str_split(
                c(
                  lha_history$
                    members_before_1,
                  lha_history$
                    members_before_2
                ),
                ";"
              )
            )
          )
        )
      },
    
    Maximum_LHA_merge_step =
      if (nrow(lha_history) == 0L) {
        0L
      } else {
        max(
          lha_history$merge_step
        )
      },
    
    HSDA_merge_steps =
      nrow(
        hsda_history
      ),
    
    HAs_with_HSDA_merges =
      n_distinct(
        hsda_history$parent_id
      ),
    
    Original_HSDAs_involved_in_merges =
      if (nrow(hsda_history) == 0L) {
        0L
      } else {
        length(
          unique(
            unlist(
              str_split(
                c(
                  hsda_history$
                    members_before_1,
                  hsda_history$
                    members_before_2
                ),
                ";"
              )
            )
          )
        )
      },
    
    Maximum_HSDA_merge_step =
      if (nrow(hsda_history) == 0L) {
        0L
      } else {
        max(
          hsda_history$merge_step
        )
      },
    
    LHA_merges_sharing_boundary_n =
      sum(
        lha_history$
          shares_boundary %in% TRUE
      ),
    
    LHA_merges_sharing_boundary_pct =
      if (nrow(lha_history) == 0L) {
        NA_real_
      } else {
        100 *
          mean(
            lha_history$
              shares_boundary %in% TRUE
          )
      },
    
    HSDA_merges_sharing_boundary_n =
      sum(
        hsda_history$
          shares_boundary %in% TRUE
      ),
    
    HSDA_merges_sharing_boundary_pct =
      if (nrow(hsda_history) == 0L) {
        NA_real_
      } else {
        100 *
          mean(
            hsda_history$
              shares_boundary %in% TRUE
          )
      }
  )
}


merge_process_summary <- imap_dfr(
  pooled_geography,
  summarize_merge_process
)

print(
  merge_process_summary,
  n = Inf,
  width = Inf
)



# 35. COMBINE CLUSTER ANALYSIS SUMMARY ############################

cluster_grouping_summary <- original_lha_summary |>
  left_join(
    final_group_summary,
    by = "Stratum"
  ) |>
  left_join(
    merge_process_summary,
    by = "Stratum"
  ) |>
  arrange(
    factor(
      Stratum,
      levels = names(STRATUM_LABELS)
    )
  )

print(
  cluster_grouping_summary,
  n = Inf,
  width = Inf
)

write_csv(
  cluster_grouping_summary,
  file.path(
    RES_DIR,
    "Tables",
    "cluster_grouping_summary_detailed.csv"
  )
)




# 36. CREATE COMPACT GROUPING TABLE ############################

format_five_number <- function(minimum, q1, median, q3, maximum, digits = 0L) {
  
  output <- rep(NA_character_, length(minimum))
  
  complete_values <- complete.cases(
    minimum,
    q1,
    median,
    q3,
    maximum
  )
  
  format_string <- paste0(
    "%.", digits, "f [%.", digits, "f, %.",
    digits, "f, %.", digits, "f] %.", digits, "f"
  )
  
  output[complete_values] <- sprintf(
    format_string,
    minimum[complete_values],
    q1[complete_values],
    median[complete_values],
    q3[complete_values],
    maximum[complete_values]
  )
  
  output
}


format_n_pct <- function(n, denominator, digits = 1L) {
  
  output <- rep(NA_character_, length(n))
  
  valid <- !is.na(n) &
    !is.na(denominator) &
    denominator > 0
  
  format_string <- paste0(
    "%s (%.",
    digits,
    "f%%)"
  )
  
  output[valid] <- sprintf(
    format_string,
    format(
      n[valid],
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    100 * n[valid] / denominator[valid]
  )
  
  output
}


required_summary_columns <- c(
  "Stratum",
  "Stratum_label",
  "Episodes",
  "Original_LHAs",
  "HSDAs",
  "HAs",
  "Episodes_per_LHA_min",
  "Episodes_per_LHA_Q1",
  "Episodes_per_LHA_median",
  "Episodes_per_LHA_Q3",
  "Episodes_per_LHA_max",
  "Strategies_per_LHA_min",
  "Strategies_per_LHA_Q1",
  "Strategies_per_LHA_median",
  "Strategies_per_LHA_Q3",
  "Strategies_per_LHA_max",
  "LHAs_with_4_strategies_n",
  "Episodes_in_all_four_strategy_LHAs_n",
  "Final_geographic_groups",
  "Original_LHA_groups_n",
  "LHA_pooled_groups_n",
  "HSDA_pooled_groups_n",
  "Original_LHAs_per_group_min",
  "Original_LHAs_per_group_Q1",
  "Original_LHAs_per_group_median",
  "Original_LHAs_per_group_Q3",
  "Original_LHAs_per_group_max",
  "Episodes_per_final_group_min",
  "Episodes_per_final_group_Q1",
  "Episodes_per_final_group_median",
  "Episodes_per_final_group_Q3",
  "Episodes_per_final_group_max",
  "Minimum_strategy_cell_n",
  "Groups_with_any_cell_lt2_n",
  "Groups_with_any_cell_lt5_n",
  "Groups_with_any_cell_lt10_n",
  "LHA_merge_steps",
  "HSDA_merge_steps",
  "All_groups_have_4_strategies"
)

missing_summary_columns <- setdiff(
  required_summary_columns,
  names(cluster_grouping_summary)
)

if (length(missing_summary_columns) > 0L) {
  stop(
    "Missing variables in cluster_grouping_summary: ",
    paste(missing_summary_columns, collapse = ", ")
  )
}


cluster_grouping_compact <- cluster_grouping_summary |>
  mutate(
    Episodes_per_LHA_summary = format_five_number(
      Episodes_per_LHA_min,
      Episodes_per_LHA_Q1,
      Episodes_per_LHA_median,
      Episodes_per_LHA_Q3,
      Episodes_per_LHA_max
    ),
    
    Strategies_per_LHA_summary = format_five_number(
      Strategies_per_LHA_min,
      Strategies_per_LHA_Q1,
      Strategies_per_LHA_median,
      Strategies_per_LHA_Q3,
      Strategies_per_LHA_max
    ),
    
    LHAs_with_4_strategies = format_n_pct(
      LHAs_with_4_strategies_n,
      Original_LHAs
    ),
    
    Episodes_in_4_strategy_LHAs = format_n_pct(
      Episodes_in_all_four_strategy_LHAs_n,
      Episodes
    ),
    
    Final_group_types = sprintf(
      "%s / %s / %s",
      Original_LHA_groups_n,
      LHA_pooled_groups_n,
      HSDA_pooled_groups_n
    ),
    
    LHAs_per_final_group_summary = format_five_number(
      Original_LHAs_per_group_min,
      Original_LHAs_per_group_Q1,
      Original_LHAs_per_group_median,
      Original_LHAs_per_group_Q3,
      Original_LHAs_per_group_max
    ),
    
    Episodes_per_final_group_summary = format_five_number(
      Episodes_per_final_group_min,
      Episodes_per_final_group_Q1,
      Episodes_per_final_group_median,
      Episodes_per_final_group_Q3,
      Episodes_per_final_group_max
    ),
    
    Groups_with_cell_lt2 = format_n_pct(
      Groups_with_any_cell_lt2_n,
      Final_geographic_groups
    ),
    
    Groups_with_cell_lt5 = format_n_pct(
      Groups_with_any_cell_lt5_n,
      Final_geographic_groups
    ),
    
    Groups_with_cell_lt10 = format_n_pct(
      Groups_with_any_cell_lt10_n,
      Final_geographic_groups
    ),
    
    Sequential_merge_steps = sprintf(
      "%s LHA / %s HSDA",
      LHA_merge_steps,
      HSDA_merge_steps
    ),
    
    Final_support_status = if_else(
      All_groups_have_4_strategies,
      "Yes",
      "No",
      missing = NA_character_
    )
  ) |>
  transmute(
    Stratum = Stratum_label,
    
    Episodes = format(
      Episodes,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    
    `Original LHAs` =
      Original_LHAs,
    
    HSDAs =
      HSDAs,
    
    HAs =
      HAs,
    
    `Episodes per original LHA, min [Q1, median, Q3] max` =
      Episodes_per_LHA_summary,
    
    `Strategies per original LHA, min [Q1, median, Q3] max` =
      Strategies_per_LHA_summary,
    
    `Original LHAs with all 4 strategies, n (%)` =
      LHAs_with_4_strategies,
    
    `Episodes in all-4-strategy LHAs, n (%)` =
      Episodes_in_4_strategy_LHAs,
    
    `Final geographic groups` =
      Final_geographic_groups,
    
    `Final groups: original LHA / LHA pool / HSDA pool` =
      Final_group_types,
    
    `Original LHAs per final group, min [Q1, median, Q3] max` =
      LHAs_per_final_group_summary,
    
    `Episodes per final group, min [Q1, median, Q3] max` =
      Episodes_per_final_group_summary,
    
    `Minimum strategy-cell count` =
      Minimum_strategy_cell_n,
    
    `Final groups with a strategy cell <2, n (%)` =
      Groups_with_cell_lt2,
    
    `Final groups with a strategy cell <5, n (%)` =
      Groups_with_cell_lt5,
    
    `Final groups with a strategy cell <10, n (%)` =
      Groups_with_cell_lt10,
    
    `Sequential merge steps: LHA / HSDA` =
      Sequential_merge_steps,
    
    `All final groups contained 4 strategies` =
      Final_support_status
  )


print(
  cluster_grouping_compact,
  n = Inf,
  width = Inf
)

cluster_grouping_compact$`Final groups: original LHA / LHA pool / HSDA pool`
write_csv(
  cluster_grouping_compact,
  file.path(
    RES_DIR,
    "Tables",
    "cluster_grouping_summary_compact.csv"
  )
)




# 37. FINAL CLUSTER TREATMENT DISTRIBUTION ############################

final_cluster_treatment_distribution <- imap_dfr(
  pooled_geography,
  function(grouped_result, stratum) {
    
    grouped_result$data |>
      count(
        final_geo_group,
        final_group_type,
        trt,
        name = "Episodes"
      ) |>
      complete(
        nesting(
          final_geo_group,
          final_group_type
        ),
        
        trt = factor(
          TRT_LEVELS,
          levels = TRT_LEVELS
        ),
        
        fill = list(
          Episodes = 0L
        )
      ) |>
      group_by(
        final_geo_group,
        final_group_type
      ) |>
      mutate(
        Total_cluster_episodes =
          sum(Episodes),
        
        Treatment_pct =
          100 *
          Episodes /
          Total_cluster_episodes
      ) |>
      ungroup() |>
      mutate(
        Stratum = stratum,
        
        Stratum_label =
          unname(
            STRATUM_LABELS[stratum]
          ),
        
        Treatment_label = recode(
          as.character(trt),
          "no reversion" = "No change",
          "1-<50%" = "<50% reduction",
          ">=50%/<=40mg" = ">=50% reduction",
          "increase" = "Dose increase"
        ),
        
        .before = 1
      ) |>
      arrange(
        final_geo_group,
        trt
      )
  }
)

print(
  final_cluster_treatment_distribution,
  n = 100,
  width = Inf
)

write_csv(
  final_cluster_treatment_distribution,
  file.path(
    RES_DIR,
    "Tables",
    "final_cluster_treatment_distribution_long.csv"
  )
)






# 41. FINAL CLUSTER TREATMENT DISTRIBUTION WIDE ############################

final_cluster_treatment_distribution_wide <-
  final_cluster_treatment_distribution |>
  mutate(
    Treatment_variable = recode(
      as.character(trt),
      "no reversion" = "No_change",
      "1-<50%" = "Reduction_lt50",
      ">=50%/<=40mg" = "Reduction_ge50",
      "increase" = "Dose_increase"
    )
  ) |>
  select(
    Stratum,
    Stratum_label,
    final_geo_group,
    final_group_type,
    Total_cluster_episodes,
    Treatment_variable,
    Episodes,
    Treatment_pct
  ) |>
  pivot_wider(
    names_from =
      Treatment_variable,
    
    values_from = c(
      Episodes,
      Treatment_pct
    ),
    
    names_glue =
      "{Treatment_variable}_{.value}",
    
    values_fill = list(
      Episodes = 0L,
      Treatment_pct = 0
    )
  ) |>
  arrange(
    factor(
      Stratum,
      levels = names(STRATUM_LABELS)
    ),
    final_geo_group
  ) |>
  rename(
    No_change_n =
      No_change_Episodes,
    
    No_change_pct =
      No_change_Treatment_pct,
    
    Reduction_lt50_n =
      Reduction_lt50_Episodes,
    
    Reduction_lt50_pct =
      Reduction_lt50_Treatment_pct,
    
    Reduction_ge50_n =
      Reduction_ge50_Episodes,
    
    Reduction_ge50_pct =
      Reduction_ge50_Treatment_pct,
    
    Dose_increase_n =
      Dose_increase_Episodes,
    
    Dose_increase_pct =
      Dose_increase_Treatment_pct
  )

print(
  final_cluster_treatment_distribution_wide,
  n = Inf,
  width = Inf
)

write_csv(
  final_cluster_treatment_distribution_wide,
  file.path(
    RES_DIR,
    "Tables",
    "final_cluster_treatment_distribution_wide.csv"
  )
)






# 43. DISTRIBUTION OF TREATMENT PERCENTAGES ACROSS CLUSTERS ############################

cluster_treatment_pct_summary <-
  final_cluster_treatment_distribution |>
  group_by(
    Stratum,
    Stratum_label,
    trt,
    Treatment_label
  ) |>
  summarise(
    Number_of_final_clusters =
      n(),
    
    Treatment_pct_mean =
      mean(Treatment_pct),
    
    Treatment_pct_SD =
      sd(Treatment_pct),
    
    Treatment_pct_min =
      min(Treatment_pct),
    
    Treatment_pct_Q1 =
      quantile(
        Treatment_pct,
        probs = 0.25,
        names = FALSE,
        type = 2
      ),
    
    Treatment_pct_median =
      median(Treatment_pct),
    
    Treatment_pct_Q3 =
      quantile(
        Treatment_pct,
        probs = 0.75,
        names = FALSE,
        type = 2
      ),
    
    Treatment_pct_max =
      max(Treatment_pct),
    
    Total_episodes =
      sum(Episodes),
    
    Episode_weighted_treatment_pct =
      100 *
      sum(Episodes) /
      sum(Total_cluster_episodes),
    
    .groups = "drop"
  ) |>
  arrange(
    factor(
      Stratum,
      levels = names(STRATUM_LABELS)
    ),
    trt
  )

print(
  cluster_treatment_pct_summary,
  n = Inf,
  width = Inf
)

write_csv(
  cluster_treatment_pct_summary,
  file.path(
    RES_DIR,
    "Tables",
    "distribution_treatment_percentages_across_final_clusters.csv"
  )
)
