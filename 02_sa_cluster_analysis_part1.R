# 1. PACKAGES ############################

required_packages <- c("dplyr", "readr", "purrr", "sf", "bcmaps", "stringr", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) install.packages(missing_packages)

library(dplyr)
library(readr)
library(purrr)
library(sf)
library(bcmaps)
library(stringr)
library(tidyr)

DATA_DIR <- "E:/SAS Data/Missed dose protocol_v3/"
OUT_DIR <- "H:/GitHub/Missed-Dose-Protocols-PRI/Results/MET/Tables"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Number of nearby candidates saved for each geographic unit.
# The grouping algorithm can try rank 1, then rank 2, and so on.
N_LHA_NEIGHBOURS  <- 10L
N_HSDA_NEIGHBOURS <- 4L

# 2. DOWNLOAD OFFICIAL BC HEALTH BOUNDARIES ############################

lha_sf  <- bcmaps::health_lha(ask = FALSE, force = FALSE) |> st_make_valid() |> st_transform(3005)
hsda_sf <- bcmaps::health_hsda(ask = FALSE, force = FALSE) |> st_make_valid() |> st_transform(3005)


# 3. STANDARDIZE BOUNDARY VARIABLES ############################

lha_code_col   <- "LOCAL_HLTH_AREA_CODE"
lha_name_col   <- "LOCAL_HLTH_AREA_NAME"
lha_hsda_col   <- "HLTH_SERVICE_DLVR_AREA_CODE"
lha_ha_col     <- "HLTH_AUTHORITY_CODE"

hsda_code_col  <- "HLTH_SERVICE_DLVR_AREA_CODE"
hsda_name_col  <- "HLTH_SERVICE_DLVR_AREA_NAME"
hsda_ha_col    <- "HLTH_AUTHORITY_CODE"

required_lha_cols <- c(
  lha_code_col, lha_name_col, lha_hsda_col, lha_ha_col, "geometry"
)

required_hsda_cols <- c(
  hsda_code_col, hsda_name_col, hsda_ha_col, "geometry"
)

if (!all(required_lha_cols %in% names(lha_sf))) {
  stop(
    "Missing required LHA columns: ",
    paste(setdiff(required_lha_cols, names(lha_sf)), collapse = ", ")
  )
}

if (!all(required_hsda_cols %in% names(hsda_sf))) {
  stop(
    "Missing required HSDA columns: ",
    paste(setdiff(required_hsda_cols, names(hsda_sf)), collapse = ", ")
  )
}


lha_std <- lha_sf |>
  transmute(
    lha_id = str_pad(
      str_remove_all(as.character(.data[[lha_code_col]]), "[^0-9]"),
      width = 3,
      pad = "0"
    ),
    
    lha_name = as.character(.data[[lha_name_col]]),
    
    hsda_id = str_pad(
      str_remove_all(as.character(.data[[lha_hsda_col]]), "[^0-9]"),
      width = 2,
      pad = "0"
    ),
    
    ha_id = str_pad(
      str_remove_all(as.character(.data[[lha_ha_col]]), "[^0-9]"),
      width = 1,
      pad = "0"
    )
  ) |>
  filter(
    str_detect(lha_id, "^[1-5][0-9]{2}$"),
    str_detect(hsda_id, "^[1-5][0-9]$"),
    ha_id %in% as.character(1:5)
  ) |>
  group_by(lha_id, lha_name, hsda_id, ha_id) |>
  summarise(do_union = TRUE, .groups = "drop") |>
  arrange(lha_id)



hsda_std <- hsda_sf |>
  transmute(
    hsda_id = str_pad(
      str_remove_all(as.character(.data[[hsda_code_col]]), "[^0-9]"),
      width = 2,
      pad = "0"
    ),
    
    hsda_name = as.character(.data[[hsda_name_col]]),
    
    ha_id = str_pad(
      str_remove_all(as.character(.data[[hsda_ha_col]]), "[^0-9]"),
      width = 1,
      pad = "0"
    )
  ) |>
  filter(
    str_detect(hsda_id, "^[1-5][0-9]$"),
    ha_id %in% as.character(1:5)
  ) |>
  group_by(hsda_id, hsda_name, ha_id) |>
  summarise(do_union = TRUE, .groups = "drop") |>
  arrange(hsda_id)


# 4. VALIDATE STANDARDIZED GEOGRAPHY ############################

VALID_HSDA_CODES <- c(
  "11", "12", "13", "14",
  "21", "22", "23",
  "31", "32", "33",
  "41", "42", "43",
  "51", "52", "53"
)

lha_validation <- lha_std |>
  st_drop_geometry() |>
  summarise(
    Number_of_LHAs = n(),
    Unique_LHA_codes = n_distinct(lha_id),
    Missing_LHA_codes = sum(is.na(lha_id)),
    Missing_LHA_names = sum(is.na(lha_name)),
    Missing_HSDA_codes = sum(is.na(hsda_id)),
    Missing_HA_codes = sum(is.na(ha_id)),
    Invalid_HSDA_codes = sum(!hsda_id %in% VALID_HSDA_CODES),
    Invalid_HA_codes = sum(!ha_id %in% as.character(1:5))
  )

hsda_validation <- hsda_std |>
  st_drop_geometry() |>
  summarise(
    Number_of_HSDAs = n(),
    Unique_HSDA_codes = n_distinct(hsda_id),
    Missing_HSDA_codes = sum(is.na(hsda_id)),
    Missing_HSDA_names = sum(is.na(hsda_name)),
    Missing_HA_codes = sum(is.na(ha_id)),
    Invalid_HSDA_codes = sum(!hsda_id %in% VALID_HSDA_CODES),
    Invalid_HA_codes = sum(!ha_id %in% as.character(1:5)),
    Number_of_HAs = n_distinct(ha_id)
  )

print(lha_validation, width = Inf)
print(hsda_validation, width = Inf)

stopifnot(
  lha_validation$Number_of_LHAs == lha_validation$Unique_LHA_codes,
  lha_validation$Missing_LHA_codes == 0L,
  lha_validation$Missing_HSDA_codes == 0L,
  lha_validation$Missing_HA_codes == 0L,
  lha_validation$Invalid_HSDA_codes == 0L,
  lha_validation$Invalid_HA_codes == 0L,
  
  hsda_validation$Number_of_HSDAs == 16L,
  hsda_validation$Unique_HSDA_codes == 16L,
  hsda_validation$Missing_HSDA_codes == 0L,
  hsda_validation$Missing_HA_codes == 0L,
  hsda_validation$Invalid_HSDA_codes == 0L,
  hsda_validation$Invalid_HA_codes == 0L,
  hsda_validation$Number_of_HAs == 5L
)

hierarchy_disagreement <- lha_std |>
  st_drop_geometry() |>
  distinct(hsda_id, ha_id) |>
  left_join(
    hsda_std |>
      st_drop_geometry() |>
      distinct(hsda_id, hsda_ha_id = ha_id),
    by = "hsda_id"
  ) |>
  filter(is.na(hsda_ha_id) | ha_id != hsda_ha_id)

if (nrow(hierarchy_disagreement) > 0L) {
  print(hierarchy_disagreement, n = Inf)
  stop("The LHA and HSDA layers disagree about the HSDA-to-HA hierarchy.")
}


# 5. CREATE NEARBY-AREA TABLES ############################

if (!requireNamespace("units", quietly = TRUE)) install.packages("units")
make_neighbour_file <- function(x, id, name, parent, max_neighbours) {
  
  id_vec <- x[[id]]
  name_vec <- x[[name]]
  parent_vec <- x[[parent]]
  
  distance_km <- units::drop_units(sf::st_distance(x)) / 1000
  touching <- sf::st_touches(x)
  
  map_dfr(seq_len(nrow(x)), function(i) {
    
    candidates <- which(
      parent_vec == parent_vec[i] &
        seq_len(nrow(x)) != i
    )
    
    if (length(candidates) == 0L) return(tibble())
    
    tibble(
      area_id = id_vec[i],
      area_name = name_vec[i],
      parent_id = parent_vec[i],
      nearby_id = id_vec[candidates],
      nearby_name = name_vec[candidates],
      shares_boundary = candidates %in% touching[[i]],
      boundary_distance_km = as.numeric(distance_km[i, candidates])
    ) |>
      arrange(
        desc(shares_boundary),
        boundary_distance_km,
        nearby_id
      ) |>
      mutate(
        nearby_rank = row_number(),
        .before = nearby_id
      ) |>
      slice_head(n = max_neighbours)
  })
}




# 6. GENERATE ADJACENCY FILES ############################

N_LHA_NEIGHBOURS  <- 10L
N_HSDA_NEIGHBOURS <- 4L

lha_neighbours <- make_neighbour_file(
  x = lha_std,
  id = "lha_id",
  name = "lha_name",
  parent = "hsda_id",
  max_neighbours = N_LHA_NEIGHBOURS
) |>
  rename(
    lha_id = area_id,
    lha_name = area_name,
    hsda_id = parent_id,
    neighbor_lha_id = nearby_id,
    neighbor_lha_name = nearby_name
  ) |>
  select(
    lha_id, lha_name, hsda_id, nearby_rank,
    neighbor_lha_id, neighbor_lha_name,
    shares_boundary, boundary_distance_km
  )

hsda_neighbours <- make_neighbour_file(
  x = hsda_std,
  id = "hsda_id",
  name = "hsda_name",
  parent = "ha_id",
  max_neighbours = N_HSDA_NEIGHBOURS
) |>
  rename(
    hsda_id = area_id,
    hsda_name = area_name,
    ha_id = parent_id,
    neighbor_hsda_id = nearby_id,
    neighbor_hsda_name = nearby_name
  ) |>
  select(
    hsda_id, hsda_name, ha_id, nearby_rank,
    neighbor_hsda_id, neighbor_hsda_name,
    shares_boundary, boundary_distance_km
  )



# 7. EXPORT ADJACENCY FILES ############################

write_csv(
  lha_neighbours,
  file.path(DATA_DIR, "BC_LHA_nearby_within_HSDA.csv")
)

write_csv(
  hsda_neighbours,
  file.path(DATA_DIR, "BC_HSDA_nearby_within_HA.csv")
)

write_csv(
  st_drop_geometry(lha_std),
  file.path(DATA_DIR, "BC_LHA_HSDA_HA_crosswalk.csv")
)

write_csv(
  st_drop_geometry(hsda_std),
  file.path(DATA_DIR, "BC_HSDA_HA_crosswalk.csv")
)

cat(
  "\nCreated geographic files successfully:\n",
  file.path(DATA_DIR, "BC_LHA_nearby_within_HSDA.csv"), "\n",
  file.path(DATA_DIR, "BC_HSDA_nearby_within_HA.csv"), "\n",
  file.path(DATA_DIR, "BC_LHA_HSDA_HA_crosswalk.csv"), "\n",
  file.path(DATA_DIR, "BC_HSDA_HA_crosswalk.csv"), "\n",
  sep = ""
)





# 8. SAFE MINIMUM FUNCTION ############################

safe_min <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) NA_real_ else min(x)
}

lha_crosswalk <- st_drop_geometry(lha_std)
hsda_crosswalk <- st_drop_geometry(hsda_std)

lha_validation <- lha_neighbours |>
  left_join(
    lha_crosswalk |>
      select(lha_id, source_hsda_id = hsda_id),
    by = "lha_id"
  ) |>
  left_join(
    lha_crosswalk |>
      select(
        neighbor_lha_id = lha_id,
        neighbor_hsda_id = hsda_id
      ),
    by = "neighbor_lha_id"
  ) |>
  mutate(
    same_HSDA =
      source_hsda_id ==
      neighbor_hsda_id
  )

hsda_validation <- hsda_neighbours |>
  left_join(
    hsda_crosswalk |>
      select(hsda_id, source_ha_id = ha_id),
    by = "hsda_id"
  ) |>
  left_join(
    hsda_crosswalk |>
      select(
        neighbor_hsda_id = hsda_id,
        neighbor_ha_id = ha_id
      ),
    by = "neighbor_hsda_id"
  ) |>
  mutate(
    same_HA =
      source_ha_id ==
      neighbor_ha_id
  )

stopifnot(
  all(lha_validation$same_HSDA),
  all(hsda_validation$same_HA),
  all(lha_neighbours$lha_id != lha_neighbours$neighbor_lha_id),
  all(hsda_neighbours$hsda_id != hsda_neighbours$neighbor_hsda_id),
  !anyDuplicated(
    lha_neighbours[
      c("lha_id", "neighbor_lha_id")
    ]
  ),
  !anyDuplicated(
    hsda_neighbours[
      c("hsda_id", "neighbor_hsda_id")
    ]
  )
)

cat(
  "\nGeographic validation passed:",
  "\n- Every LHA candidate is in the same HSDA.",
  "\n- Every HSDA candidate is in the same HA.",
  "\n- No area is listed as its own neighbour.",
  "\n- No duplicate directed neighbour pairs were found.\n"
)




# 9. SUMMARIZE NEIGHBOUR AVAILABILITY ############################

lha_neighbour_summary <- lha_crosswalk |>
  select(lha_id, lha_name, hsda_id) |>
  left_join(lha_neighbours, by = c("lha_id", "lha_name", "hsda_id")) |>
  group_by(lha_id, lha_name, hsda_id) |>
  summarise(
    Candidate_LHAs = sum(!is.na(neighbor_lha_id)),
    Directly_adjacent_LHAs = sum(shares_boundary %in% TRUE),
    Nearest_nonadjacent_distance_km =
      safe_min(boundary_distance_km[shares_boundary %in% FALSE]),
    .groups = "drop"
  )

hsda_neighbour_summary <- hsda_crosswalk |>
  select(hsda_id, hsda_name, ha_id) |>
  left_join(hsda_neighbours, by = c("hsda_id", "hsda_name", "ha_id")) |>
  group_by(hsda_id, hsda_name, ha_id) |>
  summarise(
    Candidate_HSDAs = sum(!is.na(neighbor_hsda_id)),
    Directly_adjacent_HSDAs = sum(shares_boundary %in% TRUE),
    Nearest_nonadjacent_distance_km =
      safe_min(boundary_distance_km[shares_boundary %in% FALSE]),
    .groups = "drop"
  )



# 10. FINAL GEOGRAPHIC FILE SUMMARY ############################

geographic_file_summary <- tibble(
  Geographic_level = c("LHA", "HSDA"),
  
  Number_of_areas = c(
    n_distinct(lha_crosswalk$lha_id),
    n_distinct(hsda_crosswalk$hsda_id)
  ),
  
  Areas_with_candidates = c(
    n_distinct(lha_neighbours$lha_id),
    n_distinct(hsda_neighbours$hsda_id)
  ),
  
  Areas_without_candidates = c(
    n_distinct(lha_crosswalk$lha_id) -
      n_distinct(lha_neighbours$lha_id),
    
    n_distinct(hsda_crosswalk$hsda_id) -
      n_distinct(hsda_neighbours$hsda_id)
  ),
  
  Directed_candidate_pairs = c(
    nrow(lha_neighbours),
    nrow(hsda_neighbours)
  ),
  
  Direct_boundary_pairs = c(
    sum(lha_neighbours$shares_boundary),
    sum(hsda_neighbours$shares_boundary)
  )
)

print(
  geographic_file_summary,
  n = Inf,
  width = Inf
)

write_csv(
  geographic_file_summary,
  file.path(DATA_DIR, "BC_health_geographic_candidate_summary.csv")
)
