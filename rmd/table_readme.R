
# library -----------------------------------------------------------------

pacman::p_load(tidyverse,
               sf)


# file overview -----------------------------------------------------------

fnm <- list.files("data_raw", pattern = "\\.rds")

df_fnm <- tibble(file = fnm) %>% 
  mutate(description = case_when(file == "data_env_flow.rds" ~ "Discharge data sourced from CaMa Flood.",
                                 file == "data_env_local.rds" ~ "Site-level environmental data.",
                                 file == "data_env_wsd.rds" ~ "Watershed-level environmental data.",
                                 file == "data_fcl_src.rds" ~ "Source data for food chain length. This data contains sites that are located in semi-lentic rivers (> 5k sq-km in watershed area). This source data was used to generate `wgs84_fcl_site.rds`.",
                                 file == "wgs84_fcl_site.rds" ~ "Point layer for site coordinates (with estimates of food chain length; class `sf`). Include temporal duplicates (either seasonal or annual) at some sites and some watersheds that are removed through a screening process",
                                 file == "wgs84_outlet.rds" ~ "Outlet coordinates of study watersheds (class `sf`). Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`).",
                                 file == "wgs84_region_lev01.rds" ~ "Polygon layer for HydroBASINS level-one (class `sf`).",
                                 file == "wgs84_str_sub.rds" ~ "Polyline layer for study watersheds. Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`).",
                                 file == "wgs84_subwatershed.rds" ~ "Polygon layer for upstream watershed at each sampling site (class `sf`).",
                                 file == "wgs84_wsd_sub.rds" ~ "Polygon layer for study watersheds (class `sf`). Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`).")) %>% 
  mutate(file = paste0("`", file, "`"))

# column information ------------------------------------------------------

## file names
fid <- !str_detect(fnm, "wgs84")
fid[str_detect(fnm, "fcl_site")] <- TRUE
fnmx <- fnm[fid]

## column names
cnm <- paste0("data_raw/", fnmx) %>% 
  lapply(readRDS) %>% 
  lapply(colnames)

## n columns for each data
lcnm <- sapply(cnm, length)

## file name duplication
x <- sapply(1:length(fnmx), function(i) rep(fnmx[i], each = lcnm[i]))

## data frame for column specifications
df_cnm <- tibble(file = unlist(x),
                 column = unlist(cnm)) %>% 
  mutate(description = 
           case_when(
             column == "area" ~ "Watershed area.",
             column == "base_species" ~ "Lowest taxonomic identity of the baseline species.",
             column == "base_species_taxon" ~ "Taxonomic group of the baseline species.",
             column == "coord_precision" ~ "Coordinate precision. See Section **Coordinate Precision**.",
             column == "csize" ~ "Community size. Identical to the number of links within a watershed.",
             column == "data_from" ~ "Data location in the published article.",
             column == "date" ~ "Date for discharge.",
             column == "discharge" ~ "Discharge data simulated by CaMaFlood.",
             column == "elev_mean" ~ "Mean elevation within the watershed.",
             column == "elev_sd" ~ "SD elevation within the watershed.",
             column == "fcl" ~ "Food chain length (re-defined). Corrected for trophic enrichment factor (tef = 3.4 for all data points).",
             column == "fcl_mean" ~ "Food chain length (uncorrected). Reported values in original articles.",
             column == "fcl_se" ~ "Standard error of food chain length estimate. Uncorrected for trophic enrichment factor.",
             column == "flag" ~ "Whether flagged or not. Data are flagged when they do not meet inclusion criteria (note: flagged data have been excluded).",
             column == "forest_b1km" ~ "Fraction of forest within 1km buffer of the sampling site.",
             column == "frac_agri" ~ "Fraction of agriculture within the upstream watershed area.",
             column == "frac_forest" ~ "Fraction of forest within the upstream watershed area.",
             column == "frac_urban" ~ "Fraction of urban within the upstream watershed area.",
             column == "habitat" ~ "Habitat type of the sampling site.",
             column == "hfp" ~ "Human footprint index.",
             column == "huid" ~ "Hydrological unit ID. Grouped by HydroBASINS level-4. Used only for computational purposes.",
             column == "key" ~ "GeoTIFF key used in watershed analysis. Ignore.",
             column == "lambda" ~ "Branching rate.",
             column == "last_verification" ~ "Last verified.",
             column == "lat" ~ "Latitude (Geodetic, WGS84).",
             column == "local_area" ~ "Upstream watershed area at the sampling site.",
             column == "local_elev" ~ "Elevation at the sampling site.",
             column == "local_prec" ~ "Average annual precipitation amount [kg/(year m^2)] within the upstream watershed area.",
             column == "local_temp" ~ "Average air temperature [degree Celcius] within the upstream watershed area.",
             column == "lon" ~ "Longitude (Geodetic, WGS84).",
             column == "mean_deltaC_base" ~ "Average value of delta 13 C for the baseline species.",
             column == "mean_deltaC_top" ~ "Average value of delta 13 C for the top consumer.",
             column == "mean_deltaN_base" ~ "Average value of delta 15 N for the baseline species.",
             column == "mean_deltaN_top" ~ "Average value of delta 15 N for the top consumer.",
             column == "month" ~ "Sampling month.",
             column == "n_link" ~ "Number of links within the watershed.",
             column == "n_sample_bottom_si" ~ "Number of stable isotope samples for the baseline species.",
             column == "n_sample_top_si" ~ "Number of stable isotope samples for the top consumer.",
             column == "note" ~ "Note",
             column == "oid" ~ "Outlet ID.",
             column == "org_site_id" ~ "Site name used in the original article.",
             column == "p_branch" ~ "Branching probability.",
             column == "prec" ~ "Average annual precipitation amount [kg/(year m^2)] within the watershed.",
             column == "predator_collection_scheme" ~ "Top predator collection scheme. See Section **Predator Collection**.",
             column == "r_length" ~ "Total river length within the watershed.",
             column == "sd_deltaC_base" ~ "Standard deviation of delta 13 C for the baseline species.",
             column == "sd_deltaC_top" ~ "Standard deviation of delta 13 C for the top consumer.",
             column == "sd_deltaN_base" ~ "Standard deviation of delta 15 N for the baseline species.",
             column == "sd_deltaN_top" ~ "Standard deviation of delta 15 N for the top consumer.",
             column == "sid" ~ "Site ID. Re-defined for our analysis.",
             column == "study_id" ~ "Study ID for the published article.",
             column == "tef" ~ "Trophic enrichment factor used in the original article.",
             column == "temp" ~ "Average air temperature [degree Celcius] within the watershed.",
             column == "tifid" ~ "GeoTIFF ID produced in the process of watershed analysis. Ignore.",
             column == "top_predator_collected" ~ "Whether top predator(s) is collected or not. See Section **Predator Collection**.",
             column == "top_predator_species" ~ "Lowest taxonomic identity of the top consumer.",
             column == "top_predator_taxon" ~ "Taxonomic group of the top consumer.",
             column == "tp_base" ~ "Trophic position of the baseline species.",
             column == "uid" ~ "Unique ID for the watershed.",
             column == "unit_area" ~ "Unit of watershed area.",
             column == "unit_lambda" ~ "Unit of branching rate.",
             column == "unit_r_length" ~ "Unit of river length.",
             column == "utm" ~ "UTM code.",
             column == "who_is_extracting" ~ "Person in charge of data extraction.",
             column == "wid" ~ "Watershed ID. Possible reduntant ID numbers between huid units. Use `uid` for grouping, which is a combination of `huid` and `wid`.)",
             column == "year_hfp" ~ "Year of human footprint index used.",
             column == "year_sampled" ~ "Sampling year of stable isotope data.")
  ) %>% 
  filter(column != "geometry") %>% 
  mutate(file = paste0("`", file, "`"),
         file = ifelse(duplicated(file), NA, file))


