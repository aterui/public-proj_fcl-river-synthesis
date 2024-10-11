README
================

# Article Information

**Tentative title**: Geometric ecosystem complexity regulates food
chains

**Authors**: Terui A, Shibasaki S, Pomeranz JPF, XXX, Yamazaki D, Finlay
JC

# Analysis Flow

## Theoretical Analysis

Script `run_theory.R` runs all analysis, including analytical and
numerical explorations.

For source codes, see:

- `analysis_theory_m_heatmap.R` (analytical)

- `analysis_theory_m_line.R` (analytical, for visualization)

- `analysis_theory_si.R` (numerical)

For visualization, see:

- `figure_theo_main.R` (main figures)

- `figure_theo_si.R` (supporting figures)

## Empirical Analysis

**Bayesian data analysis**: Script `run_emp_model_hier_nimble.R` runs
the Bayesian hierarchical models. Model source codes are
`model_nimble_hi.R` (random intercept model) and `model_nimble_hsl.R`
(random intercept and slope model).

Data are formatted in and sourced from:

- `format_emp_data4nimble.R` (generate `data_fmt/data_fcl_reg.rds`,
  which was used for statistical analysis)

- `gis_emp_wsd_weight.R` (distance ratio; `data_fmt/data_weight.rds`)

For visualization, see:

- `figure_emp_main.R`

**Precipitation-flow analysis**: Script `format_emp_flow.R` runs the GAM
analysis. Script `analysis_emp_cor.R` analyze the relationship between
precipitation and discharge anomalies.

**Predictor correlation**: Script `analysis_emp_cor.R` analyzes
correlations between predictors.

# Data Sources

## File Description

| file                     | description                                                                                                                                                                                                                |
|:-------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `data_env_flow.rds`      | Discharge data sourced from CaMa Flood.                                                                                                                                                                                    |
| `data_env_local.rds`     | Site-level environmental data. Averaged for upstream watershed area.                                                                                                                                                       |
| `data_env_wsd.rds`       | Watershed-level environmental data.                                                                                                                                                                                        |
| `data_fcl_src.rds`       | Source data for food chain length. This data contains sites that are located in semi-lentic rivers (\> 5k sq-km in watershed area). This source data was used to generate `wgs84_fcl_site.rds`.                            |
| `wgs84_fcl_site.rds`     | Point layer for site coordinates (with estimates of food chain length; class `sf`). Include temporal duplicates (either seasonal or annual) at some sites and some watersheds that are removed through a screening process |
| `wgs84_outlet.rds`       | Outlet coordinates of study watersheds (class `sf`). Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`).                                                                |
| `wgs84_region_lev01.rds` | Polygon layer for HydroBASINS level-one (class `sf`).                                                                                                                                                                      |
| `wgs84_str_sub.rds`      | Polyline layer for study watersheds. Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`).                                                                                |
| `wgs84_subwatershed.rds` | Polygon layer for upstream watershed at each sampling site (class `sf`).                                                                                                                                                   |
| `wgs84_wsd_sub.rds`      | Polygon layer for study watersheds (class `sf`). Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`).                                                                    |

## Column Description

| file                 | column                     | description                                                                                                                          |
|:---------------------|:---------------------------|:-------------------------------------------------------------------------------------------------------------------------------------|
| `data_env_flow.rds`  | sid                        | Site ID. Re-defined for our analysis.                                                                                                |
|                      | org_site_id                | Site name used in the original article.                                                                                              |
|                      | date                       | Date for discharge.                                                                                                                  |
|                      | discharge                  | Discharge data simulated by CaMaFlood.                                                                                               |
| `data_env_local.rds` | sid                        | Site ID. Re-defined for our analysis.                                                                                                |
|                      | local_area                 | Upstream watershed area at the sampling site.                                                                                        |
|                      | unit_area                  | Unit of watershed area.                                                                                                              |
|                      | year_sampled               | Sampling year of stable isotope data.                                                                                                |
|                      | year_hfp                   | Year of human footprint index used.                                                                                                  |
|                      | hfp                        | Human footprint index.                                                                                                               |
|                      | local_temp                 | Average air temperature \[degree Celcius\] within the upstream watershed area.                                                       |
|                      | local_prec                 | Average annual precipitation amount \[kg/(year m^2)\] within the upstream watershed area.                                            |
|                      | frac_agri                  | Fraction of agriculture within the upstream watershed area.                                                                          |
|                      | frac_urban                 | Fraction of urban within the upstream watershed area.                                                                                |
|                      | frac_forest                | Fraction of forest within the upstream watershed area.                                                                               |
|                      | forest_b1km                | Fraction of forest within 1km buffer of the sampling site.                                                                           |
|                      | utm                        | UTM code.                                                                                                                            |
|                      | local_elev                 | Elevation at the sampling site.                                                                                                      |
| `data_env_wsd.rds`   | uid                        | Unique ID for the watershed.                                                                                                         |
|                      | wid                        | Watershed ID. Possible reduntant ID numbers between huid units. Use `uid` for grouping, which is a combination of `huid` and `wid`.) |
|                      | tifid                      | GeoTIFF id produced in the process of watershed analysis. Ignore.                                                                    |
|                      | oid                        | Outlet ID.                                                                                                                           |
|                      | huid                       | Hydrological unit ID. Grouped by HydroBASINS level-4. Used only for computational purposes.                                          |
|                      | area                       | Watershed area.                                                                                                                      |
|                      | unit_area                  | Unit of watershed area.                                                                                                              |
|                      | r_length                   | Total river length within the watershed.                                                                                             |
|                      | lambda                     | Branching rate.                                                                                                                      |
|                      | csize                      | Community size. Identical to the number of links within a watershed.                                                                 |
|                      | n_link                     | Number of links within the watershed.                                                                                                |
|                      | unit_lambda                | Unit of branching rate.                                                                                                              |
|                      | unit_r_length              | Unit of river length.                                                                                                                |
|                      | p_branch                   | Branching probability.                                                                                                               |
|                      | year_sampled               | Sampling year of stable isotope data.                                                                                                |
|                      | year_hfp                   | Year of human footprint index used.                                                                                                  |
|                      | hfp                        | Human footprint index.                                                                                                               |
|                      | temp                       | Average air temperature \[degree Celcius\] within the watershed.                                                                     |
|                      | prec                       | Average annual precipitation amount \[kg/(year m^2)\] within the watershed.                                                          |
|                      | key                        | GeoTIFF key used in watershed analysis. Ignore.                                                                                      |
|                      | elev_mean                  | Mean elevation within the watershed.                                                                                                 |
|                      | elev_sd                    | SD elevation within the watershed.                                                                                                   |
| `data_fcl_src.rds`   | lat                        | Latitude (Geodetic, WGS84).                                                                                                          |
|                      | lon                        | Longitude (Geodetic, WGS84).                                                                                                         |
|                      | study_id                   | Study identification code for published articles.                                                                                    |
|                      | sid                        | Site ID. Re-defined for our analysis.                                                                                                |
|                      | who_is_extracting          | Person in charge of data extraction.                                                                                                 |
|                      | org_site_id                | Site name used in the original article.                                                                                              |
|                      | month                      | Sampling month.                                                                                                                      |
|                      | habitat                    | Habitat type of the sampling site.                                                                                                   |
|                      | year_sampled               | Sampling year of stable isotope data.                                                                                                |
|                      | top_predator_collected     | Whether top predator(s) is collected or not. See Section **Predator Collection**                                                     |
|                      | predator_collection_scheme | Top predator collection scheme. See Section **Predator Collection**                                                                  |
|                      | fcl_mean                   | Food chain length (uncorrected). Reported values in original articles.                                                               |
|                      | fcl_se                     | Standard error of food chain length estimate. Uncorrected for trophic enrichment factor.                                             |
|                      | mean_deltaN_top            | Average value of delta 15 N for the top consumer.                                                                                    |
|                      | sd_deltaN_top              | Standard deviation of delta 15 N for the top consumer.                                                                               |
|                      | mean_deltaC_top            | Average value of delta 13 C for the top consumer.                                                                                    |
|                      | sd_deltaC_top              | Standard deviation of delta 13 C for the top consumer.                                                                               |
|                      | mean_deltaN_base           | Average value of delta 15 N for the baseline species.                                                                                |
|                      | sd_deltaN_base             | Standard deviation of delta 15 N for the baseline species.                                                                           |
|                      | mean_deltaC_base           | Average value of delta 13 C for the baseline species.                                                                                |
|                      | sd_deltaC_base             | Standard deviation of delta 13 C for the baseline species.                                                                           |
|                      | n_sample_top_si            | Number of stable isotope samples for the top consumer.                                                                               |
|                      | n_sample_bottom_si         | Number of stable isotope samples for the baseline species.                                                                           |
|                      | tef                        | Trophic enrichment factor used in the original article.                                                                              |
|                      | top_predator_taxon         | Taxonomic group of the top consumer.                                                                                                 |
|                      | top_predator_species       | Lowest taxonomic identity of the top consumer.                                                                                       |
|                      | base_species_taxon         | Taxonomic group of the baseline species.                                                                                             |
|                      | base_species               | Lowest taxonomic identity of the baseline species.                                                                                   |
|                      | tp_base                    | Trophic position of the baseline species.                                                                                            |
|                      | coord_precision            | Coordinate precision. See Section **Coordinate Precision**                                                                           |
|                      | note                       | Note                                                                                                                                 |
|                      | data_from                  | Data location in the published article.                                                                                              |
|                      | last_verification          | Last verified.                                                                                                                       |
|                      | flag                       | Whether flagged or not. Data are flagged when they do not meet inclusion criteria (note: flagged data have been excluded).           |
|                      | fcl                        | Food chain length (re-defined). Corrected for trophic enrichment factor (tef = 3.4 for all data points).                             |
| `wgs84_fcl_site.rds` | uid                        | Unique ID for the watershed.                                                                                                         |
|                      | huid                       | Hydrological unit ID. Grouped by HydroBASINS level-4. Used only for computational purposes.                                          |
|                      | wid                        | Watershed ID. Possible reduntant ID numbers between huid units. Use `uid` for grouping, which is a combination of `huid` and `wid`.) |
|                      | sid                        | Site ID. Re-defined for our analysis.                                                                                                |
|                      | study_id                   | Study identification code for published articles.                                                                                    |
|                      | org_site_id                | Site name used in the original article.                                                                                              |
|                      | year_sampled               | Sampling year of stable isotope data.                                                                                                |
|                      | month                      | Sampling month.                                                                                                                      |
|                      | fcl                        | Food chain length (re-defined). Corrected for trophic enrichment factor (tef = 3.4 for all data points).                             |
|                      | top_predator_collected     | Whether top predator(s) is collected or not. See Section **Predator Collection**                                                     |
|                      | top_predator_taxon         | Taxonomic group of the top consumer.                                                                                                 |
|                      | top_predator_species       | Lowest taxonomic identity of the top consumer.                                                                                       |
|                      | base_species_taxon         | Taxonomic group of the baseline species.                                                                                             |
|                      | base_species               | Lowest taxonomic identity of the baseline species.                                                                                   |
