README
================

# Article Information

**Title**: Geometric ecosystem complexity regulates riverine food chains

**Authors**: Terui A, Shibasaki S, Pomeranz JPF, Ibrahim M, Isaac E,
LaRoque A, Yamazaki D, Finlay JC

# Analysis Flow

## Theoretical Analysis

Script `run_theory.R` runs all analysis, including analytical and
numerical explorations.

For source codes, see:

- `analysis_theo_m_2sp.R` (analytical, predator-prey, maintext)

- `analysis_theo_m_fw.R` (analytical, food web, maintext)

- `analysis_theo_si_2sp_sbn.R` (numerical, predator-prey: dispersal
  symmetry, SI Appendix)

- `analysis_theo_si_fw.R` (numerical, food web: no disturbance scenario,
  SI Appendix)

For visualization, see:

- `figure_theo_main.R` (main figures)

- `figure_theo_si.R` (SI Appendix figures)

This theoretical analysis uses the user-defined R Packages
[`rpom`](https://github.com/aterui/rpom) and
[`ecotools`](https://github.com/aterui/ecotools).

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

| file | description |
|:---|:---|
| `data_env_flow.rds` | Discharge data sourced from CaMa Flood. |
| `data_env_local.rds` | Site-level environmental data. |
| `data_env_wsd.rds` | Watershed-level environmental data. |
| `data_fcl_src.rds` | Source data for food chain length. This data contains sites that are located in semi-lentic rivers (\> 5k sq-km in watershed area). This source data was used to generate `sf_fcl_site.rds`. |
| `sf_channel.rds` | Polyline layer for study watersheds. Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`). |
| `sf_fcl_site.rds` | Point layer for site coordinates (with estimates of food chain length; class `sf`). Include temporal duplicates (either seasonal or annual) at some sites. Some watersheds are removed before analysis through a screening process (see `format_emp_data4nimble.R`) |
| `sf_hybas_lev01.rds` | Polygon layer for HydroBASINS level-one (class `sf`). |
| `sf_hybas_lev02.rds` | Polygon layer for HydroBASINS level-two (class `sf`). |
| `sf_outlet.rds` | Outlet coordinates of study watersheds (class `sf`). Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`). |
| `sf_subwsd.rds` | Polygon layer for upstream watershed at each sampling site (class `sf`). |
| `sf_wsd.rds` | Polygon layer for study watersheds (class `sf`). Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`). |

## Column Description

| file | column | description |
|:---|:---|:---|
| `data_env_flow.rds` | sid | Site ID. Re-defined for our analysis. |
|  | org_site_id | Site name used in the original article. |
|  | date | Date for discharge. |
|  | discharge | Discharge data simulated by CaMaFlood. |
| `data_env_local.rds` | huid | Hydrological unit ID. Grouped by HydroBASINS level-4. Used only for computational purposes. |
|  | hybas_id | HYBAS ID from HydroBASINS. |
|  | oid | Watershed/Outlet ID. The unit of watershed-level analysis. |
|  | sid | Site ID. Re-defined for our analysis. |
|  | local_area | Upstream watershed area at the sampling site. |
|  | unit_area | Unit of watershed area. |
|  | year_sampled | Sampling year of stable isotope data. |
|  | year_hfp | Year of human footprint index used. |
|  | local_hfp | Human footprint index within the upstream watershed area of the sampling site. |
|  | local_temp | Average air temperature \[degree Celcius\] within the upstream watershed area of the sampling site. |
|  | local_prec | Average annual precipitation amount \[kg/(year m^2)\] within the upstream watershed area of the sampling site. |
|  | frac_agri | Fraction of agriculture within the upstream watershed area. |
|  | frac_urban | Fraction of urban within the upstream watershed area. |
|  | frac_forest | Fraction of forest within the upstream watershed area. |
|  | forest_b1km | Fraction of forest within 1km buffer of the sampling site. |
|  | local_elev | Elevation at the sampling site. |
| `data_env_wsd.rds` | huid | Hydrological unit ID. Grouped by HydroBASINS level-4. Used only for computational purposes. |
|  | hybas_id | HYBAS ID from HydroBASINS. |
|  | oid | Watershed/Outlet ID. The unit of watershed-level analysis. |
|  | area | Watershed area. |
|  | unit_area | Unit of watershed area. |
|  | unit_elev | Unit of elevation. |
|  | r_length | Total river length within the watershed. |
|  | lambda | Branching rate. |
|  | csize | Community size. Identical to the number of links within a watershed. |
|  | n_link | Number of links within the watershed. |
|  | unit_lambda | Unit of branching rate. |
|  | unit_r_length | Unit of river length. |
|  | year_hfp | Year of human footprint index used. |
|  | hfp | Human footprint index at the watershed scale. |
|  | temp | Average air temperature \[degree Celcius\] within the watershed. |
|  | prec | Average annual precipitation amount \[kg/(year m^2)\] within the watershed. |
|  | elev_mean | Mean elevation within the watershed. |
|  | elev_sd | SD elevation within the watershed. |
| `data_fcl_src.rds` | lat | Latitude (Geodetic, WGS84). |
|  | lon | Longitude (Geodetic, WGS84). |
|  | study_id | Study ID for the published article. |
|  | sid | Site ID. Re-defined for our analysis. |
|  | who_is_extracting | Person in charge of data extraction. |
|  | org_site_id | Site name used in the original article. |
|  | month | Sampling month. |
|  | habitat | Habitat type of the sampling site. |
|  | year_sampled | Sampling year of stable isotope data. |
|  | top_predator_collected | Whether top predator(s) is collected or not. See Section **Predator Collection**. |
|  | predator_collection_scheme | Top predator collection scheme. See Section **Predator Collection**. |
|  | fcl_mean | Food chain length (uncorrected). Reported values in original articles. |
|  | fcl_se | Standard error of food chain length estimate. Uncorrected for trophic enrichment factor. |
|  | mean_deltaN_top | Average value of delta 15 N for the top consumer. |
|  | sd_deltaN_top | Standard deviation of delta 15 N for the top consumer. |
|  | mean_deltaC_top | Average value of delta 13 C for the top consumer. |
|  | sd_deltaC_top | Standard deviation of delta 13 C for the top consumer. |
|  | mean_deltaN_base | Average value of delta 15 N for the baseline species. |
|  | sd_deltaN_base | Standard deviation of delta 15 N for the baseline species. |
|  | mean_deltaC_base | Average value of delta 13 C for the baseline species. |
|  | sd_deltaC_base | Standard deviation of delta 13 C for the baseline species. |
|  | n_sample_top_si | Number of stable isotope samples for the top consumer. |
|  | n_sample_bottom_si | Number of stable isotope samples for the baseline species. |
|  | tef | Trophic enrichment factor used in the original article. |
|  | top_predator_taxon | Taxonomic group of the top consumer. |
|  | top_predator_species | Lowest taxonomic identity of the top consumer. |
|  | base_species_taxon | Taxonomic group of the baseline species. |
|  | base_species | Lowest taxonomic identity of the baseline species. |
|  | tp_base | Trophic position of the baseline species. |
|  | coord_precision | Coordinate precision. See Section **Coordinate Precision**. |
|  | note | Note |
|  | data_from | Data location in the published article. |
|  | last_verification | Last verified. |
|  | flag | Whether flagged or not. Data are flagged when they do not meet inclusion criteria (note: flagged data have been excluded). |
|  | fcl | Food chain length (re-defined). Corrected for trophic enrichment factor (tef = 3.4 for all data points). |
| `sf_fcl_site.rds` | huid | Hydrological unit ID. Grouped by HydroBASINS level-4. Used only for computational purposes. |
|  | hybas_id | HYBAS ID from HydroBASINS. |
|  | oid | Watershed/Outlet ID. The unit of watershed-level analysis. |
|  | sid | Site ID. Re-defined for our analysis. |
|  | study_id | Study ID for the published article. |
|  | org_site_id | Site name used in the original article. |
|  | year_sampled | Sampling year of stable isotope data. |
|  | month | Sampling month. |
|  | fcl | Food chain length (re-defined). Corrected for trophic enrichment factor (tef = 3.4 for all data points). |
|  | top_predator_collected | Whether top predator(s) is collected or not. See Section **Predator Collection**. |
|  | top_predator_taxon | Taxonomic group of the top consumer. |
|  | top_predator_species | Lowest taxonomic identity of the top consumer. |
|  | base_species_taxon | Taxonomic group of the baseline species. |
|  | base_species | Lowest taxonomic identity of the baseline species. |

## Predator Collection

Column `top_predator_collected` indicates whether suitable top
predator(s) are collected (`Y`) or not (`N`). We consider top
predator(s) to be suitably collected (`Y`) if they meet either of the
following criteria: (**A**) the article explicitly states that the
identified top predators were collected for stable isotope analysis,
and/or (**B**) the study targets the entire aquatic or consumer
community within the system to ensure that the stable isotope signatures
of the top predators are included. If (**C**) the study targets a subset
of community members excluding potential top predators, we enter `N`.

## Coordinate Precision

Coordinate precision was categorized as follows (**4** and **5** have
been excluded in `sf_fcl_site.rds`):

**1**: Exact coordinates available (either through article or author
correspondence); **2**: Nearly exact coordinates recovered from text;
**3**: Coordinates visually recovered from a map with reference to
landmarks/river morphology; **4**: Coordinates approximated by river
name or multiple sites aggregated; **5**: Not available.

## Geospatial Data

Geospatial data (RDS files with prefix `sf_`) were generated in a
[separate
repository](https://github.com/aterui/priv-proj_gis-global-fcl). The
data package used in this repository is **GFCL version 0.3.0**.

# Session Information

    ## R version 4.5.1 (2025-06-13 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26200)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: Asia/Tokyo
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] terra_1.9-1      mapview_2.11.4   tictoc_1.1       progressr_0.11.0
    ##  [5] doFuture_1.2.1   future_1.69.0    rpom_0.2.0       mcbrnet_1.4.2   
    ##  [9] foreach_1.5.2    ape_5.8-1        ecotools_0.1.1   xtable_1.8-8    
    ## [13] nimble_1.4.2     cowplot_1.2.0    magick_2.7.3     igraph_2.1.1    
    ## [17] patchwork_1.3.2  sf_1.1-1         lubridate_1.9.5  forcats_1.0.0   
    ## [21] stringr_1.6.0    dplyr_1.2.1      purrr_1.2.2      readr_2.2.0     
    ## [25] tidyr_1.3.2      tibble_3.3.1     ggplot2_4.0.3    tidyverse_2.0.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1    farver_2.1.2        S7_0.2.2           
    ##  [4] fastmap_1.2.0       leaflet_2.2.3       pracma_2.4.6       
    ##  [7] pacman_0.5.1        digest_0.6.33       timechange_0.4.0   
    ## [10] lifecycle_1.0.5     magrittr_2.0.5      compiler_4.5.1     
    ## [13] rlang_1.2.0         tools_4.5.1         yaml_2.3.10        
    ## [16] knitr_1.50          htmlwidgets_1.6.4   sp_2.2-1           
    ## [19] classInt_0.4-11     here_1.0.1          RColorBrewer_1.1-3 
    ## [22] KernSmooth_2.23-20  withr_3.0.3         numDeriv_2016.8-1.1
    ## [25] stats4_4.5.1        grid_4.5.1          leafem_0.2.5       
    ## [28] e1071_1.7-17        globals_0.19.1      scales_1.4.0       
    ## [31] iterators_1.0.14    cli_3.6.6           rmarkdown_2.30     
    ## [34] generics_0.1.4      otel_0.2.0          rstudioapi_0.17.1  
    ## [37] future.apply_1.9.1  tzdb_0.5.0          DBI_1.3.0          
    ## [40] proxy_0.4-29        base64enc_0.1-3     vctrs_0.7.3        
    ## [43] hms_1.1.4           listenv_0.8.0       crosstalk_1.2.2    
    ## [46] units_1.0-1         glue_1.8.1          parallelly_1.46.1  
    ## [49] codetools_0.2-18    stringi_1.8.7       gtable_0.3.6       
    ## [52] raster_3.6-32       pillar_1.11.1       htmltools_0.5.8.1  
    ## [55] satellite_1.0.4     R6_2.6.1            rprojroot_2.0.3    
    ## [58] evaluate_1.0.5      lattice_0.22-7      png_0.1-7          
    ## [61] class_7.3-20        Rcpp_1.1.2          coda_0.19-4.1      
    ## [64] nlme_3.1-160        xfun_0.51           pkgconfig_2.0.3
