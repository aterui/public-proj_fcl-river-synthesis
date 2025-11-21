README
================

# Article Information

**Tentative title**: Geometric ecosystem complexity regulates food
chains

**Authors**: Terui A, Shibasaki S, Pomeranz JPF, Ibrahim M, Isaac E,
LaRoque A, Yamazaki D, Finlay JC

# Analysis Flow

## Theoretical Analysis

Script `run_theory.R` runs all analysis, including analytical and
numerical explorations.

For source codes, see:

- `analysis_theo_m_heatmap.R` (analytical)

- `analysis_theo_m_line.R` (analytical, for visualization)

- `analysis_theo_si.R` (numerical)

For visualization, see:

- `figure_theo_main.R` (main figures)

- `figure_theo_si.R` (supporting figures)

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
| `wgs84_fcl_site.rds` | Point layer for site coordinates (with estimates of food chain length; class `sf`). Include temporal duplicates (either seasonal or annual) at some sites. Some watersheds are removed before analysis through a screening process (see `format_emp_data4nimble.R`) |
| `wgs84_outlet.rds` | Outlet coordinates of study watersheds (class `sf`). Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`). |
| `wgs84_region_lev01.rds` | Polygon layer for HydroBASINS level-one (class `sf`). |
| `wgs84_str_sub.rds` | Polyline layer for study watersheds. Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`). |
| `wgs84_subwatershed.rds` | Polygon layer for upstream watershed at each sampling site (class `sf`). |
| `wgs84_wsd_sub.rds` | Polygon layer for study watersheds (class `sf`). Include some watersheds that are removed through a screening process (see `format_emp_data4nimble.R`). |

## Column Description

| file | column | description |
|:---|:---|:---|
| `data_env_flow.rds` | sid | Site ID. Re-defined for our analysis. |
|  | org_site_id | Site name used in the original article. |
|  | date | Date for discharge. |
|  | discharge | Discharge data simulated by CaMaFlood. |
| `data_env_local.rds` | sid | Site ID. Re-defined for our analysis. |
|  | local_area | Upstream watershed area at the sampling site. |
|  | unit_area | Unit of watershed area. |
|  | year_sampled | Sampling year of stable isotope data. |
|  | year_hfp | Year of human footprint index used. |
|  | hfp | Human footprint index. |
|  | local_temp | Average air temperature \[degree Celcius\] within the upstream watershed area. |
|  | local_prec | Average annual precipitation amount \[kg/(year m^2)\] within the upstream watershed area. |
|  | frac_agri | Fraction of agriculture within the upstream watershed area. |
|  | frac_urban | Fraction of urban within the upstream watershed area. |
|  | frac_forest | Fraction of forest within the upstream watershed area. |
|  | forest_b1km | Fraction of forest within 1km buffer of the sampling site. |
|  | utm | UTM code. |
|  | local_elev | Elevation at the sampling site. |
| `data_env_wsd.rds` | uid | Unique ID for the watershed. |
|  | wid | Watershed ID. Possible reduntant ID numbers between huid units. Use `uid` for grouping, which is a combination of `huid` and `wid`.) |
|  | tifid | GeoTIFF ID produced in the process of watershed analysis. Ignore. |
|  | oid | Outlet ID. |
|  | huid | Hydrological unit ID. Grouped by HydroBASINS level-4. Used only for computational purposes. |
|  | area | Watershed area. |
|  | unit_area | Unit of watershed area. |
|  | r_length | Total river length within the watershed. |
|  | lambda | Branching rate. |
|  | csize | Community size. Identical to the number of links within a watershed. |
|  | n_link | Number of links within the watershed. |
|  | unit_lambda | Unit of branching rate. |
|  | unit_r_length | Unit of river length. |
|  | p_branch | Branching probability. |
|  | year_sampled | Sampling year of stable isotope data. |
|  | year_hfp | Year of human footprint index used. |
|  | hfp | Human footprint index. |
|  | temp | Average air temperature \[degree Celcius\] within the watershed. |
|  | prec | Average annual precipitation amount \[kg/(year m^2)\] within the watershed. |
|  | key | GeoTIFF key used in watershed analysis. Ignore. |
|  | elev_mean | Mean elevation within the watershed. |
|  | elev_sd | SD elevation within the watershed. |
| `wgs84_fcl_site.rds` | uid | Unique ID for the watershed. |
|  | huid | Hydrological unit ID. Grouped by HydroBASINS level-4. Used only for computational purposes. |
|  | wid | Watershed ID. Possible reduntant ID numbers between huid units. Use `uid` for grouping, which is a combination of `huid` and `wid`.) |
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

<!-- Column `predator_collection_scheme` identifies the collection scheme as `A`, `B`, or `C` as described above. -->

<!-- In some cases, the author(s) of the original article provided detailed information to evaluate `Y` and `N`. -->

<!-- We used this information to determine the collection status of top predators, which is noted in the `note` column when applicable. -->

## Coordinate Precision

Coordinate precision was categorized as follows (**4** and **5** have
been excluded in `wgs84_fcl_site.rds`):

**1**: Exact coordinates available (either through article or author
correspondence); **2**: Nearly exact coordinates recovered from text;
**3**: Coordinates visually recovered from a map with reference to
landmarks/river morphology; **4**: Coordinates approximated by river
name or multiple sites aggregated; **5**: Not available.

## Geospatial Data

Geospatial data (RDS files with prefix `wgs84_`) were generated in a
[separate
repository](https://github.com/aterui/priv-proj_gis-global-fcl). The
data package used in this repository is **GFCL version 0.2.2**.

# Session Information

    ## R version 4.5.2 (2025-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 22631)
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
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] doSNOW_1.0.20     snow_0.4-4        doParallel_1.0.17 iterators_1.0.14 
    ##  [5] rpom_0.1.1        mcbrnet_1.4.2     foreach_1.5.2     ecotools_0.1.1   
    ##  [9] xtable_1.8-4      nimble_1.2.1      cowplot_1.2.0     magick_2.7.3     
    ## [13] igraph_2.1.1      terra_1.8-29      patchwork_1.3.2   sf_1.0-19        
    ## [17] lubridate_1.9.3   forcats_1.0.0     stringr_1.6.0     dplyr_1.1.4      
    ## [21] purrr_1.2.0       readr_2.1.5       tidyr_1.3.1       tibble_3.3.0     
    ## [25] ggplot2_4.0.1     tidyverse_2.0.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6        xfun_0.51           lattice_0.22-7     
    ##  [4] tzdb_0.5.0          numDeriv_2016.8-1.1 vctrs_0.6.5        
    ##  [7] tools_4.5.2         generics_0.1.4      proxy_0.4-27       
    ## [10] pacman_0.5.1        pkgconfig_2.0.3     KernSmooth_2.23-20 
    ## [13] RColorBrewer_1.1-3  S7_0.2.0            lifecycle_1.0.4    
    ## [16] compiler_4.5.2      farver_2.1.2        codetools_0.2-18   
    ## [19] htmltools_0.5.8.1   class_7.3-20        yaml_2.3.10        
    ## [22] pracma_2.4.4        pillar_1.11.1       classInt_0.4-11    
    ## [25] tidyselect_1.2.1    digest_0.6.33       stringi_1.8.7      
    ## [28] rprojroot_2.0.3     fastmap_1.2.0       grid_4.5.2         
    ## [31] here_1.0.1          cli_3.6.5           magrittr_2.0.4     
    ## [34] e1071_1.7-16        withr_3.0.2         scales_1.4.0       
    ## [37] timechange_0.3.0    rmarkdown_2.28      hms_1.1.3          
    ## [40] coda_0.19-4.1       evaluate_1.0.5      knitr_1.49         
    ## [43] rlang_1.1.6         Rcpp_1.1.0          glue_1.8.0         
    ## [46] DBI_1.2.3           rstudioapi_0.14     R6_2.6.1           
    ## [49] units_0.8-7
