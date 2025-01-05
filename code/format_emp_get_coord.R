#' DESCRIPTION:
#' Get coordinates for sampling sites

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# read data ---------------------------------------------------------------

## snapped site coordinates
readRDS("data_raw/wgs84_subwatershed.rds") %>% 
  as_tibble() %>% 
  left_join(readRDS("data_raw/wgs84_fcl_site.rds") %>% 
              as_tibble() %>% 
              select(sid, study_id, org_site_id) %>% 
              filter(!duplicated(sid))) %>% 
  select(sid,
         study_id,
         org_site_id,
         lon = x,
         lat = y) %>% 
  write_csv("data_raw/data_site_snap.csv")

## outlet coordinates
readRDS("data_raw/wgs84_outlet.rds") %>% 
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>% 
  as_tibble() %>% 
  select(-c(id_col, pid, geometry)) %>% 
  write_csv("data_raw/data_outlet.csv")
