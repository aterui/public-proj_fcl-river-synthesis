
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# read data ---------------------------------------------------------------

sf_subwsd <- readRDS("data_fmt/wgs84_subwatershed.rds")

df_site <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  group_by(sid) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(sid, study_id, org_site_id) %>% 
  mutate(lon0 = st_coordinates(.)[,1], lat0 = st_coordinates(.)[,2]) %>% 
  as_tibble() %>% 
  select(-geometry)

df_site_snap <- sf_subwsd %>% 
  as_tibble() %>% 
  select(sid,
         lon = x,
         lat = y) %>% 
  left_join(df_site, by = "sid") %>% 
  relocate(sid, study_id, org_site_id)

write_csv(df_site_snap,
          "data_raw/data_site_snap.csv")