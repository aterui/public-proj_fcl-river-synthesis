
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# read data ---------------------------------------------------------------

sf_subwsd <- readRDS("data_fmt/wgs84_subwatershed.rds")

sf_site <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  group_by(sid) %>% 
  slice(1) %>% 
  select(NULL)

df_site_snap <- sf_subwsd %>% 
  as_tibble() %>% 
  select(sid,
         lon = x,
         lat = y)

write_csv(df_site_snap,
          "data_raw/data_site_snap.csv")