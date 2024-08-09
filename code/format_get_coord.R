
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# read data ---------------------------------------------------------------

df_coord <- readRDS("data_fmt/wgs84_outlet.rds") %>% 
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>% 
  as_tibble() %>% 
  select(-c(id_col, pid, geometry))

write_csv(df_coord,
          "data_raw/data_outlet.csv")
