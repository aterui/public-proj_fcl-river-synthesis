
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# preformat ---------------------------------------------------------------
# update as needed
# FCL data selected
# - join level01 group id
sf_lev01 <- readRDS("data_fmt/wgs84_region_lev01.rds")

sf_fcl0 <- readRDS("data_fmt/wgs84_fcl_site.rds") %>%
  st_join(sf_lev01)

df_fcl0 <- sf_fcl0 %>%
  as_tibble() %>%
  dplyr::select(-geometry)

saveRDS(df_fcl0, "data_fmt/data_fcl.rds")


# read data ---------------------------------------------------------------

## local environment
df_fcl0 <- readRDS("data_fmt/data_fcl.rds")

## local environment
df_env_local <- readRDS("data_fmt/data_env_local.rds")

## network properties
df_env_wsd <- readRDS("data_fmt/data_env_wsd.rds")

## weight information for each watershed
df_weight <- readRDS("data_fmt/data_weight.rds")

# format ------------------------------------------------------------------

## watershed-level data frame
## - join `df_weight` for weighted regression
## - remove watersheds with less than 10 links; unreliable estimates of p_branch
df_g <- df_env_wsd %>% 
  mutate(uid = paste0(huid,
                      str_pad(wid, width = 5, pad = "0"))) %>% 
  left_join(df_weight) %>% 
  filter(n_link > 4)
  
## - uid for those included
uid_incl <- df_g %>% 
  pull(uid)

## remove within-site replicates
## - take an average across seasons or years
## - censoring = 1 if top predator not collected
df_fcl <- df_fcl0 %>% 
  group_by(sid, wid, huid, id_lev01) %>% 
  summarize(fcl = mean(fcl),
            tpc = unique(top_predator_collected)) %>% 
  ungroup() %>% 
  mutate(uid = paste0(huid,
                      str_pad(wid, width = 5, pad = "0"))) %>% 
  left_join(df_env_local) %>% 
  filter(uid %in% uid_incl) %>% 
  mutate(g = as.numeric(factor(uid)),
         h = as.numeric(factor(id_lev01)),
         censoring = ifelse(tpc == "N", 1, 0),
         fcl_min = ifelse(tpc == "N", fcl, NA)) %>% 
  mutate(fcl_obs = ifelse(tpc == "N", NA, fcl)) %>% 
  relocate(uid, h, g, sid, fcl, fcl_obs, fcl_min, tpc, censoring)

## - left join `distinct(df_fcl, uid, g)` to align group id `g` between the two data frames
## - CAUTION! Don't forget `arrange(g)`
df_g <- df_g %>% 
  left_join(distinct(df_fcl, uid, g, h)) %>%
  arrange(g) %>%
  relocate(uid, g)


# data frame for prediction -----------------------------------------------

## add scaled variables
df_g <- df_g %>% 
  mutate(scl_prec = c(scale(mean.prec)),
         scl_temp = c(scale(mean.temp)),
         scl_hfp = c(scale(hfp)))

## expanded data frame for prediction
df_x <- df_g %>% 
  group_by(h) %>% 
  reframe(x_log_area = mean(log(area)),
          x_log_pb = seq(log(min(p_branch)),
                         log(max(p_branch)),
                         length = 100),
          x_prec = mean(scl_prec),
          x_temp = mean(scl_temp),
          x_hfp = mean(scl_hfp))
