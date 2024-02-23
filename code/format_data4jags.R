
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

## FCL data selected
df_fcl0 <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  as_tibble() %>% 
  dplyr::select(-geometry)

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
  filter(n_link > 19)
  
## - uid for those included
uid_incl <- df_g %>% 
  pull(uid)

## remove within-site replicates
## - take an average across seasons or years
## - censoring = 1 if top predator not collected
df_fcl <- df_fcl0 %>% 
  group_by(sid, wid, huid) %>% 
  summarize(fcl = mean(fcl),
            tpc = unique(top_predator_collected)) %>% 
  ungroup() %>% 
  mutate(uid = paste0(huid,
                      str_pad(wid, width = 5, pad = "0"))) %>% 
  left_join(df_env_local) %>% 
  filter(uid %in% uid_incl) %>% 
  mutate(g = as.numeric(factor(uid)),
         h = as.numeric(factor(huid)),
         censoring = ifelse(tpc == "N", 1, 0),
         fcl_min = ifelse(tpc == "N", fcl, NA)) %>% 
  mutate(fcl_max = ifelse(tpc == "N", NA, fcl)) %>% 
  relocate(uid, h, g, sid, fcl, fcl_max, fcl_min, tpc, censoring)

## - `df_bias` is to account for non-random spatial sampling
## - left join `distinct(df_fcl, uid, g)` to align group id `g` between the two data frames
## - CAUTION! Don't forget `arrange(g)`

df_bias <- df_fcl %>%
  rowwise() %>%
  mutate(y = sum(fcl, fcl_min, na.rm = T)) %>%
  group_by(uid) %>%
  summarize(mu_local_area = log(local_area) %>% 
              mean() %>% 
              exp(),
            mu_hfp = log(hfp) %>% 
              mean() %>% 
              exp())

df_g <- df_g %>% 
  left_join(distinct(df_fcl, uid, g)) %>%
  left_join(df_bias) %>% 
  arrange(g) %>%
  relocate(uid, g)


df