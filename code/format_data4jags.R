
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

## FCL data selected
df_fcl0 <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  as_tibble() %>% 
  dplyr::select(-geometry)

## local environment
df_env <- readRDS("data_fmt/data_env_local.rds") %>% 
  rename(local_area = area) %>% 
  mutate(local_area = as.numeric(local_area),
         across(starts_with("mean."),
                .fns = function(x) round(x, 2)))

## network properties
df_strnet <- readRDS("data_fmt/data_strnet.rds")

## weight information for each watershed
df_weight <- readRDS("data_fmt/data_weight.rds")

# format ------------------------------------------------------------------

## watershed-level data frame
## - join `df_weight` for weighted regression
## - remove watersheds with less than 10 links; unreliable estimates of p_branch
df_g <- df_strnet %>% 
  mutate(uid = paste0(huid,
                      str_pad(wid, width = 5, pad = "0"))) %>% 
  left_join(df_weight) %>% 
  filter(n_link > 9)
  
## - uid for those included
uid_incl <- df_g %>% 
  pull(uid)

## remove within-site replicates
## - take an average across seasons or years
## - censoring = 1 if top predator not collected
df_fcl <- df_fcl0 %>% 
  filter(sid %in% pull(df_env, sid)) %>% 
  group_by(sid, wid, huid) %>% 
  summarize(fcl = mean(fcl),
            tpc = unique(top_predator_collected)) %>% 
  ungroup() %>% 
  mutate(uid = paste0(huid,
                      str_pad(wid, width = 5, pad = "0"))) %>% 
  left_join(df_env) %>% 
  filter(uid %in% uid_incl) %>% 
  mutate(g = as.numeric(factor(uid)),
         h = as.numeric(factor(huid)),
         censoring = ifelse(tpc == "N", 1, 0),
         fcl_min = ifelse(tpc == "N", fcl, NA)) %>% 
  mutate(fcl = ifelse(tpc == "N", NA, fcl)) %>% 
  relocate(uid, h, g, sid, fcl, fcl_min, tpc, censoring)

## - left join `distinct(df_fcl, uid, g)` to align group id `g` between the two data frames
## - CAUTION! Don't forget `arrange(g)`
df_g <- df_g %>% 
  left_join(distinct(df_fcl, uid, g)) %>%
  arrange(g) %>%
  relocate(uid, g)
