
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

df_fcl0 <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  as_tibble() %>% 
  dplyr::select(-geometry)

df_env <- readRDS("data_fmt/data_env.rds") %>% 
  rename(local_area = area) %>% 
  mutate(local_area = as.numeric(local_area),
         across(starts_with("mean."), .fns = function(x) round(x, 2)))

df_strnet <- readRDS("data_fmt/data_strnet.rds")

df_weight <- readRDS("data_fmt/data_weight.rds")

# format ------------------------------------------------------------------

## remove within-site replicates
## - take an average across seasons or years
## - censoring = 1 if top predator not collected
df_fcl <- df_fcl0 %>% 
  filter(sid %in% pull(df_env, sid)) %>% 
  group_by(sid, wid, huid) %>% 
  summarize(fcl = mean(fcl),
            tpc = unique(top_predator_collected)) %>% 
  ungroup() %>% 
  left_join(df_env) %>% 
  mutate(uid = paste0(huid,
                      str_pad(wid, width = 5, pad = "0")),
         g = as.numeric(factor(uid)),
         censoring = ifelse(tpc == "N", 1, 0),
         fcl_min = ifelse(tpc == "N", fcl, NA)) %>% 
  mutate(fcl = ifelse(tpc == "N", NA, fcl)) %>% 
  relocate(uid, g, sid, fcl, fcl_min, tpc, censoring)

df_g <- df_strnet %>% 
  mutate(uid = paste0(huid,
                      str_pad(wid, width = 5, pad = "0"))) %>% 
  left_join(df_weight) %>% 
  left_join(distinct(df_fcl, uid, g))
