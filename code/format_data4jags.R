
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

## FCL data selected
## - join level01 group id
sf_lev01 <- readRDS("data_fmt/wgs84_region_lev01.rds") %>% 
  st_make_valid()

sf_fcl0 <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  st_join(sf_lev01)

df_fcl0 <- sf_fcl0 %>% 
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
  filter(n_link > 5)
  
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

## - `df_bias` is to account for non-random spatial sampling
## - left join `distinct(df_fcl, uid, g)` to align group id `g` between the two data frames
## - CAUTION! Don't forget `arrange(g)`

df_bias <- df_fcl %>%
  group_by(uid) %>%
  summarize(mu_local_area = log(local_area) %>% 
              mean() %>% 
              exp(),
            mu_hfp = log(hfp) %>% 
              mean() %>% 
              exp())

df_g <- df_g %>% 
  left_join(distinct(df_fcl, uid, g, h)) %>%
  left_join(df_bias) %>% 
  arrange(g) %>%
  relocate(uid, g)


# test plot ---------------------------------------------------------------
# 
# df_y <- df_fcl %>% 
#   left_join(df_g, by = "uid") %>% 
#   group_by(id_lev01, uid) %>% 
#   summarize(fcl = mean(fcl),
#             p_branch = unique(p_branch),
#             area = unique(area),
#             w = unique(w)) %>% 
#   ungroup() %>% 
#   left_join(df_bias)  
# 
# df_y %>% 
#   ggplot(aes(x = p_branch,
#              y = fcl,
#              size = w)) +
#   geom_point() +
#   facet_wrap(facets = ~id_lev01) +
#   scale_x_continuous(trans = "log10") +
#   scale_y_continuous(trans = "log10")