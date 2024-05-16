
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# preformat ---------------------------------------------------------------
# update as needed
# FCL data selected
# - join level01 group id
sf_lev01 <- readRDS("data_fmt/wgs84_region_lev01.rds") %>% 
  st_make_valid()

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
df_env_local0 <- readRDS("data_fmt/data_env_local.rds")
df_env_flow_var <- readRDS("data_fmt/data_env_flow_var.rds")

df_env_local <- left_join(df_env_local0, df_env_flow_var)

## network properties
df_env_wsd <- readRDS("data_fmt/data_env_wsd.rds")
df_channel <- readRDS("data_fmt/wgs84_str_sub.rds") %>% 
  lapply(function(x) {
    if(!is.null(x)) {
      x %>% 
        mutate(uid = paste0(huid,
                            str_pad(wid,
                                    width = 5,
                                    pad = "0")),
               length = units::set_units(length(.), "km")) %>% 
        group_by(uid) %>% 
        summarize(r_length = sum(length)) %>% 
        as_tibble() %>% 
        transmute(uid,
                  r_length = units::drop_units(r_length),
                  unit_length = "km")
    }
  }) %>% 
  bind_rows()

df_env_wsd <- left_join(df_env_wsd,
                        df_channel)

## weight information for each watershed
df_weight <- readRDS("data_fmt/data_weight.rds")

# format ------------------------------------------------------------------

## watershed-level data frame
## - join `df_weight` for weighted regression
## - remove watersheds with less than 5 links; unreliable estimates of p_branch
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
## - censoring = 1 if top predator not collected (i.e., tpc == "N")
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

## check g has a proper vector
z <- mean(df_g$g == seq_len(nrow(df_g)))
if (z != 1) stop("Error in data formatting")
