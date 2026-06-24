#' DESCRIPTION:
#' Format data for analysis in NIMBLE

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

# preformat ---------------------------------------------------------------
# update as needed
# FCL data selected
# - join level01 group id; continental hydrological units
sf_lev01 <- readRDS("data_raw/sf_hybas_lev01.rds") %>% 
  st_make_valid() %>% 
  rename(hybas_lev01 = hybas)

sf_lev02 <- readRDS("data_raw/sf_hybas_lev02.rds") %>% 
  st_make_valid() %>% 
  filter(hybas != "lev02_5020055870") %>% # invalid polygon
  rename(hybas_lev02 = hybas)

sf_fcl0 <- readRDS("data_raw/sf_fcl_site.rds") %>%
  st_join(sf_lev01) %>% 
  st_join(sf_lev02) %>% 
  relocate(starts_with("hybas"))

df_fcl0 <- sf_fcl0 %>%
  as_tibble() %>%
  dplyr::select(-geometry)

# read data ---------------------------------------------------------------

## local environment
df_env_local0 <- readRDS("data_raw/data_env_local.rds")
df_env_fsd <- readRDS("data_fmt/data_env_fsd.rds")
df_env_local <- left_join(df_env_local0, df_env_fsd)

## network properties
df_env_wsd <- readRDS("data_raw/data_env_wsd.rds")

## weight information for each watershed
df_weight <- readRDS("data_fmt/data_weight.rds")

## flag watersheds with no flow data
df_flag <- df_env_local %>%
  left_join(df_fcl0 %>% 
              select(sid, oid)) %>% 
  group_by(oid) %>% 
  summarize(flag_flow = ifelse(all(is.na(fsd)), "Y", "N"))


# format ------------------------------------------------------------------

## watershed-level data frame
## - join `df_weight` for weighted regression
## - remove watersheds with less than 5 links; unreliable estimates of p_branch
## - remove watersheds with no human footprint
df_g <- df_env_wsd %>% 
  left_join(df_weight) %>% 
  left_join(df_flag) %>% 
  filter(n_link > 4) %>% 
  drop_na(hfp)

## - uid for those included
oid_incl <- df_g %>% 
  pull(oid)

## local-level data
## - remove within-site replicates
## - take an average across seasons or years
## - censoring = 1 if top predator not collected (i.e., tpc == "N")
df_fcl_local <- df_fcl0 %>% 
  group_by(sid, 
           oid, 
           huid, 
           hybas_lev01, 
           hybas_lev02) %>% 
  summarize(fcl = mean(fcl),
            tpc = unique(top_predator_collected),
            .groups = "drop") %>% 
  left_join(df_env_local) %>% 
  filter(oid %in% oid_incl) %>% 
  mutate(g = as.numeric(factor(oid)),
         h01 = as.numeric(factor(hybas_lev01)),
         h02 = as.numeric(factor(hybas_lev02)),
         censoring = ifelse(tpc == "N", 1, 0),
         cut = ifelse(tpc == "N", fcl, Inf)) %>% 
  mutate(fcl_obs = ifelse(tpc == "N", NA, fcl)) %>% 
  relocate(
    oid, 
    h01,
    h02,
    g,
    sid,
    fcl,
    fcl_obs, 
    cut,
    tpc,
    censoring
  )

## watershed-level data
## - left join `distinct(df_fcl, uid, g)` to align group id `g` between the two data frames
## - CAUTION! Don't forget `arrange(g)`
df_fcl_wsd <- df_g %>% 
  left_join(
    distinct(df_fcl_local, 
             oid, 
             h01,
             h02,
             g)
  ) %>%
  arrange(g) %>%
  relocate(
    oid,
    g
  )

## check g has a proper vector
z <- mean(df_fcl_wsd$g == seq_len(nrow(df_fcl_wsd)))
if (z != 1) stop("Error in data formatting")


# export ------------------------------------------------------------------

list_fcl <- list(df_fcl_local, df_fcl_wsd)
saveRDS(list_fcl, "data_fmt/data_fcl_reg.rds")


# for visual check (not for analysis) -------------------------------------
# 
# df_viz <- df_fcl_local %>%
#   left_join(df_g %>%
#               select(oid, lambda, r_length, prec, temp)) %>%
#   mutate(y = ifelse(fnu > 100, 0, 1))
# 
# 
# df_viz %>%
#   group_by(oid) %>%
#   mutate(n = n()) %>%
#   ungroup() %>%
#   filter(n >= 5) %>%
#   ggplot(aes(y = local_area,
#              x = elev,
#              color = tpc)) +
#   geom_point() +
#   facet_wrap(facets =~ oid,
#              scales = "free")
