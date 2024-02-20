
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

df_fcl0 <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  as_tibble() %>% 
  dplyr::select(-geometry)

df_env <- readRDS("data_fmt/data_env.rds") %>% 
  rename(local_area = area) %>% 
  mutate(local_area = as.numeric(local_area))

df_strnet <- readRDS("data_fmt/data_strnet.rds")


# format ------------------------------------------------------------------

df_fcl <- df_fcl0 %>% 
  filter(sid %in% pull(df_env, sid)) %>% 
  group_by(sid, wid, huid) %>% 
  summarize(fcl = mean(fcl),
            tpc = unique(top_predator_collected)) %>% 
  ungroup() %>% 
  left_join(df_env)

ggplot(df_fcl,
       aes(x = mean.agri,
           y = fcl,
           color = huid)) +
  geom_point() +
  facet_wrap(facets = ~huid) +
  #scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")
