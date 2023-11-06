
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")

df_fcl <- readRDS("data_fmt/sim_fcl_two_patch.rds")


# figure ------------------------------------------------------------------

df_fcl %>% 
  group_by(phi, m, threshold, theta) %>% 
  summarize(fcl = mean(fcl)) %>% 
  ungroup() %>% 
  ggplot(aes(x = phi,
             y = m,
             fill = fcl)) +
  geom_raster() +
  facet_grid(rows = vars(theta),
             cols = vars(threshold),
             labeller = label_both) +
  MetBrewer::scale_fill_met_c("Hiroshige", direction = -1)

df_fcl %>% 
  ggplot(aes(x = factor(phi),
             y = fcl,
             color = factor(m))) +
  geom_jitter(alpha = 0.4) +
  geom_boxplot(alpha = 0.25, outlier.color = NA) +
  facet_grid(rows = vars(m),
             cols = vars(theta, threshold))
