
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

df_fcl <- readRDS("data_fmt/sim_fcl_two_patch.rds")


# figure ------------------------------------------------------------------

g_heat <- df_fcl %>% 
  group_by(phi, k, threshold, theta) %>% 
  summarize(fcl = mean(fcl)) %>% 
  ungroup() %>% 
  ggplot(aes(x = phi,
             y = k,
             fill = fcl)) +
  geom_raster() +
  facet_grid(rows = vars(theta),
             cols = vars(threshold),
             labeller = label_both) +
  MetBrewer::scale_fill_met_c("Hiroshige", direction = -1) +
  theme_classic()

ggsave(g_heat,
       filename = "output/figure_two_patch.pdf",
       height = 7,
       width = 8)
