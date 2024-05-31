#' DESCRIPTION
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")


# data --------------------------------------------------------------------

## read data
df_fcl_sim <- readRDS("data_fmt/sim_fcl_si.rds") %>% 
  as_tibble()

## take average for food web replicates
df_plot <- df_fcl_sim %>% 
  group_by(rl, lambda, rsrc, mu0, mu_p, mu_c, rho, g, delta0) %>% 
  summarize(fcl = mean(fcl)) %>% 
  ungroup()


# figure ------------------------------------------------------------------

## w/o top-down effects
df_plot %>% 
  filter(mu_c == 0,
         g == 10) %>% 
  ggplot(aes(y = factor(rl),
             x = factor(lambda),
             fill = fcl)) +
  geom_raster(alpha = 0.8) +
  facet_grid(rows = vars(mu_p, mu0),
             cols = vars(rho, rsrc)) +
  MetBrewer::scale_fill_met_c("Hiroshige",
                              direction = -1)

## w top-down effects
df_plot %>% 
  filter(mu_c == 1,
         g == 100) %>% 
  ggplot(aes(y = factor(rl),
             x = factor(lambda),
             fill = fcl)) +
  geom_raster(alpha = 0.8) +
  facet_grid(rows = vars(mu_p, mu0),
             cols = vars(rho, rsrc)) +
  MetBrewer::scale_fill_met_c("Hiroshige",
                              direction = -1)
