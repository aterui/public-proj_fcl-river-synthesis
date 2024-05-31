#' DESCRIPTION
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

# data --------------------------------------------------------------------

## read data
df_fcl_sim <- readRDS("data_fmt/sim_fcl_si.rds")

## take average for food web replicates
df_plot <- df_fcl_sim %>% 
  group_by(rl, lambda, h, delta0, rsrc, g, mu0, mu_p, mu_c, rho, theta) %>% 
  summarize(fcl = mean(fcl)) %>% 
  ungroup()

df_parms <- df_plot %>% 
  distinct(rho, g, theta)

# figure ------------------------------------------------------------------

list_g <- foreach(i = 1:nrow(df_parms)) %do% {
  
  df_plot %>% 
    filter(rho == df_parms$rho[i],
           g == df_parms$g[i],
           theta == df_parms$theta[i]) %>% 
    ggplot(aes(y = rl,
               x = lambda,
               fill = fcl)) +
    geom_raster(alpha = 1) +
    facet_grid(rows = vars(mu_p, mu0),
               cols = vars(mu_c, rsrc)) +
    MetBrewer::scale_fill_met_c("Hiroshige",
                                direction = -1)
  
}
