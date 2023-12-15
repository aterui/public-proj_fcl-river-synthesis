
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")

df_fcl <- readRDS("data_fmt/sim_fcl_main.rds") %>% 
  mutate(foodweb = factor(foodweb),
         rate = factor(rate),
         phi = factor(phi),
         theta = factor(theta),
         k_base = factor(k_base))


# figure ------------------------------------------------------------------

df_fcl %>% 
  filter(k_base == 10) %>% 
  ggplot(aes(x = p_branch,
             y = fcl)) +
  geom_point(alpha = 0.1) +
  geom_smooth(se = TRUE, method = "loess") +
  facet_grid(rows = vars(rate), cols = vars(phi)) +
  guides(color = "none")# +
  #scale_x_continuous(trans = "log10") +
  #scale_y_continuous(trans = "log10")

m <- lme4::lmer(log(fcl) ~ 
                  log(n_patch) +
                  log(p_branch) + 
                  k_base + 
                  theta + 
                  rate + 
                  phi +
                  (1 | foodweb),
                data = df_fcl %>% dplyr::filter(phi != "0"))
