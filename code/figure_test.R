
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

df_fcl <- readRDS("data_fmt/sim_fcl_main.rds") %>% 
  mutate(foodweb = factor(foodweb),
         rate = factor(rate),
         phi = factor(phi),
         theta = factor(theta),
         k_base = factor(k_base))


# figure ------------------------------------------------------------------

df_fcl %>% 
  ggplot(aes(x = n_patch,
             y = fcl)) +
  #geom_point(alpha = 0.1) +
  geom_smooth(se = T, method = "gam") +
  facet_grid(rows = vars(rate, theta),
             cols = vars(phi),
             scales = "free",
             labeller = label_both) +
  guides(color = "none") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme_bw()

m <- lme4::lmer(log(fcl) ~ 
                  log(n_patch) * rate +
                  log(p_branch) * rate + 
                  log(n_patch) * phi +
                  log(p_branch) * phi + 
                  theta + 
                  phi +
                  (1 | foodweb),
                data = df_fcl,
                REML = FALSE)

library(MuMIn)
options(na.action = "na.fail")

ms <- dredge(m)
