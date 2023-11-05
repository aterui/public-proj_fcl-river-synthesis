
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")

df_fcl <- readRDS("data_fmt/sim_fcl_main.rds") %>% 
  mutate(foodweb = factor(foodweb),
         rate = factor(rate),
         phi = factor(phi),
         theta = factor(theta))


# figure ------------------------------------------------------------------

df_fcl %>% 
  filter(n_patch > 30) %>% 
  ggplot(aes(x = p_branch,
             y = fcl,
             size = n_patch)) +
  geom_point(alpha = 0.1) +
  geom_smooth(se = FALSE, method = "lm") +
  facet_grid(rows = vars(rate), cols = vars(phi)) +
  guides(color = "none") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")

m <- lme4::lmer(log(fcl) ~ log(n_patch) + log(p_branch) + theta + rate + phi + (1 | foodweb),
                data = df_fcl)
