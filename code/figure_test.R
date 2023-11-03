
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")

df_fcl <- readRDS("data_fmt/sim_fcl_main.rds") %>% 
  mutate(foodweb = factor(foodweb),
         rate = factor(rate),
         theta = factor(theta))


# figure ------------------------------------------------------------------

df_fcl %>% 
  ggplot(aes(x = p_branch,
             y = fcl,
             color = foodweb)) +
  geom_point(alpha = 0.1) +
  geom_smooth(se = FALSE, method = "lm") +
  facet_grid(rows = vars(rate), cols = vars(theta)) +
  guides(color = "none")

m <- lme4::lmer(log(fcl) ~ log(n_patch) + log(p_branch) + theta + rate + (1 | foodweb),
                data = df_fcl)
