#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/format_data4jags.R")

## read data
df_est <- readRDS("data_fmt/output_model_summary.rds")
df_mu_est <- readRDS("data_fmt/output_model_mu_est.rds")
df_pred <- readRDS("data_fmt/output_model_pred.rds")
mcmc <- readRDS("data_fmt/output_model_mcmc.rds")


# format data -------------------------------------------------------------

## weight values for each observation
z <- df_est %>%
  filter(parms == "z") %>%
  pull(median)

df_mu_est <- df_mu_est %>%
  mutate(w = score^z / max(score^z))

# prediction plot ---------------------------------------------------------

df_mu_est %>% 
  ggplot(aes(x = p_branch,
             y = fcl_est,
             color = factor(h))) +
  geom_point(#aes(size = w),
             alpha = 0.4) + 
  geom_line(data = df_pred,
            aes(x = pb,
                y = y),
            alpha = 1) + 
  facet_wrap(facets = ~ factor(h)) +
  scale_x_continuous(trans = "log") +
  scale_y_continuous(trans = "log") +
  labs(y = "Food chain length",
       x = "Branching probability") +
  guides(size = "none")
