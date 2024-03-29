#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# fcl vs size or p_branch plot --------------------------------------------

## read data
df_fcl_mu <- readRDS("data_fmt/sim_fcl_analytical.rds") %>% 
  group_by(focus, lambda, rl, rsrc, mu) %>% 
  summarize(mu_fcl = mean(fcl),
            min_fcl = min(fcl),
            max_fcl = max(fcl)) %>% 
  ungroup() %>% 
  mutate(p_branch = 1 - exp(-lambda))

## set plot theme
source("code/set_theme.R")
ggplot2::theme_set(default_theme)

## label mapper
lab_mu <- c(`0.25` = "Stable",
            `2.5` = "Disturbed")

lab_r <- c(`0.25` = "Unproductive",
           `2.5` = "Productive")

## food chain length vs. branching
g1 <- df_fcl_mu %>% 
  filter(focus == "branch") %>% 
  ggplot(aes(x = p_branch,
             y = mu_fcl)) + 
  geom_line(linewidth = 1,
            alpha = 0.8,
            color = "steelblue") +
  geom_ribbon(aes(ymax = max_fcl,
                  ymin = min_fcl),
              fill = "steelblue",
              alpha = 0.1,
              col = NA) +
  facet_grid(rows = vars(mu),
             cols = vars(rsrc),
             labeller = labeller(mu = lab_mu,
                                 rsrc = lab_r)) +
  labs(x = "Ecosystem complexity (branching prob.)",
       y = "Food chain length",
       color = "Ecosystem size") +
  MetBrewer::scale_color_met_d("Hiroshige") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16)) +
  guides(color = "none")

## food chain length vs. size
g2 <- df_fcl_mu %>% 
  filter(focus == "size") %>% 
  ggplot(aes(x = rl,
             y = mu_fcl)) + 
  geom_line(linewidth = 1,
            alpha = 0.8,
            color = "steelblue") +
  geom_ribbon(aes(ymax = max_fcl,
                  ymin = min_fcl),
              fill = "steelblue",
              alpha = 0.1,
              col = NA) +
  facet_grid(rows = vars(mu),
             cols = vars(rsrc),
             labeller = labeller(mu = lab_mu,
                                 rsrc = lab_r)) +
  labs(x = "Ecosystem size (river length)",
       y = "Food chain length") +
  MetBrewer::scale_color_met_d("Hiroshige") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16)) +
  guides(color = "none")


## export 

ggsave(g1, filename = "output/fig_theory_branch.pdf",
       width = 7.25,
       height = 6)

ggsave(g2, filename = "output/fig_theory_size.pdf",
       width = 7.25,
       height = 6)


# upstream distance function ----------------------------------------------

## set parameters
rl <- 1:100
lambda <- -log(1 - seq(0.2, 0.8, by = 0.2))
x_values <- expand.grid(rl = rl, lambda = lambda)

u_dist <- sapply(seq_len(nrow(x_values)), function(i) {
  with(x_values,
       u_length(L = rl[i], lambda = lambda[i]))
})

## get data frame for upstream distance
df_u <- bind_cols(x_values, u_dist) %>% 
  as_tibble() %>% 
  setNames(c("rl", "lambda", "u_dist")) %>% 
  mutate(p_branch = round(1 - exp(-lambda), 1))

ug0 <- df_u %>% 
  filter(p_branch == 0.2) %>% 
  ggplot(aes(x = rl,
             y = u_dist)) +
  geom_line(linewidth = 1.5,
            color = grey(0.3)) +
  labs(y = "E(Upstream distance)",
       x = "River length") +
  guides(color = "none")

ug <- df_u %>% 
  ggplot(aes(x = rl,
             y = u_dist,
             color = factor(p_branch))) +
  geom_line(linewidth = 1.5) +
  MetBrewer::scale_color_met_d("Hiroshige",
                               direction = -1) +
  labs(y = "E(Upstream distance)",
       x = "River length") +
  guides(color = "none")

## export
ggsave(ug0, 
       filename = "output/figure_ud0.pdf",
       width = 6, height = 5)

ggsave(ug, 
       filename = "output/figure_ud.pdf",
       width = 6, height = 5)