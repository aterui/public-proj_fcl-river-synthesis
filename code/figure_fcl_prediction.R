#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

df_fcl_mu <- readRDS("data_fmt/sim_fcl_analytical.rds") %>% 
  group_by(focus, p_branch, L, rsrc, mu) %>% 
  summarize(mu_fcl = mean(fcl),
            min_fcl = min(fcl),
            max_fcl = max(fcl)) %>% 
  ungroup()

df_fcl_raw <- readRDS("data_fmt/sim_fcl_analytical.rds")
  
# a -----------------------------------------------------------------------

source("code/set_theme.R")
ggplot2::theme_set(default_theme)

lab_mu <- c(`0.25` = "Stable",
            `2.5` = "Disturbed")

lab_r <- c(`0.25` = "Unproductive",
           `2.5` = "Productive")

g1 <- df_fcl_mu %>% 
  filter(focus == "p_branch") %>% 
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

g2 <- df_fcl_mu %>% 
  filter(focus == "size") %>% 
  ggplot(aes(x = L,
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


# export ------------------------------------------------------------------

ggsave(g1, filename = "output/fig_theory_branch.pdf",
       width = 7.25,
       height = 6)

ggsave(g2, filename = "output/fig_theory_size.pdf",
       width = 7.25,
       height = 6)

