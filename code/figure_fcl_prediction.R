#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

df_fcl <- readRDS("data_fmt/sim_fcl_analytical.rds")


# a -----------------------------------------------------------------------

source("code/set_theme.R")
ggplot2::theme_set(default_theme)

lab_mu <- c(`0.25` = "Stable",
            `2.5` = "Disturbed")

lab_r <- c(`0.25` = "Unproductive",
           `2.5` = "Productive")

g1 <- df_fcl %>% 
  mutate(L = factor(L)) %>% 
  filter(focus == "p_branch") %>% 
  ggplot(aes(x = p_branch,
             y = fcl,
             color = factor(L))) + 
  geom_line(linewidth = 2,
            alpha = 0.8) +
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

g2 <- df_fcl %>% 
  filter(focus == "size") %>% 
  ggplot(aes(x = L,
             y = fcl,
             color = factor(p_branch))) + 
  geom_line(linewidth = 2,
            alpha = 0.8) +
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
       width = 7.5,
       height = 6)

ggsave(g2, filename = "output/fig_theory_size.pdf",
       width = 7,
       height = 6)

