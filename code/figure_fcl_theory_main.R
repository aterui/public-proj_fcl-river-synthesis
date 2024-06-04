#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

# data --------------------------------------------------------------------

## read data
df_fcl_mu <- readRDS("data_fmt/sim_fcl_main.rds") %>% 
  group_by(focus, lambda, rl, rsrc, mu0) %>% 
  summarize(mu_fcl = mean(fcl),
            min_fcl = min(fcl),
            max_fcl = max(fcl)) %>% 
  ungroup() %>% 
  mutate(lab_r = ifelse(rsrc == min(rsrc),
                        sprintf("italic(r)[0]==%.2f~~Resource~low", rsrc),
                        sprintf("italic(r)[0]==%.2f~~Resource~high", rsrc)),
         lab_mu = ifelse(mu0 == min(mu0),
                         sprintf("mu^{(0)}==%.2f~~Disturbance~infreq.", mu0),
                         sprintf("mu^{(0)}==%.2f~~Disturbance~freq.", mu0)))

# fcl vs size or p_branch plot --------------------------------------------

## set plot theme
source("code/set_theme.R")
ggplot2::theme_set(default_theme)

## food chain length vs. branching
g1 <- df_fcl_mu %>% 
  filter(focus == "branch") %>% 
  ggplot(aes(x = lambda,
             y = mu_fcl,
             color = factor(rl),
             fill = factor(rl))) + 
  geom_line(linewidth = 1,
            alpha = 0.8) +
  geom_ribbon(aes(ymax = max_fcl,
                  ymin = min_fcl),
              alpha = 0.1,
              col = NA) +
  facet_grid(rows = vars(lab_mu),
             cols = vars(lab_r),
             labeller = label_parsed) +
  labs(x = expression("Branching rate"~lambda[b]),
       y = "Food chain length",
       color = expression("Ecosystem size"~italic(L)),
       fill = expression("Ecosystem size"~italic(L))) +
  MetBrewer::scale_color_met_d("Hiroshige") +
  MetBrewer::scale_fill_met_d("Hiroshige") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

## food chain length vs. size
g2 <- df_fcl_mu %>% 
  filter(focus == "size") %>% 
  ggplot(aes(x = rl,
             y = mu_fcl,
             color = factor(lambda),
             fill = factor(lambda))) + 
  geom_line(linewidth = 1,
            alpha = 0.8) +
  geom_ribbon(aes(ymax = max_fcl,
                  ymin = min_fcl),
              alpha = 0.1,
              col = NA) +
  facet_grid(rows = vars(lab_mu),
             cols = vars(lab_r),
             labeller = label_parsed) +
  labs(x = expression("River length"~italic(L)),
       y = "Food chain length",
       color = expression("Branching rate"~lambda[b]),
       fill = expression("Branching rate"~lambda[b])) +
  MetBrewer::scale_color_met_d("Hiroshige") +
  MetBrewer::scale_fill_met_d("Hiroshige") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))


## export 
ggsave(g1, filename = "output/fig_theory_branch.pdf",
       width = 10,
       height = 7.5)

ggsave(g2, filename = "output/fig_theory_size.pdf",
       width = 10,
       height = 7.5)


# upstream distance function ----------------------------------------------

## set parameters
rl <- 1:100
lambda <- seq(0.1, 0.9, by = 0.2)
x_values <- expand.grid(rl = rl, lambda = lambda)

u <- sapply(seq_len(nrow(x_values)), function(i) {
  with(x_values,
       u_length(size = rl[i], lambda = lambda[i]))
})

## get data frame for upstream distance
df_u <- bind_cols(x_values, y = u) %>% 
  as_tibble() %>% 
  setNames(c("rl", "lambda", "u"))

ug0 <- df_u %>% 
  filter(lambda == 0.1) %>% 
  ggplot(aes(x = rl,
             y = u)) +
  geom_line(linewidth = 1.5,
            color = grey(0.3)) +
  labs(y = "E(Upstream length)",
       x = "River length") +
  guides(color = "none")

ug <- df_u %>% 
  ggplot(aes(x = rl,
             y = u,
             color = factor(lambda))) +
  geom_line(linewidth = 1.5) +
  MetBrewer::scale_color_met_d("Hiroshige",
                               direction = -1) +
  labs(y = "E(Upstream length)",
       x = "River length") +
  guides(color = "none")

## export
ggsave(ug0, 
       filename = "output/fig_ud0.pdf",
       width = 6, height = 5)

ggsave(ug, 
       filename = "output/fig_ud.pdf",
       width = 6, height = 5)