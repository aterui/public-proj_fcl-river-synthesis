#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

## set plot theme
source("code/set_theme.R")
ggplot2::theme_set(default_theme)

# figure 1 ----------------------------------------------------------------

## draw scheme
#img <- image_read("data_fmt/subfig_scheme.png")
img <- image_read_pdf("data_fmt/figsub_scheme.pdf")

g_subsch <- ggdraw() +
  draw_image(img)

## upstream length function
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

# ug0 <- df_u %>% 
#   filter(lambda == 0.1) %>% 
#   ggplot(aes(x = rl,
#              y = u)) +
#   geom_line(linewidth = 1.5,
#             color = grey(0.3)) +
#   labs(y = "E(Upstream length)",
#        x = "River length") +
#   guides(color = "none")

g_uhat <- df_u %>% 
  ggplot(aes(x = rl,
             y = u,
             color = factor(lambda))) +
  geom_line(linewidth = 1.5) +
  MetBrewer::scale_color_met_d("Hiroshige",
                               direction = -1) +
  labs(y = expression("Upstream river length"~~hat(italic(u))),
       x = expression("River length"~italic(L)),
       color = expression("Branching rate"~lambda[b])) + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.2, 0.7),
        legend.background = element_blank())

## layout ####
(g_scheme <- g_subsch + g_uhat +
  plot_layout(width = c(1.5, 1)) +
  plot_annotation(tag_levels = "A"))

ggsave(g_scheme,
       filename = "data_fmt/fig_scheme.pdf",
       width = 9,
       height = 4)

# figure 2 ----------------------------------------------------------------

## heatmap ####

## read data for heatmap
## take average for food web replicates
df_heat <- readRDS("data_fmt/sim_fcl_m_heat.rds") %>% 
  group_by(rl, lambda, h, delta0, rsrc, g, mu0, mu_p, mu_c, rho, theta) %>%
  summarize(fcl = max(fcl)) %>%
  ungroup() %>%
  filter(theta == min(theta)) %>% 
  mutate(lab_mu0 = ifelse(mu0 > min(mu0),
                          sprintf("mu^{(0)}==%.2f~(freq.)", mu0),
                          sprintf("mu^{(0)}==%.2f~(infreq.)", mu0)),
         lab_rsrc = ifelse(rsrc > min(rsrc),
                           sprintf("italic(r[0])==%.2f~(high)", rsrc),
                           sprintf("italic(r[0])==%.2f~(low)", rsrc)),
         lab_rho = ifelse(rho > min(rho),
                          sprintf("rho==%.2f~(disturbance~cascade)", rho),
                          sprintf("rho==%.2f~(no~disturbance~cascade)", rho)),
         lab_d = "Disturbance~rate",
         lab_r = "Resource~supply",
         lambda_mid = ifelse(rsrc == max(rsrc) & mu0 == max(mu0),
                             mean(range(lambda)),
                             NA),
         rl_mid = ifelse(rsrc == max(rsrc) & mu0 == max(mu0),
                         mean(range(rl)),
                         NA)
  )

heatmap <- function(data) {
  data %>% 
    ggplot(aes(y = rl,
               x = lambda,
               fill = fcl)) +
    geom_raster(alpha = 1) +
    ggh4x::facet_nested(rows = vars(lab_d, lab_mu0),
                        cols = vars(lab_r, lab_rsrc),
                        labeller = label_parsed,
                        nest_line = element_line(linetype = 2)) +
    MetBrewer::scale_fill_met_c("Hiroshige",
                                direction = -1) +
    geom_vline(aes(xintercept = lambda_mid),
               color = "white",
               linetype = "dashed") +
    geom_hline(aes(yintercept = rl_mid),
               color = "white",
               linetype = "dashed") +
    labs(x = expression("Branching rate"~lambda[b]),
         y = expression("River length"~italic(L)),
         fill = "FCL") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 12))
}

(g_heat05 <- heatmap(data = df_heat %>% filter(rho == 0.5)))
(g_heat0 <- heatmap(data = df_heat %>% filter(rho == 0)))

## lineart ####

## read data
df_fcl_line <- readRDS("data_fmt/sim_fcl_m_line.rds") %>% 
  filter(rsrc == max(rsrc),
         mu0 == max(mu0))

df_fcl_mu <- df_fcl_line %>% 
  group_by(focus, lambda, rl, rsrc, mu0, rho) %>% 
  summarize(mu_fcl = mean(fcl),
            min_fcl = min(fcl),
            max_fcl = max(fcl)) %>% 
  ungroup()

lineart <- function(data1, data2,
                    x,
                    x_axis) {
  data1 %>% 
    filter(focus == x) %>%
    rename(x0 = all_of(x_axis)) %>% 
    ggplot(aes(x = x0)) + 
    geom_line(aes(y = fcl,
                  color = factor(rho),
                  linetype = factor(fw)),
              linewidth = 0.1,
              alpha = 0.3) +
    scale_linetype_manual(values = rep("solid",
                                       n_distinct(data1$fw))) +
    geom_line(data = data2 %>% 
                rename(x0 = all_of(x_axis)) %>% 
                filter(focus == x),
              aes(y = mu_fcl,
                  color = factor(rho)),
              linetype = "solid",
              linewidth = 1) +
    guides(linetype = "none") +
    scale_color_manual(values = c(`0` = "grey",
                                  `0.5` = "tomato"),
                       labels = c(`0` = "No cascade",
                                  `0.5` = "Cascade"))
}

## food chain length vs. branching
g_br <- lineart(df_fcl_line,
                df_fcl_mu,
                x = "branch",
                x_axis = "lambda") +
  labs(x = expression("Branching rate"~lambda[b]),
       y = "FCL") +
  guides(color = "none")

## food chain length vs. size
g_size <- lineart(df_fcl_line,
                  df_fcl_mu,
                  x = "size",
                  x_axis = "rl") +
  labs(x = expression("River length"~italic(L)),
       y = "",
       color = "Disturbance type")

## layout ####

layout <- "
AA
AA
BC
"
g_sim_main <- g_heat05 + g_br + g_size +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")

ggsave(g_sim_main,
       filename = "data_fmt/fig_sim_main.pdf",
       width = 9,
       height = 10)

