#' DESCRIOTION:
#' Produce scheme figure and main theoretical predictions
#' Main theoretical predictions include heatmap and line arts

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

## set plot theme
source("code/set_theme.R")
ggplot2::theme_set(default_theme)

# figure 1 ----------------------------------------------------------------
## schematic diagram

## draw scheme
img <- image_read_pdf("tex/figsub_scheme.pdf")

g_subsch <- ggdraw() +
  draw_image(img) +
  labs(tag = "A")

## colonization-extinction
df_ce <- tibble(rl = seq(100, 1000, length = 100)) %>% 
  expand_grid(
    lambda = seq(0.2, 1, by = 0.2),
    delta0 = 0.1,
    h = 1,
    nu = 0.1
  ) %>% 
  rowwise() %>% 
  mutate(
    diam = diameter(lambda = lambda, size = rl, exact = TRUE),
    d = pdist(lambda = lambda, size = rl, exact = TRUE),
    u = u_length(lambda = lambda, size = rl, exact = TRUE),
    phi = laplace_rayleigh(delta = delta0, mu = d),
    rho = laplace_rt(nu = nu, mu = diam, exact = TRUE)
  ) %>% 
  ungroup() %>% 
  mutate(
    pgle = h * rl * phi,
    extn = 1 + rho * u
  )

g_pgle <- df_ce %>% 
  ggplot(
    aes(
      x = rl,
      y = pgle,
      color = lambda,
      group = lambda
    )
  ) +
  geom_line() +
  labs(
    y = "Effective propagule",
  ) +
  MetBrewer::scale_color_met_c(
    "Hiroshige",
    direction = -1
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  guides(color = "none")

g_extn <- df_ce %>% 
  ggplot(
    aes(
      x = rl,
      y = extn,
      color = lambda,
      group = lambda
    )
  ) +
  geom_line() +
  labs(
    x = expression("Total river length"~italic(L)),
    y = "Disturbance rate",
    color = expression(lambda[b])
  ) +
  MetBrewer::scale_color_met_c(
    "Hiroshige",
    direction = -1
  )

g_pe <- (g_pgle + labs(tag = "B")) / g_extn

## layout
layout <- "
AAAB
AAAB
"

(g_scheme <- g_subsch + g_pe + 
   plot_layout(design = layout, 
               guides = "collect"))

## export
ggsave(g_scheme,
       filename = "tex/fig_theo_scheme.pdf",
       width = 8,
       height = 4)

# figure 2 ----------------------------------------------------------------
## occupancy example

df_2sp <- readRDS("data_fmt/sim_fcl_2sp_line.rds") %>% 
  pivot_longer(
    cols = c("lambda", "rl"), 
    values_to = "x",
    names_to = "xnm"
  ) %>% 
  filter(focus == xnm) %>% 
  mutate(
    lab_tl = str_to_sentence(tl) %>% 
      fct_rev(),
    lab_delta = ifelse(
      delta0 == min(delta0),
      sprintf("Long~dispersal~(delta==%.2f)", delta0),
      sprintf("Short~dispersal~(delta==%.2f)", delta0)
    ),
    lab_focus = case_when(
      focus == "lambda" ~ "Branching~rate~lambda[b]",
      focus == "rl" ~ "Total~river~length~italic(L)",
    ) %>% 
      fct_rev()
  )
  
(
  g_pp <- df_2sp %>%
    ggplot(
      aes(
        x = x,
        y = o,
        linetype = lab_tl
      )
    ) +
    geom_line() +
    ggh4x::facet_grid2(
      rows = vars(lab_delta),
      cols = vars(lab_focus),
      scales = "free",
      labeller = label_parsed,
      switch = "x",
      axes = "x"
    ) +
    labs(y = "Occupancy") +
    scale_y_continuous(lim = c(0, 1)) +
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
      strip.placement = "outside",
      legend.text = element_text(size = 10)
    )
)

## export
ggsave(g_pp,
       filename = "tex/fig_theo_2sp.pdf",
       width = 6.5,
       height = 4.5)

# figure 2 ----------------------------------------------------------------

## heatmap ####

## read data for heatmap
## take average for food web replicates
df_heat <- readRDS("data_fmt/sim_fcl_m_heat.rds") %>% 
  group_by(rl, lambda, h, delta0, r0, g, mu0, mu_p, mu_c, rho0, theta) %>%
  summarize(
    fcl = mean(fcl), 
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(theta == min(theta)) %>% 
  mutate(lab_mu0 = ifelse(mu0 > min(mu0),
                          sprintf("mu^{(0)}==%.2f~(freq.)", mu0),
                          sprintf("mu^{(0)}==%.2f~(infreq.)", mu0)),
         lab_r0 = ifelse(r0 > min(r0),
                         sprintf("italic(r[0])==%.2f~(high)", r0),
                         sprintf("italic(r[0])==%.2f~(low)", r0)),
         lab_rho = ifelse(rho0 > min(rho0),
                          sprintf("rho[0]==%.2f~(disturbance~cascade)", rho0),
                          sprintf("rho[0]==%.2f~(no~disturbance~cascade)", rho0)),
         lab_d = "Disturbance~rate",
         lab_r = "Resource~supply",
         lambda_mid = ifelse(r0 == max(r0) & mu0 == max(mu0),
                             mean(range(lambda)),
                             NA),
         rl_mid = ifelse(r0 == max(r0) & mu0 == max(mu0),
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
                        cols = vars(lab_r, lab_r0),
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
         y = expression("Total river length"~italic(L)),
         fill = "FCL") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
}

(g_heat05 <- heatmap(data = df_heat %>% filter(rho0 == 0.5)))
(g_heat0 <- heatmap(data = df_heat %>% filter(rho0 == 0)))

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
                                  `0.5` = "Cascade")) +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
}

## food chain length vs. branching
g_br <- lineart(df_fcl_line,
                df_fcl_mu,
                x = "branch",
                x_axis = "lambda") +
  labs(x = expression("Branching rate"~lambda[b]),
       y = "Food chain length") +
  guides(color = "none")

## food chain length vs. size
g_size <- lineart(df_fcl_line,
                  df_fcl_mu,
                  x = "size",
                  x_axis = "rl") +
  labs(x = expression("Total river length"~italic(L)),
       y = "",
       color = "Disturbance")

## layout ####

layout <- "
AA
AA
BC
"
g_sim_main <- g_heat05 + g_br + g_size +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")


# export ------------------------------------------------------------------

ggsave(g_sim_main,
       filename = "tex/fig_theo_main.pdf",
       width = 9,
       height = 10)

ggsave(g_heat0,
       filename = "tex/fig_theo_rho0.pdf",
       width = 6,
       height = 5)
