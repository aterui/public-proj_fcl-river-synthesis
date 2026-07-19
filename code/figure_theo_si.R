#' DESCRIPTION
#' Produce predictions with numerical analysis
#' Heatmaps for a number of parameter combinations

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

# patch occupancy ---------------------------------------------------------

df_2sp <- readRDS("data_fmt/sim_fcl_2sp_heat.rds")

## take average for food web replicates
df_plot_2sp <- df_2sp %>% 
  mutate(
    lab_tl = str_to_sentence(tl) %>% 
      factor(levels = c("Prey", "Predator")),
    lab_mu0 = 
      ifelse(mu0 == max(mu0),
             sprintf("High~disturbance~(mu^{(0)}==%.2f)", mu0),
             sprintf("Low~disturbance~(mu^{(0)}==%.2f)", mu0)
      ) %>% 
      fct_rev(),
    lab_nu = 
      ifelse(
        nu == max(nu),
        sprintf("atop(Fast~decay,nu==%.2f)", nu),
        sprintf("atop(Slow~decay,nu==%.2f)", nu)
      ) %>% 
      fct_rev(),
    lab_delta0 = 
      ifelse(
        delta0 == max(delta0),
        sprintf("atop(Short~range~dispersal,delta==%.2f)", delta0),
        sprintf("atop(Long~range~dispersal,delta==%.2f)", delta0)
      )
  )

## figure
(
  g2sp <- df_plot_2sp %>% 
    ggplot(
      aes(y = rl,
          x = lambda,
          fill = o)
    ) +
    geom_raster(alpha = 1) +
    ggh4x::facet_nested(
      rows = vars(lab_mu0, lab_nu),
      cols = vars(lab_tl, lab_delta0),
      labeller = label_parsed,
      nest_line = element_line(linetype = 3)
    ) +
    scale_fill_viridis_c() +
    labs(
      x = expression("Branching rate"~lambda[b]),
      y = expression("Total river length"~italic(L)),
      fill = "Occupancy"
    ) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      axis.text =  element_text(size = 7)
    )
)

ggsave(g2sp, 
       filename = "tex/fig_theo_2sp_heat.pdf",
       height = 6,
       width = 8)


# spatially explicit model ------------------------------------------------

## set theme
source("code/set_theme.R")
set_theme(default_theme)

## read data
df_sbn <- readRDS("data_fmt/sim_fcl_2sp_nspom.rds") %>% 
  mutate(
    lab_tl = str_to_sentence(tl) %>% 
      factor(levels = c("Prey", "Predator")),
    lab_mu0 = 
      ifelse(mu0 == max(mu0),
             sprintf("High~disturbance~(mu^{(0)}==%.2f)", mu0),
             sprintf("Low~disturbance~(mu^{(0)}==%.2f)", mu0)
      ) %>% 
      fct_rev(),
    lab_nu = 
      ifelse(
        nu == max(nu),
        sprintf("atop(Fast~decay,nu==%.2f)", nu),
        sprintf("atop(Slow~decay,nu==%.2f)", nu)
      ) %>% 
      fct_rev(),
    lab_theta = 
      ifelse(
        theta == max(theta),
        sprintf("High~dispersal~(theta[d]==%.2f)", theta),
        sprintf("Low~dispersal~(theta[d]==%.2f)", theta)
      ) %>% 
      fct_rev(),
    lab_xi = 
      ifelse(
        kappa == max(kappa),
        sprintf("atop(Downstream~biased,xi==%.2f)", kappa),
        sprintf("atop(Symmetric, xi==%.2f)", kappa)
      ) %>% 
      fct_rev()
  )

## ecosystem size
(
  g_se_n <- df_sbn %>% 
    ggplot(
      aes(x = n,
          y = p,
          color = lab_tl,
          fill = lab_tl)
    ) +
    geom_point(
      size = 0.5,
      alpha = 0.1
    ) +
    geom_smooth(linewidth = 0.5)  +
    ggh4x::facet_nested(
      rows = vars(lab_mu0, lab_nu),
      cols = vars(lab_theta, lab_xi),
      labeller = label_parsed,
      nest_line = element_line(linetype = 3)
    ) +
    theme(
      strip.background = element_blank(),
      axis.text =  element_text(size = 7)
    ) +
    labs(
      x = "Number of habitats N",
      y = "Occupancy",
      color = "Trophic level",
      fill = "Trophic level"
    )
)

## ecosystem complexity
(
  g_se_br <- df_sbn %>% 
    ggplot(
      aes(x = lambda,
          y = p,
          color = lab_tl,
          fill = lab_tl)
    ) +
    geom_point(
      size = 0.5,
      alpha = 0.1
    ) +
    geom_smooth(linewidth = 0.5)  +
    ggh4x::facet_nested(
      rows = vars(lab_mu0, lab_nu),
      cols = vars(lab_theta, lab_xi),
      labeller = label_parsed,
      nest_line = element_line(linetype = 3)
    ) +
    theme(
      strip.background = element_blank(),
      axis.text =  element_text(size = 7)
    ) +
    labs(
      x = expression("Branching rate"~lambda[b]),
      y = "Occupancy",
      color = "Trophic level",
      fill = "Trophic level"
    )
)

ggsave(g_se_n, 
       filename = "tex/fig_theo_2sp_sem_n.pdf",
       height = 6,
       width = 8)

ggsave(g_se_br, 
       filename = "tex/fig_theo_2sp_sem_br.pdf",
       height = 6,
       width = 8)


# food chain length -------------------------------------------------------

## read data
df_fcl_sim <- readRDS("data_fmt/sim_fcl_si.rds")

## take average for food web replicates
df_plot <- df_fcl_sim %>% 
  group_by(
    rl, 
    lambda, 
    r0, 
    mu0, 
    mu_p, 
    mu_c, 
    delta0, 
    nu,
    rho0, 
    theta
  ) %>% 
  summarize(
    fcl = mean(fcl),
    .groups = "drop"
  ) %>% 
  mutate(
    lab_mu0 = sprintf("mu^{(0)}==%.2f", mu0),
    lab_r0 = sprintf("italic(r[0])==%.2f", r0),
    lab_mu_p = ifelse(
      mu_p == max(mu_p),
      sprintf("Strong~prey~effect~(mu^{(1)}==%.2f)", mu_p),
      sprintf("Weak~prey~effect~(mu^{(1)}==%.2f)", mu_p)
    ) %>% 
      fct_rev(),
    lab_mu_c = ifelse(
      mu_c == max(mu_c),
      sprintf("Strong~predation~(mu^{(2)}==%.2f)", mu_c),
      sprintf("Weak~predation~(mu^{(2)}==%.2f)", mu_c)
    ) %>% 
      fct_rev()
  )

df_case <- df_plot %>% 
  distinct(delta0, nu, theta) %>% 
  arrange(desc(delta0), desc(nu))

## figure
list_g <- foreach(i = seq_len(nrow(df_case))) %do% {
  
  g <- df_plot %>% 
    filter(
      nu == df_case$nu[i],
      delta0 == df_case$delta0[i],
      theta == df_case$theta[i]
    ) %>% 
    ggplot(
      aes(
        y = rl,
        x = lambda,
        fill = fcl
      )
    ) +
    geom_raster(alpha = 1) +
    ggh4x::facet_nested(
      rows = vars(lab_mu_p, lab_mu0),
      cols = vars(lab_mu_c, lab_r0),
      labeller = label_parsed,
      nest_line = element_line(linetype = 3)
    ) +
    MetBrewer::scale_fill_met_c(
      "Hiroshige",
      direction = -1
    ) +
    labs(
      x = expression("Branching rate"~lambda[b]),
      y = expression("Total river length"~italic(L)),
      fill = "FCL"
    ) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      axis.text =  element_text(size = 7)
    )
  
  filename <- with(df_case[i, ],
                   paste0("tex/fig_theo",
                          "_delta", delta0,
                          "_nu", nu,
                          "_theta", theta)) %>% 
    str_remove_all("\\.") %>% 
    paste0(".pdf")
  
  ggsave(g, filename = filename,
         height = 6,
         width = 8)
  
  return(g)
}
