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
    lab_mu0 = sprintf("mu^{(0)}==%.2f", mu0),
    lab_nu = sprintf("nu==%.2f", nu),
    lab_tl = str_to_sentence(tl) %>% 
      factor(levels = c("Prey", "Predator")),
    lab_delta0 = sprintf("delta==%.2f", delta0)
  )

## figure
(
  g2sp <- df_plot_2sp %>% 
    ggplot(aes(y = rl,
               x = lambda,
               fill = o)) +
    geom_raster(alpha = 1) +
    ggh4x::facet_nested(rows = vars(lab_mu0, lab_nu),
                        cols = vars(lab_tl, lab_delta0),
                        labeller = label_parsed,
                        nest_line = element_line(linetype = 3)) +
    # geom_vline(xintercept = c(0.4, 0.8),
    #            color = "white",
    #            linetype = "dashed") +
    scale_fill_viridis_c() +
    labs(x = expression("Branching rate"~lambda[b]),
         y = expression("Total river length"~italic(L)),
         fill = "Occupancy") +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text =  element_text(size = 7))
)

ggsave(g2sp, 
       filename = "tex/fig_theo_2sp_heat.pdf",
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
  mutate(lab_mu0 = sprintf("mu^{(0)}==%.2f", mu0),
         lab_r0 = sprintf("italic(r[0])==%.2f", r0),
         lab_mu_p = ifelse(mu_p > min(mu_p),
                           sprintf("mu^{(1)}==%.2f~(strong~prey~effect)", mu_p),
                           sprintf("mu^{(1)}==%.2f~(weak~prey~effect)", mu_p)),
         lab_mu_c = ifelse(mu_c == min(mu_c),
                           sprintf("mu^{(2)}==%.2f~(weak~predation)", mu_c),
                           sprintf("mu^{(2)}==%.2f~(strong~predation)", mu_c)))

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
      y = expression("River length"~italic(L)),
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
         height = 4.5,
         width = 6)
  
  return(g)
}
