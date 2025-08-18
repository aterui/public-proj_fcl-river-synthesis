library(tidyverse)

## parameter setup
df_grt <- readRDS(here::here("data_fmt/sim_fcl_si.rds")) %>%
  distinct(g, rho, theta, h, delta0, z)

## text setup
cap <- with(df_grt, {
  title0 <- paste0("Numerical prediction with ",
                   ifelse(g == max(g),
                          paste0("high propagule $g_0 = ", g, "$, "),
                          paste0("low propagule $g_0 = ", g, "$, ")),
                   ifelse(rho == max(rho),
                          paste0("high synchrony $\\rho = ", rho, "$, "),
                          paste0("low synchrony $\\rho = ", rho, "$, ")),
                   ifelse(theta == max(theta),
                          paste0("and strong omnivory $\\theta = ", theta, "$."),
                          paste0("and weak omnivory $\\theta = ", theta, "$. ")))
  
  c_text0 <- "Heatmaps of FCL are expressed as a function of ecosystem size (river length, $L$)
and complexity (branching rate, $\\lambda_b$), with rows and columns displaying
different combinations of resource supply ($r_0$), disturbance regime
($\\mu^{(0)}$), prey effect ($\\mu^{(1)}$), and predation effect ($\\mu^{(2)}$).
Each cell represents the average FCL of five food webs.
Additional parameter values are: "
  
  cap0 <- paste0(title0,
                 c_text0,
                 "habitat density $h=", h, "$, ",
                 "dispersal capability $\\delta_0=", delta0, "$, and ",
                 "scaling exponent $\\psi=", z, "$.")
  
  return(cap0)
})

## graphics
pdf <- with(df_grt,
            paste0("tex/",
                   "fig_theo_rho", rho,
                   "_g", g,
                   "_theta", theta) %>%
              str_remove_all("\\.")) %>%
  paste0(".pdf")

## assemble chunks for numerical solutions
fig_num_chunk <- knit_chunk <- NULL

for(i in seq_len(nrow(df_grt))) {
  
  knit_chunk <- paste0(
    "\\begin{figure}\n",
    "\\centering\n",
    "\\includegraphics[keepaspectratio]{", pdf[i], "}\n",
    "\\caption{", cap[i], "}\n",
    "\\label{fig:fig-num", i, "}\n",
    "\\end{figure}\n",
    "\\newpage\n"
  )
  
  fig_num_chunk <- c(fig_num_chunk, knit_chunk)
  
}

## export
writeLines(paste(fig_num_chunk, collapse = "\n"),
           con = here::here("tex/figure_si.tex"))
