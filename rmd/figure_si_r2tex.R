library(tidyverse)

## parameter setup
df_grt <- readRDS(here::here("data_fmt/sim_fcl_si.rds")) %>%
  distinct(
    # fixed parameters
    g, 
    h, 
    z,
    # varied parameters
    delta0,
    nu, 
    theta
  ) %>% 
  arrange(
    desc(delta0), 
    desc(nu)
  )

## text setup
cap <- with(df_grt, {
  
  title0 <- paste0("Numerical prediction with ",
                   ifelse(
                     delta0 == max(delta0),
                     paste0("dispersal decay rate $\\delta_0 = ", delta0, "$ (for producers), "),
                     paste0("dispersal decay rate $\\delta_0 = ", delta0, "$ (for producers), ")
                   ),
                   ifelse(
                     nu == max(nu),
                     paste0("disturbance decay rate $\\nu = ", nu, "$, "),
                     paste0("disturbance decay rate $\\nu = ", nu, "$, ")
                   ),
                   ifelse(
                     theta == max(theta),
                     paste0("and strong omnivory $\\theta = ", theta, "$."),
                     paste0("and weak omnivory $\\theta = ", theta, "$. ")
                   )
  )
  
  c_text0 <- "Heatmaps of FCL are expressed as a function of ecosystem size (river length, $L$)
and complexity (branching rate, $\\lambda_b$), with rows and columns displaying
different combinations of resource supply ($r_0$), disturbance regime
($\\mu^{(0)}$), prey effect ($\\mu^{(1)}$), and predation effect ($\\mu^{(2)}$).
Each cell represents the average FCL of five food webs.
Other parameter values are: "
  
  cap0 <- paste0(
    title0,
    c_text0,
    "producer propagule $g_0 = ", g, "$, ",
    "habitat density $h = ", h, "$, ",
    "scaling exponent $\\psi = ", z, "$."
  )
  
  return(cap0)
})

## graphics
pdf <- with(df_grt,
            paste0("tex/fig_theo",
                   "_delta", delta0,
                   "_nu", nu,
                   "_theta", theta)) %>% 
  str_remove_all("\\.") %>% 
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
