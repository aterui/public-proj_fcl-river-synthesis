#' DESCRIPTION:
#' Create a table for parameter values

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")


# table -------------------------------------------------------------------

## read parameters
list_parms <- list.files("data_fmt",
                         full.names = TRUE,
                         pattern = "sim_fcl_m_heat|sim_fcl_si") %>% 
  lapply(FUN = function(x) {
    readRDS(x) %>% 
      dplyr::select(rsrc, 
                    g,
                    delta0,
                    h,
                    mu0,
                    mu_p,
                    mu_c,
                    rho,
                    z,
                    theta) %>% 
      rename(r0 = rsrc,
             g0 = g,
             psi = z) %>% 
      filter(rho != 0)
  })

## transform them into tables

df_parms <- lapply(list_parms, FUN = function(x) {
  
  uv <- sapply(x, function(x) paste(sprintf("%.2f", unique(x)), collapse = ", "))
  uv[which(colnames(x) == "g0")] <- paste(sprintf("%.0f", unique(x$g0)), collapse = ", ")
    
  parms <- names(uv)
  names(uv) <- NULL
  
  cout <- tibble(parms = parms,
                 value = uv)
  
  return(cout)
}) %>% 
  reduce(left_join, by = "parms") %>% 
  mutate(symbol = case_when(parms == "r0" ~ "$r_0$",
                            parms == "g0" ~ "$g_0$",
                            parms == "h" ~ "$h$",
                            parms == "delta0" ~ "$\\delta_0$",
                            parms == "mu0" ~ "$\\mu^{(0)}$",
                            parms == "mu_p" ~ "$\\mu^{(1)}$",
                            parms == "mu_c" ~ "$\\mu^{(2)}$",
                            parms == "rho" ~ "$\\rho$",
                            parms == "psi" ~ "$\\psi$",
                            parms == "theta" ~ "$\\theta$"),
         description = case_when(parms == "r0" ~ "Resource supply [-]",
                                 parms == "g0" ~ "Propagule size for producers [-]",
                                 parms == "h" ~ "Habitat density [per unit river distance]",
                                 parms == "delta0" ~ "Dispersal capability for producers [per unit river distance]",
                                 parms == "mu0" ~ "Disturbance rate [per unit time]",
                                 parms == "mu_p" ~ "Maximum prey-induced extinction rate [per unit time]",
                                 parms == "mu_c" ~ "Predator-induced extinction rate [per unit time]",
                                 parms == "rho" ~ "Synchrony probability [-]",
                                 parms == "psi" ~ "Scaling exponent for propagule and dispersal parameters [per unit trophic position]",
                                 parms == "theta" ~ "Degree of omnivory [per unit trophic position]")) %>% 
  rename("Value (analytical)" = value.x,
         "Value (numerical)" = value.y) %>% 
  relocate(symbol, description) %>% 
  dplyr::select(-parms) %>% 
  rename_with(.fn = str_to_sentence)


# export ------------------------------------------------------------------

options(xtable.comment = FALSE)

print(xtable(df_parms %>% dplyr::select(-`Value (numerical)`),
             caption = "Parameter descriptions and values (analytical).\\label{tab:parms}"),
      tabular.environment = "tabularx", # use \begin{tabularx}
      width = "\\textwidth", # scale table with \textwidth
      sanitize.text.function = function(x) x, # for math mode
      include.rownames = FALSE,
      caption.placement = "top",
      size = "\\small",
      file = "rmd/table.tex")

print(xtable(df_parms %>% dplyr::select(-`Value (analytical)`),
             caption = "Parameter values in numerical analysis.\\label{tab:parms-num}"),
      tabular.environment = "tabularx", # use \begin{tabularx}
      width = "\\textwidth", # scale table with \textwidth
      sanitize.text.function = function(x) x, # for math mode
      include.rownames = FALSE,
      caption.placement = "top",
      size = "\\small",
      file = "rmd/table_si.tex")
