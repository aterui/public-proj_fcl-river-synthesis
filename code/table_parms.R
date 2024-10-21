#' DESCRIPTION:
#' Create tables for parameter values

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
options(xtable.comment = FALSE)


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


## export
print(xtable(df_parms %>% dplyr::select(-`Value (numerical)`),
             caption = "Parameter descriptions and values (analytical).\\label{tab:parms}"),
      tabular.environment = "tabularx", # use \begin{tabularx}
      width = "\\textwidth", # scale table with \textwidth
      sanitize.text.function = function(x) x, # for math mode
      include.rownames = FALSE,
      caption.placement = "top",
      size = "\\small",
      file = "rmd/table_main.tex")

print(xtable(df_parms %>% dplyr::select(-`Value (analytical)`),
             caption = "Parameter values in numerical analysis.\\label{tab:parms-num}"),
      tabular.environment = "tabularx", # use \begin{tabularx}
      width = "\\textwidth", # scale table with \textwidth
      sanitize.text.function = function(x) x, # for math mode
      include.rownames = FALSE,
      caption.placement = "top",
      size = "\\small",
      file = "rmd/table_si_theo.tex")


# Bayesian model estimate -------------------------------------------------

list_est <- readRDS("data_fmt/output_model_h0.rds")

df_sum <- list_est[[1]] %>% 
  filter(!str_detect(parms, "b0")) %>% 
  mutate(varname = case_when(parms == "a[1]" ~ "$\\alpha_1$",
                             parms == "b[1]" ~ "$\\beta_0$",
                             parms == "b[2]" ~ "$\\beta_1$",
                             parms == "b[3]" ~ "$\\beta_2$",
                             parms == "b[4]" ~ "$\\beta_3$",
                             parms == "b[5]" ~ "$\\beta_4$",
                             parms == "b[6]" ~ "$\\beta_5$",
                             parms == "nu" ~ "$\\nu$",
                             parms == "sigma[1]" ~ "$\\sigma$",
                             parms == "sigma[2]" ~ "$\\sigma_{\\varepsilon}$",
                             parms == "sigma[3]" ~ "$\\sigma_{\\eta}$",
                             parms == "z[1]" ~ "$\\xi_{1}$",
                             parms == "z[2]" ~ "$\\xi_{2}$"),
         description = case_when(parms == "a[1]" ~ "Elevation",
                                 parms == "b[1]" ~ "Intercept",
                                 parms == "b[2]" ~ "log River length",
                                 parms == "b[3]" ~ "log Branching rate",
                                 parms == "b[4]" ~ "Air temperature",
                                 parms == "b[5]" ~ "Precipitation",
                                 parms == "b[6]" ~ "Human footprint",
                                 parms == "nu" ~ "Degrees of freedom",
                                 parms == "sigma[1]" ~ "Site-level standard deviation",
                                 parms == "sigma[2]" ~ "Watershed-level standard deviation",
                                 parms == "sigma[3]" ~ "Region-level standard deviation",
                                 parms == "z[1]" ~ "Scaling exponent for the number of sites",
                                 parms == "z[2]" ~ "Scaling parameter for spatial sampling randomness"),
         estimate = sprintf("$%.2f$", median),
         ci = paste0("[",
                     sprintf("$%.2f$", low),
                     ", ",
                     sprintf("$%.2f$", high),
                     "]"),
         pr_neg = sprintf("$%.2f$", pr_neg),
         pr_pos = sprintf("$%.2f$", pr_pos)) %>% 
  dplyr::select(Symbol = varname,
                Description = description,
                Estimate = estimate,
                `95\\% CI` = ci,
                `$\\Pr(< 0)$` = pr_neg,
                `$\\Pr(> 0)$` = pr_pos)

## export
print(xtable(df_sum,
             caption = "Parameter estimates of the hierarchical Bayesian model 
             with corresponding 95\\% credible intervals (CIs) and 
             posterior probabilities ($\\Pr(\\cdot)$), 
             representing the uncertainty around each parameter estimate.
             \\label{tab:parms-est}"),
      tabular.environment = "tabularx", # use \begin{tabularx}
      width = "\\textwidth", # scale table with \textwidth
      sanitize.text.function = function(x) x, # for math mode
      include.rownames = FALSE,
      caption.placement = "top",
      size = "\\small",
      file = "rmd/table_si_emp.tex")


# table for meta-analysis sources -----------------------------------------

## download source file; update as needed
# googledrive::drive_download("data_lit_list_v_0_1_0",
#                             type = "csv",
#                             path = "data_raw/data_lit_list.csv",
#                             overwrite = TRUE)

## literature list and study_id for those used in our analysis
df_lit <- read_csv("data_raw/data_lit_list.csv")
u_study <- readRDS("data_fmt/data_fcl_reg.rds") %>% 
  .[[1]] %>% 
  pull(sid) %>% 
  str_remove_all("_s\\d{1,}") %>% 
  unique()

## format
df_src <- df_lit %>% 
  filter(study_id %in% u_study) %>% 
  mutate(journal = paste0("\\textit{", journal, "}"),
         pub0 = case_when(more_than_two == "N" & is.na(second_author) ~ 
                            paste0(first_author, " ",
                                   year, " ",
                                   title, ". ",
                                   journal),
                          more_than_two == "N" & !is.na(second_author) ~
                            paste0(first_author, " and ", second_author, " ",
                                   year, " ",
                                   title, ". ",
                                   journal),
                          more_than_two == "Y" ~ 
                            paste0(first_author, " \\textit{et al}. ",
                                   year, " ",
                                   title, ". ",
                                   journal)),
         Publication = ifelse(is.na(page_st)|is.na(page_end),
                              paste0(pub0, " ", volume, ": ", art_no),
                              paste0(pub0, " ", volume, ": ", page_st, "-", page_end)),
         Code = str_remove_all(study_id, "_.{1,}")
  ) %>% 
  select(Code,
         Publication)

## export
print(xtable(df_src,
             caption = "List of publications included in the meta-analysis.
             'Code' refers to the unique identifier assigned to each study for use in the analysis.
             \\label{tab:meta-list}",
             align = "p{0}p{0.1\\textwidth}p{0.9\\textwidth}"),
      tabular.environment = "longtable", # use \begin{longtable}
      sanitize.text.function = function(x) x, # for math mode
      include.rownames = FALSE,
      caption.placement = "top",
      size = "\\small",
      file = "rmd/table_si_list.tex",
      floating = FALSE)
