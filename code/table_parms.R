#' DESCRIPTION:
#' Create tables for parameter values

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
      dplyr::select(r0, 
                    b,
                    g,
                    delta0,
                    h,
                    mu0,
                    mu_p,
                    mu_c,
                    nu,
                    z,
                    theta)
  })

## transform them into tables
df_parms <- lapply(list_parms, FUN = function(x) {
  
  uv <- sapply(x, function(x) {
    if (any(abs(x) < 1e-3) && any(x != 0)) {
      u <- unique(x)
      e <- floor(log10(abs(u)))
      m <- u / 10^e
      paste(
        sprintf("$%.2f \\times 10^{%d}$", m, e),
        collapse = ","
      )
    } else {
      paste(
        sprintf("%.2f", unique(x)),
        collapse = ", "
      )
    }
  })

  uv[which(colnames(x) == "g")] <- paste(sprintf("%.0f", unique(x$g)), collapse = ", ")
    
  parms <- names(uv)
  names(uv) <- NULL
  
  cout <- tibble(parms = parms,
                 value = uv)
  
  return(cout)
}) %>% 
  reduce(left_join, by = "parms") %>% 
  mutate(symbol = case_when(parms == "r0" ~ "$r_0$",
                            parms == "b" ~ "$b$",
                            parms == "g" ~ "$g_0$",
                            parms == "h" ~ "$h$",
                            parms == "delta0" ~ "$\\delta_0$",
                            parms == "mu0" ~ "$\\mu^{(0)}$",
                            parms == "mu_p" ~ "$\\mu^{(1)}$",
                            parms == "mu_c" ~ "$\\mu^{(2)}$",
                            parms == "nu" ~ "$\\nu$",
                            parms == "z" ~ "$\\psi$",
                            parms == "theta" ~ "$\\theta$"),
         description = case_when(parms == "r0" ~ "Resource supply [-]",
                                 parms == "b" ~ "Downstream accumulation rate of resources [per unit river length]",
                                 parms == "g" ~ "Propagule production for producers [per habitat]",
                                 parms == "h" ~ "Habitat density [per unit river length]",
                                 parms == "delta0" ~ "Dispersal decay rate for producers [unit river length]",
                                 parms == "mu0" ~ "Disturbance-induced extinction rate [per unit time]",
                                 parms == "mu_p" ~ "Maximum prey-induced extinction rate [per unit time]",
                                 parms == "mu_c" ~ "Predator-induced extinction rate [per unit time]",
                                 parms == "nu" ~ "Distance decay rate of disturbance synchrony [per unit river length]",
                                 parms == "z" ~ "Scaling exponent for propagule and dispersal parameters [per unit ln trophic position]",
                                 parms == "theta" ~ "Degree of omnivory [per unit trophic position]")) %>% 
  relocate(symbol, description) %>% 
  dplyr::select(-parms) %>% 
  rename_with(.fn = str_to_sentence)


## export
print(xtable(df_parms %>%
               dplyr::select(-Value.y) %>% 
               rename(Value = Value.x),
             caption = "Parameter descriptions and values used for analytical predictions.",
             label = "tab:parms"),
      tabular.environment = "tabularx", # use \begin{tabularx}
      width = "\\textwidth", # scale table with \textwidth
      sanitize.text.function = function(x) x, # for math mode
      include.rownames = FALSE,
      caption.placement = "top",
      size = "\\small",
      file = "tex/table_main.tex")

print(xtable(df_parms %>%
               dplyr::select(-Value.x) %>% 
               rename(Value = Value.y),
             caption = "Parameter descriptions and values used for numerical predictions.",
             label = "tab:parms-num"),
      tabular.environment = "tabularx", # use \begin{tabularx}
      width = "\\textwidth", # scale table with \textwidth
      sanitize.text.function = function(x) x, # for math mode
      include.rownames = FALSE,
      caption.placement = "top",
      size = "\\small",
      file = "tex/table_si_theo.tex")


# Bayesian model estimate -------------------------------------------------

h0 <- readRDS("data_fmt/output_model_h0.rds")
h1 <- readRDS("data_fmt/output_model_h1.rds")

list_est <- 
  list(
    h0[[1]], 
    h1[[1]]
  )

list_table <- lapply(list_est, function(data) {
  
  data %>% 
    filter(!str_detect(parms, "b0")) %>% 
    mutate(symbol = case_when(varname == "local_elev" ~ "$\\alpha_1$",
                              varname == "(Intercept)" ~ "$\\beta_0$",
                              varname == "log_rl" ~ "$\\beta_1$",
                              varname == "log_lambda" ~ "$\\beta_2$",
                              varname == "temp" ~ "$\\beta_3$",
                              varname == "prec" ~ "$\\beta_4$",
                              varname == "hfp" ~ "$\\beta_5$",
                              parms == "nu" ~ "$\\nu$",
                              parms == "sigma[1]" ~ "$\\sigma$",
                              parms == "sigma[2]" ~ "$\\sigma_{\\varepsilon}$",
                              parms == "sigma[3]" ~ "$\\sigma_{\\eta}$",
                              parms == "sigma_b[1]" ~ "$\\sigma_{\beta1}$",
                              parms == "sigma_b[2]" ~ "$\\sigma_{\beta2}$",
                              parms == "sigma_b[3]" ~ "$\\sigma_{\beta3}$",
                              parms == "z[1]" ~ "$\\xi_{1}$",
                              parms == "z[2]" ~ "$\\xi_{2}$"),
           description = case_when(varname == "local_elev" ~ "Elevation",
                                   varname == "(Intercept)" ~ "Intercept",
                                   varname == "log_rl" ~ "ln Total river length",
                                   varname == "log_lambda" ~ "ln Branching rate",
                                   varname == "temp" ~ "Air temperature",
                                   varname == "prec" ~ "Precipitation",
                                   varname == "hfp" ~ "Human footprint",
                                   parms == "nu" ~ "Degrees of freedom",
                                   parms == "sigma[1]" ~ "Site-level residual SD",
                                   parms == "sigma[2]" ~ "Watershed-level SD (random effect)",
                                   parms == "sigma[3]" ~ "Region-level SD (random effect)",
                                   parms == "sigma_b[1]" ~ "Region-level SD (intercept)",
                                   parms == "sigma_b[2]" ~ "Region-level SD (ln Total river length)",
                                   parms == "sigma_b[3]" ~ "Region-level SD (ln Branching rate)",
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
    filter(!is.na(symbol)) %>% 
    dplyr::select(Symbol = symbol,
                  Description = description,
                  Estimate = estimate,
                  `95\\% CI` = ci,
                  `$\\Pr(< 0)$` = pr_neg,
                  `$\\Pr(> 0)$` = pr_pos)
  
})

## export
m <- c("random intercept", "random intercept and slope")

lapply(1:2, 
       function(i) {
         
         print(
           xtable(list_table[[i]],
                  caption = 
                    paste0(
                      "Parameter estimates of the hierarchical Bayesian model",
                      " (", m, ") ",
                      "with corresponding 95\\% credible intervals (CIs) and posterior probabilities ($\\Pr(\\cdot)$), representing the uncertainty around each parameter estimate." 
                    ),
                  label = "tab:parms-est"),
           tabular.environment = "tabularx", # use \begin{tabularx}
           width = "\\textwidth", # scale table with \textwidth
           sanitize.text.function = function(x) x, # for math mode
           include.rownames = FALSE,
           caption.placement = "top",
           size = "\\small",
           file = paste0("tex/table_si_emp", i, ".tex")
         ) # print
         
       })

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
             `Code' refers to the unique identifier assigned to each study for use in the analysis.",
             align = "p{0}p{0.1\\textwidth}p{0.9\\textwidth}",
             label = "tab:meta-list"),
      tabular.environment = "longtable", # use \begin{longtable}
      sanitize.text.function = function(x) x, # for math mode
      include.rownames = FALSE,
      caption.placement = "top",
      size = "\\small",
      file = "tex/table_si_list.tex",
      floating = FALSE)
