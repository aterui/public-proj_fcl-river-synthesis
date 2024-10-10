#' DESCRIPTION:
#' Get sample sizes

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

## - sourced `format_data4jags.R` in the following script
source("code/format_emp_est2figure.R")


# sample size -------------------------------------------------------------

## n sites
n_distinct(df_fcl_local$sid)

## n watersheds
n_distinct(df_fcl_local$g)

## n studies
df_fcl_local %>% 
  pull(sid) %>% 
  str_remove_all("_s\\d{3,3}") %>% 
  n_distinct()
