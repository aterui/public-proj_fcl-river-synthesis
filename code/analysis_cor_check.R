
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

# set data ----------------------------------------------------------------

## read data
## - df_fcl_local, local level data
## - df_fcl_wsd, watershed level data
list_fcl <- readRDS("data_fmt/data_fcl_reg.rds")

df_x <- list_fcl[[1]] %>% 
  rename(local_hfp = hfp) %>% 
  left_join(list_fcl[[2]] %>% 
              dplyr::select(uid,
                            r_length,
                            lambda,
                            prec,
                            temp,
                            hfp)) %>% 
  mutate(resid_elev = resid(lm(log(local_elev) ~ prec, .))) %>% 
  dplyr::select(local_area,
                local_elev,
                forest_b1km,
                r_length,
                lambda,
                prec,
                temp,
                hfp) %>% 
  mutate(log_area = log(local_area),
         log_rl = log(r_length),
         log_lambda = log(lambda)) %>% 
  dplyr::select(-c(local_area,
                   r_length,
                   lambda)) %>% 
  relocate(log_area,
           local_elev,
           forest_b1km,
           log_rl,
           log_lambda)

plot(df_x)
