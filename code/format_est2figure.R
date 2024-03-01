#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")


# read data ---------------------------------------------------------------

## data for predictors: df_g
source("code/format_data4jags.R")

df_g <- df_g %>% 
  mutate(log_area = log(area),
         log_pb = log(p_branch))

## read data
df_est <- readRDS("data_fmt/output_model_summary.rds")
mcmc <- readRDS("data_fmt/output_model_mcmc.rds")


# region specific prediction ----------------------------------------------

## parameters
## - region-specific intercept
v_r <- df_est %>% 
  filter(parms_gr == "r") %>% 
  pull(median)

## - fixed effects
m_b <- df_est %>% 
  filter(parms_gr == "b", parms_num > 0) %>% 
  pull(median) %>% 
  matrix(nrow = length(.), ncol = length(v_r))

B <- rbind(v_r, m_b)

## prepare predictors
## - predictor columns
cnm <- c("log_area", "log_pb", "scl_prec", "scl_temp", "scl_hfp")

df_yh0 <- foreach(j = seq_len(length(cnm)),
                  .combine = bind_rows) %do% {
                    
                    ## set mean values for each predictor
                    df_x <- df_g %>% 
                      dplyr::select(all_of(c("h", cnm))) %>% 
                      group_by(h) %>% 
                      reframe(across(cnm,
                                     .fns = function(x) rep(mean(x), 100)))
                    
                    ## get min-max range for a focus predictor
                    cid <- which(colnames(df_g) %in% c("h", cnm[j]))
                    
                    x <- df_g %>% 
                      select(all_of(cid)) %>% 
                      group_by(h) %>% 
                      reframe(x = seq(min(.[, 2]),
                                      max(.[, 2]),
                                      length = 100)
                      ) %>% 
                      pull(x)
                    
                    ## replace the mean with the range of the focus variable
                    df_x[, which(colnames(df_x) == cnm[j])] <- x
                    
                    ## predicted values of y = fcl
                    df_x <- df_x %>% 
                      mutate(i = row_number(),
                             focus = cnm[j]) %>% 
                      relocate(focus)
                    
                    X0 <- df_x %>% 
                      dplyr::select(cnm)
                    
                    X <- model.matrix(~ ., X0)
                    
                    Y <- X %*% B
                    
                    ## Y's column corresponds to each region
                    ## Y[i, h] switch columns based on `h` column
                    cout <- df_x %>% 
                      rowwise() %>% 
                      mutate(log_y = Y[i, h],
                             y = exp(log_y))
                    
                    return(cout)
                  }

## get unscaled values
df_yh <- df_yh0 %>% 
  mutate(area = exp(log_area),
         pb = exp(log_pb),
         prec = scl_prec * sd(df_g$mean.prec) + mean(df_g$mean.prec),
         temp = scl_temp * sd(df_g$mean.temp) + mean(df_g$mean.temp),
         hfp = scl_hfp * sd(df_g$hfp) + mean(df_g$hfp))


# overall prediction ------------------------------------------------------

x <- with(df_g,
          seq(min(log(p_branch)),
              max(log(p_branch)),
              length = 100))

X <- with(df_g,
          cbind(1,
                mean(log(area)),
                x))

mcmc_b <- mcmc[, str_detect(colnames(mcmc), "b0|b\\[.{1,}\\]")]

df_y <- foreach(j = seq_len(length(cnm)),
                .combine = bind_rows) %do% {
                  ## set mean values for each predictor
                  df_x <- df_g %>% 
                    dplyr::select(all_of(c("h", cnm))) %>% 
                    reframe(across(cnm,
                                   .fns = function(x) rep(mean(x), 100)))
                  
                  ## get min-max range for a focus predictor
                  cid <- which(colnames(df_g) %in% c("h", cnm[j]))
                  
                  x <- df_g %>% 
                    select(all_of(cid)) %>% 
                    reframe(x = seq(min(.[, 2]),
                                    max(.[, 2]),
                                    length = 100)
                    ) %>% 
                    pull(x)
                  
                  ## replace the mean with the range of the focus variable
                  df_x[, which(colnames(df_x) == cnm[j])] <- x
                  
                  ## predicted values of y = fcl
                  df_x <- df_x %>% 
                    mutate(focus = cnm[j]) %>% 
                    relocate(focus)
                  
                  X0 <- df_x %>% 
                    dplyr::select(cnm)
                  
                  X <- model.matrix(~ ., X0)
                  
                  mcmc_y <- X %*% t(mcmc_b)
                  
                  ## get quantiles for predicted values
                  cout <- mcmc_y %>% 
                    apply(MARGIN = 1, quantile, c(0.5, 0.025, 0.975)) %>% 
                    t() %>% 
                    as_tibble() %>% 
                    setNames(c("y", "y_low", "y_high")) %>% 
                    bind_cols(df_x) %>% 
                    mutate(prec = scl_prec * sd(df_g$mean.prec) + mean(df_g$mean.prec),
                           temp = scl_temp * sd(df_g$mean.temp) + mean(df_g$mean.temp),
                           hfp = scl_hfp * sd(df_g$hfp) + mean(df_g$hfp)) %>% 
                    relocate(focus)
                }
