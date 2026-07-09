#' DESCRIPTION:
#' Format model estimates for visualization

rm(list = ls())
source("code/set_library.R")
source("code/format_emp_data4nimble.R")

# read data ---------------------------------------------------------------

df_fcl_wsd <- df_fcl_wsd %>%
  mutate(
    log_rl = log(r_length),
    log_lambda = log(lambda)
  ) %>%
  mutate(
    across(
      c(log_rl, log_lambda, prec, temp, hfp),
      ~ as.numeric(scale(.x)),
      .names = "scl_{.col}"
    )
  )

## read data
list_est <- list.files("data_fmt",
                       pattern = "output_.{1,}_h0|output_.{1,}_h1",
                       full.names = TRUE) %>% 
  lapply(FUN = function(x) {
    readRDS(x)
  })

## select the model with the lowest WAIC
id_best <- sapply(list_est, FUN = function(data) data$waic$WAIC) %>% 
  which.min()

## get estimates
df_est <- list_est[[id_best]]$est %>% 
  mutate(parms_gr = str_remove_all(parms, "\\[.\\]"),
         parms_num = str_extract(parms, "\\[\\d{1,}\\]") %>% 
           str_remove_all("\\[|\\]") %>% 
           as.numeric()) %>% 
  relocate(parms, parms_gr, parms_num)

# region specific prediction ----------------------------------------------

## parameters
## - region-specific intercept
v_r <- df_est %>% 
  filter(parms_gr == "b0") %>% 
  pull(median)

## - fixed effects
m_b <- df_est %>% 
  filter(parms_gr == "b", parms_num > 1) %>% 
  pull(median) %>% 
  matrix(nrow = length(.), ncol = length(v_r))

B <- unname(rbind(v_r, m_b))

## prepare predictors
## - predictor columns
## - NOTE: sensitive the order of the columns
cnm <- c("scl_log_rl", "scl_log_lambda", "scl_temp", "scl_prec", "scl_hfp")

## set group-specific means for each predictor
df_x0g <- df_fcl_wsd %>% 
  group_by(ecor) %>% 
  reframe(
    across(
      all_of(cnm),
      .fns = function(x) rep(mean(x), 100)
    )
  )

df_yh0 <- foreach(j = seq_len(length(cnm)),
                  .combine = bind_rows) %do% {
                    
                    ## reset df_x
                    df_x <- df_x0g
                    
                    ## get min-max range for a focus predictor
                    cid <- which(colnames(df_fcl_wsd) %in% c("ecor", cnm[j]))
                    
                    x <- df_fcl_wsd %>% 
                      dplyr::select(all_of(cid)) %>% 
                      rename_with(
                        .fn = function(z) ifelse(z %in% cnm, "x", z)
                      ) %>% 
                      group_by(ecor) %>% 
                      reframe(
                        x = seq(min(x),
                                max(x),
                                length = 100)
                      ) %>% 
                      pull(x)
                    
                    ## replace the mean with the range of the focus variable
                    df_x[[cnm[j]]] <- x
                    
                    ## predicted values of y = fcl
                    df_x <- df_x %>% 
                      mutate(
                        i = row_number(),
                        focus = cnm[j]
                      ) %>% 
                      relocate(focus)
                    
                    X0 <- df_x %>% 
                      dplyr::select(all_of(cnm))
                    
                    X <- model.matrix(~ ., X0)
                    
                    Y <- X %*% B
                    
                    ## Y's column corresponds to each region
                    ## Y[i, h] switch columns based on `h` column
                    cout <- df_x %>% 
                      rowwise() %>% 
                      mutate(log_y = Y[i, ecor],
                             y = exp(log_y)) %>% 
                      ungroup()
                  }

## get unscaled values
df_yh <- df_yh0 %>% 
  mutate(r_length = 
           exp(
             scl_log_rl * sd(df_fcl_wsd$log_rl) + mean(df_fcl_wsd$log_rl)
           ),
         lambda = exp(scl_log_lambda * sd(df_fcl_wsd$log_lambda) + 
                        mean(df_fcl_wsd$log_lambda)),
         prec = scl_prec * sd(df_fcl_wsd$prec) + mean(df_fcl_wsd$prec),
         temp = scl_temp * sd(df_fcl_wsd$temp) + mean(df_fcl_wsd$temp),
         hfp = scl_hfp * sd(df_fcl_wsd$hfp) + mean(df_fcl_wsd$hfp))


# overall prediction ------------------------------------------------------

x <- with(df_fcl_wsd,
          seq(min(log(lambda)),
              max(log(lambda)),
              length = 100))

X <- with(df_fcl_wsd,
          cbind(1,
                mean(log(r_length)),
                x))

## set mean values for each predictor
df_x0 <- df_fcl_wsd %>% 
  dplyr::select(all_of(c("ecor", cnm))) %>% 
  reframe(
    across(all_of(cnm),
           .fns = function(x) rep(mean(x), 100)
    )
  )

v_b <- df_est %>% 
  filter(parms_gr == "b") %>% 
  pull(median)

mcmc_b <- with(list_est[[id_best]],
               MCMCvis::MCMCchains(mcmc) %>% 
                 .[, str_detect(colnames(.), "b\\[.{1,}\\]")])

df_y <- foreach(j = seq_len(length(cnm)),
                .combine = bind_rows) %do% {
                  
                  ## reset df_x
                  df_x <- df_x0
                  
                  ## get min-max range for a focus predictor
                  cid <- which(colnames(df_fcl_wsd) %in% c("ecor", cnm[j]))
                  
                  x <- df_fcl_wsd %>% 
                    select(all_of(cid)) %>% 
                    rename_with(.fn = function(z) ifelse(z %in% cnm, "x", z)) %>% 
                    reframe(x = seq(min(x),
                                    max(x),
                                    length = 100)) %>% 
                    pull(x)
                  
                  ## replace the mean with the range of the focus variable
                  df_x[[cnm[j]]] <- x
                  
                  ## predicted values of y = fcl
                  df_x <- df_x %>% 
                    mutate(focus = cnm[j]) %>% 
                    relocate(focus)
                  
                  X0 <- df_x %>% 
                    dplyr::select(all_of(cnm))
                  
                  X <- model.matrix(~ ., X0)
                  
                  mu_y <- drop(X %*% v_b)
                  mcmc_y <- X %*% t(mcmc_b)
                  
                  ## get quantiles for predicted values
                  cout <- mcmc_y %>% 
                    apply(MARGIN = 1, quantile, c(0.025, 0.975)) %>% 
                    t() %>% 
                    as_tibble() %>% 
                    setNames(c("log_y_low", "log_y_high")) %>% 
                    mutate(log_y = mu_y) %>% 
                    mutate(y = exp(log_y),
                           y_low = exp(log_y_low),
                           y_high = exp(log_y_high)) %>% 
                    bind_cols(df_x) %>% 
                    mutate(r_length = exp(scl_log_rl * sd(df_fcl_wsd$log_rl) + 
                                            mean(df_fcl_wsd$log_rl)),
                           lambda = exp(scl_log_lambda * sd(df_fcl_wsd$log_lambda) + 
                                          mean(df_fcl_wsd$log_lambda)),
                           prec = scl_prec * sd(df_fcl_wsd$prec) + mean(df_fcl_wsd$prec),
                           temp = scl_temp * sd(df_fcl_wsd$temp) + mean(df_fcl_wsd$temp),
                           hfp = scl_hfp * sd(df_fcl_wsd$hfp) + mean(df_fcl_wsd$hfp)) %>% 
                    relocate(focus)
                }

