
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")
source("code/format_emp_data4jags.R")


# set data ----------------------------------------------------------------

## read data
## - df_fcl_local, local level data
## - df_fcl_wsd, watershed level data
list_fcl <- readRDS("data_fmt/data_fcl_reg.rds")
df_fcl_local <- list_fcl[[1]] %>% 
  rename(local_hfp = hfp) %>% 
  left_join(list_fcl[[2]] %>% 
              dplyr::select(uid, r_length, lambda, prec, temp, hfp)) %>% 
  mutate(resid_elev = resid(lm(log(local_elev) ~ prec, .)),
         resid_temp = resid(lm(temp ~ prec, .)))

df_fcl_wsd <- list_fcl[[2]]

## set data for JAGS
## - local level data
X1 <- df_fcl_local %>% 
  mutate(log_area = log(local_area) - mean(log(local_area))) %>% 
  dplyr::select(local_elev) %>% 
  mutate(across(.cols = -starts_with("log"),
                .fns = function(x) c(scale(x)))) %>% 
  data.matrix()

list_const_local <- with(df_fcl_local,
                         list(logCut = log(cut),
                              X1 = X1,
                              G = g,
                              K1 = ncol(X1),
                              Ns = length(fcl),
                              Nw = n_distinct(g)))

## - watershed level data
df_x2 <- df_fcl_wsd %>% 
  mutate(log_rl = log(r_length) - mean(log(r_length)),
         log_lambda = log(lambda) - mean(log(lambda))) %>% 
  dplyr::select(log_rl,
                log_lambda,
                temp,
                prec,
                hfp) %>% 
  mutate(across(.cols = -starts_with("log"),
                .fns = function(x) c(scale(x))))


X2 <- model.matrix(~ ., df_x2)

list_const_wsd <- with(df_fcl_wsd,
                       list(X2 = X2,
                            K2 = ncol(X2),
                            H = h,
                            Nh = n_distinct(h),
                            Ratio = d_ratio,
                            N_site = n_site))

## combine
list_data <- with(df_fcl_local,
                  list(logY = log(fcl_obs),
                       C = censoring))

list_const <- c(list_const_local, list_const_wsd)


# mcmc setup --------------------------------------------------------------

n_iter <- 3e+4
n_sample <- 1000
n_burn <- floor(n_iter / 2)
n_thin <- (n_iter - n_burn) / n_sample
n_chain <- 3
s0 <- 0.01


# model run ---------------------------------------------------------------

## random-intercept model ####

## model object `m0`
source("code/model_nimble_h0.R")

## initial value
list_inits_h0 <- with(df_fcl_local,
                      list(z = runif(2, min = 0, max = 1),
                           nu = runif(1, min = 4, max = 5),
                           a = rnorm(1, sd = s0),
                           b = rnorm(ncol(X2), sd = s0),
                           sigma = runif(3, min = 0.1, max = 0.5),
                           r = rnorm(n_distinct(h), sd = s0),
                           a0 = rnorm(n_distinct(uid),
                                      mean = mean(log(fcl_obs),
                                                  na.rm = TRUE),
                                      sd = s0),
                           logY = ifelse(censoring == 1,
                                         log(cut + 1),
                                         NA))
)

## model setup as nimbleModel
post_h0 <- nimbleMCMC(code = m0,
                      constants = list_const,
                      data = list_data,
                      inits = list_inits_h0,
                      niter = n_iter, 
                      nburnin = n_burn,
                      thin = n_thin,
                      nchains = n_chain,
                      progressBar = TRUE,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)

pr_neg <- MCMCvis::MCMCpstr(post_h0$samples,
                            func = function(x) mean(x < 0)) %>% 
  unlist()

(df_est_h0 <- MCMCvis::MCMCsummary(post_h0$samples) %>% 
    as_tibble(rownames = "parms") %>% 
    transmute(parms,
              median = `50%`,
              low = `2.5%`,
              high = `97.5%`,
              rhat = Rhat,
              pr_neg = pr_neg,
              pr_pos = 1 - pr_neg) %>% 
    relocate(pr_neg, pr_pos,
             .before = rhat))


## random-slope model ####

## model object `m1`
source("code/model_nimble_h1.R")

## initial value
list_inits_h1 <- with(df_fcl_local,
                      list(z = runif(2, min = 0, max = 1),
                           nu = runif(1, min = 4, max = 5),
                           a = rnorm(1, sd = s0),
                           b = matrix(rnorm(n_distinct(h) * ncol(X2),
                                            sd = s0),
                                      nrow = n_distinct(h),
                                      ncol = ncol(X2)),
                           sigma = runif(2),
                           v_mu_b = rnorm(ncol(X2), sd = s0),
                           v_sigma_b = runif(ncol(X2), min = 0, max = 1),
                           Ustar = diag(ncol(X2)),
                           a0 = rnorm(n_distinct(uid),
                                      mean = mean(log(fcl_obs),
                                                  na.rm = TRUE),
                                      sd = s0),
                           logY = ifelse(censoring == 1,
                                         log(cut + 1),
                                         NA))
)

## model setup as nimbleModel
post_h1 <- nimbleMCMC(code = m1,
                      constants = list_const,
                      data = list_data,
                      inits = list_inits_h1,
                      niter = n_iter, 
                      nburnin = n_burn,
                      thin = n_thin,
                      nchains = n_chain,
                      progressBar = TRUE,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE)

pr_neg <- MCMCvis::MCMCpstr(post_h1$samples,
                            func = function(x) mean(x < 0)) %>% 
  unlist()

(df_est_h1 <- MCMCvis::MCMCsummary(post_h1$samples) %>% 
    as_tibble(rownames = "parms") %>% 
    transmute(parms,
              median = `50%`,
              low = `2.5%`,
              high = `97.5%`,
              rhat = Rhat,
              pr_neg = pr_neg,
              pr_pos = 1 - pr_neg) %>% 
    relocate(pr_neg, pr_pos,
             .before = rhat))
