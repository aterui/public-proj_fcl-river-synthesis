#' DESCRIPTION:
#' Run nimble models
#' Model code: `model_nimble_hi.R` or `model_nimble_hsl.R`

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")
source("code/format_emp_data4nimble.R")


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

n_iter <- 6e+4
n_sample <- 500
n_burn <- floor(n_iter / 2)
n_thin <- floor((n_iter - n_burn) / n_sample)
n_chain <- 4
s0 <- 0.01


# model run ---------------------------------------------------------------

## h0: random-intercept ###################################################

## model object `m0`
source("code/model_nimble_hi.R")

## initial value
f_inits_h0 <- function() {
  with(df_fcl_local,
       list(z = runif(2, min = 0.5, max = 0.55),
            nu = runif(1, min = 4, max = 5),
            a = rnorm(ncol(X1), sd = s0),
            b = rnorm(ncol(X2), sd = s0),
            sigma = runif(3, min = 0.05, max = 0.1),
            r = rnorm(n_distinct(h), sd = s0),
            a0 = rnorm(n_distinct(uid),
                       mean = mean(log(fcl_obs),
                                   na.rm = TRUE),
                       sd = s0),
            logY = ifelse(censoring == 1,
                          log(cut + 1),
                          NA))
  ) 
}

parms <- c("a",
           "a0",
           "b",
           "b0",
           "nu",
           "sigma",
           "z")

## model setup as nimbleModel
post_h0 <- nimbleMCMC(code = m0,
                      constants = list_const,
                      data = list_data,
                      inits = f_inits_h0,
                      niter = n_iter, 
                      nburnin = n_burn,
                      thin = n_thin,
                      nchains = n_chain,
                      monitors = parms,
                      progressBar = TRUE,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE,
                      setSeed = TRUE)

(df_est_h0 <- MCMCvis::MCMCsummary(post_h0$samples,
                                   func = function(x) mean(x < 0),
                                   func_name = "pr_neg") %>%  
    as_tibble(rownames = "parms") %>% 
    filter(!str_detect(parms, "a0\\[.{1,}\\]")) %>% 
    transmute(parms,
              median = `50%`,
              low = `2.5%`,
              high = `97.5%`,
              rhat = Rhat,
              pr_neg = pr_neg,
              pr_pos = 1 - pr_neg) %>% 
    relocate(pr_neg, pr_pos,
             .before = rhat))

## append variable names
vn <- with(list_const, c(colnames(X1), colnames(X2)))

df_est_h0 <- df_est_h0 %>% 
  mutate(varname = c(vn,
                     rep(NA, nrow(df_est_h0) - length(vn))))

## export
if(max(df_est_h0$rhat, na.rm = TRUE) < 1.1) {
  print("yes, converged!")
  
  saveRDS(list(est = df_est_h0,
               waic = post_h0$WAIC,
               mcmc = post_h0$samples,
               setup = c(iter = n_iter,
                         burn = n_burn,
                         thin = n_thin,
                         nchain = n_chain)),
          file = "data_fmt/output_model_h0.rds")
  
} else {
  print(paste("f*ck", max(df_est_h0$rhat, na.rm = TRUE)))
}


## h1: random-slope, length and branch ####################################

## model object `m1`
source("code/model_nimble_hsl.R")
R <- 3

## initial value
f_inits_h1 <- function() {
  with(df_fcl_local,
       list(z = runif(2, min = 0.5, max = 0.55),
            nu = runif(1, min = 4, max = 5),
            a = rnorm(ncol(X1), sd = s0),
            b = matrix(rnorm(n_distinct(h) * R,
                             sd = s0),
                       nrow = n_distinct(h),
                       ncol = ncol(X2)),
            b_mu = rnorm(R, sd = s0),
            b_prime = rnorm(ncol(X2) - R, sd = s0),
            sigma = runif(2, min = 0.05, max = 0.1),
            sigma_b = runif(R, min = 0.05, max = 0.1),
            Ustar = diag(ncol(X2)),
            a0 = rnorm(n_distinct(uid),
                       mean = mean(log(fcl_obs),
                                   na.rm = TRUE),
                       sd = s0),
            logY = ifelse(censoring == 1,
                          log(cut + 1),
                          NA))
  )
}

parms <- c("a",
           "a0",
           "b_mu",
           "b_prime", 
           "z",
           "nu",
           "sigma",       
           "sigma_b",
           "rho")

## model setup as nimbleModel
list_const$R <- R

post_h1 <- nimbleMCMC(code = m1,
                      constants = list_const,
                      data = list_data,
                      inits = f_inits_h1,
                      niter = n_iter, 
                      nburnin = n_burn,
                      thin = n_thin,
                      nchains = n_chain,
                      monitors = parms,
                      progressBar = TRUE,
                      samplesAsCodaMCMC = TRUE,
                      WAIC = TRUE,
                      setSeed = TRUE)

(df_est_h1 <- MCMCvis::MCMCsummary(post_h1$samples,
                                   func = function(x) mean(x < 0),
                                   func_name = "pr_neg") %>% 
    as_tibble(rownames = "parms") %>% 
    filter(!str_detect(parms, "a0\\[.{1,}\\]")) %>% 
    transmute(parms,
              median = `50%`,
              low = `2.5%`,
              high = `97.5%`,
              rhat = Rhat,
              pr_neg = pr_neg,
              pr_pos = 1 - pr_neg) %>% 
    relocate(pr_neg, pr_pos,
             .before = rhat))

## append variable names
vn <- with(list_const, c(colnames(X1), colnames(X2)))

df_est_h1 <- df_est_h1 %>% 
  mutate(varname = c(vn,
                     rep(NA, nrow(df_est_h1) - length(vn))))

## export
if(max(df_est_h1$rhat, na.rm = TRUE) < 1.1) {
  print("yes, converged!")
  
  saveRDS(list(est = df_est_h1,
               waic = post_h1$WAIC,
               mcmc = post_h1$samples,
               setup = c(iter = n_iter,
                         burn = n_burn,
                         thin = n_thin,
                         nchain = n_chain)),
          file = "data_fmt/output_model_h1.rds")
} else {
  print(paste("f*ck", max(df_est_h1$rhat, na.rm = TRUE)))
}

# ## h2: random-slope, length ###############################################
# 
# ## model object `m1`
# source("code/model_nimble_hsl.R")
# R <- 2
# list_const$X2 <- X2
# 
# ## initial value
# f_inits_h2 <- function() {
#   with(df_fcl_local,
#        list(z = runif(2, min = 0.5, max = 0.55),
#             nu = runif(1, min = 4.5, max = 5),
#             a = rnorm(1, sd = s0),
#             b = matrix(rnorm(n_distinct(h) * R,
#                              sd = s0),
#                        nrow = n_distinct(h),
#                        ncol = ncol(X2)),
#             b_mu = rnorm(R, sd = s0),
#             b_prime = rnorm(ncol(X2) - R, sd = s0),
#             sigma = runif(2, min = 0.05, max = 0.1),
#             sigma_b = runif(R, min = 0.05, max = 0.1),
#             Ustar = diag(ncol(X2)),
#             a0 = rnorm(n_distinct(uid),
#                        mean = mean(log(fcl_obs),
#                                    na.rm = TRUE),
#                        sd = s0),
#             logY = ifelse(censoring == 1,
#                           log(cut + 1),
#                           NA))
#   )
# }
# 
# parms <- c("a",
#            "a0",
#            "b_mu",
#            "b_prime",
#            "z",
#            "nu",
#            "sigma",
#            "sigma_b",
#            "rho")
# 
# ## model setup as nimbleModel
# list_const$R <- R
# 
# post_h2 <- nimbleMCMC(code = m1,
#                       constants = list_const,
#                       data = list_data,
#                       inits = f_inits_h2,
#                       niter = n_iter,
#                       nburnin = n_burn,
#                       thin = n_thin,
#                       nchains = n_chain,
#                       monitors = parms,
#                       progressBar = TRUE,
#                       samplesAsCodaMCMC = TRUE,
#                       WAIC = TRUE,
#                       setSeed = TRUE)
# 
# (df_est_h2 <- MCMCvis::MCMCsummary(post_h2$samples,
#                                    func = function(x) mean(x < 0),
#                                    func_name = "pr_neg") %>%
#     as_tibble(rownames = "parms") %>%
#     filter(!str_detect(parms, "a0\\[.{1,}\\]")) %>% 
#     transmute(parms,
#               median = `50%`,
#               low = `2.5%`,
#               high = `97.5%`,
#               rhat = Rhat,
#               pr_neg = pr_neg,
#               pr_pos = 1 - pr_neg) %>%
#     relocate(pr_neg, pr_pos,
#              .before = rhat))
# 
# ## append variable names
# vn <- with(list_const, c(colnames(X1), colnames(X2)))
# 
# df_est_h2 <- df_est_h2 %>%
#   mutate(varname = c(vn,
#                      rep(NA, nrow(df_est_h2) - length(vn))))
# 
# ## export
# if(max(df_est_h2$rhat, na.rm = TRUE) < 1.1) {
#   print("yes, converged!")
#   
#   saveRDS(list(est = df_est_h2,
#                waic = post_h2$WAIC,
#                mcmc = post_h2$samples,
#                setup = c(iter = n_iter,
#                          burn = n_burn,
#                          thin = n_thin,
#                          nchain = n_chain)),
#           file = "data_fmt/output_model_h2.rds")
# } else {
#   print(paste("f*ck", max(df_est_h2$rhat, na.rm = TRUE)))
# }
# 
# 
# ## h3: random-slope, branch ###############################################
# 
# ## model object `m1`
# source("code/model_nimble_hsl.R")
# R <- 2
# list_const$X2 <- X2 %>%
#   as_tibble() %>%
#   relocate(log_lambda, .before = log_rl) %>%
#   data.matrix()
# 
# ## initial value
# f_inits_h3 <- function() {
#   with(df_fcl_local,
#        list(z = runif(2, min = 0.5, max = 0.55),
#             nu = runif(1, min = 4, max = 5),
#             a = rnorm(1, sd = s0),
#             b = matrix(rnorm(n_distinct(h) * R,
#                              sd = s0),
#                        nrow = n_distinct(h),
#                        ncol = ncol(X2)),
#             b_mu = rnorm(R, sd = s0),
#             b_prime = rnorm(ncol(X2) - R, sd = s0),
#             sigma = runif(2, min = 0.05, max = 0.1),
#             sigma_b = runif(R, min = 0.05, max = 0.1),
#             Ustar = diag(ncol(X2)),
#             a0 = rnorm(n_distinct(uid),
#                        mean = mean(log(fcl_obs),
#                                    na.rm = TRUE),
#                        sd = s0),
#             logY = ifelse(censoring == 1,
#                           log(cut + 1),
#                           NA))
#   )
# }
# 
# parms <- c("a",
#            "a0",
#            "b_mu",
#            "b_prime",
#            "z",
#            "nu",
#            "sigma",
#            "sigma_b",
#            "rho")
# 
# ## model setup as nimbleModel
# list_const$R <- R
# 
# post_h3 <- nimbleMCMC(code = m1,
#                       constants = list_const,
#                       data = list_data,
#                       inits = f_inits_h3,
#                       niter = n_iter,
#                       nburnin = n_burn,
#                       thin = n_thin,
#                       nchains = n_chain,
#                       monitors = parms,
#                       progressBar = TRUE,
#                       samplesAsCodaMCMC = TRUE,
#                       WAIC = TRUE,
#                       setSeed = TRUE)
# 
# (df_est_h3 <- MCMCvis::MCMCsummary(post_h3$samples,
#                                    func = function(x) mean(x < 0),
#                                    func_name = "pr_neg") %>%
#     as_tibble(rownames = "parms") %>%
#     filter(!str_detect(parms, "a0\\[.{1,}\\]")) %>% 
#     transmute(parms,
#               median = `50%`,
#               low = `2.5%`,
#               high = `97.5%`,
#               rhat = Rhat,
#               pr_neg = pr_neg,
#               pr_pos = 1 - pr_neg) %>%
#     relocate(pr_neg, pr_pos,
#              .before = rhat))
# 
# ## append variable names
# vn <- with(list_const, c(colnames(X1), colnames(X2)))
# 
# df_est_h3 <- df_est_h3 %>%
#   mutate(varname = c(vn,
#                      rep(NA, nrow(df_est_h3) - length(vn))))
# 
# ## export
# if(max(df_est_h3$rhat, na.rm = TRUE) < 1.1) {
#   print("yes, converged!")
#   
#   saveRDS(list(est = df_est_h3,
#                waic = post_h3$WAIC,
#                mcmc = post_h3$samples,
#                setup = c(iter = n_iter,
#                          burn = n_burn,
#                          thin = n_thin,
#                          nchain = n_chain)),
#           file = "data_fmt/output_model_h3.rds")
# } else {
#   print(paste("f*ck", max(df_est_h3$rhat, na.rm = TRUE)))
# }
# 
# 
# # WAIC compare ------------------------------------------------------------
# 
# print(c(h0 = post_h0$WAIC$WAIC,
#         h1 = post_h1$WAIC$WAIC,
#         h2 = post_h2$WAIC$WAIC,
#         h3 = post_h3$WAIC$WAIC) %>%
#         sort())
