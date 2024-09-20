
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


# model setup -------------------------------------------------------------

## model object `m`
source("code/model_hier_nimble.R")

## initial value
list_inits <- with(df_fcl_local,
                   list(z = rep(0.5, 2),
                        nu = 5,
                        a = 0,
                        b = rep(0, 6),
                        sigma = rep(0.25, 3),
                        r = rep(0, n_distinct(h)),
                        a0 = rep(mean(fcl_obs, na.rm = TRUE),
                                 n_distinct(uid)),
                        logY = ifelse(censoring == 1,
                                      log(cut) + 1,
                                      NA)
                   )
)

## model setup as nimbleModel
nm <- nimbleModel(code = m,
                  constants = list_const,
                  data = list_data,
                  inits = list_inits)

## compile model
cnm <- compileNimble(nm) # make compiled version from generated C++

## build MCMC function
f_mcmc <- buildMCMC(nm,
                    monitors = c("a", "b", "sigma", "z", "nu"))
# f_mcmc <- buildMCMC(nm,
#                     monitors = nm$getNodeNames(stochOnly = TRUE,
#                                                includeData = FALSE))

cf_mcmc <- compileNimble(f_mcmc, project = nm)

## run MCMC
post <- runMCMC(cf_mcmc,
                niter = 2e+4,
                nburnin = 1e+4,
                thin = 20,
                nchains = 4,
                progressBar = TRUE,
                samplesAsCodaMCMC = TRUE,
                summary = TRUE)

pr_neg <- MCMCvis::MCMCpstr(post$samples,
                            func = function(x) mean(x < 0)) %>% 
  unlist()

(df_est <- MCMCvis::MCMCsummary(post$samples) %>% 
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
