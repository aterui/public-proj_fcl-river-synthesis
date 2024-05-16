
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")
source("code/format_data4jags.R")


# set data ----------------------------------------------------------------

## set data for JAGS
## - local level data
X1 <- df_fcl %>% 
  dplyr::select(local_area,
                forest_b1km,
                fsd,
                elev) %>% 
  mutate(across(.cols = -local_area,
                .fns = function(x) c(scale(x))),
         log_area = log(local_area)) %>% 
  dplyr::select(-local_area) %>% 
  relocate(log_area, fsd, elev, forest_b1km) %>% 
  data.matrix()

list_local <- with(df_fcl,
                   list(logY = log(fcl_obs),
                        C = censoring,
                        logY_min = ifelse(is.na(fcl_min),
                                          yes = 0,
                                          no = log(fcl_min)),
                        X1 = X1,
                        G = g,
                        K1 = ncol(X1),
                        Ns = length(fcl),
                        Nw = n_distinct(g)))

## - watershed level data
X2 <- df_g %>% 
  dplyr::select(r_length,
                lambda,
                mean.prec,
                mean.temp,
                hfp) %>% 
  mutate(across(.cols = -c(r_length, lambda),
                .fns = function(x) c(scale(x))),
         log_r_length = log(r_length),
         log_lambda = log(lambda)) %>% 
  dplyr::select(-c(r_length, lambda)) %>% 
  relocate(log_r_length, log_lambda) %>% 
  data.matrix()
  
list_wsd <- with(df_g,
                 list(X2 = X2,
                      K2 = ncol(X2),
                      H = h,
                      Nh = n_distinct(h),
                      Ratio = d_ratio,
                      N_site = n_site))

list_jags <- c(list_local, list_wsd)


# jags fit ----------------------------------------------------------------

## run.jags arguments
## - mcmc setup
n_ad <- 1000
n_iter <- 1.0E+4
n_sample <- 1000
n_thin <- max(3, ceiling(n_iter / n_sample))
n_burn <- ceiling(max(10, n_iter * 0.5))
n_chain <- 3

## - initial values
inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- j * 100

## - parameters to be monitored
parms <- c("a", "b", "b0", "sigma", "z", "nu", "r", "a0")

## model files
m <- runjags::read.jagsfile("code/model_hier.R")

## run model
post <- runjags::run.jags(model = m$model,
                          data = list_jags,
                          monitor = parms,
                          adapt = n_ad,
                          burnin = n_burn,
                          sample = n_sample,
                          thin = n_thin,
                          n.chains = n_chain,
                          inits = inits,
                          method = "parallel",
                          module = "glm")


# export ------------------------------------------------------------------

## summary
pr_neg <- MCMCvis::MCMCpstr(post$mcmc,
                            func = function(x) mean(x < 0)) %>% 
  unlist()

df_est <- MCMCvis::MCMCsummary(post$mcmc) %>% 
  as_tibble(rownames = "parms") %>% 
  transmute(parms,
            median = `50%`,
            low = `2.5%`,
            high = `97.5%`,
            rhat = Rhat,
            pr_neg = pr_neg,
            pr_pos = 1 - pr_neg) %>% 
  relocate(pr_neg, pr_pos, .before = rhat)

## - check convergence: rhat < 1.1
print(max(df_est$rhat))
  
df_summary <- df_est %>% 
  filter(!str_detect(parms, "a0|y_pred")) %>% 
  mutate(parms_gr = str_remove_all(parms, "\\[.\\]|\\d{1,}"),
         parms_num = str_extract(parms, "\\d{1,}")) %>% 
  relocate(parms, parms_gr, parms_num)

saveRDS(df_summary, "data_fmt/output_model_summary.rds")

## estimated watershed means
## - observed mean, not accounting for censoring
df_fcl_g <- df_fcl %>% 
  group_by(g) %>% 
  summarize(mu_fcl_obs = log(fcl) %>% 
              mean() %>% 
              exp(),
            tpc = ifelse(any(tpc == "N"), "N", "Y"))

## - estimated mean, after accounting for weighting and censoring
df_a0 <- df_est %>% 
  filter(str_detect(parms, "a0")) %>% 
  transmute(g = row_number(),
            fcl_est = exp(median),
            fcl_low = exp(low),
            fcl_high = exp(high))

df_wsd <- reduce(list(df_g, df_fcl_g, df_a0), left_join)

saveRDS(df_wsd, "data_fmt/output_model_mu_est.rds")

## mcmc samples
mcmc <- MCMCvis::MCMCchains(post$mcmc)
saveRDS(mcmc, "data_fmt/output_model_mcmc.rds")
