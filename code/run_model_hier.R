
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")
source("code/format_data4jags.R")


# set data ----------------------------------------------------------------

## set data for JAGS
## - local level data
list_local <- with(df_fcl,
                   list(logY = log(fcl_obs),
                        C = censoring,
                        logY_min = ifelse(is.na(fcl_min), 0, log(fcl_min)),
                        Hsize = local_area,
                        Forest = forest_b1km,
                        G = g,
                        Ns = length(fcl),
                        Nw = n_distinct(g)))

## - watershed level data
list_wsd <- with(df_g,
                 list(Esize = r_length,
                      Pbranch = p_branch,
                      Prec = mean.prec,
                      Temp = mean.temp,
                      Hfp = hfp,
                      H = h,
                      Nh = n_distinct(h),
                      Ratio = d_ratio,
                      N_site = n_site))

list_jags <- c(list_local, list_wsd)


# jags fit ----------------------------------------------------------------

## run.jags arguments
## - initial values
n_chain <- 3

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- 10 * j

## - parameters to be monitored
parms <- c("a", "b0", "b", "sigma", "z", "nu", "r", "a0")

## model files
m <- runjags::read.jagsfile("code/model_hier.R")

## run model
post <- runjags::run.jags(model = m$model,
                          data = list_jags,
                          monitor = parms,
                          burnin = 10000,
                          sample = 1000,
                          thin = 15,
                          n.chains = n_chain,
                          inits = inits,
                          method = "parallel",
                          module = "glm")


# export ------------------------------------------------------------------

## summary
pr_neg <- MCMCvis::MCMCpstr(post$mcmc, func = function(x) mean(x < 0)) %>% 
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
