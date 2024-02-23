
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")
source("code/format_data4jags.R")


# set data ----------------------------------------------------------------

## set data for JAGS
## - local level data
list_local <- with(df_fcl,
                   list(logY = log(fcl_obs),
                        C = censoring,
                        logY_min = ifelse(is.na(fcl_min), 0, log(fcl_min)),
                        Hsize = local_area,
                        Agri = frac_agri,
                        G = g,
                        Ns = length(fcl),
                        Nw = n_distinct(g))
)

## - watershed level data
list_wsd <- with(df_g,
                 list(Esize = area,
                      Pbranch = p_branch,
                      Prec = mean.prec,
                      Temp = mean.temp,
                      Hfp = hfp,
                      H = h,
                      Nh = n_distinct(h),
                      Score = score))

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
parms <- c("a", "b0", "b", "sigma", "z", "nu")

## model files
m <- runjags::read.jagsfile("code/model_hier.R")

## run model
post <- runjags::run.jags(model = m$model,
                          data = list_jags,
                          monitor = parms,
                          burnin = 10000,
                          sample = 1000,
                          thin = 20,
                          n.chains = n_chain,
                          inits = inits,
                          method = "parallel",
                          module = "glm")

## summary
MCMCvis::MCMCsummary(post$mcmc)
