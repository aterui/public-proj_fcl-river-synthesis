
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
inits <- replicate(3,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:3) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1

## - parameters to be monitored
parms <- c("a", "b", "sigma", "z")

## model files
m <- runjags::read.jagsfile("code/model_hier.R")

## run model
post <- runjags::run.jags(model = m$model,
                          data = list_jags,
                          monitor = parms,
                          burnin = 10000,
                          sample = 1000,
                          thin = 10,
                          n.chains = 3,
                          inits = inits,
                          method = "parallel",
                          module = "glm")

## summary
MCMCvis::MCMCsummary(post$mcmc)
