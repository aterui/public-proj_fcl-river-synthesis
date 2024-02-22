
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")
source("code/format_data4jags.R")


# set data ----------------------------------------------------------------

## set data for JAGS
## - local level data
list_jags <- with(df_fcl,
                  list(logY = log(fcl),
                       C = censoring,
                       logY_min = ifelse(is.na(fcl_min), 0, log(fcl_min)),
                       Hsize = local_area,
                       Forest = mean.forest,
                       Prec = mean.prec,
                       Temp = mean.temp,
                       Ns = length(fcl),
                       Nw = n_distinct(uid),
                       G = g)
)

## - watershed level data
list_jags$Esize <- as.numeric(df_g$area)
list_jags$Pbranch <- df_g$p_branch
list_jags$W <- df_g$w


# jags fit ----------------------------------------------------------------

## run.jags arguments
## - initial values
inits <- replicate(3,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:3) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1

## - parameters to be monitored
parms <- c("a", "b", "sigma")

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
