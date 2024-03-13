
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_function.R")
source("code/set_library.R")


# prediction --------------------------------------------------------------

set.seed(123)
S <- 32
fw <- ppm(n_species = S,
          n_basal = round(S * 0.18),
          l = round(S^2 * 0.11),
          theta = 0.25)

rsrc <- c(0.25, 2.5)
mu <- c(0.25, 2.5)
rho <- 0.5

## focus, branching
parms <- expand.grid(L = c(20, 80),
                     rsrc = rsrc,
                     mu = mu,
                     rho = rho)

lambda <- seq(-log(1 - 0.2), -log(1 - 0.8), length = 100)

df_fcl_pb <- foreach(i = seq_len(nrow(parms)),
                     .combine = bind_rows) %do% {
                       
                       y <- sapply(lambda,
                                   function(x) {
                                     with(parms, 
                                          fcl(foodweb = fw,
                                              lambda = x,
                                              L = L[i],
                                              rsrc = rsrc[i],
                                              mu_base = mu[i],
                                              mu_cnsm = mu[i],
                                              rho = rep(rho[i], 2))
                                     )
                                   })
                       
                       cout <- with(parms,
                                    tibble(fcl = y,
                                           lambda = lambda,
                                           p_branch = 1 - exp(-lambda),
                                           L = L[i],
                                           rsrc = rsrc[i],
                                           mu = mu[i],
                                           rho = rho[i],
                                           focus = "p_branch")
                       )
                       
                       return(cout)
                     }


## focus, size
parms <- expand.grid(p_branch = c(0.2, 0.6),
                     rsrc = rsrc,
                     mu = mu,
                     rho = rho)

L <- seq(10, 100, length = 100)

df_fcl_l <- foreach(i = seq_len(nrow(parms)),
                    .combine = bind_rows) %do% {
                      
                      y <- sapply(L,
                                  function(x) {
                                    with(parms, 
                                         fcl(foodweb = fw,
                                             lambda = -log(1 - p_branch[i]),
                                             L = x,
                                             rsrc = rsrc[i],
                                             mu_base = mu[i],
                                             mu_cnsm = mu[i],
                                             rho = rep(rho[i], 2))
                                    )
                                  })
                      
                      cout <- with(parms,
                                   tibble(fcl = y,
                                          lambda = -log(1 - p_branch[i]),
                                          p_branch = p_branch[i],
                                          L = L,
                                          rsrc = rsrc[i],
                                          mu = mu[i],
                                          rho = rho[i],
                                          focus = "size")
                      )
                      
                      return(cout)
                    }


# export ------------------------------------------------------------------

df_fcl <- bind_rows(df_fcl_pb, df_fcl_l)
saveRDS(df_fcl, file = "data_fmt/sim_fcl_analytical.rds")
