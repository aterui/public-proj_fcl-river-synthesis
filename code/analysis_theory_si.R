#' DESCRIPTION:
#' This script explores possible outcomes of food chain length
#' in a variety of ecological scenarios

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

# set parameters ----------------------------------------------------------

## generate food webs
## - S: number of nodes/species in a community
## - theta: degree of omnivory
S <- 32

set.seed(10)
v_theta <- rep(c(0.25, 0.50), each = 5)
v_n_base <- rpois(length(v_theta), lambda = S * 0.18)
v_l <- rpois(length(v_theta), lambda = S^2 * 0.11)

list_fw <- lapply(seq_len(length(v_theta)), function(i) {
  ## - set seed for reproducibility
  set.seed(i * 10)
  
  mcbrnet::ppm(n_species = S,
               n_basal = v_n_base[i],
               l = v_l[i],
               theta = v_theta[i])
})

## other parameters
## - rl, ecoystem size
## - lambda, branching rate
## - h, habitat density
## - delta0, dispersal capability for basal
## - rsrc, resource supply rate
## - g, propagule number
## - mu0, disturbance rate
## - mu_p, prey-induced extinction rate
## - mu_c, predator-induced extinction rate
## - rho, synchrony prob.
## - fw, foodweb index
parms <- expand.grid(rl = seq(10, 100, length = 20),
                     lambda = seq(0.1, 1, length = 20),
                     h = 10,
                     delta0 = 0.05,
                     rsrc = c(0.25, 0.5),
                     g = c(150, 300),
                     mu0 = c(0.25, 2.5),
                     mu_p = c(2.5, 5),
                     mu_c = c(0, 1),
                     rho = c(0, 0.5),
                     z = 0.5,
                     fw = seq_len(length(list_fw))) %>% 
  mutate(theta = v_theta[fw],
         n_base = v_n_base[fw],
         link = v_l[fw]) %>% 
  as_tibble()

## - nt0, # time steps for numerical simulations
## - tol, tolerance value for convergence check
## - p_th, absorbing condition for p
nt0 <- 50
max_nt <- 1000
tol <- 1e-5
p_th <- 1e-5

# run ---------------------------------------------------------------------

ncore <- floor(0.8 * detectCores())
cl <- makeCluster(ncore)
registerDoSNOW(cl)

pb <- txtProgressBar(min = 0,
                     max = nrow(parms),
                     initial = 0,
                     style = 3) 
fun_progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = fun_progress)

tictoc::tic()
df_y <- foreach(i = seq_len(nrow(parms)),
                .combine = bind_rows,
                .options.snow = opts) %dopar% {
                  
                  y <- with(parms, {
                    
                    ## vector of binary trophic positions if all present
                    v_tp <- attr(list_fw[[parms$fw[i]]], "tp")
                    
                    if (mu_c[i] == 0) {
                      ## w/o predation effect
                      
                      ## - analytical equilibrium with rpom::fcl()
                      y0 <- rpom::fcl(foodweb = list_fw[[fw[i]]],
                                      lambda = lambda[i],
                                      size = rl[i],
                                      h = h[i],
                                      delta = delta0[i] * v_tp^z[i],
                                      rsrc = rsrc[i],
                                      g = g[i] * v_tp^z[i],
                                      mu0 = mu0[i],
                                      mu_p = mu_p[i],
                                      rho = rho[i],
                                      weight = TRUE)
                    } else {
                      ## w/ predation effect
                      
                      ## - initialize
                      y0 <- NA
                      attr(y0, "p_hat") <- 0.5
                      attr(y0, "convergence") <- 1
                      nt <- 0
                      
                      ## - numerical equilibrium with rpom::nfcl()
                      while (attr(y0, "convergence") > 0) {
                        y0 <- rpom::nfcl(foodweb = list_fw[[fw[i]]],
                                         lambda = lambda[i],
                                         size = rl[i],
                                         h = h[i],
                                         delta = delta0[i] * v_tp^z[i],
                                         rsrc = rsrc[i],
                                         g = g[i] * v_tp^z[i],
                                         mu0 = mu0[i],
                                         mu_p = mu_p[i],
                                         mu_c = mu_c[i],
                                         rho = rho[i],
                                         x0 = attr(y0, "p_hat"),
                                         n_timestep = nt0,
                                         interval = 0.01,
                                         threshold = p_th,
                                         n_plus = 10,
                                         weight = TRUE,
                                         tol = tol)
                        
                        ## - # time steps updated
                        nt <- nt + nt0
                        attr(y0, "nt") <- nt
                        if (nt >= max_nt) break
                      }
                      
                    }
                    
                    return(y0)
                  })
                  
                  nt <- ifelse(is.null(attr(y, "nt")),
                               NA,
                               attr(y, "nt")) 
                  
                  zo <- ifelse(is.null(attr(y, "convergence")),
                               NA,
                               attr(y, "convergence"))
                  
                  return(dplyr::tibble(fcl = y,
                                       nt = nt,
                                       converge = zo))
                }
tictoc::toc()

stopCluster(cl)

# export ------------------------------------------------------------------

df_fcl <- parms %>% 
  bind_cols(df_y)

saveRDS(df_fcl, file = "data_fmt/sim_fcl_si.rds")
