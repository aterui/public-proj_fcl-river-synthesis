#' DESCRIPTION:
#' Explores possible outcomes of food chain length
#' No predation and analytical

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

# set parameters ----------------------------------------------------------

## generate food webs
## - S: number of nodes/species in a community
## - theta: degree of omnivory
S <- 32
n_rep <- 20

set.seed(10)
v_theta <- rep(0.25, each = n_rep)
v_n_base <- rpois(length(v_theta), lambda = S * 0.19)
v_l <- rpois(length(v_theta), lambda = S^2 * 0.11)

list_fw <- lapply(seq_len(length(v_theta)), function(i) {
  ## - set seed for reproducibility
  set.seed(i * 10)
  
  ecotools::ppm(n_species = S,
                n_basal = v_n_base[i],
                l = v_l[i],
                theta = v_theta[i])
})

## other parameters
## - rl, ecosystem size
## - lambda, branching rate
## - h, habitat density
## - delta0, dispersal capability for basal
## - rsrc, resource supply rate
## - zeta, stream size influence on establishment prob of basal species
##  -- expressed in a form of p = rsrc + zeta * u (truncate at zero and one)
## - g, propagule number
## - mu0, disturbance rate
## - mu_p, prey-induced extinction rate
## - mu_c, predator-induced extinction rate
## - rho, synchrony prob.
## - z, scaling exponent for delta and g
## - fw, foodweb index
v_rl <- seq(10, 100, length = 50)

parms <- expand.grid(rl = v_rl,
                     lambda = seq(0.1, 1, length = 50),
                     h = 2.5,
                     delta0 = 0.5,
                     rsrc = c(0.25, 0.5),
                     zeta = 1 / max(v_rl),
                     g = 150,
                     mu0 = c(2.5, 5),
                     mu_p = 5,
                     mu_c = 0, # must be mu_c = 0 for analytical prediction
                     rho = c(0, 0.5),
                     z = 0.5,
                     fw = seq_len(length(list_fw))) %>% 
  mutate(theta = v_theta[fw],
         n_base = v_n_base[fw],
         link = v_l[fw]) %>% 
  as_tibble()

## - v, scaling factor for non-spatial disturbance
## - mu0 is multiplied by v to make it comparable with spatial case
## - mu0_scl = mu0 * (1 + rho * u_hat(lambda_mid, rl_mid))
v <- with(parms, 
          1 + max(rho) * rpom::u_length(mean(range(lambda)), mean(range(rl))))

parms <- parms %>% 
  mutate(mu0_scl = ifelse(rho == 0,
                          mu0 * v,
                          mu0))

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
                    ## trophic position
                    v_tp <- attr(list_fw[[fw[i]]], "tp")
                    
                    ## analytical output
                    y0 <- rpom::fcl(w = list_fw[[fw[i]]],
                                    lambda = lambda[i],
                                    size = rl[i],
                                    h = h[i],
                                    delta = delta0[i] * v_tp^z[i],
                                    rsrc = rsrc[i],
                                    zeta = zeta[i],
                                    g = g[i] * v_tp^(-z[i]),
                                    mu0 = mu0_scl[i],
                                    mu_p = mu_p[i],
                                    rho = rho[i],
                                    weight = TRUE)
                    
                    return(y0)
                  })
                  
                  return(dplyr::tibble(fcl = y,
                                       s_alpha = sum(attr(y, "p_hat")),
                                       s_gamma = sum(attr(y, "p_hat") > 0)
                  )
                  )
                }
tictoc::toc()

stopCluster(cl)

# export ------------------------------------------------------------------

df_fcl <- parms %>% 
  bind_cols(df_y)

saveRDS(df_fcl, file = "data_fmt/sim_fcl_m_heat.rds")
