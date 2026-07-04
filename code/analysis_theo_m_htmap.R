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

withr::with_seed(10, {
  v_theta <- rep(0.25, each = n_rep)
  v_n_base <- rpois(length(v_theta), lambda = S * 0.19)
  v_l <- rpois(length(v_theta), lambda = S^2 * 0.11)
})

list_fw <- lapply(seq_len(length(v_theta)), function(i) {
  ## - set seed for reproducibility
  withr::with_seed(i * 10, {
    ecotools::ppm(n_species = S,
                  n_basal = v_n_base[i],
                  l = v_l[i],
                  theta = v_theta[i])
  })
})

## other parameters
## - rl, ecosystem size
## - lambda, branching rate
## - h, habitat density
## - delta0, dispersal capability for basal
## - r0, resource supply rate
## - b, stream size influence on establishment prob of basal species
##  -- expressed in a form of p = r0 + b * u
## - g, propagule number
## - mu0, disturbance rate
## - mu_p, prey-induced extinction rate
## - mu_c, predator-induced extinction rate
## - rho0, baseline synchrony prob.
## - z, scaling exponent for delta and g
## - fw, foodweb index
v_rl <- seq(100, 1000, length = 10)
v_lambda <- seq(0.1, 1, length = 10)
v_r0 <- c(0.4, 0.8)

parms <- expand.grid(rl = v_rl,
                     lambda = v_lambda,
                     h = 1,
                     delta0 = 0.1,
                     r0 = v_r0,
                     b = (1 - max(v_r0)) / max(v_rl),
                     g = 5,
                     mu0 = c(2.5, 5),
                     mu_p = 5,
                     mu_c = 0, # must be mu_c = 0 for analytical prediction
                     rho0 = c(0, 0.5),
                     nu = 1 / max(v_rl),
                     z = 0.5,
                     fw = seq_len(length(list_fw))) %>% 
  mutate(theta = v_theta[fw],
         n_base = v_n_base[fw],
         link = v_l[fw]) %>% 
  as_tibble()

## - v, scaling factor for non-spatial disturbance
## - mu0 is multiplied by v to make it comparable with spatial case
## - mu0_scl = mu0 * (1 + rho * u_hat(lambda_mid, rl_mid))
v <- with(parms, {
  
  diam <- rpom::diameter(
    lambda = mean(range(lambda)),
    size = mean(range(rl))
  )
  
  u <- rpom::u_length(
    lambda = mean(range(lambda)),
    size = mean(range(rl))
  )
  
  rho <- (1 - unique(nu) * (diam/3))
  
  1 + rho * u
})

parms <- parms %>% 
  mutate(
    mu0_scl = 
      ifelse(rho0 == 0,
             mu0 * v,
             mu0)
  )

# run ---------------------------------------------------------------------

# parallel setup
n_workers <- floor(0.8 * detectCores())
plan(multisession, workers = n_workers)
registerDoFuture()
handlers(global = TRUE)

tic(msg = "create stream vector")

df_y <- 
  with_progress({
    
    p <- progressor(steps = nrow(parms))
    
    foreach(i = seq_len(nrow(parms)),
            .combine = bind_rows,
            .options.future = list(
              packages = c("rpom"),
              seed = TRUE
            )) %dofuture% {
              
              y <- with(parms, {
                ## trophic position
                v_tp <- attr(list_fw[[fw[i]]], "tp")
                
                ## analytical output
                y0 <- rpom::fcl(w = list_fw[[fw[i]]],
                                lambda = lambda[i],
                                size = rl[i],
                                h = h[i],
                                delta = delta0[i] * v_tp^z[i],
                                r0 = r0[i],
                                b = b[i],
                                g = g[i] * v_tp^(-z[i]),
                                mu0 = mu0_scl[i],
                                mu_p = mu_p[i],
                                rho0 = rho0[i],
                                weight = TRUE,
                                exact = TRUE)
                
                return(y0)
              })
              
              ## progress
              p()
              
              ## output
              dplyr::tibble(
                fcl = y,
                s_alpha = sum(attr(y, "p_hat")),
                s_gamma = sum(attr(y, "p_hat") > 0)
              )
            } # foreach
    
  }) # with_progress

toc()

## free RAM + reset workers
plan(sequential)
gc()

# export ------------------------------------------------------------------

df_fcl <- parms %>% 
  bind_cols(df_y)

saveRDS(df_fcl, file = "data_fmt/sim_fcl_m_heat.rds")
