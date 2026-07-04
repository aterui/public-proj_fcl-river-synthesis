#' DESCRIPTION:
#' Explores possible outcomes of food chain length
#' Expanded analysis using numerical methods
#' 
#' PARAMETERS:
#' - rl, ecosystem size
#' - lambda, branching rate
#' - h, habitat density
#' - delta0, dispersal capability for basal
#' - r0, resource supply rate
#' - b, stream size influence on establishment prob of basal species
#'    -- expressed in a form of p = r0 + b * u
#' - g, propagule number
#' - mu0, disturbance rate
#' - mu_p, prey-induced extinction rate
#' - mu_c, predator-induced extinction rate
#' - rho0, baseline synchrony prob.
#' - z, scaling exponent for delta and g
#' - fw, foodweb index

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

# parameters --------------------------------------------------------------
## generate food webs
## - S: number of nodes/species in a community
## - theta: degree of omnivory
S <- 32
n_rep <- 5

withr::with_seed(10, {
  v_theta <- rep(c(0.25, 0.50), each = n_rep)
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

## fixed parameters
v_rl <- seq(100, 1000, length = 20)
v_lambda <- seq(0.2, 1, length = 20)

df_parms <- 
  expand_grid(
    rl = v_rl,
    lambda = v_lambda,
    h = 1,
    g = 5,
    delta0 = 0.1,
    r0 = c(0.4, 0.8),
    rho0 = 1,
    nu = c(0.05, 0.1),
    kernel_type = "exp",
    mu0 = c(2.5, 5),
    mu_p = c(2.5, 5),
    mu_c = c(1.25, 2.5),
    estb_type = "prey",
    z = 0.5,
    svar = 0.1, # species-level variation in g and delta
    fw = seq_len(length(list_fw))
  ) %>% 
  mutate(
    b = (1 - max(r0)) / max(rl),
    theta = v_theta[fw],
    n_base = v_n_base[fw],
    link = v_l[fw]
  ) %>% 
  as_tibble()

## - nt0, # time steps for numerical simulations
## - tol, tolerance value for convergence check
## - p_th, absorbing condition for p
nt0 <- 50
max_nt <- 1000
tol <- 1e-5
p_th <- 1e-5

# run ---------------------------------------------------------------------

## RUN
n_workers <- floor(0.8 * detectCores())
plan(multisession, workers = n_workers)
registerDoFuture()
handlers(global = TRUE)

tic(msg = "fcl heatmap run")

df_y <- 
  with_progress({
    
    p <- progressor(steps = nrow(df_parms))
    
    foreach(i = seq_len(nrow(df_parms)),
            .combine = bind_rows,
            .options.future = list(
              packages = c("rpom"),
              seed = TRUE
            )) %dofuture% {
              
              y <- with(df_parms[i, ], {
                
                ## trophic position
                v_tp <- attr(list_fw[[fw]], "tp")
                
                ## propagule
                v_g <- withr::with_seed(fw, {
                  rnorm(
                    n = length(v_tp),
                    mean = log(g * v_tp^(-z)),
                    sd = svar
                  ) %>% 
                    exp()
                })
                
                ## dispersal
                v_delta <- withr::with_seed(fw, {
                  rnorm(
                    n = length(v_tp),
                    mean = log(delta0 * v_tp^(-z)),
                    sd = svar
                  ) %>% 
                    exp()
                })
                
                ## - initialize
                y0 <- NA
                attr(y0, "p_hat") <- 0.5
                attr(y0, "convergence") <- 1
                nt <- 0
                
                ## - numerical equilibrium with rpom::nfcl()
                while (attr(y0, "convergence") > 0) {
                  y0 <- rpom::nfcl(w = list_fw[[fw]],
                                   size = rl,
                                   lambda = lambda,
                                   h = h,
                                   g = v_g,
                                   delta = v_delta,
                                   r0 = r0,
                                   b = b,
                                   mu0 = mu0,
                                   mu_p = mu_p,
                                   mu_c = mu_c,
                                   rho0 = rho0,
                                   nu = nu,
                                   kernel = as.character(kernel_type),
                                   x0 = attr(y0, "p_hat"),
                                   n_timestep = nt0,
                                   intv = 0.01,
                                   threshold = p_th,
                                   n_plus = 10,
                                   weight = TRUE,
                                   tol = tol,
                                   exact = TRUE)
                  
                  ## - # time steps updated
                  nt <- nt + nt0
                  attr(y0, "nt") <- nt
                  if (nt >= max_nt) break
                }
                
                return(y0)
              })
              
              nt <- ifelse(is.null(attr(y, "nt")),
                           NA,
                           attr(y, "nt")) 
              
              zo <- ifelse(is.null(attr(y, "convergence")),
                           NA,
                           attr(y, "convergence"))
              
              ## progress
              p()
              
              ## output
              dplyr::tibble(fcl = y,
                            s_alpha = sum(attr(y, "p_hat")),
                            s_gamma = sum(attr(y, "p_hat") > 0),
                            nt = nt,
                            converge = zo)
            } #foreach
    
  }) #with_progress

tictoc::toc()

## free RAM + reset workers
plan(sequential)
gc()

# export ------------------------------------------------------------------

df_fcl <- parms %>% 
  bind_cols(df_y)

saveRDS(df_fcl, file = "data_fmt/sim_fcl_si.rds")
