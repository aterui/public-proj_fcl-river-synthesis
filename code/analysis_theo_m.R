#' DESCRIPTION:
#' Explores possible outcomes of food chain length
#' No predation and analytical
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

## fixed parameters
v_rl <- seq(100, 1000, length = 10)
v_lambda <- seq(0.2, 1, length = 10)

df_parms <- 
  expand_grid(
    h = 1,
    g = 5,
    delta0 = 0.1,
    rho0 = c(0, 1),
    nu = 0.1,
    kernel_type = "exp",
    mu_p = 5,
    mu_c = 0, # must be mu_c = 0 for analytical prediction
    estb_type = "prey",
    z = 0.5,
    svar = 0.1, # species-level variation in g and delta
    fw = seq_len(length(list_fw))
  ) %>% 
  mutate(
    theta = v_theta[fw],
    n_base = v_n_base[fw],
    link = v_l[fw]
  ) %>% 
  as_tibble()

## - v, scaling factor for non-spatial disturbance
## - mu0 is multiplied by v to make it comparable with spatial case
## - mu0_scl = mu0 * (1 + rho * u_hat(lambda_mid, rl_mid))
mu_lambda <- mean(range(v_lambda))
mu_rl <- mean(range(v_rl))

v <- with(df_parms, {
  
  diam <- rpom::diameter(
    lambda = mu_lambda,
    size = mu_rl
  )
  
  u <- rpom::u_length(
    lambda = mu_lambda,
    size = mu_rl
  )
  
  rho <- laplace_rt(
    nu = unique(nu), 
    mu = diam, 
    exact = TRUE
  )
  
  1 + rho * u
})

# heatmap -----------------------------------------------------------------

## for heatmap
df_heat <- df_parms %>% 
  expand_grid(
    rl = v_rl,
    lambda = v_lambda,
    r0 = c(0.4, 0.8),
    mu0 = c(2.5, 5)
  ) %>% 
  mutate(
    b = (1 - max(r0)) / max(rl),
    mu0_scl = ifelse(
      rho0 == 0,
      mu0 * v,
      mu0
    )
  )

## RUN
n_workers <- floor(0.8 * detectCores())
plan(multisession, workers = n_workers)
registerDoFuture()
handlers(global = TRUE)

tic(msg = "fcl heatmap run")

df_y <- 
  with_progress({
    
    p <- progressor(steps = nrow(df_heat))
    
    foreach(i = seq_len(nrow(df_heat)),
            .combine = bind_rows,
            .options.future = list(
              packages = c("rpom"),
              seed = TRUE
            )) %dofuture% {
              
              y <- with(df_heat[i, ], {
                ## trophic position
                v_tp <- attr(list_fw[[fw]], "tp")
                
                v_g <- withr::with_seed(fw, {
                  rnorm(
                    n = length(v_tp),
                    mean = log(g * v_tp^(-z)),
                    sd = svar
                  ) %>% 
                    exp()
                })
                
                v_delta <- withr::with_seed(fw, {
                  rnorm(
                    n = length(v_tp),
                    mean = log(delta0 * v_tp^(-z)),
                    sd = svar
                  ) %>% 
                    exp()
                })
                
                ## analytical output
                y0 <- rpom::fcl(
                  w = list_fw[[fw]],
                  lambda = lambda,
                  size = rl,
                  h = h,
                  g = v_g,
                  delta = v_delta,
                  r0 = r0,
                  b = b,
                  rho0 = rho0,
                  nu = nu,
                  kernel = as.character(kernel_type),
                  mu0 = mu0_scl,
                  mu_p = mu_p,
                  estb_type = as.character(estb_type),
                  weight = TRUE,
                  exact = TRUE
                )
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

## EXPORT
df_fcl_heat <- df_heat %>% 
  bind_cols(df_y)

saveRDS(df_fcl_heat, file = "data_fmt/sim_fcl_m_heat.rds")

# line art ----------------------------------------------------------------

## for line art
df_line <- 
  bind_rows(
    focus("lambda", v_lambda, "rl", v_rl, n = 1),
    focus("rl", v_rl, "lambda", v_lambda, n = 1)
  ) %>% 
  mutate(
    r0 = 0.8,
    mu0 = 5
  ) %>% 
  expand_grid(df_parms) %>% 
  mutate(
    b = (1 - max(r0)) / max(rl),
    mu0_scl = ifelse(
      rho0 == 0,
      mu0 * v,
      mu0
    )
  )

## RUN
n_workers <- floor(0.8 * detectCores())
plan(multisession, workers = n_workers)
registerDoFuture()
handlers(global = TRUE)

tic(msg = "fcl line-art run")

df_y <- 
  with_progress({
    
    p <- progressor(steps = nrow(df_line))
    
    foreach(i = seq_len(nrow(df_line)),
            .combine = bind_rows,
            .options.future = list(
              packages = c("rpom"),
              seed = TRUE
            )) %dofuture% {
              
              y <- with(df_line[i, ], {
                
                ## trophic position
                v_tp <- attr(list_fw[[fw]], "tp")
                
                ## propagule
                ## - fixed by foodweb id, fw
                v_g <- withr::with_seed(fw, {
                  rnorm(
                    n = length(v_tp),
                    mean = log(g * v_tp^(-z)),
                    sd = 0.1
                  ) %>% 
                    exp()
                }) %>% print()
                
                ## dispersal
                ## - fixed by foodweb id, fw
                v_delta <- withr::with_seed(fw, {
                  rnorm(
                    n = length(v_tp),
                    mean = log(delta0 * v_tp^(-z)),
                    sd = 0.1
                  ) %>% 
                    exp()
                })
                
                ## analytical output
                y0 <- rpom::fcl(
                  w = list_fw[[fw]],
                  lambda = lambda,
                  size = rl,
                  h = h,
                  g = v_g,
                  delta = v_delta,
                  r0 = r0,
                  b = b,
                  rho0 = rho0,
                  nu = nu,
                  kernel = as.character(kernel_type),
                  mu0 = mu0_scl,
                  mu_p = mu_p,
                  estb_type = as.character(estb_type),
                  weight = TRUE,
                  exact = TRUE
                )
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

## EXPORT
df_fcl_line <- df_line %>% 
  bind_cols(df_y)

saveRDS(df_fcl_line, file = "data_fmt/sim_fcl_m_line.rds")
