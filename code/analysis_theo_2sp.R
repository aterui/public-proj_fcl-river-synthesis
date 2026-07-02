#' DESCRIPTION:
#' Explores predator-prey system

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

## food web: simple predator-prey system
w <- matrix(0, 2, 2)
w[1, 2] <- w[2, 1] <- 1
v_tp <- c(1, 2)

# for line figures --------------------------------------------------------

## other parameters
v_rl <- seq(100, 1000, length = 100)
v_lambda <- seq(0.2, 0.8, length = 100)

df_parms <- bind_rows(
  focus("lambda", v_lambda, "rl", v_rl, n = 1),
  focus("rl", v_rl, "lambda", v_lambda, n = 1)
) %>% 
  expand_grid(
    h = 1,
    delta0 = c(0.01, 0.1),
    r0 = 0.4,
    rho0 = 1,
    nu = 0.1,
    g0 = 5,
    mu0 = 2.5,
    mu_p = 5,
    z = 0,
    kernel = "exp"
  ) %>% 
  mutate(b = (1 - max(r0)) / max(v_rl),
         cascade = ifelse(rho0 == 0, "No cascade", "Cascade"))

## set parameters when no disturbance cascade exists
v <- with(df_parms, {
  
  diam <- rpom::diameter(
    lambda = mean(range(lambda)),
    size = mean(range(rl))
  )
  
  u <- rpom::u_length(
    lambda = mean(range(lambda)),
    size = mean(range(rl))
  )
  
  rho <- laplace_rt(nu = unique(nu), mu = diam)#(1 - unique(nu) * (diam/3))
  
  1 + rho * u
})

df_parms <- df_parms %>% 
  mutate(
    mu0_scl = 
      ifelse(rho0 == 0,
             mu0 * v,
             mu0)
  )

## occupancy prediction
df_fcl <- foreach(i = seq_len(nrow(df_parms)),
                  .combine = bind_rows) %do% {
                    
                    df_i <- df_parms[i, ]        
                    
                    y <- 
                      with(df_i, {
                        rpom::fcl(w = w, 
                                  lambda = lambda, 
                                  size = rl,
                                  h = h,
                                  delta = delta0 * v_tp^(-z),
                                  r0 = r0,
                                  b = b,
                                  rho0 = rho0,
                                  nu = nu,
                                  g = g0 * v_tp^(-z),
                                  mu0 = mu0_scl,
                                  mu_p = mu_p,
                                  kernel = kernel,
                                  exact = TRUE) %>% 
                          attr("p_hat")
                      })
                    
                    tibble(o = y,
                           tl = c("prey", "predator")) %>% 
                      bind_cols(df_i)
                  }

## export
saveRDS(df_fcl, "data_fmt/sim_fcl_2sp_line.rds")

# df_fcl %>%
#   filter(focus == "rl") %>%
#   ggplot(
#     aes(
#       x = rl,
#       y = o,
#       color = tl,
#       linetype = factor(delta0)
#     )
#   ) +
#   geom_line() +
#   #facet_wrap(facets = ~delta0) +
#   scale_y_continuous(lim = c(0, 1))# +
# #scale_x_continuous(trans = "log10")
# 
# df_fcl %>%
#   filter(focus == "lambda") %>%
#   ggplot(
#     aes(x = lambda,
#         y = o,
#         color = tl,
#         linetype = factor(delta0)
#     )
#   ) +
#   geom_line() +
#   #facet_wrap(facets = ~delta0) +
#   scale_y_continuous(lim = c(0, 1))

# for heatmap -------------------------------------------------------------

## other parameters
df_parms <- expand_grid(
    rl = seq(100, 1000, length = 20),
    lambda = seq(0.2, 1, length = 20),
    h = 1,
    delta0 = c(0.05, 0.1),
    r0 = c(0.4, 0.8),
    rho0 = 1,
    nu = 0.1,
    g0 = 5,
    mu0 = c(2.5, 5),
    mu_p = 5,
    z = 0,
    kernel = "exp"
  ) %>% 
  mutate(b = (1 - max(r0)) / max(v_rl),
         cascade = ifelse(rho0 == 0, "No cascade", "Cascade"))

## occupancy prediction
df_fcl <- foreach(i = seq_len(nrow(df_parms)),
                  .combine = bind_rows) %do% {
                    
                    df_i <- df_parms[i, ]        
                    
                    y <- 
                      with(df_i, {
                        rpom::fcl(w = w, 
                                  lambda = lambda, 
                                  size = rl,
                                  h = h,
                                  delta = delta0 * v_tp^(-z),
                                  r0 = r0,
                                  b = b,
                                  rho0 = rho0,
                                  nu = nu,
                                  g = g0 * v_tp^(-z),
                                  mu0 = mu0,
                                  mu_p = mu_p,
                                  kernel = kernel,
                                  exact = TRUE) %>% 
                          attr("p_hat")
                      })
                    
                    tibble(o = y,
                           tl = c("prey", "predator")) %>% 
                      bind_cols(df_i)
                  }

## export
saveRDS(df_fcl, "data_fmt/sim_fcl_2sp_heat.rds")
