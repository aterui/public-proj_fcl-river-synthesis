#' DESCRIPTION:
#' This script draw representative prediction for some parameter combos
#' Dotted lines in the heatmap, `figure_fcl_theory_main.R`
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

## parms data frame
## - get vectors for lambda and river length
n0 <- 200
rho0 <- c(0, 0.5)
x_lambda <- rep(seq(0.1, 1, length = n0),
                times = length(rho0))
x_rl <- rep(seq(10, 100, length = n0), 
            times = length(rho0))
c_lambda <- rep(mean(range(x_lambda)),
                times = length(x_lambda))
c_rl <- rep(mean(range(x_rl)),
            times = length(x_rl))

## - fix values when varying another
lambda <- c(x_lambda, c_lambda)
rl <- c(c_rl, x_rl)

## - assemble x values along with `focus` and `id`
df_x <- tibble(lambda,
               rl,
               focus = rep(c("branch", "size"),
                           each = length(lambda) * 0.5)) %>% 
  mutate(id = row_number())

## - v, scaling factor for non-spatial disturbance
## - mu0 is multiplied by v to make it comparable with spatial case
## - mu0_scl = mu0 * (1 + rho * u_hat(lambda_mid, rl_mid))
v <- 1 + max(rho0) * rpom::u_length(unique(c_lambda), unique(c_rl))

## - combine with other parameters
parms <- expand.grid(id = seq_len(length(lambda)),
                     rsrc = c(0.4, 0.8),
                     mu0 = c(2.5, 5),
                     mu_p = 5,
                     delta = 0.5,
                     h = 2.5,
                     g = 150,
                     rho = rho0,
                     z = 0.5,
                     fw = seq_len(length(list_fw))) %>% 
  left_join(df_x) %>% 
  mutate(mu0_scl = ifelse(rho == 0,
                          mu0 * v,
                          mu0))

# prediction --------------------------------------------------------------

## parallel setup
ncore <- floor(0.8 * detectCores())
cl <- makeCluster(ncore)
registerDoSNOW(cl)

pb <- txtProgressBar(min = 0,
                     max = nrow(parms),
                     initial = 0,
                     style = 3) 
fun_progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = fun_progress)

## run
tictoc::tic()
y <- foreach(i = seq_len(nrow(parms)),
             .combine = c,
             .options.snow = opts) %dopar% {
               
               y <- with(parms, {
                 
                 ## trophic position
                 v_tp <- attr(list_fw[[fw[i]]], "tp")
                 
                 ## analytical calculation
                 y0 <- rpom::fcl(foodweb = list_fw[[fw[i]]],
                                 lambda = lambda[i],
                                 size = rl[i],
                                 rsrc = rsrc[i],
                                 mu0 = mu0_scl[i],
                                 mu_p = mu_p[i],
                                 delta = delta[i] * v_tp^z[i],
                                 h = h[i],
                                 g = g[i] * v_tp^(-z[i]),
                                 rho = rho[i])
                 
                 return(y0)
               })
               
               return(y)
             }
tictoc::toc()

stopCluster(cl)

# export ------------------------------------------------------------------

df_fcl <- parms %>% 
  mutate(fcl = y)

saveRDS(df_fcl, file = "data_fmt/sim_fcl_m_line.rds")