
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")


# set parameters ----------------------------------------------------------

## generate food webs
## - S: number of nodes/species in a community
## - theta: degree of omnivory
S <- 32
n_rep <- 10

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
x_lambda <- seq(0.1, 1, length = 100)
x_rl <- seq(1, 100, length = 100)

## - fix values when varying another
lambda <- c(x_lambda, rep(0.5, length(x_lambda)))
rl <- c(rep(50, length(x_rl)), x_rl)

## - assemble x values along with `focus` and `id`
df_x <- tibble(lambda,
               rl,
               focus = rep(c("branch", "size"),
                           each = length(x_lambda))) %>% 
  mutate(id = row_number())

## - combine with other parameters
parms <- expand.grid(id = seq_len(length(lambda)),
                     rsrc = c(0.25, 0.5),
                     mu0 = c(0.25, 2.5),
                     mu_p = 2.5,
                     delta = 0.05,
                     h = 10,
                     g = 150,
                     rho = 0.5,
                     z = 0.5,
                     fw = seq_len(length(list_fw))) %>% 
  left_join(df_x)


# prediction --------------------------------------------------------------

pb <- txtProgressBar(min = 0,
                     max = nrow(parms),
                     initial = 0,
                     style = 3) 

y <- foreach(i = seq_len(nrow(parms)),
             .combine = c) %do% {
               
               setTxtProgressBar(pb, i)
               
               y <- with(parms, {
                 
                 v_tp <- attr(list_fw[[fw[i]]], "tp")
                 
                 y0 <- fcl(foodweb = list_fw[[fw[i]]],
                           lambda = lambda[i],
                           size = rl[i],
                           rsrc = rsrc[i],
                           mu0 = mu0[i],
                           mu_p = mu_p[i],
                           delta = delta[i] * v_tp^z[i],
                           h = h[i],
                           g = g[i] * v_tp^(-z[i]),
                           rho = rho[i])
                 
                 return(y0)
               })
               
               return(y)
             }

close(pb)

# export ------------------------------------------------------------------

df_fcl <- parms %>% 
  mutate(fcl = y)

saveRDS(df_fcl, file = "data_fmt/sim_fcl_main.rds")
