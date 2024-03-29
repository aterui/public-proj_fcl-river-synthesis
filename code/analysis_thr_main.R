
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_function.R")
source("code/set_library.R")


# set parameters ----------------------------------------------------------

## generate food webs
## - set seed for reproducibility
set.seed(123)

## - generate food web matrix
## - S: number of nodes/species in a community
S <- 32
list_fw <- replicate(10, {
  fw <- ppm(n_species = S,
            n_basal = round(S * 0.18),
            l = round(S^2 * 0.11),
            theta = 0.25)},
  simplify = FALSE)

## other parameters
rsrc <- c(0.25, 2.5)
mu <- c(0.25, 2.5)
rho <- 0.5
delta <- 1.5
fw <- seq_len(length(list_fw))

## parms data frame
## - get vectors for lambda and river length
x_lambda <- seq(-log(1 - 0.2), -log(1 - 0.8), length = 100)
x_rl <- seq(1, 100, length = 100)

## - fix values when varying another
lambda <- c(x_lambda, rep(-log(0.5), length(x_lambda)))
rl <- c(rep(50, length(x_rl)), x_rl)

## - assemble x values along with `focus` and `id`
df_x <- tibble(lambda,
               rl,
               focus = rep(c("branch", "size"),
                           each = length(x_lambda))) %>% 
  mutate(id = row_number())

## - combine with other parameters
parms <- expand.grid(id = seq_len(length(lambda)),
                     rsrc = rsrc,
                     mu = mu,
                     rho = rho,
                     delta = delta,
                     fw = fw) %>% 
  left_join(df_x)


# prediction --------------------------------------------------------------

y <- foreach(i = seq_len(nrow(parms)),
             .combine = c) %do% {
               
               print(i)
               y <- with(parms, 
                         fcl(foodweb = list_fw[[fw[i]]],
                             lambda = lambda[i],
                             L = rl[i],
                             rsrc = rsrc[i],
                             mu_base = mu[i],
                             mu_cnsm = mu[i],
                             delta = rep(delta[i], 2),
                             rho = rep(rho[i], 2)))
               
             }


# export ------------------------------------------------------------------

df_fcl <- parms %>% 
  mutate(fcl = y)

saveRDS(df_fcl, file = "data_fmt/sim_fcl_analytical.rds")
