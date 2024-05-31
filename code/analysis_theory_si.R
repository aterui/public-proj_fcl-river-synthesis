
# setup -------------------------------------------------------------------

rm(list = ls())
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
            l = round(S ^ 2 * 0.11),
            theta = 0.25)},
  simplify = FALSE)

## other parameters
rsrc <- c(0.25, 2.5)
mu0 <- c(0.25, 2.5)
rho <- c(0, 0.5)
mu_p <- c(0, 1)
mu_c <- c(0, 1)
g <- c(10, 100)
delta0 <- 1
fw <- seq_len(length(list_fw))

## parms data frame
parms <- expand.grid(rl = seq(1, 100, length = 10),
                     lambda = seq(0.1, 1, length = 10),
                     rsrc = rsrc,
                     mu0 = mu0,
                     mu_p = mu_p,
                     mu_c = mu_c,
                     rho = rho,
                     g = g,
                     delta0 = delta0,
                     fw = fw)


# prediction --------------------------------------------------------------

parms <- filter(parms, fw == 1)

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
y <- foreach(i = seq_len(nrow(parms)),
             .combine = c,
             .options.snow = opts) %dopar% {
               
               y <- with(parms, {
                 
                 v_tp <- attr(list_fw[[parms$fw[i]]], "tp")
                 
                 if (mu_c[i] == 0) {
                   y0 <- rpom::fcl(foodweb = list_fw[[fw[i]]],
                                   lambda = lambda[i],
                                   size = rl[i],
                                   rsrc = rsrc[i],
                                   mu0 = mu0[i],
                                   mu_p = mu_p[i],
                                   delta = delta0[i] * sqrt(v_tp),
                                   g = g[i],
                                   rho = rho[i],
                                   weight = TRUE)
                 } else {
                   y0 <- rpom::nfcl(foodweb = list_fw[[fw[i]]],
                                    lambda = lambda[i],
                                    size = rl[i],
                                    rsrc = rsrc[i],
                                    mu0 = mu0[i],
                                    mu_p = mu_p[i],
                                    mu_c = mu_c[i],
                                    delta = delta0[i] * sqrt(v_tp),
                                    g = g[i],
                                    rho = rho[i],
                                    weight = TRUE)
                 }
                 
                 return(y0)
               })
               
               return(y)
             }
tictoc::toc()

stopCluster(cl)

# export ------------------------------------------------------------------

df_fcl <- parms %>% 
  mutate(fcl = y)

saveRDS(df_fcl, file = "data_fmt/sim_fcl_si.rds")
