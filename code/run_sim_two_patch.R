
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/run_parms.R")

cl <- makeCluster(10)
registerDoSNOW(cl)

# parameter setup ---------------------------------------------------------

## food web replicate
list_fw <- readRDS("data_fmt/parms_foodweb.rds")

## number of replications within parameter combo
n_rep <- 10

## parameter combinations
# xi = search interval for findr()
parms <- expand.grid(n_timestep = 200,
                     n_species = ncol(list_fw[[1]]),
                     n_patch = 2,
                     phi = seq(0, 0.01, length = 10),
                     m = 0,
                     rate = 0.1,
                     e = 0.9,
                     k = seq(1, 50, length = 10),
                     s = 0.5,
                     threshold = c(1E-4, 1E-3),
                     foodweb = seq_len(length(list_fw)),
                     xi = 0.05) %>%
  mutate(theta = sapply(list_fw, function(x) attr(x, "theta"))[foodweb],
         i = row_number()) %>% 
  as_tibble()

pb <- txtProgressBar(max = nrow(parms), style = 3)
fun_progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = fun_progress)

tictoc::tic()
df_fcl <- foreach(x = iterators::iter(parms, by = "row"),
                  .combine = rbind,
                  .packages = c("mcbrnet",
                                "foreach"),
                  .options.snow = opts) %dopar% {
                    
                    # index for seeding
                    i <- with(x, i)
                    
                    # define food web
                    # a = interaction matrix
                    # e = disturbance mortality rate
                    # k = carrying capacity
                    a <- with(x, list_fw[[foodweb]])
                    xi <- with(x, xi)
                    k <- with(x, k)
                    adj <- with(x, matrix(1, n_patch, n_patch))
                    diag(adj) <- 0
                    
                    dt_j <- foreach(j = seq_len(n_rep),
                                    .combine = rbind) %do% {
                                      
                                      # m_r = species x patch matrix of r
                                      r <- findr(alpha = a,
                                                 k0 = k,
                                                 interval = xi)[, 1]
                                      
                                      v_r <- with(x, rep(r, n_patch))
                                      m_r <- with(x,
                                                  matrix(v_r,
                                                         nrow = n_species,
                                                         ncol = n_patch))
                                      
                                      seed <- i * 1000 + j
                                      set.seed(seed)
                                      n <- with(x,
                                                sglv(n_species = n_species,
                                                     n_patch = n_patch,
                                                     n_timestep = n_timestep,
                                                     r = m_r,
                                                     alpha = a,
                                                     dispersal = list(adj = adj,
                                                                      phi = phi,
                                                                      m = m),
                                                     disturb = list(int = -log(1 - e),
                                                                    rate = rate,
                                                                    s = s),
                                                     threshold = threshold,
                                                     n0 = list(min = threshold,
                                                               max = 1))
                                      )
                                      
                                      fcl <- with(x, foodchain(n,
                                                               n_species = n_species,
                                                               n_patch = n_patch,
                                                               alpha = a))
                                      
                                      return(data.table::data.table(rep = j,
                                                                    sglv_seed = seed,
                                                                    fcl = fcl))
                                    }
                    
                    out <- data.table::data.table(x, dt_j)
                    return(out)
                  }
tictoc::toc()

stopCluster(cl); gc(); gc()


# export ------------------------------------------------------------------

saveRDS(df_fcl, "data_fmt/sim_fcl_two_patch.rds")
