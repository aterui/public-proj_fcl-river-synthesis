
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
ncore <- 6

cl <- makeCluster(ncore - 1)
registerDoSNOW(cl)

# parameter setup ---------------------------------------------------------

## food web replicate
list_fw <- readRDS("data_fmt/parms_foodweb.rds")

## number of replications within parameter combo
n_rep <- 10

## parameter combinations
parms <- expand.grid(n_timestep = 200,
                     n_species = ncol(list_fw[[1]]),
                     p_branch = seq(0.1, 0.9, by = 0.2),
                     n_patch = seq(10, 50, by = 10),
                     phi = 1E-3,
                     m = 1E-4,
                     rate = 0.1,
                     s = 0.25,
                     threshold = 1E-3,
                     k_base = 0.1,
                     foodweb = seq_len(length(list_fw))) %>% 
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
                    
                    # produce random river networks
                    adj <- with(x, brnet(n_patch = n_patch,
                                         p_branch = p_branch))
                    
                    # e = disturbance mortality rate
                    # k = carrying capacity for basal
                    # a = interaction matrix
                    e <- with(adj, -log(1 - df_patch$disturbance))
                    k <- with(adj, x$k_base * sqrt(df_patch$n_patch_upstream))
                    a <- with(x, list_fw[[foodweb]])
                    
                    # m_r = species x patch matrix of r
                    m_r <- sapply(k, function(u) findr(alpha = a,
                                                       k0 = u,
                                                       interval = 0.05)[, 1])
                    
                    dt_j <- foreach(j = seq_len(n_rep),
                                    .combine = rbind) %do% {
                                      n <- with(x,
                                                sglv(n_species = n_species,
                                                     n_patch = n_patch,
                                                     n_timestep = n_timestep,
                                                     r = m_r,
                                                     alpha = a,
                                                     dispersal = list(adj = with(adj, adjacency_matrix),
                                                                      phi = phi,
                                                                      m = m),
                                                     disturb = list(int = e,
                                                                    rate = rate,
                                                                    s = s),
                                                     threshold = threshold)
                                      )
                                      
                                      fcl <- with(x, foodchain(n,
                                                               n_species = n_species,
                                                               n_patch = n_patch,
                                                               alpha = a))
                                      
                                      return(data.table::data.table(rep = j, fcl = fcl))
                                    }
                    
                    out <- data.table::data.table(x, dt_j)
                    return(out)
                  }
tictoc::toc()

stopCluster(cl)


# export ------------------------------------------------------------------

saveRDS(df_fcl, "data_fmt/sim_fcl_main.rds")
