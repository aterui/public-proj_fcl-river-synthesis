
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/run_parms.R")

cl <- makeCluster(10)
registerDoSNOW(cl)

# parameter setup ---------------------------------------------------------

## food web replicate
list_fw <- readRDS("data_fmt/parms_foodweb.rds")

## network replicate
list_net <- readRDS("data_fmt/parms_net_obj.rds")
parms_net <- readRDS("data_fmt/parms_net_value.rds")

## number of replications within parameter combo
n_rep <- length(list_net)

## parameter combinations
# z = from Finlay 2011, Ecosphere
# xi = search interval for findr()
parms <- expand.grid(n_timestep = 200,
                     n_species = ncol(list_fw[[1]]),
                     phi = 1E-3,
                     rate = c(0.01, 0.1),
                     s = 0.25,
                     threshold = 1E-3,
                     k_base = 1,
                     z = 0.54,
                     foodweb = seq_len(length(list_fw)),
                     xi = 0.025) %>%
  mutate(m = phi * 0.1,
         theta = sapply(list_fw, function(x) attr(x, "theta"))[foodweb],
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
                    # k_base = carrying capacity at stream sources
                    # z = scaling exponent for K
                    a <- with(x, list_fw[[foodweb]])
                    k_base <- with(x, k_base)
                    z <- with(x, z)
                    xi <- with(x, xi)
                    
                    dt_j <- foreach(j = seq_len(n_rep),
                                    .combine = rbind) %do% {
                                      
                                      # network list
                                      graph <- list_net[[j]]
                                      
                                      # e = disturbance mortality rate
                                      # k = carrying capacity
                                      e <- with(graph, -log(1 - attr$death))
                                      k <- with(graph, k_base * attr$wa^z)
                                      
                                      # m_r = species x patch matrix of r
                                      m_r <- sapply(k, function(u) findr(alpha = a,
                                                                         k0 = u,
                                                                         interval = xi)[, 1])
                                      
                                      seed <- i * 1000 + j
                                      set.seed(seed)
                                      n <- with(x,
                                                sglv(n_species = n_species,
                                                     n_patch = with(graph, ncol(adj)),
                                                     n_timestep = n_timestep,
                                                     r = m_r,
                                                     alpha = a,
                                                     dispersal = list(adj = with(graph, adj),
                                                                      phi = phi,
                                                                      m = m),
                                                     disturb = list(int = e,
                                                                    rate = rate,
                                                                    s = s),
                                                     threshold = threshold)
                                      )
                                      
                                      fcl <- with(x, foodchain(n,
                                                               n_species = n_species,
                                                               n_patch = with(graph, ncol(adj)),
                                                               alpha = a))
                                      
                                      return(data.table::data.table(rep = j,
                                                                    sglv_seed = seed,
                                                                    fcl = fcl, 
                                                                    parms_net[j,]))
                                    }
                    
                    out <- data.table::data.table(x, dt_j)
                    return(out)
                  }
tictoc::toc()

stopCluster(cl); gc(); gc()


# export ------------------------------------------------------------------

saveRDS(df_fcl, "data_fmt/sim_fcl_main.rds")
