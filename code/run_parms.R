
# library -----------------------------------------------------------------

rm(list = ls())
source("code/library.R")

# foodweb -----------------------------------------------------------------

## replicate number
n_fw <- 5

## food web setup
## 0.18 proportion of basal species - see Briand and Cohen 1987 Nature
## 0.11 connectance - see Dunne et al. 2002 PNAS
## 0.25 theta - see Johnson et al. 2014 PNAS
s <- 16
b <- floor(s * 0.18)
l <- floor(0.11 * s ^ 2)
theta <- 0.25

## generate food webs
set.seed(123)
list_fw <- replicate(n = n_fw,
                     foodweb(n_species = s,
                             n_basal = b,
                             l = l,
                             theta = theta),
                     simplify = FALSE)

saveRDS(list_fw, "data_fmt/parms_foodweb.rds")


# network -----------------------------------------------------------------

n_rep <- 100
sigma_src <- 1
sigma_lon <- 0.01

repeat {
  n_patch <- round(10^runif(n_rep, log(10, base = 10), log(50, base = 10)))
  p_branch <- 10^runif(n_rep, log(0.05, base = 10), log(0.95, base = 10))
  
  if(min(n_patch) < 15 & 
     max(n_patch) > 45 &
     min(p_branch) < 0.10 &
     max(p_branch) > 0.90) break    
}

parms <- data.table::data.table(id = seq_len(n_rep),
                                net_seed = seq_len(n_rep) + 1000,
                                n_patch,
                                p_branch,
                                sigma_src,
                                sigma_lon)

list_net <- lapply(seq_len(n_rep), function(i) {
  
  set.seed(with(parms, net_seed[i])) # for reproducibility
  
  out <- brnet(n_patch = n_patch[i],
               p_branch = p_branch[i],
               sd_disturb_source = sigma_src,
               sd_disturb_lon = sigma_lon)
  
  adj <- with(out, adjacency_matrix)
  death <- with(out, df_patch$disturbance)
  wa <- with(out, df_patch$n_patch_upstream)
  
  return(list(adj = adj, attr = data.table::data.table(death, wa)))
})

saveRDS(parms, "data_fmt/parms_net_value.rds")
saveRDS(list_net, "data_fmt/parms_net_obj.rds")