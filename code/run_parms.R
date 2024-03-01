
# library -----------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

# foodweb -----------------------------------------------------------------

## replicate number
n_fw <- 20

## food web setup
## 0.18 proportion of basal species - see Briand and Cohen 1987 Nature
## 0.11 connectance - see Dunne et al. 2002 PNAS
## 0.25 theta - degree of omnivory; see Johnson et al. 2014 PNAS
s <- 16
b <- round(s * 0.18)
l <- round(0.11 * s ^ 2)
theta <- c(0.25, 1)

## generate food webs
set.seed(123)

list_fw <- foreach(x = theta,
                   .combine = c) %do% {
                     replicate(n = n_fw, {
                       ## generate food webs
                       fw <- foodweb(n_species = s,
                                     n_basal = b,
                                     l = l,
                                     theta = x,
                                     convert = list(min = 0.1, max = 0.5))
                       
                       while(any(rowSums(fw[1:b, (b + 1):s]) == 0)) {
                         ## make sure that all the basal species have 
                         ## at least one consumer
                         fw <- foodweb(n_species = s,
                                       n_basal = b,
                                       l = l,
                                       theta = x,
                                       convert = list(min = 0.1, max = 0.5))
                       }
                       
                       attr(fw, "theta") <- x
                       return(fw)
                     },
                     simplify = FALSE) 
                   }

saveRDS(list_fw, "data_fmt/parms_foodweb.rds")

# network -----------------------------------------------------------------

## set network parameters
## n_rep = number of network replicates
## sigma_src = SD for disturbance values at headwaters
## sigma_lon = longitudinal SD for disturbance values
n_rep <- 50
sigma_src <- 1
sigma_lon <- 0.01

## np_xx = min or max values for n_patch
## pb_xx = min or max values for p_branch
np_min <- 10
np_max <- 50
pb_min <- 0.05
pb_max <- 0.95

## sample parameter values
repeat {
  n_patch <- round(runif(n_rep, np_min, np_max))
  p_branch <- runif(n_rep, pb_min, pb_max)
  
  if(min(n_patch) < np_min + 5 & 
     max(n_patch) > np_max - 5 &
     min(p_branch) < pb_min + 0.05 &
     max(p_branch) > pb_max - 0.05) break    
}

parms <- data.table::data.table(id = seq_len(n_rep),
                                net_seed = seq_len(n_rep) + 1000,
                                n_patch,
                                p_branch,
                                sigma_src,
                                sigma_lon,
                                np_min,
                                np_max,
                                pb_min,
                                pb_max)

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
