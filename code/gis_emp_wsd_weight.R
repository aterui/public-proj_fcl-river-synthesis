#' DESCRIPTION:
#' GIS analysis to obtain distance ratio d_w

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

# read data ---------------------------------------------------------------

sf_str <- readRDS("data_raw/sf_channel.rds") %>% 
  arrange(huid)

sf_fcl_point <- readRDS("data_raw/sf_fcl_site.rds") %>% 
  group_by(sid) %>% 
  slice(1) %>% # some sites (`sid`) has seasonal replicates
  ungroup()%>% 
  dplyr::select(huid, oid, sid) %>% 
  arrange(huid)

# weight calculation ------------------------------------------------------
# calculate distance ratio for each watershed

v_oid <- unique(sf_str$oid)

tictoc::tic()
df_ratio <- foreach(i = seq_along(v_oid),
                    .combine = bind_rows) %do% {
                      
                      print(i)
                      
                      p <- sf_fcl_point %>% 
                        filter(oid == v_oid[i])
                      
                      if (n_distinct(p$sid) == 1) {
                        ## if only one point, return zero
                        cout <- 0
                        
                      } else {
                        ## if more than one site, measure distance ratio
                        ## - randomly generate random points (n = n site)
                        ## - calculate median distance between points
                        ## - repeat n times
                        ## - take a mean of median distances
                        ## - compare with median observed distance between sites
                        x <- sf_str %>% 
                          filter(oid == v_oid[i])
                        
                        cout <- d_ratio(x = x,
                                        point = p,
                                        n_rand = 100,
                                        seed = 123)
                      }
                      
                      return(tibble(oid = v_oid[i],
                                    d_ratio = cout,
                                    n_site = n_distinct(p$sid)))
                    }
tictoc::toc()


# score -------------------------------------------------------------------

## export
saveRDS(df_ratio, "data_fmt/data_weight.rds")
