
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")


# read data ---------------------------------------------------------------

sf_wsd <- readRDS("data_fmt/wgs84_wsd_sub.rds")

sf_str <- readRDS("data_fmt/wgs84_str_sub.rds") %>% 
  bind_rows() %>% 
  mutate(wid = str_pad(wid, width = 5, pad = "0"),
         uid = paste0(huid, wid)) %>% 
  dplyr::select(uid) %>% 
  arrange(uid)

sf_fcl_point <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  group_by(sid) %>% 
  slice(1) %>% # some sites (`sid`) has seasonal replicates
  ungroup() %>% 
  mutate(wid = str_pad(wid, width = 5, pad = "0"),
         uid = paste0(huid, wid)) %>% 
  dplyr::select(uid, sid, study_id) %>% 
  arrange(uid)


# weight calculation ------------------------------------------------------
# calculate distance ratio for each watershed

v_uid <- unique(sf_str$uid)

df_ratio <- foreach(i = 1:length(v_uid),
                    .combine = bind_rows) %do% {
                      
                      print(i)
                      
                      p <- sf_fcl_point %>% 
                        filter(uid == v_uid[i])
                      
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
                          filter(uid == v_uid[i])
                        
                        cout <- d_ratio(x = x,
                                        point = p,
                                        n_rand = 100,
                                        seed = 123)
                      }
                      
                      return(tibble(uid = v_uid[i],
                                    d_ratio = cout,
                                    n_site = n_distinct(p$sid)))
                    }


# score -------------------------------------------------------------------

df_weight <- left_join(sf_wsd, df_ratio) %>% 
  mutate(area = units::set_units(st_area(.), "100km^2")) %>% 
  mutate(score = n_site * (d_ratio + 1))

## export
saveRDS(df_weight, "data_fmt/data_weight.rds")
