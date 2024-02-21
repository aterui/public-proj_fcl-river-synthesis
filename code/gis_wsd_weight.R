
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")


# read data ---------------------------------------------------------------

df_wsd <- readRDS("data_fmt/wgs84_wsd_sub.rds") %>% 
  mutate(wid = str_pad(wid, width = 5, pad = "0"),
         uid = paste0(huid, wid)) %>% 
  st_cast("LINESTRING") %>% 
  mutate(length = units::set_units(st_length(.), "km")) %>% 
  group_by(uid) %>% 
  summarize(perimeter = round(sum(length), 2)) %>% 
  as_tibble() %>% 
  dplyr::select(-geometry)

sf_fcl_point <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  group_by(sid) %>% 
  slice(1) %>% # some sites (`sid`) has seasonal replicates
  ungroup() %>% 
  mutate(wid = str_pad(wid, width = 5, pad = "0"),
         uid = paste0(huid, wid)) %>% 
  dplyr::select(uid, sid, study_id) %>% 
  arrange(uid)


# weight calculation based on perimeter -----------------------------------
## weight based on the # sites and spatial coverage
## spatial coverage was estimated as follows:
## - if only one site, take the length of one pixel (longitude-wise)
## - if two sites, take the distance between the two
## - if > two sites, make a convex hull and measure the perimeter
## - weight score = (spatial coverage / full perimeter of the watershed) * n_site

## get a vector of unique uid
uid <- sf_fcl_point %>% 
  pull(uid) %>% 
  unique() %>% 
  sort()

df_cov <- foreach(x = uid, .combine = bind_rows) %do% {
  
  ## calculate "spatial coverage" for weight
  ## - choose one watershed
  sf_i <- sf_fcl_point %>% 
    filter(uid == x)
  
  if (n_distinct(sf_i$sid) == 1) {
    
    ## - if only one site
    ## - get coordinates
    xy <- st_coordinates(sf_i)
    arc3sec <- (1 / 3600) * 3 * 0.5
    
    ## - measure length of diagonal line of a 3-arc second pixel
    ## - 3-arc second is a resolution of flow direction/upstream watershed area
    y <- st_linestring(rbind(c(xy[,1] - arc3sec, xy[,2]),
                             c(xy[,1] + arc3sec, xy[,2]))) %>% 
      st_sfc() %>% 
      st_set_crs(4326) %>% 
      st_length() %>% 
      units::set_units("km")
    
  } else {
    
    if (n_distinct(sf_i$sid) == 2) {
      
      ## - if two sites
      ## - measure length of a line connecting two points
      y <- st_union(sf_i) %>% 
        st_cast("LINESTRING") %>% 
        st_length() %>% 
        units::set_units("km")
      
    } else {
      if (n_distinct(sf_i$sid) > 2) {
        
        ## - if > two sites
        ## - measure perimeter of a convex hull connecting sampling sites
        y <- st_union(sf_i) %>% 
          st_convex_hull() %>% 
          st_cast("LINESTRING") %>% 
          st_length() %>% 
          units::set_units("km")
        
      } else {
        
        y <- NA
        
      } # if     
      
    } # if
  } # if
  
  return(tibble(uid = x,
                coverage = round(y, 2),
                n_site = n_distinct(sf_i$sid)))
}

## join watershed polygons
df_weight <- df_cov %>% 
  left_join(df_wsd, by = 'uid') %>% 
  mutate(score = as.numeric((coverage / perimeter) * n_site),
         w = as.numeric(score / max(score)))

## export
saveRDS(df_weight, "data_fmt/data_weight.rds")
