
pacman::p_load(tidyverse)

# split a polygon into sub-polygons with equal area -----------------------

## function inspired by the following page:
## - https://gis.stackexchange.com/questions/375345/dividing-polygon-into-parts-which-have-equal-area-using-r
split_poly <- function(x, n_div, n_point = 10000, seed) {
  
  if (!missing(seed)) set.seed(seed)
  
  df_point <- sf::st_sample(x, size = n_point) %>% 
    sf::st_coordinates() %>%
    tidyr::as_tibble() %>%
    setNames(c("lon", "lat"))
  
  kcl <- kmeans(df_point, centers = n_div)
  
  x_voronoi <- dismo::voronoi(kcl$centers, ext = x)
  terra::crs(x_voronoi) <- terra::crs(x)
  
  cout <- x_voronoi %>% 
    sf::st_as_sf() %>% 
    sf::st_intersection(x) %>%
    sf::st_make_valid() %>% 
    dplyr::mutate(area = units::set_units(st_area(.), "km^2"))
  
  return(cout)
}


# compare median distance with random points ------------------------------

d_ratio <- function(x, point, n_rand = 100, seed) {
  
  ## observed median distance
  m_dist <- sf::st_distance(point)
  d_obs <- units::set_units(median(m_dist[upper.tri(m_dist)]),
                            "km") %>% 
    units::drop_units()
  
  ## median distance random points
  if (!missing(seed)) set.seed(seed)
  
  d_rand <- replicate(n_rand,
                      {
                        p_rand <- suppressMessages(sf::st_sample(x, size = nrow(point)))
                        p_rand <- p_rand[!sf::st_is_empty(p_rand)] %>% 
                          sf::st_as_sf() %>% 
                          sf::st_cast("POINT")
                        
                        while(nrow(p_rand) != nrow(point)) {
                          p_rand <- suppressMessages(sf::st_sample(x, size = nrow(point)))
                          p_rand <- p_rand[!sf::st_is_empty(p_rand)]
                        } 
                          
                        m_dist_rand <- sf::st_distance(p_rand)
                        d <- units::set_units(median(m_dist_rand[upper.tri(m_dist_rand)]),
                                              "km") %>% 
                          units::drop_units()
                        
                        return(d)
                      },
                      simplify = TRUE)
  
  return(d_obs / mean(d_rand))
}



