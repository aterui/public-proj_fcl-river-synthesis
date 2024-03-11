
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


# get a value of within-branch combinations -------------------------------

## expected number of within-branch combos
get_mu_phi <- function(n, p) {
  
  ## calculate summed probability of getting odd numbers from binom dist
  pr <- dbinom(0:n, size = n, prob = p)
  odd <- which(0:n %% 2 == 1)
  pr_odd <- sum(pr[odd])
  
  phi_b <- sapply(X = 1:n, function(b) {
    
    # prob # branches
    ## b: possible values of branch numbers
    ## exclude even numbers
    if (b %% 2 == 0) {
      pr_b <- 0
    } else {
      ## use the truncated probability distribution so that sum(pr_b) = 1.0
      pr_b <- dbinom(x = b, size = n, prob = p) / pr_odd
    }
    
    # expected # of within-branch combo, conditional on B = b
    if (b == 1) {
      ## if only one branch present
      ## E(phi) is "pick 2 out of N"
      phi0 <- choose(n, 2)
      
    } else {
      
      if (b == n) {
        ## if #branches = n patch, no within-branch combo
        phi0 <- 0
      } else {
        ## if 1 < b < N
        ## use beta-binomial distribution
        w <- sapply(X = 0:(n - b),
                    function(k) {
                      
                      ## with at least one non-branching node
                      w0 <- choose(k + 1, 2) *
                        choose(n - b, k) *
                        beta(k + 1, n - k - 1) *
                        (b - 1)
                      
                      return(w0)
                    }) # sapply()
        
        phi0 <- sum(w)
      }
      
    } # ifelse for b
    
    return(pr_b * b * phi0)
  }) #sapply
  
  return(sum(phi_b))
}


## simulation-based realized number of within-branch combos
get_phi <- function(n, p) {
  
  ## produce # of branches
  b <- rbinom(n = 1, size = n, prob = p)
  while(b %% 2 == 0) {
    b <- rbinom(n = 1, size = n, prob = p)
  }  
  
  #print(paste("b =", b))
  
  ## produce # nodes for each branch
  repeat {
    repeat {
      k <- rgeom(n = b, prob = p) + 1
      if (sum(k) >= n) break
    }
    if (sum(k) == n) break
  }
  
  #print(paste("k =", k))
  phi <- sum(sapply(k, function(x) choose(x, 2)))
  
  return(phi)
}


# flow-connected fraction -------------------------------------------------

## get probability of drawing a link with m magnitude from a network with M magnitude
## @m m magnitude
## @M maximum magnitude of a network
p_mag <- function(m, M) {
  
  if (any(m > M)) stop("M must be >= m")
  if (any(m <= 0) | M <= 0) stop("m and M must be positive integer")
  
  if (M <= 500) {
    ## if M <= 500, exact calculation
    pr_m <- (2 * m - 1)^(-1) * 
      choose(2 * m, m) * 
      choose(2 * (M - m), M - m) / choose(2 * M, M) 
  } else {
    ## if M > 500, approximation as if M -> infinity
    pr_m <- 2 ^ (-(2 * m - 1)) * 
      (2 * m - 1) ^ (-1) * 
      choose(2 * m - 1, m)
  }
  
  return(pr_m)
}

## get expected upstream length
## @lambda branching rate (per unit distance)
## @q patch rate parameter (per unit distance)
## @L total stream length

u_length <- function(lambda, q, L) {
  ## z: (number of links/branches) minus 1
  ## pr_z: probability of b - 1 (= z) links/branches
  pois_max <- qpois(1 - 1e-10, lambda = lambda * L)
  v_z <- 0:pois_max
  pr_z <- dpois(v_z, lambda = lambda * L)
  
  ## pr_z_tr: truncate probabilities for z taking odd numbers
  ## (= n links/branches becomes even)
  odd_id <- which(v_z %% 2 == 1)
  pr_odd <- sum(pr_z[odd_id])
  pr_z_tr <- pr_z / pr_odd
  pr_z_tr[odd_id] <- 0
  
  ## expected total length of upstream links/branches, conditional on z
  u_z <- sapply(v_z, function(z) {
    if (z %% 2 == 0) {
      ## when z is even = n branches is odd
      
      ## b: number of links/branches
      ## l_hat: expected length of a link/branch, conditional on b
      ## - l_hat derived from a Beta distribution Beta(1, z)
      b <- (z + 1)
      l_hat <- L / b
      
      ## maximum magnitude in a network
      M <- 0.5 * (b + 1)
      
      ## weighted values for the number of upstream links/branches
      ## - `2m - 2` is the number of upstream links/branches, ub
      ## - weighted by probability of drawing a link with m magnitude, w_ub
      m <- 1:M
      w_ub <- (2 * (1:M) - 2) * p_mag(1:M, M)
      
      ## expected value for the number of upstream links ul
      ub_hat <- sum(w_ub)
      
      ## expected value for upstream stream length
      ## - add within-branch upstream portion (l_hat * (1 + q * l_hat) ^ (-1))
      ## - q * l_hat = patch intensity (rate) parameter given the link length
      ## - (1 + q * l_hat) ^ (-1) = expected fraction of link length to a focal patch
      ## - this fraction derived from a beta distribution Beta(1, q * l_hat)
      u <- ub_hat * l_hat + l_hat * (1 + q * l_hat) ^ (-1)
    } else {
      ## when z is odd = n branches is even
      u <- -1
    }
    
    return(u)
  })
  
  ## sum over z to get an expected value of upstream link length
  u_hat <- sum(pr_z_tr * u_z)
  
  return(u_hat)
}
