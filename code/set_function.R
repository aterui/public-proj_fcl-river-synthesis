
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
## @L total stream length

u_length <- function(lambda, L) {
  ## check input
  if (lambda < 0) 
    stop("invalid input: lambda must be >= 0")
  
  if (L <= 0)
    stop("invalid input: L must be > 0")
  
  ## z: (number of links/branches) minus 1
  ## pr_z: probability of b - 1 (= z) links/branches
  pois_max <- qpois(1 - 1e-10, lambda = lambda * L)
  v_z <- 0:pois_max
  pr_z <- dpois(v_z, lambda = lambda * L)
  
  ## pr_z_tr: truncate probabilities for z taking odd numbers
  ## (z odd = n links/branches even)
  even_id <- which(v_z %% 2 == 0)
  odd_id <- which(v_z %% 2 == 1)
  pr_even <- sum(pr_z[even_id])
  pr_z_tr <- pr_z / pr_even
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
      
      if (M > 514)
        stop(paste0("Stream magnitude M exceeds 514, which will return Inf in choose(2 * M, M); consider smaller values of lambda and/or L"))
      
      ## weighted values for the number of upstream links/branches
      ## - `2m - 2` is the number of upstream links/branches, ub
      ## - weighted by probability of drawing a link with m magnitude, w_ub
      m <- 1:M
      w_ub <- (2 * (1:M) - 2) * p_mag(1:M, M)
      
      ## expected value for the number of upstream links ul
      ub_hat <- sum(w_ub)
      
      ## expected value for upstream stream length
      ## - if poisson distributed,
      ## - the arrival time (distance to a given patch) will be a uniform dist
      ## - thus, expectation is (max - min) / 2 = l / 2
      u <- ub_hat * l_hat + 0.5 * l_hat
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



# metapopulation model ----------------------------------------------------

## for basal species
p_base <- function(lambda,
                   L,
                   qr = 1,
                   delta = 1,
                   rsrc = 1,
                   mu = c(1, 1),
                   rho = 0.5,
                   pg = 100) {
  
  ## n_patch: scalar, # habitat patches
  ## - qr: scalar, rate parameter for patch density (unit distance)
  ## - L: scalar, ecosystem size or total river length
  n_patch <- qr * L
  
  ## clnz: colonization rate
  ## - s: scalar, survival probability during migration
  ## - delta: vector, rate parameter(s) for species i's movement ability
  ## - rsrc: resource availability
  s <- 1 - exp(- delta * qr)
  pgle <- ifelse(s * pg < n_patch, 
                 yes = s * pg,
                 no = n_patch)
  clnz <- rsrc * pgle
  
  ## extn: extinction rate
  if (length(mu) == 1) {
    mu <- rep(mu, 2)
  } else {
    if (length(mu) != 2) stop("error in mu")
  }
  
  extn <- mu[1] + mu[2] * rho * u_length(lambda = lambda, L = L)
  
  p_hat <- 1 - (extn / clnz)
  p_hat <- ifelse(p_hat > 0, p_hat, 0)
  
  return(p_hat)
}

## for consumers
p_cnsm <- function(lambda,
                   L,
                   qr = 1,
                   delta = 1,
                   prey,
                   max_prey,
                   mu = c(1, 1, 1),
                   rho = 0.5,
                   pg = 100) {
  
  ## n_patch: scalar, # habitat patches
  ## - qr: scalar, rate parameter for patch density (unit distance)
  ## - L: scalar, ecosystem size or total river length
  n_patch <- qr * L
  
  ## clnz: colonization rate
  ## - s: scalar, survival probability during migration
  ## - delta: vector, rate parameter(s) for species i's movement ability
  ## - rsrc: resource availability
  s <- 1 - exp(- delta * qr)
  
  pgle <- ifelse(s * pg < n_patch, 
                 yes = s * pg,
                 no = n_patch)
  
  clnz <- prey * pgle
  
  ## extn: extinction rate
  if (length(mu) == 1) {
    mu <- rep(mu, 3)
  } else {
    if (length(mu) != 3) stop("error in mu")
  }
  
  extn <- mu[1] * 
    mu[2] * rho * u_length(lambda = lambda, L = L) +
    mu[3] * (1 - (prey / max_prey))
  
  p_hat <- 1 - (extn / clnz)
  p_hat <- ifelse(p_hat > 0, p_hat, 0)
  
  return(p_hat)
}

fcl <- function(foodweb,
                lambda,
                L,
                qr = 1,
                delta = c(1, 1),
                rsrc = 1,
                pg = c(100, 100),
                mu_base = 1,
                mu_cnsm = 1,
                rho = 0.5) {
  
  ## foodweb: matrix, consumer-resource matrix. produce with ppm()
  foodweb <- abs(foodweb)
  foodweb[lower.tri(foodweb)] <- 0
  
  ## get binary trophic position
  tp <- attr(foodweb, "tp")
  
  ## p_hat: vector, equilibrium occpancy
  ## max_prey: vector, maximum number of prey items for consumer j
  p_hat <- rep(-1, ncol(foodweb))
  max_prey <- colSums(foodweb)
  
  ## seqencial determination of equilibrium occupancies
  for (j in seq_len(ncol(foodweb))) {
    
    if (max_prey[j] == 0) {
      ## basal species
      p_hat[j] <- p_base(lambda = lambda,
                         L = L,
                         qr = qr,
                         delta = delta[1],
                         rsrc = rsrc,
                         mu = mu_base,
                         rho = rho,
                         pg = pg[1])
      
    } else {
      ## consumers
      
      ## index of prey species for consumer j
      index_prey <- which(foodweb[, j] == 1)
      
      ## mean-field prey richness
      prey <- sum(p_hat[index_prey])
      
      ## possible maximum of prey richness
      n_prey <- max_prey[j]
      
      p_hat[j] <- p_cnsm(lambda = lambda,
                         L = L,
                         qr = qr,
                         delta = delta[2],
                         prey = prey,
                         max_prey = n_prey,
                         mu = mu_cnsm,
                         rho = rho,
                         pg = pg[2])
      
    } # ifelse
  } # for j
  
  ## report fcl as the maximum binary trophic position in the landscape
  fcl <- ifelse(any(p_hat > 0), max(tp[p_hat > 0]), 0)
  attr(fcl, "p_hat") <- p_hat
  
  return(fcl)
}
