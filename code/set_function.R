
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
## arguments
## - m, magnitude of a given link
## - M, maximum magnitude of a network

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
## arguments
## - lambda, branching rate (per unit distance)
## - L, total stream length

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
## arguments
## - lambda, branching rate. branching probability = 1 - exp(-lambda)
## - L, total river length
## - h, habitat patch density (or rate) per unit distance
## - delta, species' dispersal capability
## - rsrc, resource availability for basal species
## - mu, disturbance rate: mu[1], base rate; mu[2], spatial rate
## - rho, sychrony probability of spatial disturbance
## - pg, propagule size

p_base <- function(lambda,
                   L,
                   h = 1,
                   delta = 1,
                   rsrc = 1,
                   mu = c(1, 1),
                   rho = 0.5,
                   pg = 100) {
  
  ## n_patch: scalar, # habitat patches
  n_patch <- h * L
  
  ## clnz: colonization rate
  ## - s: scalar, survival probability during migration
  s <- 1 - exp(- delta * h)
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
  
  ## equilibrium patch occupancy
  if (extn == 0 & clnz == 0) 
    stop("both colonization and extinction rates are zero; equilibrium undefined.")
  
  p_hat <- 1 - (extn / clnz)
  p_hat <- ifelse(p_hat > 0, p_hat, 0)
  
  return(p_hat)
}

## for consumers
## arguments
## - lambda, branching rate. branching probability = 1 - exp(-lambda)
## - L, total river length
## - h, habitat patch density (or rate) per unit distance
## - delta, species' dispersal capability
## - prey, expected prey richness (i.e., sum of equilibrium occupancies of prey)
## - max_prey, the number of prey species the consumer can prey on
## - mu, disturbance rate: mu[1] = base rate, mu[2] = spatial rate, mu[3] = prey
## - rho, synchronous extinction probability due to spatial disturbance
## - pg, propagule size

p_cnsm <- function(lambda,
                   L,
                   h = 1,
                   delta = 1,
                   prey,
                   max_prey,
                   mu = c(1, 1, 1),
                   rho = 0.5,
                   pg = 100) {
  
  ## n_patch: scalar, # habitat patches
  n_patch <- h * L
  
  ## clnz: colonization rate
  ## - s: scalar, survival probability during migration
  s <- 1 - exp(- delta * h)
  
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
  
  extn <- mu[1] +
    mu[2] * rho * u_length(lambda = lambda, L = L) +
    mu[3] * (1 - (prey / max_prey))
  
  ## equilibrium patch occupancy
  if (extn == 0 & clnz == 0) 
    stop("both colonization and extinction rates are zero; equilibrium undefined.")
  
  p_hat <- 1 - (extn / clnz)
  p_hat <- ifelse(p_hat > 0, p_hat, 0)
  
  return(p_hat)
}

## for food chain length
## arguments
## - foodweb, binary food web matrix from ppm()
## - lambda, branching rate. branching probability = 1 - exp(-lambda)
## - L, total river length
## - h, habitat patch density (or rate) per unit distance
## - delta, species' dispersal capability
## - rsrc, resource availability for basal species
## - pg, propagule size
## - mu_base, disturbance rate for basal species
## - mu_cnsm, disturbance rate for consumers
## - rho, synchronous extinction probability due to spatial disturbance

fcl <- function(foodweb,
                lambda,
                L,
                h = 1,
                delta = c(1, 1),
                rsrc = 1,
                pg = c(10, 10),
                mu_base = 1,
                mu_cnsm = 1,
                rho = c(0.5, 0.5)) {
  
  ## foodweb: matrix, consumer-resource matrix. produce with ppm()
  foodweb <- abs(foodweb)
  foodweb[lower.tri(foodweb)] <- 0
  
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
                         h = h,
                         delta = delta[1],
                         rsrc = rsrc,
                         mu = mu_base,
                         rho = rho[1],
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
                         h = h,
                         delta = delta[2],
                         prey = prey,
                         max_prey = n_prey,
                         mu = mu_cnsm,
                         rho = rho[2],
                         pg = pg[2])
      
    } # ifelse
  } # for j
  
  ## report fcl as the maximum binary trophic position in the landscape
  if (any(p_hat > 0)) {
    ## at least one species persist
    
    ## index of persistent species
    ## subset the foodweb by persistent species
    index_p <- which(p_hat > 0)
    sub_fw <- foodweb[index_p, index_p]
    tp <- rep(-1, ncol(sub_fw))
    
    ## index of basal species, number of basal, number of prey
    index_b <- which(colSums(sub_fw) == 0)
    n_b <- sum(colSums(sub_fw) == 0)
    n_p <- colSums(sub_fw)
    
    if (any(n_p > 0)) {
      ## tp = 1 if basal
      tp[index_b] <- 1  
      
      for (i in (n_b + 1):length(tp)) {
        ## update tp recursively if consumers
        tp[i] <-  (drop(tp %*% sub_fw)[i]) / n_p[i] + 1
      }
      
    } else {
      
      tp[index_b] <- 1
      
    }
    
    fcl <- max(tp)
    attr(fcl, "tp") <- tp
  } else {
    ## no species persist
    fcl <- 0
  }
  
  attr(fcl, "p_hat") <- p_hat
  
  return(fcl)
}


# numerical metapopulation model ------------------------------------------

## utility function
to_v <- function(x, n) {
  
  if (length(x) == 1) {
    
    v_x <- rep(x, n)
    
  } else {
    
    if (n != length(x))
      stop("incorrect input in one or more of the parameters")
    
    v_x <- x
  }
  
  return(v_x) 
}

## numerical solver with top down effects
p_fw <- function(foodweb,
                 h = 1,
                 delta = 1,
                 r = 1,
                 g = 1,
                 mu0 = 1,
                 mu_p = 1,
                 mu_c = 1,
                 mu_s = 1,
                 rho = 0.5,
                 L,
                 lambda,
                 x0 = 0.5,
                 n_timestep = 100,
                 interval = 0.01,
                 threshold = 1E-5) {
  
  n_species <- nrow(foodweb)
  Mc <- Mm <- abs(foodweb)
  
  Mc[upper.tri(Mc)] <- 0
  Mm[lower.tri(Mm)] <- 0
  
  s_prey <- rowSums(Mc)
  inv_s_prey <- ifelse(s_prey > 0,
                       yes = 1 / s_prey,
                       no = -1)
  
  id_b <- which(s_prey == 0)
  n_b <- length(id_b)
  
  id_c <- which(s_prey > 0)
  n_c <- length(id_c)
  
  ## colonization rate
  ## - propagule survival
  v_phi <- to_v(h, n = n_species)
  
  ## - resource availability
  r0 <- to_v(r, n = n_b)
  v_r <- c(r0, rep(0, n_c))
  
  ## - propagule
  n_patch <- h * L
  v_g <- to_v(g, n = n_species)
  v_phi_x_g <- v_phi * v_g
  
  id_phi_x_g <- which(v_phi_x_g < n_patch)
  id_n_patch <- which(v_phi_x_g >= n_patch)
  v_g[id_phi_x_g] <- v_phi_x_g[id_phi_x_g]
  v_g[id_n_patch] <- n_patch
  
  ## extinction rate
  ## - base rate
  v_mu0 <- to_v(mu0, n = n_species)
  
  ## - prey availability
  v_mu_p <- to_v(mu_p, n = n_species)
  v_mu_p[id_b] <- 0
  
  ## - predation influence
  v_mu_c <- to_v(mu_c, n = n_species)
  
  ## - spatial
  v_mu_s <- to_v(mu_s, n = n_species)
  v_rho <- to_v(rho, n = n_species)
  u <- u_length(lambda = lambda, L = L)
  
  ## derivative
  derivr <- function(t, x, parms) {
    with(parms, {
      
      dx <- phi * g * (Mc %*% x + r) * x * (1 - x) -
        (mu0 + 
           mu_p * (1 - (Mc %*% x) * inv_s_prey) + 
           mu_c * Mm %*% x + 
           mu_s * rho * u) * x
      
      list(dx)
    })
  }
  
  ## set parameters for ode()
  parms <- list(phi = v_phi,
                g = v_g,
                Mc = Mc,
                r = v_r,
                mu0 = v_mu0,
                mu_p = v_mu_p,
                mu_c = v_mu_c,
                mu_s = v_mu_s,
                inv_s_prey = inv_s_prey,
                u = u)
  
  x_init <- rep(x0, n_species)
  times <- seq(0, n_timestep, by = interval)
  
  ## define absorbing condition
  ## - root function
  rootfun <- function(t, x, parms) {
    return(x - threshold)
  }
  
  ## - extinction: triggered when "x - threshold = 0"
  eventfun <- function(t, x, parms) {
    x <- ifelse(x <= threshold, 0, x)
    return(x)
  }
  
  # run ode solver
  cout <- deSolve::ode(y = x_init,
                       times = times,
                       func = derivr,
                       parms = parms,
                       events = list(func = eventfun,
                                     root = TRUE),
                       rootfun = rootfun)
  
  return(cout)
}
