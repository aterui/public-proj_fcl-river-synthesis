
pacman::p_load(tidyverse,
               foreach,
               deSolve)


# extra prey function -----------------------------------------------------

extra_prey <- function(A, j, i0, tp, theta, kappa, cannibal = FALSE) {
  # A:      adjacency matrix defining trophic interactions
  # j:      consumer's index (j - 1 is the number of possible prey)
  # i0:     index of the first prey chosen (1 =< i0 < j)
  # tp:     vector of trophic positions for possible prey nodes
  #         (length = S, j - 1 elements are non-NA)
  # kappa:  number of extra prey nodes in addition to the first prey i0
  # theta:  scale parameter for an exponential decay of prey preference
  
  # verify inputs
  if (length(tp) != ncol(A)) stop(paste("'tp' length seems incorrect;",
                                        "the length must be",
                                        ncol(A)))
  
  if (!(i0 < j && 1 <= i0)) stop("i0 must be a non-zero intger smaller than j")
  
  if (kappa >= j) stop("kappa must be an integer smaller than j")
  
  # probability of consumer j picks prey i
  lambda <- 1 / theta
  tp[j] <- tp[i0] + 1
  p_ij <- exp(- lambda * abs(tp[i0] - tp))
  
  # pick kappa prey species out of j - 1 nodes (without cannibalism) or j nodes
  if (cannibal) {
    i_index <- 1:j
    p_ij <- p_ij[1:j]
  } else {
    i_index <- 1:(j - 1)
    p_ij <- p_ij[1:(j - 1)]
  }
  
  # exclude i0 from resample() because i0 is the first pick
  # pick 'kappa' samples from index[-i0]
  # weighted by p_ij
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  if (length(i_index) > 1) {
    i_pick <- resample(i_index[-i0],
                       size = kappa,
                       prob = p_ij[-i0])
    
    i <- c(i0, i_pick)
  } else {
    # when j = 2 with no cannibalism
    # no choice but i0 available as prey
    i <- i0  
  }
  
  A[i, j] <- 1
  
  return(A)
}


# preferential prey model function ----------------------------------------

ppm <- function(n_species,
                n_basal,
                l,
                theta,
                cannibal = F) {
  # n_species: number of species in the pool
  # n_basal: number of basal species
  # l: expected number of links
  # theta: scale parameter
  
  # verify inputs
  if (n_basal < 1) stop("At least one basal species needed to construct a food web")
  if (n_basal >= n_species) stop("n_basal must be smaller than n_species")
  
  # declare new parameter
  n_c <- n_species - n_basal
  
  # tp: trophic position
  # assign trophic position for basal species
  # assign -1 for consumers as initial values (to be updated)
  tp <- rep(-1, n_species)
  tp[seq_len(n_basal)] <- 1
  
  # beta parameter for beta distribution
  # determined so that E(L) = L
  if (l < n_c) stop("l must be at least equal to the number of consumers (n_species - n_basal)")
  
  if (cannibal) {
    max_l <- sum((n_basal + 1):n_species)
    if (l > max_l) stop(paste("maximum l is", max_l))
    if (n_species + n_basal <= 1) stop("n_species + n_basal must be greater than 1")
    b <- ((n_species + n_basal - 1) * n_c) / (2 * (l - n_c)) - 1
  } else {
    max_l <- sum((n_basal + 1):n_species - 1)
    if (l > max_l) stop(paste("maximum l is", max_l))
    if (n_species + n_basal <= 3) stop("n_species + n_basal must be greater than 3")
    b <- ((n_species + n_basal - 3) * n_c) / (2 * (l - n_c)) - 1
  }
  
  # vector for link proportions for all consumers
  v_xi <- rbeta(n_c,
                shape1 = 1,
                shape2 = b)
  
  # vector for initial prey choice for all consumers
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  v_i0 <- c(rep(-1, n_basal),
            sapply(X = (n_basal + 1):n_species,
                   FUN = function(i0) resample(seq_len(i0 - 1),
                                               size = 1)))
  
  # realized number of prey nodes for all consumers
  # kappa follows a beta-binomial distribution
  # size parameter is the possible maximum number of prey nodes present
  if (cannibal) {
    v_kappa <- c(rep(-1, n_basal),
                 rbinom(n = n_c,
                        size = n_basal:(n_species - 1),
                        prob = v_xi))
  } else {
    v_kappa <- c(rep(-1, n_basal),
                 rbinom(n = n_c,
                        size = (n_basal - 1):(n_species - 2),
                        prob = v_xi))
  }
  
  # A: S x S interaction binary matrix
  # initialized with all zero
  # update the initial consumer's first prey
  A <- matrix(0, n_species, n_species)
  A[v_i0[n_basal + 1], n_basal + 1] <- 1
  
  for (j in (n_basal + 1):n_species) {
    A <- extra_prey(A = A,
                    j = j,
                    i0 = v_i0[j],
                    tp = tp,
                    kappa = v_kappa[j],
                    theta = theta,
                    cannibal = cannibal)
    
    v_n_prey <- colSums(A)
    v_n_prey[v_n_prey == 0] <- 1
    
    tp_new <- (tp %*% A) / v_n_prey + 1
    tp[j] <- tp_new[j]
  }
  
  return(list(A = A, tp = tp))
}


# convert A ---------------------------------------------------------------

to_alpha <- function(A,
                     attack = list(min = 0,
                                   max = 1),
                     convert = list(min = 0,
                                    max = 1),
                     mortal = list(min = 0,
                                   max = 1)) {
  
  # identify basal species
  uA <- A
  basal <- which(colSums(A) == 0)
  
  # scale by # of resources
  suA <- t(t(A[, -basal]) / colSums(A)[-basal])
  
  # generate random parameters
  ## matrix for attack rates
  a <- with(attack, matrix(runif(ncol(A)^2, min = min, max = max),
                           nrow = nrow(A),
                           ncol = ncol(A)))
  
  ## matrix for conversion efficiency
  b <- with(convert, matrix(runif(ncol(A)^2, min = min, max = max),
                            nrow = nrow(A),
                            ncol = ncol(A)))
  
  ## vector for intraspecific competition
  m <- with(mortal, runif(ncol(A), min = min, max = max))
  
  # interaction matrix
  uA[, -basal] <- a[, -basal] * suA
  alpha <- b * uA - t(uA)
  
  # intraspecific interaction
  diag(alpha) <- - m
  
  return(alpha)
}


# deSolve wrapper ---------------------------------------------------------

sglv <- function(n_species,
                 n_patch,
                 n_timestep = 100,
                 interval = 0.01,
                 r,
                 alpha,
                 phi,
                 connectivity,
                 disturb = list(int = 1,
                                rate = (1 / n_timestep) * 10,
                                s = interval * 10),
                 n0 = list(min = 0,
                           max = 1),
                 threshold = 1E-4,
                 ...) {
  
  # verify inputs -----------------------------------------------------------
  
  # r input
  if (length(r) != 1 && length(r) != n_species)
    stop("r must have length one or n_species")
  
  # alpha input
  if (!all(dim(alpha) == n_species))
    stop("alpha's dimensions must be n_species x n_species")
  
  # connectivity input
  if (!all(dim(connectivity) == n_patch))
    stop("connectivity's dimensions must be n_patch x n_patch")
  
  # disturbance intensity
  if (length(disturb$int) != 1 && length(disturb$int) != n_patch)
    stop("disturbance intensity (int) must have length one or n_patch")
    
  # disturbance setup -------------------------------------------------------
  
  # n possible disturbance points
  # + 100 to ensure cumsum(x) > n_timestep
  # remove points > n_timestep
  psi <- with(disturb,
              cumsum(rexp(n = ceiling((1 / rate) + 100), rate)))
  psi <- psi[psi <= n_timestep]
  
  # pseudo time steps
  time <- seq(0, n_timestep, by = disturb$s)
  t_psi <- sapply(psi, FUN = function(x) which.min(abs(time - x)))
  
  # create dummy time-series for signals
  signal <- data.frame(time = time,
                       psi = rep(0, length(time)))
  signal$psi[t_psi] <- 1
  
  # linear interpolation function
  input <- approxfun(signal, rule = 2)
  
  # n_patch x n_species disturbance intensity matrix
  ## fun_to_m() returns n_species x n_patch matrix
  ## transpose to n_patch x n_species
  Et <- with(disturb, mcbrnet::fun_to_m(int, 
                                        n_species = n_species,
                                        n_patch = n_patch,
                                        param_attr = "patch"))
  
  E <- t(Et$m_x)
  
  # intrinsic growth --------------------------------------------------------
  
  Rt <- mcbrnet::fun_to_m(r,
                          n_species = n_species,
                          n_patch = n_patch,
                          param_attr = "species")  
  
  R <- t(Rt$m_x)
  
  # ode ---------------------------------------------------------------------
  
  # n: n_patch x n_species abundance vector
  # N: n_patch x n_species abundance matrix
  # R: n_patch x n_species growth rate matrix
  # psi: indicator parameter scalar
  # E: n_patch x n_species disturbance intensity matrix
  # A: n_species x n_species interaction matrix
  # phi: migration rate
  # C: n_patch x n_patch connectivity matrix
  
  deriv <- function(t, n, parms) {
    with(parms, {
      psi <- input(t)
      N <- matrix(n, nrow = n_patch, ncol = n_species)
      dN <- N * (R + psi * E + N %*% A) + phi * C %*% N
      list(dN)
    })
  }
  
  # define absorbing condition
  ## root function
  rootfun <- function (t, n, parms) {
    return(n - threshold)
  }
  
  ## extinction: triggered when "n - threshold = 0"
  eventfun <- function(t, n, pars) {
    n <- ifelse(n <= threshold, 0, n)
    return(n)
  }
  
  # parameter list
  parms <- list(n_species = n_species,
                n_patch = n_patch,
                R = R,
                E = E,
                A = alpha,
                phi = phi,
                C = connectivity)
  
  # initial values for a state variable
  n_init <- with(n0, runif(n_species * n_patch,
                           min = min,
                           max = max))
  
  # time-series
  times <- seq(0, n_timestep, by = interval)
  
  # run ode solver
  cout <- deSolve::ode(y = n_init,
                       times = times,
                       func = deriv,
                       parms = parms,
                       events = list(func = eventfun, root = TRUE),
                       rootfun = rootfun,
                       ...)
  
  return(cout)
}

