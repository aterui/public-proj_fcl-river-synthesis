
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
  i_pick <- resample(i_index[-i0],
                     size = kappa,
                     prob = p_ij[-i0])
  
  i <- c(i0, i_pick)
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

