#' DESCRIPTION:
#' model code for a random intercept-slope model

## define function
uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2), vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p, ncol = p, init = FALSE)
    for(i in 1:p)
      out[ , i] <- mat[ , i] * vec[i]
    return(out)
  })

m1 <- nimble::nimbleCode({
  
  ## sd for coef prior
  sigma0 <- 5
  
  ## sd and degrees of freedom for residual/random effects
  sigma_r <- 1
  df_sigma <- 10
  
  # prior -------------------------------------------------------------------
  
  ## variance parameters
  ## - residual SD, sigma[1]
  ## - watershed random SD, sigma[2]
  for (k in 1:2) {
    sigma[k] ~ T(dt(0, sigma = sigma_r, df = df_sigma), 0, )
  }
  
  ## local level
  ## - coefficients
  for (l in 1:K1) {
    a[l] ~ dnorm(0, sd = sigma0)
  }
  
  ## watershed level
  ## - intercept/coefficients
  ## - intercept, length and/or branch
  for (h in 1:Nh) {
    b[h, 1:R] ~ dmnorm(b_mu[1:R], cholesky = U[1:R, 1:R], prec_param = 0)
  }
  
  for (k in 1:R) {
    b_mu[k] ~ dnorm(0, sd = sigma0)
    sigma_b[k] ~ T(dt(0, sigma = sigma_r, df = df_sigma), 0, )
  }
  
  Ustar[1:R, 1:R] ~ dlkj_corr_cholesky(2, R)
  U[1:R,1:R] <- uppertri_mult_diag(Ustar[1:R, 1:R], sigma_b[1:R]) 
  rho[1:R, 1:R] <- t(Ustar[1:R, 1:R]) %*% Ustar[1:R, 1:R]
  
  ## - other covariates
  for (k in 1:(K2 - R)) {
    b_prime[k] ~ dnorm(0, sd = sigma0)
  }
  
  ## weight scaling exponent
  z[1] ~ T(dnorm(0, sd = 1), 0, )
  z[2] ~ T(dnorm(0, sd = log(10)), 0, )
  
  ## degree of freedom
  nu ~ T(dexp(0.1), 2, )
  
  # likelihood --------------------------------------------------------------
  
  ## local level
  for (i in 1:Ns) {
    ## - censoring
    C[i] ~ dinterval(logY[i], logCut[i])
    logY[i] ~ dt(mu[i], sigma = sigma[1], df = nu)
    
    mu[i] <- a0[G[i]] + a[1] * X1[i, 1]
  }
  
  ## watershed level
  for (j in 1:Nw) {
    
    ## - watershed regression
    ## - 1:R, random slope
    ## - R:K2, fixed slope
    a0[j] ~ dnorm(a0_hat[j], tau = tau_w[j])
    a0_hat[j] <- 
      sum(b[H[j], 1:R] * X2[j, 1:R]) + 
      sum(b_prime[1:(K2 - R)] * X2[j, (R + 1):K2])
    
    ## - "scl_w[j]" scaled weight
    ## - "Ratio[j]" is the distance ratio to randomly-generated sites
    ## - "xi[j]", squared deviation from the random samples
    ## - "N_site[j]" is the number of sites within a watershed
    ## - note: w = N_site^z[1] * exp(-z[2] * xi)
    tau_w[j] <- pow(sigma[2], -2) * scl_w[j]
    scl_w[j] <- w[j] / max(w[1:Nw])
    log(w[j]) <- z[1] * log(N_site[j]) - z[2] * xi[j]
    xi[j] <- pow((Ratio[j] - 1), 2)
  }
  
})