#' DESCRIPTION:
#' model code for a random intercept model

m0 <- nimble::nimbleCode({
  
  ## sd for coef prior
  sigma0 <- 5
  
  ## sd and degrees of freedom for residual/random effects
  sigma_r <- 1
  df_sigma <- 10
  
  # prior -------------------------------------------------------------------
  
  ## variance parameters
  ## - residual SD, sigma[1]
  ## - watershed random SD, sigma[2]
  ## - region random SD, sigma[3]
  for (k in 1:3) {
    sigma[k] ~ T(dt(0, sigma = sigma_r, df = df_sigma), 0, )
  }
  
  ## local level
  ## - coefficients
  for (l in 1:K1) {
    a[l] ~ dnorm(0, sd = sigma0)
  }
  
  ## watershed level
  ## - intercept/coefficients
  for (m in 1:K2) {
    b[m] ~ dnorm(0, sd = sigma0)
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
    a0[j] ~ dnorm(a0_hat[j], tau = tau_w[j])
    a0_hat[j] <- sum(b[1:K2] * X2[j, 1:K2]) + r[H[j]]
    
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
  
  for (h in 1:Nh) {
    r[h] ~ dnorm(0, sd = sigma[3])
    b0[h] <- b[1] + r[h]
  }
  
})