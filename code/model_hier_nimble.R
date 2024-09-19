nimble::nimbleCode({
  
  sigma0 <- 10
  df_sigma <- 15
  
  # prior -------------------------------------------------------------------
  
  ## variance parameters
  ## - residual SD, sigma[1]
  ## - watershed random SD, sigma[2]
  ## - region random SD, sigma[3]
  for (k in 1:3) {
    sigma[k] ~ dt(0, sigma = 2.5, df = df_sigma)
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
    logY[i] ~ dt(mu[i], sd = sigma[1], nu)
    
    mu[i] <- a0[G[i]] + inprod(a[], X1[i, ])
  }
  
  ## watershed level
  for (j in 1:Nw) {
    
    ## - watershed regression
    a0[j] ~ dnorm(a0_hat[j], var = pow(sigma[2], 2) / scl_w[j])
    a0_hat[j] <- inprod(b[], X2[j, ]) + r[H[j]]
    
    ## - "scl_w[j]" scaled weight
    ## - "Ratio[j]" is the distance ratio to randomly-generated sites
    ## - "xi[j]", squared deviation from the random samples
    ## - "N_site[j]" is the number of sites within a watershed
    ## - note: w = N_site^z[1] * exp(-z[2] * xi)
    scl_w[j] <- w[j] / max(w[])
    log(w[j]) <- z[1] * log(N_site[j]) - z[2] * xi[j]
    xi[j] <- pow((Ratio[j] - 1), 2)
  }
  
  for (h in 1:Nh) {
    r[h] ~ dnorm(0, sd = sigma[3])
    b0[h] <- b[1] + r[h]
  }
  
})