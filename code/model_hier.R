model {
  
  tau0 <- 0.1
    
  # prior -------------------------------------------------------------------
  
  ## variance parameters
  ## - residual SD, sigma[1]
  ## - random SD, sigma[2]
  for (k in 1:2) {
    tau[k] ~ dscaled.gamma(2.5, 6)
    sigma[k] <- pow(sqrt(tau[k]), -1)
  }
  
  ## regression coefficients
  ## - local level
  for (l in 1:5) {
    a[l] ~ dnorm(0, tau0)
  }
    
  ## - watershed level
  for (m in 1:2) {
    b[m] ~ dnorm(0, tau0)
  }
      
  # likelihood --------------------------------------------------------------
  
  ## local level regression
  for (i in 1:Ns) {
    
    C[i] ~ dinterval(logY[i], logY_min[i])
    logY[i] ~ dnorm(mu[i], tau[1])
    
    mu[i] <- 
      a[1] + r[G[i]] +
      a[2] * log(Hsize[i]) +
      a[3] * Forest[i] +
      a[4] * log(Prec[i]) +
      a[5] * Temp[i]
    
  }
  
  ## watershed level regression
  for (j in 1:Nw) {
    r[j] <- 
      b[1] * log(Esize[j]) +
      b[2] * log(Pbranch[j]) +
      eps[j]
    
    eps[j] ~ dnorm(0, tau[2] * W[j])
  }
  
}