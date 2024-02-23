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
  for (l in 1:3) {
    a[l] ~ dnorm(0, tau0)
  }
    
  ## - watershed level
  for (m in 1:4) {
    b[m] ~ dnorm(0, tau0)
  }
      
  # likelihood --------------------------------------------------------------
  
  ## local level regression
  for (i in 1:Ns) {
    
    C[i] ~ dinterval(logY[i], logY_min[i])
    logY[i] ~ dnorm(mu[i], tau[1])
    
    mu[i] <- 
      a[1] + 
      a[2] * log(Hsize[i]) +
      a[3] * scl_hfp[i] +
      r[G[i]]
  
    scl_hfp[i] <- (Hfp[i] - mean(Hfp[])) / sd(Hfp[])
    
  }
  
  ## watershed level regression
  for (j in 1:Nw) {
    r[j] <- 
      b[1] * log(Esize[j]) +
      b[2] * log(Pbranch[j]) +
      b[3] * scl_prec[j] +
      b[4] * scl_temp[j] +
      eps[j]
    
    eps[j] ~ dnorm(0, tau[2] * W[j])
    
    scl_prec[j] <- (Prec[j] - mean(Prec[])) / sd(Prec[])
    scl_temp[j] <- (Temp[j] - mean(Temp[])) / sd(Temp[])
  }
  
}