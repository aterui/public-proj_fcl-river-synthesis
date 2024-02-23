model {
  
  tau0 <- 0.1
  df0 <- 6  
  
  # prior -------------------------------------------------------------------
  
  ## variance parameters
  ## - residual SD, sigma[1]
  ## - random SD, sigma[2]
  for (k in 1:3) {
    tau[k] ~ dscaled.gamma(2.5, df0)
    sigma[k] <- pow(sqrt(tau[k]), -1)
  }

  ## regression coefficients
  ## - local level
  for (l in 1:1) {
    a[l] ~ dt(0, tau0, df0)
  }
  
  ## regression coefficients
  ## - watershed level
  for (m in 1:6) {
    b[m] ~ dt(0, tau0, df0)
  }
   
  z ~ dunif(0, 1)
     
  # likelihood --------------------------------------------------------------
  
  ## local level
  for (i in 1:Ns) {
    ## - censoring
    C[i] ~ dinterval(logY[i], logY_min[i])
    logY[i] ~ dt(mu[i], tau[1], df0)
    
    mu[i] <- 
      a0[G[i]] +
      a[1] * log(Hsize[i])
  }
  
  ## watershed level
  for (j in 1:Nw) {
    
    ## - watershed regression
    a0[j] <- 
      b[1] +
      b[2] * log(Esize[j]) +
      b[3] * log(Pbranch[j]) +
      b[4] * scl_prec[j] +
      b[5] * scl_temp[j] +
      b[6] * scl_hfp[j] +
      r[H[j]] +
      eps[j]
    
    eps[j] ~ dnorm(0, tau[2] * scl_w[j])
    scl_w[j] <- w[j] / max(w[])
    log(w[j]) <- z * log(Score[j])
    
    ## - standardization
    scl_prec[j] <- (Prec[j] - mean(Prec[])) / sd(Prec[])
    scl_temp[j] <- (Temp[j] - mean(Temp[])) / sd(Temp[])
    scl_hfp[j] <- (Hfp[j] - mean(Hfp[])) / sd(Hfp[])
  }
  
  for (h in 1:Nh) {
    r[h] ~ dnorm(0, tau[3])
  }
}