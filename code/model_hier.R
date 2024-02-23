model {
  
  tau0 <- 0.1
  df_tau <- 2 
  
  # prior -------------------------------------------------------------------
  
  ## variance parameters
  ## - residual SD, sigma[1]
  ## - watershed random SD, sigma[2]
  ## - region random SD, sigma[3]
  for (k in 1:3) {
    tau[k] ~ dscaled.gamma(2.5, df_tau)
    sigma[k] <- pow(sqrt(tau[k]), -1)
  }

  ## local level
  ## - coefficients
  for (l in 1:2) {
    a[l] ~ dnorm(0, tau0)
  }
  
  ## watershed level
  ## - intercept
  b0 ~ dnorm(0, tau0)
  
  ## - coefficients
  for (m in 1:5) {
    b[m] ~ dnorm(0, tau0)
  }
  
  ## weight scaling exponent
  z ~ dunif(0, 1)
   
  ## degree of freedom
  nu ~ dexp(0.01)T(2, )
    
  # likelihood --------------------------------------------------------------
  
  ## local level
  for (i in 1:Ns) {
    ## - censoring
    C[i] ~ dinterval(logY[i], logY_min[i])
    logY[i] ~ dt(mu[i], tau[1], nu)
    
    mu[i] <- a0[G[i]] + 
      a[1] * log(Hsize[i]) +
      a[2] * scl_forest[i]
    
    scl_forest[i] <- (Forest[i] - mean(Forest[])) / sd(Forest[])
  }
  
  ## watershed level
  for (j in 1:Nw) {
    
    ## - watershed regression
    a0[j] <- 
      r[H[j]] +
      b[1] * log(Esize[j]) +
      b[2] * log(Pbranch[j]) +
      b[3] * scl_prec[j] +
      b[4] * scl_temp[j] +
      b[5] * scl_hfp[j] +
      eps[j]
    
    eps[j] ~ dnorm(0, tau[2] * scl_w[j])
    
    ## - scaled weight
    scl_w[j] <- w[j] / max(w[])
    log(w[j]) <- z * log(Score[j])
    
    ## - standardization
    scl_prec[j] <- (Prec[j] - mean(Prec[])) / sd(Prec[])
    scl_temp[j] <- (Temp[j] - mean(Temp[])) / sd(Temp[])
    scl_hfp[j] <- (Hfp[j] - mean(Hfp[])) / sd(Hfp[])
  }
  
  for (h in 1:Nh) {
    r[h] ~ dnorm(b0, tau[3])
  }
}