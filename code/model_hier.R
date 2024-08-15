model {
  
  sd0 <- 1
  tau0 <- pow(sd0, -2)
  df_tau <- 10
  
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
  for (l in 1:K1) {
    a[l] ~ dnorm(0, tau0)
  }
  
  ## watershed level
  ## - intercept
  b0 ~ dnorm(0, tau0)

  ## - coefficients
  for (m in 1:K2) {
    b[m] ~ dnorm(0, tau0)
  }
  
  ## weight scaling exponent
  z ~ dnorm(0, tau0)T(0, )
  
  ## degree of freedom
  nu ~ dexp(0.1)T(2, )
  
  # likelihood --------------------------------------------------------------
  
  ## local level
  for (i in 1:Ns) {
    ## - censoring
    C[i] ~ dinterval(logY[i], logY_min[i])
    logY[i] ~ dt(mu[i], tau[1], nu)
    
    mu[i] <- a0[G[i]] + inprod(a[], X1[i, ])
  }
  
  ## watershed level
  for (j in 1:Nw) {
    
    ## - watershed regression
    a0[j] ~ dnorm(a0_hat[j], tau[2] * scl_w[j])
    a0_hat[j] <- r[H[j]] + inprod(b[], X2[j, ])
    
    ## - "scl_w[j]" scaled weight
    ## - "Ratio[j]" is the distance ratio to randomly-generated sites
    ## - "xi[j]", squared deviation from the random samples
    ## - "N_site[j]" is the number of sites within a watershed
    scl_w[j] <- w[j] / max(w[])
    log(w[j]) <- z * (log(N_site[j]) - xi[j])
    xi[j] <- pow((Ratio[j] - 1), 2)
  }
  
  for (h in 1:Nh) {
    r[h] ~ dnorm(b0, tau[3])
  }

}