
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
  
  sigma0 <- 10
  df_sigma <- 15
  
  # prior -------------------------------------------------------------------
  
  ## variance parameters
  ## - residual SD, sigma[1]
  ## - watershed random SD, sigma[2]
  for (k in 1:2) {
    sigma[k] ~ T(dt(0, sigma = 2.5, df = df_sigma), 0, )
  }
  
  ## local level
  ## - coefficients
  for (l in 1:K1) {
    a[l] ~ dnorm(0, sd = sigma0)
  }
  
  ## watershed level
  ## - intercept/coefficients
  for (h in 1:Nh) {
    b[h, 1:K2] ~ dmnorm(v_mu_b[1:K2], cholesky = U[1:K2, 1:K2], prec_param = 0)
  }
  
  for (k in 1:K2) {
    v_mu_b[k] ~ dnorm(0, sd = sigma0)
    v_sigma_b[k] ~ T(dt(0, sigma = 2.5, df = df_sigma), 0, )
  }
  
  Ustar[1:K2, 1:K2] ~ dlkj_corr_cholesky(2, K2)
  U[1:K2,1:K2] <- uppertri_mult_diag(Ustar[1:K2, 1:K2], v_sigma_b[1:K2]) 
  rho[1:K2, 1:K2] <- t(Ustar[1:K2, 1:K2]) %*% Ustar[1:K2, 1:K2]
  
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
    a0_hat[j] <- sum(b[H[j], 1:K2] * X2[j, 1:K2])
    
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