#' DESCRIPTION:
#' Produce non-random subset of river networks

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

b <- replicate(30, {
  n <- 200
  net <- mcbrnet::brnet(
    n_patch = n,
    p_branch = 0.40
  )
  
  frac <- 0.1
  u <- net$df_patch[["n_patch_upstream"]]
  f <- stats::ecdf(u)
  a <- 1:ceiling(n * frac)
  pr <- 1 - f(a)
  m <- lm(log(pr) ~ log(a))
  
  a <- 1:(n - 1)
  pr <- 1 - f(a)
  plot(pr ~ a, log = "xy")
  abline(m)
  abline(v = ceiling(n*frac), lty = 2)
  
  coef(m)[2]
})
