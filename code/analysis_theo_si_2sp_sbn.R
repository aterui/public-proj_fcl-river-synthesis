#' DESCRIPTION:
#' Patch-level simulation using OCN

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

# read sbn ----------------------------------------------------------------

list_sbn <- readRDS("data_fmt/sim_sbn.rds")

df_parms <- expand_grid(
  theta = c(40, 80),
  kappa = c(0.5, 0.8),
  nu = c(0.1, 0.05),
  mu0 = c(2.5, 5),
  cp = 0.5,
  mu_p = 5,
  nid = seq_len(length(list_sbn)) # network id
) %>% 
  mutate(
    theta_d = kappa * theta,
    theta_u = (1 - kappa) * theta,
    .after = theta
  )


# simulation run ----------------------------------------------------------

n_workers <- floor(0.8 * detectCores())
plan(multisession, workers = n_workers)
registerDoFuture()
handlers(global = TRUE)

tic(msg = "spatially explicit model run")

df_sim <- 
  with_progress({
    
    p <- progressor(steps = nrow(df_parms))
    
    foreach(
      i = seq_len(nrow(df_parms)),
      .combine = bind_rows,
      .options.future = 
        list(
          packages = c("rpom",
                       "tidyverse",
                       "igraph"),
          seed = TRUE
        )
    ) %dofuture% {
      ## parameter set
      df_i <- df_parms[i, ]
      nid <- df_i$nid
      sbn <- list_sbn[[nid]]
      
      ## read sbn (igraph object)
      m <- as_adjacency_matrix(sbn) %>% 
        as.matrix()
      n <- vcount(sbn)
      
      ## upstream connectivity matrix "xi"
      d <- distances(sbn)
      phi <- d[1, ] # distance to the outlet
      o <- outer(phi, phi, FUN = "-")
      up <- (d - o) / 2
      down <- (d + o) / 2
      xi <- as.matrix(m * up)
      
      ## upstream river length
      ## original value included the focal patch itself (-1 removes it)
      u <- V(sbn)$u - 1
      
      ## run simulation
      m_p <- with(df_i, {
        ## expected synchrony for each patch
        ## - length-bias correction "w", then take expectations
        ## - '(1 - exp(-nu * D2t)) / (nu * D2t)' is E(exp(-nu d)) when the disturbance origin is uniform over the flow path
        idx <- which(V(sbn)$u == 1) # tip index
        D2t <- up[, idx, drop = FALSE]
        D2t[down[, idx] > 0] <- NA
        w <- D2t / rowSums(D2t, na.rm = TRUE)
        rho <- rowSums(w * ((1 - exp(-nu * D2t)) / (nu * D2t)), na.rm = TRUE)
        
        nspom(
          m = m, # adj matrix
          xi = xi, # upstream adj matrix
          theta = c(theta_u, theta_d), # c(up-move, down-move)
          cp = cp, # establishment prob.
          mu = c(mu0, mu_p), # c(disturbance, prey-induced)
          u = u, # upstream river length
          rho = rho,
          intv = 1,
          threshold = 1E-5,
          n_timestep = 100
        )
        
      })
      
      v_p <- m_p[nrow(m_p), -1]
      p_prey <- mean(v_p[seq_len(n)])
      p_pred <- mean(v_p[seq(n + 1, 2 * n, by = 1)])
      
      ## progress
      p()
      
      ## output
      df_i %>% 
        bind_cols(
          tibble(
            p = c(p_prey, p_pred),
            tl = c("prey", "predator"),
            n = n,
            lambda = attr(sbn, "lambda")
          )
        )
      
    }# foreach
  })# with_progress

toc()

## release RAM
plan(sequential)
gc()

## export
saveRDS(df_sim, "data_fmt/sim_fcl_2sp_nspom.rds")


# scaling exponent test ---------------------------------------------------

n <- 500
lambda <- c(0.4, 0.55, 0.7)
p <- 1 - exp(-lambda)

n_workers <- length(p)
plan(multisession, workers = n_workers)
registerDoFuture()
handlers(global = TRUE)

df_pr <- foreach(
  i = seq_along(p),
  .combine = bind_rows,
  .options.future = 
    list(seed = TRUE)
) %dofuture% {
  
  a <- replicate(10, {
    net <- mcbrnet::brnet(n_patch = n, p_branch = p[i])
    net$df_patch$n_patch_upstream    
  }) |>
    as.vector()
  
  f <- ecdf(a)
  x <- seq(1, max(a), by = 1)
  
  tibble(
    pr = 1 - f(x),
    a = seq_len(max(a)),
    p = p[i],
    lambda = lambda[i]
  )
}

plan(sequential)
gc()

split(df_pr, df_pr$lambda) %>% 
  lapply(function(data) {
    
    lm(log(pr) ~ log(a), filter(data, a <= ceiling(n * 0.025)))
    
  })

df_pr %>% 
  filter(a < max(a)) %>% 
  ggplot() +
  geom_point(
    aes(
      x = a,
      y = pr,
      color = factor(lambda)
    )
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    y = "Pr(A > a)",
    x = "a",
    color = expression(lambda[b])
  ) +
  geom_vline(xintercept = ceiling(n * 0.02))

