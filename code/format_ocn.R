#' DESCRIPTION:
#' Generate OCN landscapes across different landscape sizes and cooling rates.
#' For each landscape, estimate a scaling coefficient (lambda) from the
#' relationship between habitat density (inverse mean link length) and
#' aggregation threshold.

source("code/set_library.R")

# parameter ---------------------------------------------------------------

# Landscape size
dim <- seq(10, 60, by = 5)

# Cooling rate (higher values produce more dendritic networks)
cr <- seq(0.5, 10, length = 10)

# Aggregation thresholds
res <- seq(1, 15, by = 1)

# Parameter combinations
df_dc <- expand_grid(dim, cr)

# generate OCN ------------------------------------------------------------

## RUN
n_workers <- floor(0.8 * detectCores())
plan(multisession, workers = n_workers)
registerDoFuture()
handlers(global = TRUE)

tic(msg = "generate ocn")

list_ocn <- 
  with_progress({
    
    p <- progressor(steps = nrow(df_dc))
    
    foreach(
      i = seq_len(nrow(df_dc)),
      .options.future = 
        list(
          packages = c("OCNet",
                       "tidyverse"),
          seed = TRUE
        )
    ) %dofuture% {
      
      # Generate OCN landscape
      df_i <- df_dc[i, ]
      ocn <- with(df_i, {
        create_OCN(
          dimX = dim,
          dimY = dim,
          coolingRate = cr
        ) %>%
          landscape_OCN()
      })
      
      # Branching rate (= inverse mean link length)
      y <- sapply(res, function(x) {
        v_link <- ocn %>%
          aggregate_OCN(thrA = x) %>%
          paths_OCN(includePaths = TRUE) %>%
          { .@AG$leng[.@AG$leng > 0] }
        
        if (length(v_link) < 5) NA else 1 / mean(v_link)
      })
      
      # Aggregate network
      ocn <- ocn %>%
        aggregate_OCN(thrA = 3) %>%
        paths_OCN(includePaths = TRUE)
      
      # Scaling coefficient (intercept of log-log regression)
      lambda <- MASS::rlm(log(y) ~ log(res), maxit = 100) %>%
        { coef(.)[1] } %>%
        exp()
      
      # Store lambda
      attr(ocn, "lambda") <- unname(lambda)
      
      # Progress
      p()
      
      ocn
    } #foreach
    
  })# with progress

toc()

## export
saveRDS(list_ocn, "data_fmt/sim_ocn.rds")
