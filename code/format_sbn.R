#' DESCRIPTION:
#' Produce non-random subset of river networks

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

# generate networks -------------------------------------------------------

## parameters
v_n <- seq(50, 250, by = 50)
v_lambda <- seq(0.2, 1, by = 0.2)
v_rep <- seq_len(4)

df_parms <- expand_grid(
  n = v_n,
  lambda = v_lambda,
  rep = v_rep
) %>% 
  mutate(p = 1 - exp(-lambda))


# generate networks -------------------------------------------------------

n_workers <- floor(0.8 * detectCores())
plan(multisession, workers = n_workers)
registerDoFuture()
handlers(global = TRUE)

tic(msg = "generate SBNs")

list_sbn <- with_progress({
  
  p <- progressor(steps = nrow(df_parms))
  
  foreach(
    i = seq_len(nrow(df_parms)),
    .options.future = 
      list(
        packages = c("mcbrnet",
                     "tidyverse",
                     "igraph"),
        seed = TRUE
      )
  ) %dofuture% {
    
    ## progress
    p()
    
    df_i <- df_parms[i, ]
    
    with(df_i, {
      
      net <- mcbrnet::brnet(
        n_patch = n,
        p_branch = p
      )
      
      g <- net$adjacency_matrix %>% 
        graph_from_adjacency_matrix(mode = "undirected")
      
      V(g)$u <- net$df_patch$n_patch_upstream
      V(g)$b <- net$df_patch$branch_id
      
      attr(g, "n") <- df_i$n
      attr(g, "lambda") <- df_i$lambda
      
      g
    })
    
  } # foreach
}) #with progress

toc()

## release RAM
plan(sequential)
gc()

## export
saveRDS(list_sbn, "data_fmt/sim_sbn.rds")
