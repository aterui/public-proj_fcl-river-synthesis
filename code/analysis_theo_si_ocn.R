#' DESCRIPTION:
#' Patch-level simulation using OCN

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

# read ocn ----------------------------------------------------------------

list_ocn <- readRDS("data_fmt/sim_ocn.rds")

## read OCN
ocn <- list_ocn[[1]]

## convert to igraph
g <- OCN_to_igraph(ocn, level = "RN") 

## full adjacency matrix
(m <- g %>% 
    as_undirected() %>% 
    as_adjacency_matrix())

## downstream connectivity
(xi <- g %>% 
  as_adjacency_matrix())

D <- igraph::distances(g, mode = "in")

rn <- ocn@RN
tips <- which(rn$nUpstream == 1)
(D2t <- D[, tips])

rho <- rowMeans(D2t * ((1 - exp(-nu * D2t)) / (nu * D2t)), na.rm = TRUE)
rho[tips] <- 0


# V(g)$u <- rn$nUpstream
# ggraph::ggraph(g) +
#   ggraph::geom_node_point(aes(color = factor(u))) +
#   ggraph::geom_edge_arc()
