#' DESCRIPTION:
#' Set packages
 
pacman::p_load(
  ## general use
  tidyverse,
  patchwork,
  sf,
  terra,
  igraph,
  magick,
  cowplot,
  nimble,
  xtable,
  ecotools,
  ## simulation
  foreach,
  mcbrnet,
  rpom,
  doParallel,
  doSNOW)
