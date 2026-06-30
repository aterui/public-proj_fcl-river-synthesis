#' DESCRIPTION:
#' Set packages
 
pacman::p_load(
  ## general use
  tidyverse,
  patchwork,
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
  parallel,
  doFuture,
  progressr,
  tictoc,
  ## gis
  mapview,
  sf,
  terra)
