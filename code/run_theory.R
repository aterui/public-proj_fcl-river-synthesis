#' DESCRIPTION:
#' Run a set of theoretical analysis

## main: predator-prey analysis
source("code/analysis_theo_m_2sp.R")
rm(list = ls()); gc()

## main: food web analysis
source("code/analysis_theo_m_fw.R")
rm(list = ls()); gc()

## SI analysis for heatmap
source("code/analysis_theo_si_fw.R")
rm(list = ls()); gc()

## SI analysis for spatially-explicit simulations
source("code/analysis_theo_si_2sp_sbn.R")
rm(list = ls()); gc()