#' DESCRIPTION:
#' Run a set of theoretical analysis

## main analysis for heatmap
source("code/analysis_theo_m_htmap.R")
rm(list = ls()); gc()

## main analysis for lineart
source("code/analysis_theo_m_line.R")
rm(list = ls()); gc()

## SI analysis for heatmap
source("code/analysis_theo_si.R")
rm(list = ls()); gc()