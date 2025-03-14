#' DESCRIPTION:
#' Run a set of theoretical analysis

## main analysis for heatmap
source("code/analysis_theory_m_heatmap.R")
rm(list = ls()); gc()

## main analysis for lineart
source("code/analysis_theory_m_line.R")
rm(list = ls()); gc()

## SI analysis for heatmap
source("code/analysis_theory_si.R")
rm(list = ls()); gc()