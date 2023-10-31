
# library -----------------------------------------------------------------

source("code/library.R")

# foodweb -----------------------------------------------------------------

## replicate number
n_fw <- 5

## food web setup
## 0.18 proportion of basal species - see Briand and Cohen 1987 Nature
## 0.11 connectance - see Dunne et al. 2002 PNAS
## 0.25 theta - see Johnson et al. 2014 PNAS
s <- 16
b <- floor(s * 0.18)
l <- floor(0.11 * s ^ 2)
theta <- 0.25

## generate food webs
set.seed(123)
list_fw <- replicate(n = n_fw,
                     foodweb(n_species = s,
                             n_basal = b,
                             l = l,
                             theta = theta),
                     simplify = FALSE)

saveRDS(list_fw, "data_fmt/parms_foodweb.rds")
