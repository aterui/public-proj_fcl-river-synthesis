#' DESCRIPTION:
#' Analyze the relationship  between river discharge and regional precipitation

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")


# data --------------------------------------------------------------------

## read data used for main regression analysis
list_reg <- readRDS("data_fmt/data_fcl_reg.rds")
df_fsd <- readRDS("data_fmt/data_env_fsd.rds") %>% 
  left_join(list_reg[[1]] %>% 
              select(sid, oid, h01, local_area))

## data formatting
df_m <- left_join(df_fsd,
                  list_reg[[2]] %>% 
                    select(oid, prec, area, lambda)) %>% 
  mutate(scl_fnu = scale(log(fnu)),
         scl_prec = scale(log(prec))) %>% 
  drop_na(fnu, oid, h01)

# analysis ----------------------------------------------------------------

fit0 <- lme4::lmer(log(fnu) ~ log(prec) + (1 | oid) + (1 | h01),
                   df_m,
                   REML = FALSE)

summary(fit0)

# df_m %>%
#   ggplot(aes(x = lambda,
#              y = fmu)) +
#   geom_point() +
#   facet_wrap(facets =~ h01,
#              scales = "free") +
#   scale_x_continuous(trans = "log10") +
#   scale_y_continuous(trans = "log10")
