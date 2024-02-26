#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/format_data4jags.R")

## read data
df_est <- readRDS("data_fmt/output_model_summary.rds")
df_mu_est <- readRDS("data_fmt/output_model_mu_est.rds")
df_pred <- readRDS("data_fmt/output_model_pred.rds")
mcmc <- readRDS("data_fmt/output_model_mcmc.rds")


# format data -------------------------------------------------------------

## weight values for each observation
z <- df_est %>%
  filter(parms == "z") %>%
  pull(median)

df_mu_est <- df_mu_est %>%
  mutate(w = score^z / max(score^z))

x <- with(df_g,
          seq(min(log(p_branch)),
              max(log(p_branch)),
              length = 100))

X <- with(df_g,
          cbind(1,
                mean(log(area)),
                x))

b <- rbind(mcmc[,"b0"], mcmc[,"b[1]"], mcmc[,"b[2]"])

df_y <- (X %*% b) %>% 
  exp() %>% 
  apply(MARGIN = 1, quantile, c(0.025, 0.5, 0.975)) %>% 
  t() %>% 
  as_tibble() %>% 
  setNames(c("low", "median", "high")) %>% 
  mutate(x_log_pb = x,
         x_pb = exp(x))


# prediction plot ---------------------------------------------------------

source("code/set_theme.R")
ggplot2::theme_set(default_theme)

## fcl vs. p_branch
g_pb <- df_mu_est %>% 
  ggplot(aes(x = p_branch,
             y = fcl_est)) +
  geom_point(aes(size = w,
                 color = factor(h)),
             alpha = 0.4) + 
  geom_line(data = df_pred,
            aes(x = pb,
                y = y,
                color = factor(h)),
            alpha = 1,
            linetype = "dashed") + 
  geom_line(data = df_y,
            aes(x = x_pb,
                y = median)) +
  geom_ribbon(data = df_y,
              aes(y = median,
                  ymin = low,
                  ymax = high,
                  x = x_pb),
              alpha = 0.1) +
  #facet_wrap(facets = ~ factor(h)) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(y = "Food chain length",
       x = "Branching probability") +
  guides(size = "none",
         color = "none")

## map
sf_lev01 <- readRDS("data_fmt/wgs84_region_lev01.rds")

sf_site <- readRDS("data_fmt/wgs84_fcl_site.rds") %>% 
  group_by(sid) %>% 
  slice(1) %>% 
  dplyr::select(NULL) %>% 
  st_join(sf_lev01) %>% 
  mutate(h = as.numeric(factor(id_lev01)))

df_site_ref <- sf_site %>% 
  as_tibble() %>% 
  dplyr::select(-geometry) %>% 
  distinct(id_lev01, h) %>% 
  arrange(id_lev01)

sf_region <- sf_lev01  %>% 
  rmapshaper::ms_simplify() %>%
  left_join(df_site_ref) %>% 
  st_make_valid()

g_map <- ggplot(sf_region) +
  geom_sf(aes(fill = factor(h)),
          alpha = 0.4) +
  geom_sf(data = sf_site,
          aes(color = factor(h))) +
  guides(fill = "none",
         color = "none")

g_comb <- g_pb / g_map


ggsave(g_comb,
       filename = "output/figure_fcl_pb.pdf",
       width = 7,
       height = 9)
