#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

## read data
## - sourced `format_data4jags.R` in the following script
source("code/format_est2figure.R")

## latent variables
z <- df_est %>% 
  filter(parms == "z") %>% 
  pull(median)

df_mu_est <- readRDS("data_fmt/output_model_mu_est.rds") %>% 
  mutate(w = score^z) %>% 
  rename(pb = p_branch)


# prediction plot ---------------------------------------------------------

source("code/set_theme.R")
ggplot2::theme_set(default_theme)

x0 <- c("area", "p_branch")
x1 <- c("log_area", "log_pb")

## fcl vs. watershed area
g_area <- df_mu_est %>% 
  ggplot(aes(x = area,
             y = fcl_est)) +
  geom_point(aes(color = factor(h),
                 size = w),
             alpha = 0.4) + 
  geom_line(data = filter(df_yh, focus == "log_area"),
            aes(x = area,
                y = y,
                color = factor(h)),
            alpha = 1,
            linetype = "dashed") + 
  geom_line(data = filter(df_y, focus == "log_area"),
            aes(x = area,
                y = y)) +
  geom_ribbon(data = filter(df_y, focus == "log_area"),
              aes(y = y,
                  ymin = y_low,
                  ymax = y_high,
                  x = area),
              alpha = 0.1) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(y = "Food chain length",
       x = expression("Watershed area ("*km^2*")")) +
  guides(size = "none",
         color = "none")

## fcl vs. branching prob
g_pb <- df_mu_est %>% 
  ggplot(aes(x = pb,
             y = fcl_est)) +
  geom_point(aes(color = factor(h),
                 size = w),
             alpha = 0.4) + 
  geom_line(data = filter(df_yh, focus == "log_pb"),
            aes(x = pb,
                y = y,
                color = factor(h)),
            alpha = 1,
            linetype = "dashed") + 
  geom_line(data = filter(df_y, focus == "log_pb"),
            aes(x = pb,
                y = y)) +
  geom_ribbon(data = filter(df_y, focus == "log_pb"),
              aes(y = y,
                  ymin = y_low,
                  ymax = y_high,
                  x = pb),
              alpha = 0.1) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(y = "Food chain length",
       x = "Branching probability") +
  guides(size = "none",
         color = "none")


# map ---------------------------------------------------------------------

## map layer
sf_lev01 <- readRDS("data_fmt/wgs84_region_lev01.rds")

## site layer
sf_site <- readRDS("data_fmt/wgs84_fcl_site.rds") %>%
  filter(uid %in% uid_incl) %>% 
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
  guides(color = "none",
         fill = "none")


# arrange -----------------------------------------------------------------

layout <- "
AABB
CCCC
"

g_reg <- g_pb + (g_area + theme(axis.title.y = element_blank()))

g_comb <- g_reg + g_map + plot_layout(nrow = 2, 
                                      widths = c(1, 3),
                                      design = layout)

ggsave(g_comb,
       filename = "output/figure_fcl.pdf",
       width = 9,
       height = 8)
