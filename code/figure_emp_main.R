#' DESCRIOTION:
#' xxx

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")

## - sourced `format_data4jags.R` in the following script
source("code/format_emp_est2figure.R")

# format data -------------------------------------------------------------

## latent variables
z <- df_est %>% 
  filter(str_detect(parms, "z\\[.\\]")) %>% 
  pull(median)

df_mu_est <- list_est[[id_best]]$mcmc %>% 
  ggmcmc::ggs() %>% 
  filter(str_detect(Parameter, "a0\\[.{1,}\\]")) %>% 
  group_by(Parameter) %>% 
  summarize(fcl_est = exp(median(value))) %>% 
  mutate(g = str_extract(Parameter, "\\[.{1,}\\]") %>% 
           str_remove_all("\\[|\\]") %>% 
           as.numeric()) %>% 
  right_join(df_fcl_wsd) %>% 
  mutate(w = (n_site^z[1] * exp(-z[2] * (d_ratio - 1)^2)),
         scl_w = w / max(w))


# map ---------------------------------------------------------------------

## map layer
sf_lev01 <- readRDS("data_raw/wgs84_region_lev01.rds") %>% 
  st_make_valid()

## site layer
sf_site <- readRDS("data_raw/wgs84_outlet.rds") %>%
  left_join(readRDS("data_raw/wgs84_wsd_sub.rds") %>% 
              as_tibble() %>% 
              dplyr::select(uid, oid),
            by = "oid") %>%
  filter(uid %in% uid_incl) %>% 
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
          alpha = 0.5) +
  geom_sf(data = sf_site,
          aes(color = factor(h)),
          size = 1) +
  guides(color = "none",
         fill = "none") +
  coord_sf(xlim = c(-163.8, 163.8)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.ticks = element_blank(),
        panel.border = element_blank())


# prediction plot ---------------------------------------------------------

source("code/set_theme.R")
ggplot2::theme_set(default_theme)

## fcl vs. watershed area
(g_size <- df_mu_est %>% 
  ggplot(aes(x = r_length,
             y = fcl_est)) +
  geom_point(aes(color = factor(h),
                 size = scl_w),
             alpha = 0.5) + 
  # geom_line(data = filter(df_yh, focus == "log_r_length"),
  #           aes(x = r_length,
  #               y = y,
  #               color = factor(h)),
  #           alpha = 1,
  #           linetype = "dashed") + 
  # geom_line(data = filter(df_y, focus == "log_r_length"),
  #           aes(x = r_length,
  #               y = y)) +
  # geom_ribbon(data = filter(df_y, focus == "log_r_length"),
  #             aes(y = y,
  #                 ymin = y_low,
  #                 ymax = y_high,
  #                 x = r_length),
  #             alpha = 0.1) +
  scale_x_log10(labels = scales::label_log(digits = 2)) +
  scale_y_log10() +
  labs(y = "Food chain length",
       x = "River length (km)",
       size = "Weight") +
    guides(color = "none"))

## fcl vs. branching prob
(g_b <- df_mu_est %>% 
  ggplot(aes(x = lambda,
             y = fcl_est)) +
  geom_point(aes(color = factor(h),
                 size = scl_w),
             alpha = 0.4) +
  geom_line(data = filter(df_yh, focus == "log_lambda"),
            aes(x = lambda,
                y = y,
                color = factor(h)),
            alpha = 1,
            linetype = "dashed") + 
  geom_line(data = filter(df_y, focus == "log_lambda"),
            aes(x = lambda,
                y = y)) +
  geom_ribbon(data = filter(df_y, focus == "log_lambda"),
              aes(y = y,
                  ymin = y_low,
                  ymax = y_high,
                  x = lambda),
              alpha = 0.1) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(y = "Food chain length",
       x = expression("Branching rate ("*km^-1*")")) +
  guides(size = "none",
         color = "none"))

# arrange -----------------------------------------------------------------

layout <- "
AAAA
BBCC
"

g_comb <- g_map + 
  g_b + 
  (g_size + theme(axis.title.y = element_blank())) +
  plot_layout(nrow = 2, 
              widths = c(3, 1),
              design = layout) +
  plot_annotation(tag_levels = "A")

ggsave(g_comb,
       filename = "output/fig_emp_fcl.pdf",
       width = 9.5,
       height = 8)
