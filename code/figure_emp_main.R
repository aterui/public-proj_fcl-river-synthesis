#' DESCRIOTION:
#' Produce empirical figures for maintext
#' Map/prediction plots, and posterior distributions

rm(list = ls())
source("code/set_library.R")
source("code/format_emp_est2figure.R")

# region code mapping
ecor_labels <- c(
  `1` = "Africa",
  `2` = "Europe and Middle East",
  `3` = "Asia",
  `4` = "Australia and Pacific",
  `5` = "South America",
  `6` = "North and Central America",
  `7` = "Arctic"
)

# pull posterior median estimates of the latent z parameters (z[1], z[2], ...)
z <- df_est %>% 
  filter(str_detect(parms, "z\\[.\\]")) %>% 
  pull(median)

# summarize a0 (log-scale FCL intercepts) from the best-fitting model's MCMC samples,
# back-transform to natural scale, and join with observed FCL/weighting data
df_mu_est <- list_est[[id_best]]$mcmc %>% 
  ggmcmc::ggs() %>%                                    # convert coda/mcmc object to tidy long format
  filter(
    str_detect(Parameter, "a0\\[.{1,}\\]")             # keep only a0[g] parameters
  ) %>%   
  group_by(Parameter) %>% 
  summarize(fcl_est = exp(median(value))) %>%          # posterior median, exponentiated back to response scale
  mutate(
    g = str_extract(Parameter, "\\[.{1,}\\]") %>%      # extract group index from parameter name, e.g. "a0[3]" -> 3
      str_remove_all("\\[|\\]") %>% 
      as.numeric()
  ) %>% 
  right_join(df_fcl_wsd) %>%                           # join to FCL/site summary data by group index
  mutate(
    w = (n_site^z[1] * exp(-z[2] * (d_ratio - 1)^2)),  # site-weighting function using latent z params
    scl_w = w / max(w),
    region = ecor_labels[as.character(ecor)]
  )                                                    # rescale weights to max = 1

df_yh <- df_yh %>% 
  mutate(region = ecor_labels[as.character(ecor)])

# first figure ------------------------------------------------------------

## fcl data
# ecoregion / basin lookup table, keyed by outlet id (oid)
df_ecor <- readRDS("data_fmt/data_fcl_reg.rds")[[2]] %>% 
  dplyr::select(
    oid, 
    ecor, 
    hybas_id
  ) %>% 
  mutate(
    region = ecor_labels[as.character(ecor)]
  )

## map layer
# level-01 HydroBASINS polygons, repaired for invalid geometries
sf_lev01 <- readRDS("data_raw/sf_hybas_lev01.rds") %>% 
  st_make_valid()

## site layer
# outlet points joined to ecoregion info, restricted to sites included in the analysis
sf_site <- readRDS("data_raw/sf_outlet.rds") %>%
  left_join(df_ecor,
            by = "oid") %>%
  filter(oid %in% oid_incl)

# simplified + validity-checked basin polygons used as the background map layer
sf_region <- sf_lev01  %>% 
  rmapshaper::ms_simplify() %>% 
  st_make_valid()

# extract x/y from sf geometry since geom_mark_ellipse isn't sf-aware
site_coords <- sf_site %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  )

# build the map: basin polygons + site points colored by ecoregion +
# ellipses summarizing the spatial extent of each ecoregion's sites
g_map <- ggplot(sf_region) +
  geom_sf(
    alpha = 0.5                        # background basin polygons, semi-transparent
  ) +
  geom_sf(
    data = sf_site,
    aes(color = region),         # site points colored by ecoregion
    size = 1
  ) +
  ggforce::geom_mark_ellipse(
    data = site_coords,
    aes(
      x = x, 
      y = y, 
      color = region,
      fill = region
    ),
    linewidth = 0.1,
    alpha = 0.15,
    expand = unit(2, "mm")             # padding around enclosed points; smaller = tighter ellipse
  ) + 
  labs(
    color = "Region",
    fill = "Region"
  ) +
  coord_sf(xlim = c(-163.8, 163.8)) +  # restrict longitude extent of the map
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.ticks = element_blank(),
    panel.border = element_blank()
  )

## prediction plot 
# source("code/set_theme.R")
# ggplot2::theme_set(default_theme)

## fcl vs. watershed area
(g_size <- df_mu_est %>% 
    ggplot(
      aes(x = r_length,
          y = fcl_est)
    ) +
    geom_point(
      aes(
        color = region,
        size = scl_w
      ),
      alpha = 0.5) + 
    # geom_line(data = filter(df_yh, focus == "scl_log_rl"),
    #           aes(x = r_length,
    #               y = y,
    #               color = region),
    #           alpha = 1,
    #           linetype = "dashed") +
    # geom_line(data = filter(df_y, focus == "scl_log_rl"),
    #           aes(x = r_length,
    #               y = y)) +
    # geom_ribbon(data = filter(df_y, focus == "scl_log_rl"),
    #             aes(y = y,
    #                 ymin = y_low,
    #                 ymax = y_high,
    #                 x = r_length),
    #             alpha = 0.1) +
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10() +
    labs(y = "Food chain length",
         x = "Total river length (km)",
         size = "Weight") +
    guides(color = "none") +
    theme_classic() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.title.x = element_text(margin = margin(t = 3))
    )
)

## fcl vs. branching prob
(g_b <- df_mu_est %>% 
    ggplot(aes(x = lambda,
               y = fcl_est)) +
    geom_point(aes(color = region,
                   size = scl_w),
               alpha = 0.4) +
    geom_line(data = filter(df_yh, focus == "scl_log_lambda"),
              aes(x = lambda,
                  y = y,
                  color = region),
              alpha = 1,
              linetype = "dashed") + 
    geom_line(data = filter(df_y, focus == "scl_log_lambda"),
              aes(x = lambda,
                  y = y)) +
    geom_ribbon(data = filter(df_y, focus == "scl_log_lambda"),
                aes(y = y,
                    ymin = y_low,
                    ymax = y_high,
                    x = lambda),
                alpha = 0.1) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(
      y = "Food chain length",
      x = expression("Branching rate (per km)")
    ) +
    guides(size = "none",
           color = "none") +
    theme_classic() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.title.x = element_text(margin = margin(t = 3))
    )
)

## arrange 
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
       filename = "tex/fig_emp_fcl.pdf",
       width = 10,
       height = 7)


# second figure -----------------------------------------------------------

var_name <- c("b[2]" = "ln Total river length",
              "b[3]" = "ln Branching rate",
              "b[4]" = "Air temperature",
              "b[5]" = "Precipitation",
              "b[6]" = "Human footprint",
              "a[1]" = "Elevation")

df_ridge <- list_est[[id_best]]$mcmc %>% 
  ggmcmc::ggs() %>% 
  rename_with(.fn = str_to_lower) %>% 
  filter(str_detect(parameter, "b\\[\\d{1,}\\]|^a\\[\\d{1,}\\]"),
         parameter != "b[1]") %>% 
  mutate(var = var_name[as.character(parameter)]) 

df_stats <- df_ridge %>% 
  group_by(var) %>% 
  summarize(
    b = median(value),
    pp = sprintf(fmt = "%.2f", mean(value > 0)),
    pn = sprintf(fmt = "%.2f", mean(value < 0))
  ) 

var_level <- df_stats %>% 
  arrange(b) %>% 
  pull(var)

(g_ridge <- df_ridge %>% 
    mutate(var = factor(var, levels = var_level)) %>% 
    ggplot(
      aes(
        x = value,
        y = var,
        fill = 0.5 - abs(0.5 - after_stat(ecdf))
      )
    ) +
    ggridges::stat_density_ridges(
      quantile_lines = TRUE,
      calc_ecdf = TRUE,
      geom = "density_ridges_gradient",
      quantiles = 0.5,
      color = grey(0, 0.2), 
      size = 0.25,
      scale = 0.95
    ) +
    scale_fill_gradient(
      low = "white",
      high = "salmon",
      name = "Tail prob."
    ) +
    geom_vline(
      xintercept = 0, 
      linetype = "dashed",
      color = grey(0.5, 0.5),
      linewidth = 0.25
    ) +
    ggridges::theme_ridges() +
    geom_text(
      data = df_stats %>%
        mutate(var = factor(var, levels = var_level)),
      aes(x = Inf, y = var, label = pp),
      inherit.aes = FALSE,
      hjust = 1.0,
      vjust = -0.3,
      size = 3
    ) +
    geom_text(
      data = df_stats %>%
        mutate(var = factor(var, levels = var_level)),
      aes(x = -Inf, y = var, label = pn),
      inherit.aes = FALSE,
      hjust = -1.0,
      vjust = -0.3,
      size = 3
    ) +
    theme(axis.title.x = element_text(hjust = 0.5),  # center x-axis label
          axis.title.y = element_text(hjust = 0.5)   # center y-axis label
    ) +
    labs(x = "Posterior estimate",
         y = "Predictor"))


## export
ggsave(g_ridge,
       filename = "tex/fig_emp_ridge.pdf",
       width = 7,
       height = 4)
