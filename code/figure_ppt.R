#' DESCRIOTION:
#' Produce scheme figure and main theoretical predictions
#' Main theoretical predictions include heatmap and line arts

# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")

## set plot theme
source("code/set_theme.R")
ggplot2::theme_set(default_theme)

# figure theory -----------------------------------------------------------

## read data ####
df_fcl_line <- readRDS("data_fmt/sim_fcl_m_line.rds") %>% 
  filter(rsrc == max(rsrc),
         mu0 == max(mu0))

df_fcl_mu <- df_fcl_line %>% 
  group_by(focus, lambda, rl, rsrc, mu0, rho) %>% 
  summarize(mu_fcl = mean(fcl),
            min_fcl = min(fcl),
            max_fcl = max(fcl)) %>% 
  ungroup() %>% 
  filter(rho > 0)

## food chain length vs. branching ####
g_br <- df_fcl_mu %>% 
  filter(focus == "branch") %>% 
  ggplot(aes(x = lambda,
             y = mu_fcl)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = min_fcl,
                  ymax = max_fcl),
              color = NA,
              fill = grey(0, 0.2)) +
  labs(x = expression("Branching rate"~lambda[b]),
       y = "Food chain length")

## food chain length vs. size ####
g_size <- df_fcl_mu %>% 
  filter(focus == "size") %>% 
  ggplot(aes(x = rl,
             y = mu_fcl)) +
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = min_fcl,
                  ymax = max_fcl),
              color = NA,
              fill = grey(0, 0.2)) +
  labs(x = expression("Total river length"~italic(L)),
       y = "Food chain length")

## export ####

ggsave(g_br,
       filename = "output/fig_ppt_theo_br.pdf",
       width = 5,
       height = 4.5)

ggsave(g_size,
       filename = "output/fig_ppt_theo_size.pdf",
       width = 5,
       height = 4.5)

# figure empirical --------------------------------------------------------

## - sourced `format_data4jags.R` in the following script
source("code/format_emp_est2figure.R")

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

## fcl vs. watershed area
(g_size <- df_mu_est %>% 
    ggplot(aes(x = r_length,
               y = fcl_est)) +
    geom_point(aes(color = factor(h),
                   size = scl_w),
               alpha = 0.5) + 
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10() +
    labs(y = "Food chain length",
         x = "Total river length (km)") +
    guides(color = "none",
           size = "none"))

## fcl vs. branching prob
(g_b <- df_mu_est %>% 
    ggplot(aes(x = lambda,
               y = fcl_est)) +
    geom_point(aes(color = factor(h),
                   size = scl_w),
               alpha = 0.4) +
    geom_line(data = filter(df_yh, focus == "scl_lambda"),
              aes(x = lambda,
                  y = y,
                  color = factor(h)),
              alpha = 1,
              linetype = "dashed") + 
    geom_line(data = filter(df_y, focus == "scl_lambda"),
              aes(x = lambda,
                  y = y)) +
    geom_ribbon(data = filter(df_y, focus == "scl_lambda"),
                aes(y = y,
                    ymin = y_low,
                    ymax = y_high,
                    x = lambda),
                alpha = 0.1) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(y = "",
         x = expression("Branching rate ("*km^-1*")")) +
    guides(size = "none",
           color = "none"))

## arrange ################################################################

g_comb <- g_size + g_b

ggsave(g_comb,
       filename = "output/fig_ppt_emp_fcl.pdf",
       width = 8,
       height = 4)

