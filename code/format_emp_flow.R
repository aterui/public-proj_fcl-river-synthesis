
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# discharge variation -----------------------------------------------------

## - discharge values contain negative - this represents down to upstream flow
## - take absolute() to represent "habitat size" or similar

df_flow <- readRDS("data_raw/data_env_flow.rds") %>% 
  filter(discharge != 0) %>% 
  mutate(month = month(date),
         flow = abs(discharge),
         log_flow = log(flow, 10),
         julian = format(date, "%j") %>% 
           as.numeric()) %>% 
  mutate(log_flow_s = if (min(log_flow) < 0) log_flow + abs(min(log_flow))) %>% 
  group_by(sid) %>% 
  mutate(norm_log_flow = log_flow_s / mean(log_flow_s)) %>% 
  ungroup()

v_sid <- sort(unique(df_flow$sid))

cl <- makeCluster(detectCores() - 2)
registerDoSNOW(cl)

df_fsd <- foreach(i = seq_len(length(v_sid)),
                  .combine = bind_rows,
                  .packages = c("tidyverse",
                                "mgcv")) %dopar% {
                    
                    ## subset data by site
                    df_i <- df_flow %>%
                      filter(sid == v_sid[i]) %>%
                      arrange(date) %>%
                      mutate(year = format(date, "%Y"))

                    # ## gam fitting for seasonality
                    # ## error normal
                    # fit <- mgcv::gam(norm_log_flow ~ s(julian),
                    #                  data = df_i)
                    # 
                    # ## get residual SD
                    # sigma0 <- sqrt(fit$sig2)
                    # 
                    # ## output set
                    # cout <- tibble(sid = v_sid[i],
                    #                fsd = sigma0,
                    #                fmu = mean(df_i$log_flow))
                    
                    ## gam fitting for seasonality
                    ## error structure student-t, with minimum d.f. = 2
                    fit <- mgcv::gam(norm_log_flow ~ s(julian),
                                     family = mgcv::scat(min.df = 2),
                                     data = df_i)

                    ## degree of freedom in student-t
                    nu <- with(family(fit), getTheta(trans = TRUE))[1]

                    ## SD in student-t
                    sigma0 <- with(family(fit), getTheta(trans = TRUE))[2]

                    ## output set
                    cout <- tibble(sid = v_sid[i],
                                   fsd = sigma0,
                                   fnu = nu,
                                   fmu = mean(df_i$log_flow))
                    
                    return(cout)
                  }

stopCluster(cl)

# visualization for error check -------------------------------------------

# for (i in seq_len(length(v_sid))) {
# 
#   print(i)
# 
#   g_i <- df_flow %>%
#     filter(sid == v_sid[i]) %>%
#     ggplot(aes(x = julian,
#                y = norm_log_flow)) +
#     geom_point(alpha = 0.1,
#                size = 1) +
#     geom_smooth(color = "black") +
#     theme_classic() +
#     labs(y = "Normalized flow (Y / mean(Y))",
#          x = "Julian date")
# 
#   ggsave(g_i, width = 5, height = 4,
#          filename = paste0("tools/figure_flow/", v_sid[i], ".pdf"))
# }


# error correction --------------------------------------------------------

## sid: 080_itakura_s005
## - this site seems misassined to an adjacent small tributary
## - clearly deviate from a general trend
## - assign values of 080_itakura_s006 (an upstream site of the same river)
val <- df_fsd[str_detect(df_fsd$sid, "080_itakura_s006"), c("fsd", "fnu", "fmu")]
df_fsd[str_detect(df_fsd$sid, "080_itakura_s005"), c("fsd", "fnu", "fmu")] <- val

# export ------------------------------------------------------------------

saveRDS(df_fsd, "data_fmt/data_env_fsd.rds")
