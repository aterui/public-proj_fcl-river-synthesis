
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# discharge variation -----------------------------------------------------

## - discharge values contain negative - this represents down to upstream flow
## - take absolute() to represent "habitat size" or similar

df_flow <- readRDS("data_raw/data_env_flow.rds") %>% 
  mutate(month = month(date),
         flow = abs(discharge) + 1,
         log_flow = log(flow, 10),
         julian = format(date, "%j") %>% 
           as.numeric()) %>% 
  group_by(sid) %>% 
  mutate(norm_log_flow = log_flow / mean(log_flow)) %>% 
  ungroup()

v_sid <- sort(unique(df_flow$sid))

cl <- makeCluster(detectCores() - 2)
registerDoSNOW(cl)

df_fsd <- foreach(i = seq_len(length(v_sid)),
                  .combine = bind_rows,
                  .packages = c("tidyverse",
                                "mgcv")) %dopar% {
                    
                    print(i)
                    
                    ## subset data by site
                    df_i <- df_flow %>% 
                      filter(sid == v_sid[i])
                    
                    fit <- mgcv::gam(norm_log_flow ~ s(julian),
                                     data = df_i)
                    
                    sigma0 <- sqrt(fit$sig2)
                    
                    # ## gam fitting for seasonality
                    # ## error structure student-t, with minimum d.f. = 2
                    # fit <- mgcv::gam(norm_log_flow ~ s(julian),
                    #                  family = mgcv::scat(min.df = 2),
                    #                  data = df_i)
                    # 
                    # ## degree of freedom in student-t
                    # nu <- with(family(fit), getTheta(trans = TRUE))[1]
                    # 
                    # ## SD in student-t
                    # sigma0 <- with(family(fit), getTheta(trans = TRUE))[2]
                    
                    ## output set
                    cout <- tibble(sid = v_sid[i],
                                   fsd = sigma0)
                    
                    return(cout)
                  }

stopCluster(cl)

# visualization for error check -------------------------------------------

# for (i in seq_len(length(v_sid))) {
# 
#   print(i)
# 
#   g_i <- df_flow %>%
#     filter(sid == "003_jackson_s007") %>%
#     ggplot(aes(x = julian,
#                y = log_flow)) +
#     geom_point(alpha = 0.1,
#                size = 1) +
#     geom_smooth(color = "black") +
#     theme_classic() +
#     labs(y = "Normalized flow (Y / mean(Y))",
#          x = "Julian date")
# 
#   ggsave(g_i, width = 5, height = 4,
#          filename = paste0("output/figure_flow/", v_sid[i], ".pdf"))
# }


# export ------------------------------------------------------------------

saveRDS(df_fsd, "data_fmt/data_env_fsd.rds")
