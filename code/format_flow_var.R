
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# discharge variation -----------------------------------------------------

## - discharge values contain negative - this represents down to upstream flow
## - take absolute() to represent "habitat size" or similar

df_flow <- readRDS("data_fmt/data_env_flow.rds") %>% 
  mutate(month = month(date),
         flow = abs(discharge),
         julian = format(date, "%j") %>% 
           as.numeric()) %>% 
  group_by(sid) %>% 
  mutate(flow_norm = flow / median(flow),
         flow_median = median(flow)) %>% 
  ungroup()

v_sid <- sort(unique(df_flow$sid))

df_fsd <- foreach(i = seq_len(length(v_sid)),
                  .combine = bind_rows) %do% {
                    
                    print(i)
                    
                    ## subset data by site
                    df_i <- df_flow %>% 
                      filter(sid == v_sid[i])
                    
                    ## gam fitting for seasonality
                    ## gaussian fit
                    fit <- mgcv::gam(flow_norm ~ s(julian), data = df_i)
                    
                    ## flow SD of residuals
                    fsd <- sqrt(fit$sig2)
                    
                    ## flow inter-quantile of residuals
                    fqr <- quantile(fit$residuals, c(0.25, 0.75)) %>% 
                      diff()
                    
                    ## output set
                    cout <- tibble(sid = v_sid[i],
                                   fsd = fsd,
                                   fqr = fqr)
                    
                    return(cout)
                  }


# visualization for error check -------------------------------------------

# for (i in seq_len(length(v_sid))) {
#   
#   print(i)
#   
#   g_i <- df_flow %>% 
#     filter(sid == v_sid[i]) %>% 
#     ggplot(aes(x = julian,
#                y = flow_norm)) +
#     geom_point(alpha = 0.1,
#                size = 1) +
#     geom_smooth(color = "black") +
#     theme_classic() +
#     labs(y = "Normalized flow (Y / median(Y))",
#          x = "Julian date")
#   
#   ggsave(g_i, width = 5, height = 4,
#          filename = paste0("output/figure_flow/", v_sid[i], ".pdf"))
# }


# export ------------------------------------------------------------------

saveRDS(df_fsd, "data_fmt/data_env_flow_var.rds")
