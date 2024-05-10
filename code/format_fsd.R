
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/set_library.R")
source("code/set_function.R")


# discharge variation -----------------------------------------------------

df_flow <- readRDS("data_fmt/data_env_flow.rds") %>% 
  mutate(month = month(date),
         discharge = ifelse(discharge <= 0, 0, discharge),
         log_flow = log(discharge + 1),
         norm_flow = log_flow / mean(log_flow),
         julian = format(date, "%j") %>% 
           as.numeric())

v_sid <- sort(unique(df_flow$sid))

df_fsd <- foreach(i = seq_len(length(v_sid)),
                  .combine = bind_rows) %do% {
                    print(i)
                    df_i <- df_flow %>% 
                      filter(sid == v_sid[i])
                    
                    fit <- mgcv::gam(norm_flow ~ s(julian), data = df_i)
                    cout <- tibble(sid = v_sid[i], fsd = sd(fit$residuals))
                    
                    return(cout)
                  }


# export ------------------------------------------------------------------

saveRDS(df_fsd, "data_fmt/data_env_flowsd.rds")
