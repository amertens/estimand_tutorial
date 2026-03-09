
library(data.table)
library(dplyr)
library(tidyr)
library(lmtp)
library(SuperLearner)
library(furrr)
library(here)
library(MatchIt)
library(survival)
library(broom)

source("DGP.R")

set.seed(12345)


df_complex_A1 <- generate_hcv_data(        N               = 1000000,
                                           np_hazard     = TRUE,
                                    dep_censor    = TRUE,
                                    complexity    = TRUE,
                                    switch_on = FALSE,
                                    censor_base = 0.00000000000000000000000000001, #need to update to allow no cens
                                    treat_override  = "all_treated",
                                    seed=1234)

df_complex_A0 <- generate_hcv_data(     N               = 1000000,
                                        np_hazard     = TRUE,
                                       dep_censor    = TRUE,
                                       complexity    = TRUE,
                                       switch_on = FALSE,
                                       censor_base = 0.00000000000000000000000000001,
                                       treat_override  = "all_control",
                                       seed=1234)

mean(df_complex_A1$event) *100

#CI at 180 days:
df_complex_A1$Y180 <- df_complex_A1$event
df_complex_A0$Y180 <- df_complex_A0$event
df_complex_A1$Y180[df_complex_A1$follow_time>180] <- 0
df_complex_A0$Y180[df_complex_A0$follow_time>180] <- 0

mean(df_complex_A1$event)
(mean(df_complex_A1$Y180) - mean(df_complex_A0$Y180))*100
(mean(df_complex_A1$Y180) / mean(df_complex_A0$Y180))

#true CID: 0.0002619194
#true CIR: 1.011929


saveRDS(all_contrasts, file = here("results","sim_results", "lmtp_sim_all.rds"))

lmtp_res= readRDS(here("results","sim_results", "lmtp_sim_all.rds"))
ps_res = read.csv(here("results","summary_ps.csv"))
mean(ps_res$cox_hr)


temp=lmtp_res[[1]]
class(temp)
temp$estimates

for(i in 1:length(lmtp_res)){
  if (i==1){
    all_contrasts <- lmtp_res[[i]]$estimates
  } else {
    all_contrasts <- rbind(all_contrasts, lmtp_res[[i]]$estimates)
  }
}

mean(all_contrasts$estimate)
