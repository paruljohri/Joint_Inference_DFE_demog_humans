#On my laptop:
#divergence rate is being calculated for 4Na_scaled+200 generations for each sim ID

library(abc)

setwd("Work/Projects/BgsDfeDemo_Human/ABC/")


#Inference:
#read in stats:
s_cat <- "high_div"
t_abc <- read.table(paste("demo_dfe_SingExon_human_v2_stats_filtered_5p_", s_cat, "_abc.txt", sep=""), h=T)
t_par <- t_abc[,c(2:7)]


#scale divergence:
#tau <- 6e6 #min
#tau <- 12e6 #max
tau <- (6e6 + 12e6)/2 #mean
gen_time <- 25
num_gens <- tau/gen_time
t_abc$func_div_scaled_m <- t_abc$func_div_m*((num_gens*t_abc$scaling_factor)/((4.0*t_abc$Na_scaled)+200))
t_abc$func_div_scaled_sd <- t_abc$func_div_sd*((num_gens*t_abc$scaling_factor)/((4.0*t_abc$Na_scaled)+200))**2
t_abc$linked1_div_scaled_m <- t_abc$linked1_div_m*((num_gens*t_abc$scaling_factor)/((4.0*t_abc$Na_scaled)+200))
t_abc$linked1_div_scaled_sd <- t_abc$linked1_div_sd*((num_gens*t_abc$scaling_factor)/((4.0*t_abc$Na_scaled)+200))**2
t_abc$linked2_div_scaled_m <- t_abc$linked2_div_m*((num_gens*t_abc$scaling_factor)/((4.0*t_abc$Na_scaled)+200))
t_abc$linked2_div_scaled_sd <- t_abc$linked2_div_sd*((num_gens*t_abc$scaling_factor)/((4.0*t_abc$Na_scaled)+200))**2

#make the final stats matrix:
t_stats_m <- cbind(t_abc[,c(14:43)], t_abc[,c(80,82,84)])
t_stats_sd <- cbind(t_abc[,c(47:76)], t_abc[,c(81,83,85)])
t_stats <- cbind(t_stats_m, t_stats_sd)

#read in human data:
t_data <- read.table(paste("../HumanStats/YRI_50/STATS_FILTERED/by_category/human_phase3_", s_cat, "_abc.txt", sep=""), h=T)

#inference with all stats:
t_obs_5p <- t_data[1,c(2:dim(t_data)[2])]
t_obs_3p <- t_data[2,c(2:dim(t_data)[2])]

nnet_5p <- abc(target=t_obs_5p, param=t_par, sumstat=t_stats, tol=0.08, method="neuralnet")
summary(nnet_5p)

nnet_3p <- abc(target=t_obs_3p, param=t_par, sumstat=t_stats, tol=0.05, method="neuralnet")
summary(nnet_3p)







