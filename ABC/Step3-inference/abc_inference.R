#On my laptop:
#divergence rate is being calculated for 4Na_scaled+200 generations for each sim ID

library(abc)
library(spatstat)
setwd("Work/Projects/BgsDfeDemo_Human/ABC/")


#Inference:
#with time:
t_abc <- read.table("demo_dfe_SingExon_human_v2_stats_filtered_5p_abc.txt", h=T)
time_gen <- t_abc$scaling_factor*200
t_par0 <- t_abc[,c(2:7)]
t_par <- cbind(t_par0, time_gen)

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
t_data <- read.table("../HumanStats/YRI_50/STATS_FILTERED/human_phase3_abc.txt", h=T)

#inference with all stats:
t_obs_5p <- t_data[1,c(2:dim(t_data)[2])]
t_obs_3p <- t_data[2,c(2:dim(t_data)[2])]

nnet_5p <- abc(target=t_obs_5p, param=t_par, sumstat=t_stats, tol=0.08, method="neuralnet")
summary(nnet_5p)

nnet_3p <- abc(target=t_obs_3p, param=t_par, sumstat=t_stats, tol=0.05, method="neuralnet")
summary(nnet_3p)


#Inference performed multiple times and then taking a mean:
num_reps <- 50
f0 <- c() 
f1 <- c()
f2 <- c()
f3 <- c()
Na <- c() 
Nc <- c()
time_gen <- c()
i = 1
while(i<=num_reps){
	print (i)
	nnet_5p <- abc(target=t_obs_5p, param=t_par, sumstat=t_stats, tol=0.1, method="neuralnet")
	#calculating errors:
	est1 <- weighted.median(nnet_5p$adj.values[,1], nnet_5p$weights, na.rm = TRUE)
	est2 <- weighted.median(nnet_5p$adj.values[,2], nnet_5p$weights, na.rm = TRUE)
	est3 <- weighted.median(nnet_5p$adj.values[,3], nnet_5p$weights, na.rm = TRUE)
	est4 <- weighted.median(nnet_5p$adj.values[,4], nnet_5p$weights, na.rm = TRUE)
	est5 <- weighted.median(nnet_5p$adj.values[,5], nnet_5p$weights, na.rm = TRUE)
	est6 <- weighted.median(nnet_5p$adj.values[,6], nnet_5p$weights, na.rm = TRUE)
	est7 <- weighted.median(nnet_5p$adj.values[,7], nnet_5p$weights, na.rm = TRUE)
	if (est1 < 0){f0 <- c(f0, 0.0)}
	else{f0 <- c(f0, est1)}
	if (est2 < 0){f1 <- c(f1, 0.0)}
	else{f1 <- c(f1, est2)}
	if (est3 < 0){f2 <- c(f2, 0.0)}
	else{f2 <- c(f2, est3)}
	if (est4 < 0){f3 <- c(f3, 0.0)}
	else{f3 <- c(f3, est4)}
	if (est5 < 0){Na <- c(Na, 0.0)}
	else{Na <- c(Na, est5)}
	if (est6 < 0){Nc <- c(Nc, 0.0)}
	else{Nc <- c(Nc, est6)}
	if (est7 < 0){time_gen <- c(time_gen, 0.0)}
	else{time_gen <- c(time_gen, est7)}
	i = i + 1
	}	
mean(f0)/(mean(f0) + mean(f1) + mean(f2) + mean(f3))
mean(f1)/(mean(f0) + mean(f1) + mean(f2) + mean(f3))
mean(f2)/(mean(f0) + mean(f1) + mean(f2) + mean(f3))
mean(f3)/(mean(f0) + mean(f1) + mean(f2) + mean(f3))
mean(Na)
mean(Nc)
mean(time_gen)

