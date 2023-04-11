#On my laptop:
#divergence rate is being calculated for 4Na_scaled+200 generations for each sim ID

library(abc)
library(spatstat)

setwd("/Users/paruljohri/Work/Projects/BgsDfeDemo_Human/ABC")


#with time:
t_abc <- read.table("demo_dfe_SingExon_human_v2_stats_filtered_5p_abc.txt", h=T)
time_gen <- t_abc$scaling_factor*200
t_par0 <- t_abc[,c(2:7)]
t_par <- cbind(t_par0, time_gen)

#make the final stats matrix:
t_stats <- t_abc[,c(14:79)]

#With and without mutation rate variation:
t_aa <- read.table("ModelViolations_aa_stats_filtered_5p_abc.txt", h=T)

#inference with all stats:
#Note that only sim2, sim4 and sim6 are relevant here as they are with constant mutation rates.
simID <- 6 #2/4/6
t_obs_5p <- t_aa[simID,c(17:dim(t_aa)[2])]

nnet_5p <- abc(target=t_obs_5p, param=t_par, sumstat=t_stats, tol=0.08, method="neuralnet")
summary(nnet_5p)




#Inference performed multiple times and then taking a mean:
#To get absolute error-
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
	nnet_5p <- abc(target=t_obs_5p, param=t_par, sumstat=t_stats, tol=0.08, method="neuralnet")
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

f0_norm <- f0/(f0 + f1 + f2 + f3)
f1_norm <- f1/(f0 + f1 + f2 + f3)
f2_norm <- f2/(f0 + f1 + f2 + f3)
f3_norm <- f3/(f0 + f1 + f2 + f3)

mean(f0_norm)
mean(f1_norm)
mean(f2_norm)
mean(f3_norm)
mean(Na)
mean(Nc)
mean(time_gen)

sd(f0_norm)
sd(f1_norm)
sd(f2_norm)
sd(f3_norm)
sd(Na)
sd(Nc)
sd(time_gen)

