#On my laptop:
#divergence rate is being calculated for 4Na_scaled+200 generations for each sim ID

library(abc)
library(spatstat)
library(weights)
setwd("Work/Projects/BgsDfeDemo_Human/ABC/")


t_abc <- read.table("demo_dfe_SingExon_human_v2_stats_filtered_5p_abc.txt", h=T)
time_gen <- t_abc$scaling_factor*200
t_par0 <- t_abc[,c(2:7)]
t_par <- cbind(t_par0, time_gen)

#scale divergence:
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
	nnet_5p <- abc(target=t_obs_5p, param=t_par, sumstat=t_stats, tol=0.08, method="neuralnet")
	if (i==1){
		m_adj_values <- nnet_5p$adj.values
		m_weights <- nnet_5p$weights
	}
	else{
		m_adj_values <- rbind(m_adj_values, nnet_5p$adj.values)
		m_weights <- c(m_weights, nnet_5p$weights)
	}
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



#Plot posteriors:
par(mfrow=c(4,2))
par(mar=c(4,5,2,1))

#f0
h_f0_pos <- wtd.hist(x=m_adj_values[,1], weight=m_weights, plot=F)
h_f0_pos$counts <- h_f0_pos$counts/sum(h_f0_pos$counts)
plot(h_f0_pos, main="", xlab="f0", ylab="Frequency", cex.lab=1.5, xlim=c(0,100), ylim=c(0,0.7))
h_f0_prior <- hist(t_par$f0, plot=F)
h_f0_prior$counts <- h_f0_prior$counts/sum(h_f0_prior$counts)
h_f0_prior$mids <- c(0.0, h_f0_prior$mids, 100)
h_f0_prior$counts <- c(0.0, h_f0_prior$counts, 0.0)
lines(h_f0_prior$mids, h_f0_prior$counts, type="l", lty=2)
abline(v=mean(f0), col="red")

#f1
h_f1_pos <- wtd.hist(x=m_adj_values[,2], weight=m_weights, plot=F)
h_f1_pos$counts <- h_f1_pos$counts/sum(h_f1_pos$counts)
plot(h_f1_pos, main="", xlab="f1", ylab="Frequency", cex.lab=1.5, xlim=c(0,100), ylim=c(0,0.7))
h_f1_prior <- hist(t_par$f1, plot=F)
h_f1_prior$counts <- h_f1_prior$counts/sum(h_f1_prior$counts)
h_f1_prior$mids <- c(0.0, h_f1_prior$mids, 100)
h_f1_prior$counts <- c(0.0, h_f1_prior$counts, 0.0)
lines(h_f1_prior$mids, h_f1_prior$counts, type="l", lty=2)
abline(v=mean(f1), col="red")

#f2
h_f2_pos <- wtd.hist(x=m_adj_values[,3], weight=m_weights, plot=F)
h_f2_pos$counts <- h_f2_pos$counts/sum(h_f2_pos$counts)
plot(h_f2_pos, main="", xlab="f2", ylab="Frequency", cex.lab=1.5, xlim=c(0,100), ylim=c(0,0.7))
h_f2_prior <- hist(t_par$f2, plot=F)
h_f2_prior$counts <- h_f2_prior$counts/sum(h_f2_prior$counts)
h_f2_prior$mids <- c(0.0, h_f2_prior$mids, 100)
h_f2_prior$counts <- c(0.0, h_f2_prior$counts, 0.0)
lines(h_f2_prior$mids, h_f2_prior$counts, type="l", lty=2)
abline(v=mean(f2), col="red")

#f3
h_f3_pos <- wtd.hist(x=m_adj_values[,4], weight=m_weights, plot=F)
h_f3_pos$counts <- h_f3_pos$counts/sum(h_f3_pos$counts)
plot(h_f3_pos, main="", xlab="f3", ylab="Frequency", cex.lab=1.5, xlim=c(0,100), ylim=c(0,0.7))
h_f3_prior <- hist(t_par$f3, plot=F)
h_f3_prior$counts <- h_f3_prior$counts/sum(h_f3_prior$counts)
h_f3_prior$mids <- c(0.0, h_f3_prior$mids, 100)
h_f3_prior$counts <- c(0.0, h_f3_prior$counts, 0.0)
lines(h_f3_prior$mids, h_f3_prior$counts, type="l", lty=2)
abline(v=mean(f3), col="red")

#Na
h_Na_pos <- wtd.hist(x=log(m_adj_values[,5]), weight=m_weights, plot=F)
h_Na_pos$counts <- h_Na_pos$counts/sum(h_Na_pos$counts)
plot(h_Na_pos, main="", xlab="log(Na)", ylab="Frequency", cex.lab=1.5, xlim=c(min(log(t_par$Na)),11), ylim=c(0,0.7))
h_Na_prior <- hist(log(t_par$Na), plot=F)
h_Na_prior$counts <- h_Na_prior$counts/sum(h_Na_prior$counts)
h_Na_prior$mids <- c(min(log(t_par$Na)), h_Na_prior$mids, max(log(t_par$Na)))
h_Na_prior$counts <- c(0.0, h_Na_prior$counts, 0.0)
lines(h_Na_prior$mids, h_Na_prior$counts, type="l", lty=2)
abline(v=mean(log(Na)), col="red")

#Nc
h_Nc_pos <- wtd.hist(x=log(m_adj_values[,6]), weight=m_weights, plot=F)
h_Nc_pos$counts <- h_Nc_pos$counts/sum(h_Nc_pos$counts)
plot(h_Nc_pos, main="", xlab="log(Nc)", ylab="Frequency", cex.lab=1.5, xlim=c(min(log(t_par$Nc)),11.5), ylim=c(0,0.8))
h_Nc_prior <- hist(log(t_par$Nc), plot=F)
h_Nc_prior$counts <- h_Nc_prior$counts/sum(h_Nc_prior$counts)
h_Nc_prior$mids <- c(min(log(t_par$Nc)), h_Nc_prior$mids, max(log(t_par$Nc)))
h_Nc_prior$counts <- c(0.0, h_Nc_prior$counts, 0.0)
lines(h_Nc_prior$mids, h_Nc_prior$counts, type="l", lty=2)
abline(v=mean(log(Nc)), col="red")

#time_gen
h_time_gen_pos <- wtd.hist(x=m_adj_values[,7], weight=m_weights, plot=F)
h_time_gen_pos$counts <- h_time_gen_pos$counts/sum(h_time_gen_pos$counts)
plot(h_time_gen_pos, main="", xlab="time (gen)", ylab="Frequency", cex.lab=1.5, xlim=c(min(t_par$time_gen),max(t_par$time_gen)), ylim=c(0,0.7))
h_time_gen_prior <- hist(t_par$time_gen, plot=F)
h_time_gen_prior$counts <- h_time_gen_prior$counts/sum(h_time_gen_prior$counts)
h_time_gen_prior$mids <- c(min(t_par$time_gen), h_time_gen_prior$mids, max(t_par$time_gen))
h_time_gen_prior$counts <- c(0.0, h_time_gen_prior$counts, 0.0)
lines(h_time_gen_prior$mids, h_time_gen_prior$counts, type="l", lty=2)
abline(v=mean(time_gen), col="red")










