#On my laptop:
#divergence rate is being calculated for 4Na_scaled+200 generations for each sim ID

library(abc)
library(spatstat)

setwd("Work/Projects/BgsDfeDemo_Human/ABC/")

#stats without filtering
t_abc <- read.table("demo_dfe_SingExon_human_v2_stats_abc.txt", h=T)

#5prime, filtered stats:
t_abc <- read.table("demo_dfe_SingExon_human_v2_stats_filtered_5p_abc.txt", h=T)

#get time and make stats:
time_gen <- t_abc$scaling_factor*200
t_par <- cbind(t_abc[,c(2:7)],time_gen)
t_stats <- t_abc[,c(14:79)]

#cross validation:

cv_ridge <- cv4abc(t_par, t_stats, nval=100, tols=0.05, statistic="median", method="ridge", transf="none")
summary(cv_ridge)
plot(cv_ridge)

cv_nnet <- cv4abc(t_par, t_stats,nval=100, tols=c(0.05, 0.08, 0.1), statistic="median", method="neuralnet", transf="none")
summary(cv_nnet)


#5prime:
t_abc <- read.table("demo_dfe_SingExon_human_v2_stats_filtered_3p_abc.txt", h=T)
time_gen <- t_abc$scaling_factor*200
t_par <- cbind(t_abc[,c(2:7)],time_gen)
t_stats <- t_abc[,c(14:79)]

#cross validation:
cv_nnet <- cv4abc(t_par, t_stats,nval=100, tols=c(0.05, 0.08, 0.1), statistic="median", method="neuralnet", transf="none")
summary(cv_nnet)

#Plotting:
par(mfrow=c(2,4))
par(mar=c(4,4,2,1))
plot(cv_nnet)
>> save as 6 x 11.69 (landscape)
