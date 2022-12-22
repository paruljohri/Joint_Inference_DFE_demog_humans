#!/bin/env Rscript
#This is the final list of stats:
#chmod +x get_final_statistics.R
#Rscript ./get_final_divergence_phase3.R YRI_50 exon

#options("winSize" = ()$win_size)
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)
s_popn <- args[1] #YRI/YRI_50
region <- args[2] #exon/5p/3p

t <- read.table(paste("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/", s_popn, "/divergence_", region, "_filtered.txt", sep=""),h=T)
t_sub <- t[, c(2:length(t))]
stats <- c()
col_names <- c()
n <- length(t_sub)
i <- 1
while (i <= n){
    stats <- cbind(stats, mean(t_sub[,i], na.rm=T), sd(t_sub[,i], na.rm=T))
    col_names <- c(col_names, paste(colnames(t_sub)[i], "_m", sep=""), paste(colnames(t_sub)[i], "_sd",sep=""))
    i <- i + 1
    }
#colMeans(t_sub, na.rm=T)
colnames(stats) <- col_names



write.table(stats, file=paste("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/", s_popn, "/STATS_FILTERED/human_phase3_bigwindow_", region, ".divergence", sep=""), sep="\t", quote=F, row.names=F)

print ("done")


