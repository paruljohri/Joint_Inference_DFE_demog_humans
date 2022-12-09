#!/bin/env Rscript
#This is the final list of stats:
#chmod +x get_final_statistics.R
#Rscript ./get_final_divergence_phase3.R exon

#options("winSize" = ()$win_size)
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)
region <- args[1]
cat <- args[2]

t <- read.table(paste("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/YRI_50/by_category/divergence_", region, "_filtered_", cat, ".txt", sep=""),h=T)
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



write.table(stats, file=paste("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/YRI_50/STATS_FILTERED/by_category/human_phase3_bigwindow_", region, "_", cat, ".divergence", sep=""), sep="\t", quote=F, row.names=F)

print ("done")


