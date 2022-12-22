#To plot the fit of the best model:

setwd("/Users/paruljohri/Work/Projects/BgsDfeDemo_Human")

t_model_5p <- read.table("model_fit_stats_filtered_5p/sim1_bigwindow.stats", h=T)

#read data:
t_5p_stats <- read.table("HumanStats/YRI_50/STATS_FILTERED/human_phase3_bigwindow_5p.stats", h=T)
t_exon_stats <- read.table("HumanStats/YRI_50/STATS_FILTERED/human_phase3_bigwindow_exon.stats", h=T)
t_3p_stats <- read.table("HumanStats/YRI_50/STATS_FILTERED/human_phase3_bigwindow_3p.stats", h=T)

t_exon_div <- read.table("HumanStats/divergence/YRI_50/divergence_exon_filtered.txt", h=T)
t_5p_div <- read.table("HumanStats/divergence/YRI_50/divergence_5p_filtered.txt", h=T)
t_3p_div <- read.table("HumanStats/divergence/YRI_50/divergence_3p_filtered.txt", h=T)

col_model <- rgb(1,0,0,0.5)

#plot for each stat separately:
#This is for thetapi, thetaw, thetah, and numSing
stat <- "numSing"
if (stat=="thetapi"){xlab_name <- expression(pi)} #don't replace these
if (stat=="thetaw"){xlab_name <- expression(theta[w])} #don't replace these
if (stat=="thetah"){xlab_name <- expression(theta[H])}#don't repalce these
if (stat=="numSing"){xlab_name <- "singleton density"}#don't replace these
pdf(file=paste("/Users/paruljohri/Work/Projects/BgsDfeDemo_Human/plots/model_fit/sim1_ten_reps/numSing.pdf", sep=""), width=8.27, height=3.0)
par(mfrow=c(1,3))
#func:
t_model <- t_model_5p$numSing[which(t_model_5p$WinType=="func")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="func")]
t_data <- t_exon_stats$numSing/t_exon_stats$WinSize_filtered
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
hist_model1 <- hist(t_model1, plot=F)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="functional", xlab=xlab_name, cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(t_model1,t_data1, na.rm=T),max(t_model1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked1
t_model <- t_model_5p$numSing[which(t_model_5p$WinType=="linked1")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="linked1")]
t_data <- t_5p_stats$numSing[which(t_5p_stats$WinType=="linked1")]/t_5p_stats$WinSize_filtered[which(t_5p_stats$WinType=="linked1")]
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
hist_model1 <- hist(t_model1, plot=F)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="linked", xlab=xlab_name, cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(t_model1,t_data1, na.rm=T),max(t_model1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked2
t_model <- t_model_5p$numSing[which(t_model_5p$WinType=="linked2")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="linked2")]
t_data <- t_5p_stats$numSing[which(t_5p_stats$WinType=="linked2")]/t_5p_stats$WinSize_filtered[which(t_5p_stats$WinType=="linked2")]
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
hist_model1 <- hist(t_model1, plot=F)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="less linked", xlab=xlab_name, cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(t_model1,t_data1, na.rm=T),max(t_model1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)
dev.off()

>> save as
A4 (landscape)
8.27 x 3

#This is for plotting log of thetapi, thetaw, thetah, and numSing
stat <- "thetapi"
if (stat=="thetapi"){xlab_name <- expression(paste("ln(", pi, ")"))} #don't replace these
if (stat=="thetaw"){xlab_name <- expression(paste("ln(", theta[w], ")"))} #don't replace these
if (stat=="thetah"){xlab_name <- expression(paste("ln(", theta[H], ")"))} #don't replace these
if (stat=="numSing"){xlab_name <- "ln(singleton density)"} #don't replace these
pdf(file=paste("/Users/paruljohri/Work/Projects/BgsDfeDemo_Human/plots/model_fit/sim1_ten_reps/log_thetapi.pdf", sep=""), width=8.27, height=3.0)
par(mfrow=c(1,3))
#func:
t_model <- t_model_5p$thetapi[which(t_model_5p$WinType=="func")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="func")]
t_data <- t_exon_stats$thetapi/t_exon_stats$WinSize_filtered
t_model1 <- log(t_model[which(is.finite(t_model))])
t_data1 <- log(t_data[which(is.finite(t_data))])
hist_model1 <- hist(t_model1, plot=F)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="functional", xlab=xlab_name, cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(hist_model1$breaks,hist_data1$breaks, na.rm=T),max(hist_model1$breaks,hist_data1$breaks, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked1
t_model <- t_model_5p$thetapi[which(t_model_5p$WinType=="linked1")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="linked1")]
t_data <- t_5p_stats$thetapi[which(t_5p_stats$WinType=="linked1")]/t_5p_stats$WinSize_filtered[which(t_5p_stats$WinType=="linked1")]
t_model1 <- log(t_model[which(is.finite(t_model))])
t_data1 <- log(t_data[which(is.finite(t_data))])
hist_model1 <- hist(t_model1, plot=F)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="linked", xlab=xlab_name, cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(hist_model1$breaks,hist_data1$breaks, na.rm=T),max(hist_model1$breaks,hist_data1$breaks, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked2
t_model <- t_model_5p$thetapi[which(t_model_5p$WinType=="linked2")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="linked2")]
t_data <- t_5p_stats$thetapi[which(t_5p_stats$WinType=="linked2")]/t_5p_stats$WinSize_filtered[which(t_5p_stats$WinType=="linked2")]
t_model1 <- log(t_model[which(is.finite(t_model))])
t_data1 <- log(t_data[which(is.finite(t_data))])
hist_model1 <- hist(t_model1, plot=F)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="less linked", xlab=xlab_name, cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(hist_model1$breaks,hist_data1$breaks, na.rm=T),max(hist_model1$breaks,hist_data1$breaks, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)
dev.off()

>> save as
A4 (landscape)
8.27 x 3


#This is for hprime, tajimasd, hapdiv, rsq, Dprime
stat <- "hprime"
if (stat=="hprime"){xlab_name <- expression(paste(italic(H), "'"))} #don't replace these
if (stat=="tajimasd"){xlab_name <- expression(paste("Tajima's ", italic(D)))} #don't replace these
if (stat=="hapdiv"){xlab_name <- "haplotype diversity"} #don't replace these
if (stat=="rsq"){xlab_name <- expression('r'^2)} #don't replace these
if (stat=="Dprime"){xlab_name <- expression(paste(italic(D), "'"))} #don't replace these
pdf(file=paste("/Users/paruljohri/Work/Projects/BgsDfeDemo_Human/plots/model_fit/sim1_ten_reps/hprime.pdf", sep=""), width=8.27, height=3.0)
par(mfrow=c(1,3))
#func:
t_model <- t_model_5p$hprime[which(t_model_5p$WinType=="func")]
t_data <- t_exon_stats$hprime
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
hist_model1 <- hist(t_model1, plot=F)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="functional", xlab=xlab_name, cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(t_model1,t_data1, na.rm=T),max(t_model1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked1
t_model <- t_model_5p$hprime[which(t_model_5p$WinType=="linked1")]
t_data <- t_5p_stats$hprime[which(t_5p_stats$WinType=="linked1")]
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
hist_model1 <- hist(t_model1, plot=F)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="linked", xlab=xlab_name, cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(t_model1,t_data1, na.rm=T),max(t_model1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked2
t_model <- t_model_5p$hprime[which(t_model_5p$WinType=="linked2")]
t_data <- t_5p_stats$hprime[which(t_5p_stats$WinType=="linked2")]
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
hist_model1 <- hist(t_model1, plot=F)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="less linked", xlab=xlab_name, cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(t_model1,t_data1, na.rm=T),max(t_model1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)
dev.off()

>> save as
A4 (landscape)
8.27 x 3


#This is for divergence per site
tau <- (6e6 + 12e6)/2 #mean
gen_time <- 25
num_gens <- tau/gen_time
Nanc <- 7509
time_of_change <- 332
pdf(file=paste("/Users/paruljohri/Work/Projects/BgsDfeDemo_Human/plots/model_fit/sim1_ten_reps/divergence.pdf", sep=""), width=8.27, height=3.0)
par(mfrow=c(1,3))
#func:
t_model <- t_model_5p$div[which(t_model_5p$WinType=="func")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="func")]
t_data <- t_exon_div$func_hg19.anc
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
#transform divergence:
t_model2 <- t_model1*((num_gens)/((4.0*Nanc)+time_of_change))
hist_model1 <- hist(t_model2, plot=F, breaks=100)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="functional", xlab="divergence per site", cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(t_model1,t_data1, na.rm=T),max(t_model1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked1
t_model <- t_model_5p$div[which(t_model_5p$WinType=="linked1")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="linked1")]
t_data <- t_5p_div$linked1_hg19.anc
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
#transform divergence:
t_model2 <- t_model1*((num_gens)/((4.0*Nanc)+time_of_change))
hist_model1 <- hist(t_model2, plot=F, breaks=100)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="linked", xlab="divergence per site", cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(t_model1,t_data1, na.rm=T),max(t_model1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked2
t_model <- t_model_5p$div[which(t_model_5p$WinType=="linked2")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="linked2")]
t_data <- t_5p_div$linked2_hg19.anc
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
#transform divergence:
t_model2 <- t_model1*((num_gens)/((4.0*Nanc)+time_of_change))
hist_model1 <- hist(t_model2, plot=F, breaks=100)
hist_data1 <- hist(t_data1, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="less linked", xlab="divergence per site", cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(t_model1,t_data1, na.rm=T),max(t_model1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)
dev.off()


#This log of divergence
tau <- (6e6 + 12e6)/2 #mean
gen_time <- 25
num_gens <- tau/gen_time
Nanc <- 7509
time_of_change <- 332
pdf(file=paste("/Users/paruljohri/Work/Projects/BgsDfeDemo_Human/plots/model_fit/sim1_ten_reps/log_divergence.pdf", sep=""), width=8.27, height=3.0)
par(mfrow=c(1,3))
#func:
t_model <- t_model_5p$div[which(t_model_5p$WinType=="func")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="func")]
t_data <- t_exon_div$func_hg19.anc
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
#transform divergence:
t_model2 <- log(t_model1*((num_gens)/((4.0*Nanc)+time_of_change)))
t_data2 <- log(t_data1)
hist_model1 <- hist(t_model2, plot=F)
hist_data1 <- hist(t_data2, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="functional", xlab="ln(divergence per site)", cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(hist_model1$breaks,hist_data1$breaks, na.rm=T),max(hist_model1$breaks,hist_data1$breaks, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked1
t_model <- t_model_5p$div[which(t_model_5p$WinType=="linked1")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="linked1")]
t_data <- t_5p_div$linked1_hg19.anc
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
#transform divergence:
t_model2 <- log(t_model1*((num_gens)/((4.0*Nanc)+time_of_change)))
t_data2 <- log(t_data1)
hist_model1 <- hist(t_model2, plot=F)
hist_data1 <- hist(t_data2, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="linked", xlab="ln(divergence per site)", cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(hist_model1$breaks,hist_data1$breaks, na.rm=T),max(hist_model1$breaks,hist_data1$breaks, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)

#linked2
t_model <- t_model_5p$div[which(t_model_5p$WinType=="linked2")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="linked2")]
t_data <- t_5p_div$linked2_hg19.anc
t_model1 <- t_model[which(is.finite(t_model))]
t_data1 <- t_data[which(is.finite(t_data))]
#transform divergence:
t_model2 <- log(t_model1*((num_gens)/((4.0*Nanc)+time_of_change)))
t_data2 <- log(t_data1)
hist_model1 <- hist(t_model2, plot=F)
hist_data1 <- hist(t_data2, plot=F)
hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
plot(hist_model1, freq=FALSE, main="less linked", xlab="ln(divergence per site)", cex.lab=1.5, border=col_model, col=col_model, xlim=c(min(hist_model1$breaks,hist_data1$breaks, na.rm=T),max(hist_model1$breaks,hist_data1$breaks, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_data1$density)))
plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)
dev.off()

>> save as
A4 (landscape)
8.27 x 3