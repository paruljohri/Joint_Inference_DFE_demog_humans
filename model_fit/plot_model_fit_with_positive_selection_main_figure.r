#To plot the fit of the best model:

setwd("/Users/paruljohri/Work/Projects/BgsDfeDemo_Human")

#read the best-fit model
t_model_5p <- read.table("model_fit_stats_filtered_5p/sim1_bigwindow.stats", h=T)

#read data:
t_exon_stats <- read.table("HumanStats/YRI_50/STATS_FILTERED/human_phase3_bigwindow_exon.stats", h=T)
t_exon_div <- read.table("HumanStats/divergence/YRI_50/divergence_exon_filtered.txt", h=T)

#red blue color scheme:
col_model <- rgb(1,0,0,0.5)
col_pos <- rgb(0.0,0.0,1.0,0.5)


#infrequent and weak positive selection (sim2)
scenarios <- c(2,6,10)
pdf(file=paste("/Users/paruljohri/Work/Projects/BgsDfeDemo_Human/plots/model_fit/sim2_6_10_td_rsq_div.pdf", sep=""), width=8.27, height=8.5)
par(mfrow=c(3,3))#fill by column
for (simID in scenarios){
	#read data:
	t_pos_5p <- read.table(paste("model_fit_stats_filtered_5p/sim", simID, "_bigwindow.stats", sep=""), h=T)
		
	#tajimasd
	xlab_name <- expression(paste("Tajima's ", italic(D))) #don't replace these
	#func:
	t_model <- t_model_5p$tajimasd[which(t_model_5p$WinType=="func")]
	t_pos <- t_pos_5p$tajimasd[which(t_pos_5p$WinType=="func")]
	t_data <- t_exon_stats$tajimasd
	t_model1 <- t_model[which(is.finite(t_model))]
	t_pos1 <- t_pos[which(is.finite(t_pos))]
	t_data1 <- t_data[which(is.finite(t_data))]
	hist_model1 <- hist(t_model1, plot=F)
	hist_pos1 <- hist(t_pos1, plot=F)
	hist_data1 <- hist(t_data1, plot=F)
	hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
	hist_pos1$density <- hist_pos1$counts/sum(hist_pos1$counts)
	hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
	plot(hist_model1, freq=FALSE, main="", xlab=xlab_name, cex.lab=1.5, cex.axis=1.2, border=col_model, col=col_model, xlim=c(min(t_model1,t_pos1,t_data1, na.rm=T),max(t_model1,t_pos1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_pos1$density, hist_data1$density)))
	plot(hist_pos1, freq=FALSE, add=T, border=col_pos, col=col_pos)
	plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)
	
	#rsq
	xlab_name <- expression('r'^2)
	#func:
	t_model <- t_model_5p$rsq[which(t_model_5p$WinType=="func")]
	t_pos <- t_pos_5p$rsq[which(t_pos_5p$WinType=="func")]
	t_data <- t_exon_stats$rsq
	t_model1 <- t_model[which(is.finite(t_model))]
	t_pos1 <- t_pos[which(is.finite(t_pos))]
	t_data1 <- t_data[which(is.finite(t_data))]
	hist_model1 <- hist(t_model1, plot=F)
	hist_pos1 <- hist(t_pos1, plot=F)
	hist_data1 <- hist(t_data1, plot=F)
	hist_model1$density <- hist_model1$counts/sum(hist_model1$counts)
	hist_pos1$density <- hist_pos1$counts/sum(hist_pos1$counts)
	hist_data1$density <- hist_data1$counts/sum(hist_data1$counts)
	plot(hist_model1, freq=FALSE, main="", xlab=xlab_name, cex.lab=1.5, cex.axis=1.2, border=col_model, col=col_model, xlim=c(min(t_model1,t_pos1,t_data1, na.rm=T),max(t_model1,t_pos1,t_data1, na.rm=T)), ylim=c(0,max(hist_model1$density, hist_pos1$density, hist_data1$density)))
	plot(hist_pos1, freq=FALSE, add=T, border=col_pos, col=col_pos)
	plot(hist_data1, freq=FALSE, add=T, border="black", col=NULL)
	
	#log of divergence per site
	tau <- (6e6 + 12e6)/2 #mean
	gen_time <- 25
	num_gens <- tau/gen_time
	Nanc <- 7509
	time_of_change <- 332
	#func:
	t_model <- t_model_5p$div[which(t_model_5p$WinType=="func")]/t_model_5p$WinSize_filtered[which(t_model_5p$WinType=="func")]
	t_pos <- t_pos_5p$div[which(t_pos_5p$WinType=="func")]/t_pos_5p$WinSize_filtered[which(t_pos_5p$WinType=="func")]
	t_data <- t_exon_div$func_hg19.anc
	t_model1 <- t_model[which(is.finite(t_model))]
	t_pos1 <- t_pos[which(is.finite(t_pos))]
	t_data1 <- t_data[which(is.finite(t_data))]
	#transform divergence:
	t_model2 <- log(t_model1*((num_gens)/((4.0*Nanc)+time_of_change)))
	t_pos2 <- log(t_pos1*((num_gens)/((4.0*Nanc)+time_of_change)))
	t_data2 <- log(t_data1)
	hist_model <- hist(t_model2, plot=F)
	hist_pos <- hist(t_pos2, plot=F)
	hist_data <- hist(t_data2, plot=F)
	hist_model$density <- hist_model$counts/sum(hist_model$counts)
	hist_pos$density <- hist_pos$counts/sum(hist_pos$counts)
	hist_data$density <- hist_data$counts/sum(hist_data$counts)
	plot(hist_model, freq=FALSE, main="", xlab="ln(divergence per site)", cex.lab=1.5, cex.axis=1.2, border=col_model, col=col_model, xlim=c(min(hist_model$breaks,hist_pos$breaks,hist_data$breaks, na.rm=T),max(hist_model$breaks,hist_pos$breaks,hist_data$breaks, na.rm=T)), ylim=c(0,max(hist_model$density, hist_pos$density, hist_data$density)))
	plot(hist_pos, freq=FALSE, add=T, border=col_pos, col=col_pos)
	plot(hist_data, freq=FALSE, add=T, border="black", col=NULL)
	}
dev.off()

