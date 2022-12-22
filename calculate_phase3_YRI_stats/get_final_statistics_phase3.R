#!/bin/env Rscript
#This is the final list of stats:
#chmod +x get_final_statistics.R
#Rscript ./get_final_statistics_phase3.R YRI_50 exon

options(scipen=999)
args = commandArgs(trailingOnly=TRUE)
yri_folder <- args[1]
region <- args[2]

t <- read.table(paste("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/", yri_folder, "/STATS_FILTERED/human_phase3_bigwindow_", region, ".stats", sep=""),h=T)


#To calculate mean and variance of all statistics per window:
t_windows <- unique(t$WinType)
for (win in t_windows) {
        wpii <- t$thetapi[which(t$WinType==win)]/t$WinSize_filtered[which(t$WinType==win)]
        wwi <- t$thetaw[which(t$WinType==win)]/t$WinSize_filtered[which(t$WinType==win)]
        whi <- t$thetah[which(t$WinType==win)]/t$WinSize_filtered[which(t$WinType==win)]
        whpi <- t$hprime[which(t$WinType==win)]
        wtdi <- t$tajimasd[which(t$WinType==win)]
        wsingi <- t$numSing[which(t$WinType==win)]/t$WinSize_filtered[which(t$WinType==win)]
        whapdivi <- t$hapdiv[which(t$WinType==win)]
        wrsqi <- as.numeric(as.character(t$rsq[which(t$WinType==win & t$rsq!="<NA>")]))
        wDi <- t$D[which(t$WinType==win)]
        wDpri <- t$Dprime[which(t$WinType==win)]
        #wdivi <- c()
        #wdivi <- t$div[which(t$WinType==win)]/t$WinSize[which(t$WinType==win)]
        v_win_m <- c(mean(wpii, na.rm=T), mean(wwi, na.rm=T), mean(whi, na.rm=T), mean(whpi, na.rm=T), mean(wtdi, na.rm=T), mean(wsingi, na.rm=T), mean(whapdivi, na.rm=T), mean(wrsqi, na.rm=T), mean(wDi, na.rm=T), mean(wDpri, na.rm=T))
        v_win_sd <- c(sd(wpii, na.rm=T), sd(wwi, na.rm=T), sd(whi, na.rm=T), sd(whpi, na.rm=T), sd(wtdi, na.rm=T), sd(wsingi, na.rm=T), sd(whapdivi, na.rm=T), sd(wrsqi, na.rm=T), sd(wDi, na.rm=T), sd(wDpri, na.rm=T))
        if (win==as.character(t_windows[1])){
                stats <- matrix(c(v_win_m, v_win_sd))
                colnames(stats) <- "win_name"
                rownames(stats) <- c("thetapi_m", "thetaw_m", "thetah_m", "hprime_m", "tajimasd_m", "numSing_m", "hapdiv_m", "rsq_m", "D_m", "Dprime_m", "thetapi_sd", "thetaw_sd", "thetah_sd", "hprime_sd", "tajimasd_sd", "numSing_sd", "hapdiv_sd", "rsq_sd", "D_sd", "Dprime_sd")
                }
        else{

                stats <- cbind(stats, c(v_win_m, v_win_sd))
                }
        }
colnames(stats, do.NULL = FALSE)
colnames(stats) <- t_windows

write.table(stats, file=paste("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/Numbp50/", yri_folder, "/STATS_FILTERED/human_phase3_bigwindow_", region, ".bigwinsummary", sep=""), sep="\t", quote=F)

print ("done")


