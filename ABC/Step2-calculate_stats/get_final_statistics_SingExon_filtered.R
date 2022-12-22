#!/bin/env Rscript
#This is the final list of stats:
#chmod +x get_final_statistics.R
#Rscript ./get_final_statistics.R demo_dfe_SingExon_human_v2_stats 1

#options("winSize" = ()$win_size)
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)
folder <- args[1]
simID <- args[2]

t <- read.table(paste("/scratch/pjohri1/BgsDfeDemo_Human/", folder, "/sim", simID, "_bigwindow.stats", sep=""),h=T)


#To calculate mean and variance of all statistics per window:
t_windows <- c("func", "linked1", "linked2")
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
        wdivi <- t$div[which(t$WinType==win)]/t$WinSize_filtered[which(t$WinType==win)]
        v_win_m <- c(mean(wpii[is.finite(wpii)], na.rm=T), mean(wwi[is.finite(wwi)], na.rm=T), mean(whi[is.finite(whi)], na.rm=T), mean(whpi[is.finite(whpi)], na.rm=T), mean(wtdi[is.finite(wtdi)], na.rm=T), mean(wsingi[is.finite(wsingi)], na.rm=T), mean(whapdivi[is.finite(whapdivi)], na.rm=T), mean(wrsqi[is.finite(wrsqi)], na.rm=T), mean(wDi[is.finite(wDi)], na.rm=T), mean(wDpri[is.finite(wDpri)], na.rm=T), mean(wdivi[is.finite(wdivi)], na.rm=T))
        v_win_sd <- c(sd(wpii[is.finite(wpii)], na.rm=T), sd(wwi[is.finite(wwi)], na.rm=T), sd(whi[is.finite(whi)], na.rm=T), sd(whpi[is.finite(whpi)], na.rm=T), sd(wtdi[is.finite(wtdi)], na.rm=T), sd(wsingi[is.finite(wsingi)], na.rm=T), sd(whapdivi[is.finite(whapdivi)], na.rm=T), sd(wrsqi[is.finite(wrsqi)], na.rm=T), sd(wDi[is.finite(wDi)], na.rm=T), sd(wDpri[is.finite(wDpri)], na.rm=T), sd(wdivi[is.finite(wdivi)], na.rm=T))
        if (win=="func"){
                stats <- matrix(c(v_win_m, v_win_sd))
                colnames(stats) <- "win_name"
                rownames(stats) <- c("thetapi_m", "thetaw_m", "thetah_m", "hprime_m", "tajimasd_m", "numSing_m", "hapdiv_m", "rsq_m", "D_m", "Dprime_m", "div_m", "thetapi_sd", "thetaw_sd", "thetah_sd", "hprime_sd", "tajimasd_sd", "numSing_sd", "hapdiv_sd", "rsq_sd", "D_sd", "Dprime_sd", "div_sd")
                }
        else{

                stats <- cbind(stats, c(v_win_m, v_win_sd))
                }
        }
colnames(stats, do.NULL = FALSE)
colnames(stats) <- t_windows

write.table(stats, file=paste("/scratch/pjohri1/BgsDfeDemo_Human/", folder, "/sim", simID, ".bigwinsummary", sep=""), sep="\t", quote=F)

print ("done")


