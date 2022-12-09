#!/bin/bash
#SBATCH -n 1                        # number of cores
#SBATCH -t 0-15:00                  # wall time (D-HH:MM)
#SBATCH -o /home/pjohri1/LOGFILES/geneconv_humans_%j.out
#SBATCH -e /home/pjohri1/LOGFILES/geneconv_humans_%j.err
#SBATCH --mail-type=ALL             # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=pjohri1@asu.edu # send-to address

#Set path to working directory
cd /home/pjohri1/BgsDfeDemo_Human/programs_gene_conv

#Set the environment
module load slim/3.1

#rec_rate="1"
#simID=2 #2 #no gene conv, no sweeps, rec=1
#simID=11 #11 #no gene conv, sweeps, rec=1

rec_rate="10"
#simID=3 #3 #no gene conv, no sweeps, rec=10
simID=12 #12 #no gene conv, sweeps, rec=10

declare -i repEnd=1001
declare -i repID=1

while [ $repID -lt $repEnd ];                                                   
do                                                                              
    echo "rep"$repID
    mkdir /scratch/pjohri1/BgsDfeDemo_Human/geneconv_humans/sim${simID}
    ##With no sweeps:
    #2/3
    #slim -d d_seed=${repID} -d d_f0=0.22 -d d_f1=0.27 -d d_f2=0.13 -d d_f3=0.38 -d d_f_pos=0.0 -d d_rec_rate=${rec_rate} -d "d_simID='${simID}'" -d "d_repID='${repID}'" eqm_disc_5_SingExon_nogeneconv.slim
    ##With sweeps:
    slim -d d_seed=${repID} -d d_f0=0.22 -d d_f1=0.27 -d d_f2=0.13 -d d_f3=0.38 -d d_f_pos=0.01 -d d_rec_rate=${rec_rate} -d "d_simID='${simID}'" -d "d_repID='${repID}'" eqm_disc_5_SingExon_nogeneconv.slim
    repID=$(( ${repID} + 1 ))
done
echo "Finished"

