#!/bin/bash
#SBATCH --mail-user=pjohri1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1 #number of tasks
#SBATCH --time=0-20:00
#SBATCH -o /home/pjohri1/LOGFILES/get_vcfs_numbp50_%j.out
#SBATCH -e /home/pjohri1/LOGFILES/get_vcfs_numbp50_%j.err


######################################
#This script is for getting VCF files 
#for each exon, for only yri individuals 
#
#To be run on AGAVE!
#######################################

yriFolder="YRI" #YRI/YRI_50
cd /home/pjohri1/BgsDfeDemo_Human/programs_stats
declare -i chrID=1                                                              
while [ $chrID -lt 23 ];                                                      
do                                                                              
    echo "starting chromosome " $chrID
    python get_exon_vcfs_numbp50.py ${yriFolder} ${chrID}
    echo "Finished chromosome " $chrID
    chrID=$(( $chrID + 1 ))                                    
done







