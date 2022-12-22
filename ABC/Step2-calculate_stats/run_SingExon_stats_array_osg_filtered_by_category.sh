#!/bin/bash
#SBATCH --mail-user=pjohri1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1 #number of tasks
#SBATCH --time=0-1:00
##SBATCH --mem=5000
#SBATCH -a 1-100%100
#SBATCH -o /scratch/pjohri1/BgsDfeDemo_Human/demo_dfe_SingExon_human_v2_logfiles/sim_%A_%a.out
#SBATCH -e /scratch/pjohri1/BgsDfeDemo_Human/demo_dfe_SingExon_human_v2_logfiles/sim_%A_%a.err


echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
######################################
#This script is for performing detailed 
#statistics for simualtions of single-exon 
#with diff recombination rate etc.
#To be run on AGAVE!
#######################################
declare -i startID=1700
declare -i simID=${startID}+$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 #This is for pylibseq to work correctly. It might change later!
module load r/2.15.3

#The region you want to use to filter:
intergenic="3p" #5p/3p

#unzip the simulation folder:
cd /scratch/pjohri1/BgsDfeDemo_Human/demo_dfe_SingExon_human_v2
unzip sim${simID}.zip

#filter the .ms and .fixed files:
cd /home/pjohri1/BgsDfeDemo_Human/programsABC
python filter_ms_phase3_data.py ${simID} ${intergenic}

#get stats:
cd /home/pjohri1/BgsDfeDemo_Human/programsABC

#categories=( "low_rec_low_div" "low_rec_high_div" "high_rec_low_div" "high_rec_high_div")
categories=( "low_div" "high_div")
for category in "${categories[@]}"
do
    echo ${category}
    mkdir /scratch/pjohri1/BgsDfeDemo_Human/demo_dfe_SingExon_human_v2_stats_filtered_${intergenic}_${category}
    
    if [ ! -f "/scratch/pjohri1/BgsDfeDemo_Human/demo_dfe_SingExon_human_v2_stats_filtered_${intergenic}_${category}/sim"$simID"_bigwindow.stats" ]
    then
        echo "Calculating stats for neutral, linked, and functional regions"
        python statistics_bigwindow_pylibseq_SingExon_osg_human_filtered_by_category.py -folder demo_dfe_SingExon_human_v2 -outFolder demo_dfe_SingExon_human_v2_stats_filtered_${intergenic}_${category} -interRegion ${intergenic} -simID ${simID} -numRep 1 -category ${category}
    fi

    Rscript ./get_final_statistics_SingExon_filtered.R demo_dfe_SingExon_human_v2_stats_filtered_${intergenic}_${category} ${simID}
done

#remove the unzipped folder:
cd /scratch/pjohri1/BgsDfeDemo_Human/demo_dfe_SingExon_human_v2
if [ -f "sim"$simID".zip" ]
then
    rm -r sim${simID}
fi
echo "Done"



