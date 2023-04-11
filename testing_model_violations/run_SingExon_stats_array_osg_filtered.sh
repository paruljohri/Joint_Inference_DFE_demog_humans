#!/bin/bash
#SBATCH --mail-user=pjohri1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1 #number of tasks
#SBATCH --time=0-1:00
##SBATCH --mem=5000
#SBATCH -a 1-1%1
#SBATCH -o /home/pjohri1/LOGFILES/model_violations_stats_%A_%a.out
#SBATCH -e /home/pjohri1/LOGFILES/model_violations_stats_%A_%a.err


echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
######################################
#This script is for performing detailed 
#statistics for simualtions of single-exon 
#with diff recombination rate etc.
#To be run on AGAVE!
#######################################
declare -i startID=4
declare -i simID=${startID}+$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 #This is for pylibseq to work correctly. It might change later!
module load r/2.15.3

#The region you want to use to filter:
intergenic="5p" #5p/3p

#unzip the simulation folder:
cd /scratch/pjohri1/BgsDfeDemo_Human/ModelViolations
unzip sim${simID}.zip

#filter the .ms and .fixed files:
cd /home/pjohri1/BgsDfeDemo_Human/programs_model_violations
python filter_ms_phase3_data.py ${simID} ${intergenic}

#get stats:
cd /home/pjohri1/BgsDfeDemo_Human/programs_model_violations

if [ ! -f "/scratch/pjohri1/BgsDfeDemo_Human/ModelViolations_stats_filtered_${intergenic}/sim"$simID"_bigwindow.stats" ]
then
    echo "Calculating stats for neutral, linked, and functional regions"
    python statistics_bigwindow_pylibseq_SingExon_osg_human_filtered.py -folder ModelViolations -outFolder ModelViolations_stats_filtered_${intergenic} -interRegion ${intergenic} -simID ${simID} -numGenes 465 -numRep 1
fi

cd /home/pjohri1/BgsDfeDemo_Human/programsABC
Rscript ./get_final_statistics_SingExon_filtered.R ModelViolations_stats_filtered_${intergenic} ${simID}

#remove the unzipped folder:
cd /scratch/pjohri1/BgsDfeDemo_Human/ModelViolations
if [ -f "sim"$simID".zip" ]
then
    rm -r sim${simID}
fi
echo "Done"



