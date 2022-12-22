#!/bin/bash
#SBATCH --mail-user=pjohri1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1 #number of tasks
#SBATCH --time=0-5:00
#SBATCH -a 1-1%1
#SBATCH -o /scratch/pjohri1/BgsDfeDemo_Human/LOGFILES/gene_%A_rep_%a.out
#SBATCH -e /scratch/pjohri1/BgsDfeDemo_Human/LOGFILES/gene_%A_rep_%a.err

#module load perl/5.22.1
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
######################################
#This script is for performing detailed 
#simulations using differnet exon size, 
#recombination rate etc.
#To be run on Agave!
#Ten simulation per submitted job!
#######################################

simID=$1
geneID=$2
declare -i repID=0+$SLURM_ARRAY_TASK_ID
seed=${geneID}${repID}

module load slim/3.1

#run the simulations
echo "Loading Slim"
cd /scratch/pjohri1/BgsDfeDemo_Human/BASHFILES/

echo "Running Slim"
if [ ! -f "/scratch/pjohri1/BgsDfeDemo_Human/model_fit/sim"${simID}"_gene"${geneID}"_rep"${repID}".ms" ]
then
    bash sim${simID}_gene${geneID}.sh ${repID} ${seed}
fi





