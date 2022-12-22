#!/bin/bash

sim_start=10 #$1
sim_end=10 #$2

declare -i simID=$sim_start
declare -i simEnd=$sim_end
simEnd=$(( ${simEnd} + 1 ))

#make a new folder for this simulation ID
cd /scratch/pjohri1/BgsDfeDemo_Human/model_fit

while [ $simID -lt $simEnd ];
do
	echo "sim"$simID
	if [ ! -f "/scratch/pjohri1/BgsDfeDemo_Human/model_fit/sim"$simID".zip" ]
	then
		mkdir /scratch/pjohri1/BgsDfeDemo_Human/model_fit/sim$simID #remove this line if you want this script to only re-run parameter combinations.
		if [ -d "/scratch/pjohri1/BgsDfeDemo_Human/model_fit/sim${simID}" ]
		then
        		python /home/pjohri1/BgsDfeDemo_Human/programs_model_fit/write_slim_scripts_SingExon_human_model_fit_agave.py -numRep 1 -simID $simID
        		bash /scratch/pjohri1/BgsDfeDemo_Human/BASHFILES/sim${simID}_submit.sh
		fi
	fi
	simID=$(( ${simID} + 1 ))
done
echo "Finished"

