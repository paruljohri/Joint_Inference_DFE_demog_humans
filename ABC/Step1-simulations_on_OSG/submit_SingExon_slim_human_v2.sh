#!/bin/sh
#
sim_start=$1
sim_end=$2

echo "sim start "$sim_start
echo "sim_end "$sim_end

declare -i simulationID=$sim_start
declare -i simulationEnd=$sim_end
simulationEnd=$(( ${simulationEnd} + 1 ))

while [ $simulationID -lt $simulationEnd ];
do
	echo "creating submit files for simulation " $simulationID
	python write_slim_scripts_SingExon_human_v2.py -numRep 1 -simID ${simulationID}
	simulationID=$(( ${simulationID} + 1 ))
done

echo "submitting files"
cd /home/pjohri/SUBMITFILES

declare -i simulationID=$sim_start
while [ $simulationID -lt ${simulationEnd} ];
do
        echo "submitting submit files for simulation " $simulationID
        condor_submit SingExon_slim_sim${simulationID}.submit
        simulationID=$(( ${simulationID} + 1 ))
done
echo "done"
