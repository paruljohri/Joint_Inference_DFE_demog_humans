#!/bin/bash

simID=$1

cd /scratch/pjohri1/BgsDfeDemo_Human/geneconv_humans
unzip sim${simID}.zip

cd /home/pjohri1/BgsDfeDemo_Human/programs_gene_conv
module load gcc/6.3.0
python statistics_slidingwindow_pylibseq.py -winSize 500 -stepSize 500 -folder geneconv_humans -simID ${simID} -noncodingLen 50000 -codingLen 4000
python statistics_slidingwindow_pylibseq.py -winSize 1000 -stepSize 1000 -folder geneconv_humans -simID ${simID} -noncodingLen 50000 -codingLen 4000

cd /scratch/pjohri1/BgsDfeDemo_Human/geneconv_humans
if [ -f "sim"$simID".zip" ]                                                     
then                                                                            
    rm -r sim${simID}                                                           
fi                                                                              
echo "Done"

