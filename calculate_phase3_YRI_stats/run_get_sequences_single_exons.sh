#!/bin/bash

declare -i chrID=1 
while [ ${chrID} -lt 23 ]; #should be 23                                                       
do                                                                              
    echo "starting chromosome " $chrID
    
    ##ancestral sequences:                                        
    python get_ancestral_sequences_single_exons.py "chr"${chrID}
    
    ##hg19 sequences:
    #cd /scratch/pjohri1/BgsDfeDemo_Human/HumanData/RefGenome/hg19_ucsc/
    #gzip -d chr${chrID}.fa.gz
    #cd /home/pjohri1/BgsDfeDemo_Human/programs_stats
    #python get_hg19_sequences_single_exons.py "chr"${chrID}
    #cd /scratch/pjohri1/BgsDfeDemo_Human/HumanData/RefGenome/hg19_ucsc/
    #gzip chr${chrID}.fa
    
    ##grch37 sequences:
    #python get_grch37_sequences_single_exons.py "chr"${chrID}

    echo "Finished chromosome " $chrID                                          
    chrID=$(( ${chrID} + 1 ))                                                     
done
