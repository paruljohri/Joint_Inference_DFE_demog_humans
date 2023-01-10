#This script is simply to take divergence values from the Refseq exons and add it to another file:
#Note this is a very rough estiamte of divergence which was simply used inititally as a check. 
#The values of divergence used for inference were carefully calculated later and the scripts are here - https://github.com/paruljohri/Joint_Inference_DFE_demog_humans/tree/main/calculate_phase3_YRI_stats

import sys

#read divergence values for exons:
d_div = {}
d_col = {}
f_div = open("/home/pjohri1/BgsDfeDemo_Human/Humans/divergence/GRCh37_exons_divergence.txt", 'r')
for line in f_div:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "chrom_ucsc":
        col = 0
        for x in line2:
            d_col[x] = col
            col = col + 1
    else:
        exon = line2[d_col["chrom_refseq"]] + "-" + line2[d_col["start"]] + "-" + line2[d_col["end"]]
        d_div[exon] = line2[d_col["soft_div"]]
f_div.close()

#Write it in this file:
#f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon1-6kb_CNE500_rec_numbp50_final.tab", 'r')
f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final.tab", 'r')
#result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon1-6kb_CNE500_rec_numbp50_final_div.tab", 'w+')
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_exons_firstelement_noTEs_500bp/GRCh37_exon2-6kb_CNE500_rec_numbp50_final_div.tab", 'w+')

d_col = {}
for line in f_exons:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] == "exon":
        result.write(line1 + '\t' + "divergence" + '\n')
        col = 0
        for x in line2:
            d_col[x] = col
            col = col + 1
    else:
        exon = line2[d_col["exon"]]
        result.write(line1 + '\t' + d_div.get(exon, "NA") + '\n')
result.close()
f_exons.close()
print ("done")



