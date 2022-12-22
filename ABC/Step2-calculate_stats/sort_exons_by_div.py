#This is to get values of divergence and recombination rates for every exon:
#Seprate recombination rates into min < rec <= 2.35 and 2.35 < rec< max
#divergence: 0 < 0.0037 < max

import sys

#read the exon names:
f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/single_exon_coordinates.txt", 'r')
l_exons = []
for line in f_exons:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] != "exon":
        if line2[0] not in l_exons:
            l_exons.append(line2[0])
f_exons.close()

#read the functional divergence per site:
d_div = {}
f_div = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/YRI_50/divergence_exon_filtered.txt", 'r')
for line in f_div:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if line2[0] != "exon":
        s_exon = line2[0].replace("_exon", "")
        d_div[s_exon] = line2[1]
f_div.close()


#write it in 2 files:
result1 = open("/home/pjohri1/BgsDfeDemo_Human/Humans/single_exon_low_div.txt", 'w+')
result2 = open("/home/pjohri1/BgsDfeDemo_Human/Humans/single_exon_high_div.txt", 'w+')
t_files = (result1, result2)
for result in t_files:
    result.write("exon" + '\t' + "gene_num" + '\t' + "func_div" + '\n')
for s_exon in l_exons:
    s_gene = s_exon.split("_")[0].replace("exon", "gene")
    if d_div[s_exon] != "NA":
        if float(d_div[s_exon]) <= 0.0037: #low div
            result1.write(s_exon + '\t' + s_gene + '\t' + d_div[s_exon] + '\n')
        else:
            result2.write(s_exon + '\t' + s_gene + '\t' + d_div[s_exon] + '\n')

for result in t_files:
    result.close()
print ("done")

