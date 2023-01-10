#This is to sort the BestRefSeq gff file by left coordinates. 

import sys

#read in chromosome names:
f_info = open("/home/pjohri1/BgsDfeDemo_Human/Humans/annotation/GRCh37_info.txt", 'r')
d_col = {} 
l_chrom = []                                                                    
for line in f_info:                                                             
    line1 = line.strip('\n')                                                    
    line2 = line1.split('\t')                                                   
    if line2[0] == "Sequence-Name":                                             
        col = 0                                                                 
        for x in line2:                                                         
            d_col[x] = col                                                      
            col = col + 1                                                       
    else:
        s_chr = line2[d_col["RefSeq-Accn"]]
        if s_chr not in l_chrom:
            l_chrom.append(s_chr)                        
f_info.close()

def get_start(val):
    return(d_start[val])
l_weird_lines = []
num_genes = 0
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/RefGenome_hg19/GRCh37_latest_genomic_sorted.gff", 'w+')
for chrom in l_chrom:
    print(chrom)
    l_lines = []
    d_start = {}
    f_gff = open("/home/pjohri1/BgsDfeDemo_Human/Humans/RefGenome_hg19/GRCh37_latest_genomic.gff", 'r')
    for line in f_gff:
        line1 = line.strip('\n')
        line2 = line1.split('\t')
        if line1[0] != "#":
            if line2[0] == chrom:
                #print ("yes")
                l_lines.append(line1)
                d_start[line1] = int(line2[3])
                if line2[2] == "gene":
                    num_genes += 1
    f_gff.close()
    l_lines.sort(key=get_start)
    for x in l_lines:
        result.write(x + '\n')
#write weird lines:
f_gff = open("/home/pjohri1/BgsDfeDemo_Human/Humans/RefGenome_hg19/GRCh37_latest_genomic.gff", 'r')
for line in f_gff:
    line1 = line.strip('\n')
    if line1[0] != "#":
        line2 = line1.split('\t')
        if line2[0] not in l_chrom:
            l_weird_lines.append(line1)
f_gff.close()
for x in l_weird_lines:
    result.write(x + '\n')
result.close()
print (num_genes)
print("done")

