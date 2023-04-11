#This is to get the average mutation rate in each region of the single exons:

import sys

def read_fasta(FILE):                                                           
    d_SEQ = {}                                                                  
    for line in FILE:                                                           
        line1 = line.strip('\n')                                                
        if line1[0] == ">":                                                     
            s_seq_name = line1.replace(">", "")                                 
            d_SEQ[s_seq_name] = ""                                              
        else:                                                                   
            d_SEQ[s_seq_name] = d_SEQ[s_seq_name] + line1                       
    return(d_SEQ)

#get the list of exons:
f_exons = open("/home/pjohri1/BgsDfeDemo_Human/Humans/single_exon_coordinates.txt", 'r')
d_col = {}                                                                      
l_exons = []                                                                    
d_chrom, d_coord_start, d_coord_end = {}, {}, {}
l_regions = ["exon", "5p", "3p"]
for region in l_regions:
    d_coord_start[region] = {}
    d_coord_end[region] = {}
for line in f_exons:                                                            
    line1 = line.strip('\n')                                                    
    line2 = line1.split('\t')                                                   
    if line2[0] == "exon":                                                      
        col_num = 0                                                             
        for x in line2:                                                         
            d_col[x] = col_num                                                  
            col_num += 1
    else:
        s_exon = line2[d_col["exon"]]
        if s_exon not in l_exons:
            l_exons.append(s_exon)
        d_chrom[s_exon] = line2[d_col["chrom"]]
        d_coord_start[line2[d_col["region"]]][s_exon] = int(line2[d_col["start"]])
        d_coord_end[line2[d_col["region"]]][s_exon] = int(line2[d_col["end"]])

#Get the nucleotide content of each exon:
l_ntds = ["A", "C", "G", "T"]
d_ntd = {}
for region in l_regions:
    d_ntd[region] = {}
    chr_num = 1
    while chr_num <= 22:
        f_fa = open("/home/pjohri1/BgsDfeDemo_Human/Humans/Exons/divergence/HG19/" + region + "/chr" + str(chr_num) + ".fa", 'r')
        d_fa = read_fasta(f_fa)
        for x in d_fa.keys():
            x1 = x.split("_")
            exon = x1[0] + "_" + x1[1]
            d_ntd[region][exon] = {}
            d_ntd[region][exon]["A"] = float((d_fa[x].count("A")+d_fa[x].count("a"))/len(d_fa[x]))
            d_ntd[region][exon]["C"] = float((d_fa[x].count("C")+d_fa[x].count("c"))/len(d_fa[x]))
            d_ntd[region][exon]["G"] = float((d_fa[x].count("G")+d_fa[x].count("g"))/len(d_fa[x]))
            d_ntd[region][exon]["T"] = float((d_fa[x].count("T")+d_fa[x].count("t"))/len(d_fa[x]))
        chr_num += 1

print (d_ntd["exon"]["exon275_STN1"])
print (d_ntd["5p"]["exon275_STN1"])
print (d_ntd["3p"]["exon275_STN1"])

#get mutation rate and write it:
#read the mutation rate map for each nucleotide:
d_rate = {}
for region in l_regions:
    d_rate[region] = {}
f_mutn = open("/home/pjohri1/BgsDfeDemo_Human/Humans/MutationRateMap/local_mutation_rate.bias_corrected.SEXAVG.bed", 'r')
d_col = {}
for line in f_mutn:
    line1 = line.strip('\n')                                                    
    line2 = line1.split('\t')                                                   
    if line2[0] == "chr":                                                      
        col_num = 0                                                             
        for x in line2:                                                         
            d_col[x] = col_num                                                  
            col_num += 1                                                        
    else:
        for exon in l_exons:
            if d_chrom[exon] == line2[d_col["chr"]]:
                for region in l_regions:
                    win_start = int(line2[d_col["start"]])
                    win_end = int(line2[d_col["end"]])
                    if win_start <= d_coord_start[region][exon] and d_coord_end[region][exon] <= win_end:
                        rateA = d_ntd[region][exon]["A"]*(float(line2[d_col["AtoC"]]) + float(line2[d_col["AtoG"]]) + float(line2[d_col["AtoT"]]))
                        rateC = d_ntd[region][exon]["C"]*(float(line2[d_col["CtoA"]]) + float(line2[d_col["CtoG"]]) + float(line2[d_col["CtoT"]]))
                        rateG = d_ntd[region][exon]["G"]*(float(line2[d_col["GtoA"]]) + float(line2[d_col["GtoC"]]) + float(line2[d_col["GtoT"]]))
                        rateT = d_ntd[region][exon]["T"]*(float(line2[d_col["TtoA"]]) + float(line2[d_col["TtoC"]]) + float(line2[d_col["TtoG"]]))
                        d_rate[region][exon] = str(rateA+rateC+rateG+rateT)
f_mutn.close()

#write it out
result = open("/home/pjohri1/BgsDfeDemo_Human/Humans/MutationRateMap/mutn_rate_single_exons.txt", 'w+')
#mean_rate = 1.25e-8
mean_rate="NA"
result.write("exon" + '\t' + "region" + '\t' + "rate" + '\n')
for exon in l_exons:
    for region in l_regions:
        print (exon + '\t' + region)
        result.write(exon + '\t' + region + '\t' + str(d_rate[region].get(exon, mean_rate)) + '\n')
result.close()

print("done")
